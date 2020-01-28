
moment_match <- function(mean, var){
  list(alpha = mean^2/var,
       beta  = mean/var) 
}

time_to_event <- function(n, mean, var){
  if (var > 0){
    parms <- moment_match(mean, var)
    return(rgamma(n, shape = parms$alpha, rate = parms$beta))
  } else{
    return(rep(mean, n))
  }
}

generate_histories <- function(input){
  
  with(input,
       data.frame(
         incu = time_to_event(n = sims, mean =  mu_inc, var =  sigma_inc),
         inf  = time_to_event(sims, mu_inf, sigma_inf)) %>%
         mutate(flight.departure = runif(sims, min = 0, max =2*(incu + inf)),
                flight.arrival   = flight.departure + dur_flight) )
  
}


calc_probs <- function(dur.flight,
                       mu_inc,
                       sigma_inc,
                       mu_inf,
                       sigma_inf,
                       sens.exit,
                       sens.entry,
                       prop.asy,
                       sims){
  
  
  #convert flight time to days
  dur_flight = dur.flight / 24
  
  infection_histories <- generate_histories(list(dur.flight = dur.flight,
                                                 mu_inc     = mu_inc,
                                                 sigma_inc  = sigma_inc,
                                                 mu_inf     = mu_inf,
                                                 sigma_inf  = sigma_inf,
                                                 sens.exit  = sens.exit,
                                                 sens.entry = sens.entry,
                                                 prop.asy   = prop.asy,
                                                 sims       = sims,
                                                 dur_flight = dur_flight))
  
  infection_histories %>%
    mutate(hospitalised_prior_to_departure = inf + incu < flight.departure) %>%
    filter(hospitalised_prior_to_departure == FALSE) %>%
    mutate(exit_screening_label  = runif(n(), 0, 1) < sens.exit /100,
           entry_screening_label = runif(n(), 0, 1) < sens.entry/100) %>%
    
    mutate(symp_at_exit   = incu < flight.departure,
           symp_at_entry  = incu < flight.arrival,
           found_at_exit  =  symp_at_exit &  exit_screening_label,
           missed_at_exit = symp_at_exit &  !exit_screening_label,
           found_at_entry = symp_at_entry & entry_screening_label ,
           sev_at_exit    = 0, # no hospitalised can exit country
           sev_from_lat   = (!symp_at_exit) & 
             (incu + inf < flight.arrival),
           sev_from_symp  = symp_at_exit & (!exit_screening_label) &
             (incu + inf < flight.arrival),
           sev_at_entry   = sev_from_lat | sev_from_symp,
           found_at_entry_only = found_at_entry & (!symp_at_exit)
    ) %>%
    summarise(
      
      prop_sev_at_entry = (1-prop.asy/100)*mean(sev_at_entry),
      
      prop_symp_at_exit = (1-prop.asy/100)*mean(found_at_exit),
      
      prop_symp_at_entry = (1-prop.asy/100)*mean(
        (missed_at_exit & found_at_entry & !sev_at_entry) |
          (found_at_entry_only & !sev_at_entry))
    ) %>%
    mutate(prop_undetected = 1 - (prop_sev_at_entry +
                                    prop_symp_at_exit +
                                    prop_symp_at_entry)) %>%
    as.list %>%
    return
  
  # return(list("prob_det_symp_exit"       = prob_det_symp_exit,
  #             "prob_det_symp_entry_only" = prob_det_symp_entry_only,
  #             "prob_det_symp_entry"      = prob_det_symp_entry,
  #             "prob_det_severe_flight"   = prob_det_severe_flight,
  #             "prob_undet"               = prob_undet))
  
}

pathogen <- list(`SARS-like` = data.frame(mu_inc    =  6.4,
                                          sigma_inc = 16.7,
                                          mu_inf    =  3.8,
                                          sigma_inf =  6.0),
                 Custom     = data.frame(mu_inc    =  5.0,
                                         sigma_inc =  5.0,
                                         mu_inf    =  5.0,
                                         sigma_inf =  5.0)) %>%
  dplyr::bind_rows(., .id = "name")


make_ci_label <- function(x){
  return(sprintf("%0.2f (%0.2f, %0.2f)", x[1], x[2], x[3]))
}


# function to generate travellers to have screening applied
generate_travellers <- function(input, i){
  tibble(i = i) %>% 
    mutate(probs=future_pmap(.f=calc_probs,list(   
      dur.flight = input$dur.flight,
      mu_inc     = input$mu_inc,
      sigma_inc  = input$sigma_inc,
      mu_inf     = input$mu_inf,
      sigma_inf  = input$sigma_inf,
      sens.exit  = input$sens.exit,
      sens.entry = input$sens.entry,
      prop.asy   = input$prop.asy,
      sims       = i))) %>% 
    unnest_wider(probs)
}

# function to take travellers and work out their detection probabilities
generate_probabilities <- function(travellers){
  travellers %>%
    pivot_longer(cols = c(prop_symp_at_exit,
                          prop_symp_at_entry,
                          prop_sev_at_entry,
                          prop_undetected),
                 names_to = "screening",
                 values_to = "prob") %>% 
    group_by(screening) %>% 
    summarise(mean_prob = mean(prob*100),
              lb_prob=quantile(probs=0.025,x=prob*100),
              ub_prob=quantile(probs=0.975,x=prob*100)) %>%  
    pivot_longer(cols=c(mean_prob,
                        lb_prob,ub_prob)) %>% 
    pivot_wider(names_from = screening, values_from = value)
}

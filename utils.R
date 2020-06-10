
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
         lazy_dt() %>% 
         mutate(flight.departure = runif(sims, min = 0, max =2*(incu + inf)),
                flight.arrival   = flight.departure + dur_flight) %>% 
         as.data.frame())
  
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
    lazy_dt() %>% 
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
    as.data.frame() %>% 
    as.list %>%
    return
}

pathogen <- list(
  `nCoV-2019` = 
    data.frame(
      # (Li et al. (2020) NEJM)
      mu_inc    =  5.2,
      sigma_inc =  4.1,
      mu_inf    =  9.1,
      sigma_inf = 14.7,
      prop.asy  = 17
    ),
  # `nCoV-2019 (Backer)` = 
  #   data.frame(
  #     # (Backer et al., 2020)
  #     mu_inc    =  5.7,
  #     sigma_inc =  round(2.6^2,1),
  #     # (Huang et al., 2020)
  #     mu_inf    = 8.0,
  #     sigma_inf = round(((13-5)/1.35)^2,1), # using Higgins (2008) Cochrane Handbook
  #     prop.asy  = 17
  #   ),
  `SARS-like (2002)` = 
    data.frame(
      mu_inc    =  6.4,
      sigma_inc = 16.7,
      mu_inf    =  3.8,
      sigma_inf =  6.0,
      prop.asy  =  0.0),
  `Flu A/H1N1-like (2009)` =
    data.frame(
      mu_inc    =  4.3,
      sigma_inc =  1.05,
      mu_inf    =  9.3,
      sigma_inf =  0.7,
      prop.asy  = 16.0 # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4586318/ 
    ),
  `MERS-like (2012)` = 
    data.frame(
      mu_inc    =  5.5,
      sigma_inc = 6.25, # https://www.sciencedirect.com/science/article/pii/S1473309913703049?via%3Dihub#sec1
      mu_inf    =  5.0,  # https://www.nejm.org/doi/10.1056/NEJMoa1306742
      sigma_inf =  7.5,
      prop.asy  = 21.0  # https://doi.org/10.1016/j.tmaid.2018.12.003 citing https://www.who.int/csr/disease/coronavirus_infections/risk-assessment-august-2018.pdf?ua=1
    ),
  Custom     = 
    data.frame(
      mu_inc    =  5.0,
      sigma_inc =  5.0,
      mu_inf    =  5.0,
      sigma_inf =  5.0,
      prop.asy  = 10.0)) %>%
  dplyr::bind_rows(., .id = "name")


make_ci_label <- function(x){
  x <- round(x)
  return(sprintf("%i (%i, %i)", x[1], x[2], x[3]))
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
    dt_pivot_longer(cols = c(prop_symp_at_exit,
                          prop_symp_at_entry,
                          prop_sev_at_entry,
                          prop_undetected),
                 names_to = "screening",
                 values_to = "prob") %>% 
    group_by(screening) %>% 
    summarise(mean_prob = mean(prob*100),
              lb_prob=quantile(probs=0.025,x=prob*100),
              ub_prob=quantile(probs=0.975,x=prob*100)) %>%  
    dt_pivot_longer(cols=c(mean_prob,
                        lb_prob,ub_prob)) %>% 
    dt_pivot_wider(names_from = screening, values_from = value) %>% 
    as.data.frame()
}


get_var <- function(lower, upper, n){
  
  # back-calculate variance from the lower and upper bands of a symmetric
  # confidence interval of the mean
  
  n*((lower - upper)/(2*qt(0.975, n-1)))^2
  
}

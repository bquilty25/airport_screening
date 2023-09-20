#' Convert mean and variance to shape and rate for a gamma distribution
#' 
#'
#' @inheritParams time_to_event
#' @keywords internal
#' @return A named list, where `alpha` is the `shape` of the gamma distribution,
#' while `beta` is the `rate`.
moment_match <- function(mean, var) {
  checkmate::assert_number(
    mean
  )
  checkmate::assert_number(
    var,
    lower = 0.0
  )
  list(
    alpha = mean^2 / var,
    beta = mean / var
  )
}

#' Draw waiting times until an event occurs
#'
#' @description Draws a waiting time from either a gamma distribution
#' ([rgamma()]), or provides a repeating sequence of the mean provided, if the
#' variance is 0.0.
#' @param n The number of waiting times to draw.
#' @param mean The mean of the gamma distribution from which to draw waiting
#' times.
#' @param var The variance of gamma distribution from which to draw waiting
#' times.
#' @keywords internal
#' @return A vector of length `n`, of the waiting times, or of the provided
#' mean repeated \eqn{n} times.
time_to_event <- function(n, mean, var) {
  if (var > 0.0) {
    parms <- moment_match(mean, var)
    return(stats::rgamma(n, shape = parms$alpha, rate = parms$beta))
  } else {
    return(rep(mean, n))
  }
}

#' Simulate travel and infection histories from user input data
#'
#' @inheritParams calc_probs
#' @importFrom rlang .data
#' @keywords internal
#' @return A data.frame of travel and infection outcomes.


generate_histories <- function(dur.flight, mu_inc, sigma_inc,
                               mu_inf, sigma_inf, sens.exit, prop_fever, prop_relevant,
                               sens.entry, prop.asy, sims, n_travellers) {
 
  # Generate infection status for each individual (0 = not infected, 1 = infected)
  data.frame(i=1:n_travellers) %>% 
    rowwise() %>% 
    mutate(fever_status = rbinom(1, size = 1, prob = prop_fever > runif(n=n())),
           relevant_infection_status = ifelse(fever_status == 1, 
                                             rbinom(1, size = 1, prob = prop_relevant > runif(n=n())), 
                                             0)) %>% 
    mutate(
      # Generate incubation times for each individual
      incu = time_to_event(n = n(), mean = mu_inc, var = sigma_inc),
      inf = time_to_event(n(), mu_inf, sigma_inf),
      
      # Generate flight departure and arrival times
      flight.departure = stats::runif(n(), min = 0, max = 2 * (.data$incu + .data$inf)),
      flight.arrival = .data$flight.departure + dur.flight
    )
}



#' Calculate the probabilities of travel and infection outcomes
#'
#' @param dur.flight The flight duration in hours.
#' @param mu_inc Mean incubation period in days.
#' @param sigma_inc Variance of incubation period in days.
#' @param mu_inf Mean time to symptom onset in days.
#' @param sigma_inf Variance in time to symptom onset, in days.
#' @param sens.exit Sensitivity of tests used upon departure.
#' @param sens.entry Sensitivity of tests used upon arrival.
#' @param prop.asy Proportion of asymptomatic infections.
#' @param sims Number of simulation runs.
#'
#' @importFrom rlang .data
#' @keywords internal
#' @return A data.frame with probabilities of different travel and infection
#' outcomes.
calc_probs <- function(dur.flight, mu_inc, sigma_inc,
                       mu_inf, sigma_inf, sens.exit, prop_fever, prop_relevant,
                       sens.entry, prop.asy, sims, n_travellers) {
  
  # simulate infection histories
  .args <- as.list(match.call())[-1] # remove fn call
  
  # convert flight time to days
  .args$dur.flight <- .args$dur.flight / 24.0
  infection_histories <- do.call(generate_histories, .args)
  #browser()
  # simulate probabilities of different infection and travel related events
  infection_histories <- infection_histories %>%
    dplyr::mutate(
      exit_screening_label = stats::runif(dplyr::n(), 0, 1) < sens.exit,
      entry_screening_label = stats::runif(dplyr::n(), 0, 1) < sens.entry
    )
  
  # simulate different outcomes related to detection during travel
  infection_histories <-
    dplyr::mutate(
      infection_histories,
      symp_at_exit = .data$incu < .data$flight.departure,
      symp_at_entry = .data$incu < .data$flight.arrival,
      
      symp_fever_relevant_at_exit = .data$symp_at_exit & .data$exit_screening_label & .data$fever_status == 1 & .data$relevant_infection_status == 1,
      symp_fever_relevant_at_entry = .data$symp_at_entry & .data$entry_screening_label & .data$fever_status == 1 & .data$relevant_infection_status == 1,
      
      symp_fever_irrelevant_at_exit = .data$symp_at_exit & .data$exit_screening_label & .data$fever_status == 1 & .data$relevant_infection_status == 0,
      symp_fever_irrelevant_at_entry = .data$symp_at_entry & .data$entry_screening_label & .data$fever_status == 1 & .data$relevant_infection_status == 0,
      
      found_at_exit_relevant = .data$symp_fever_relevant_at_exit & .data$exit_screening_label,
      found_at_exit_irrelevant = .data$symp_fever_irrelevant_at_exit & .data$exit_screening_label,
      
      missed_at_exit_relevant = .data$symp_fever_relevant_at_exit & !.data$exit_screening_label,
      missed_at_exit_irrelevant = .data$symp_fever_irrelevant_at_exit & !.data$exit_screening_label,
      
      found_at_entry_relevant = .data$symp_fever_relevant_at_entry & .data$entry_screening_label,
      found_at_entry_irrelevant = .data$symp_fever_irrelevant_at_entry & .data$entry_screening_label,
      
      found_at_entry_only_relevant = .data$found_at_entry_relevant & (!.data$symp_fever_relevant_at_exit),
      found_at_entry_only_irrelevant = .data$found_at_entry_irrelevant & (!.data$symp_fever_irrelevant_at_exit)
    )
  
  # summaries detection outcomes
  
  infection_histories_summary <-
    dplyr::summarise(
      infection_histories,
      prop_symp_at_exit_relevant = (1.0 - prop.asy/100) * mean(.data$found_at_exit_relevant),
      prop_symp_at_exit_irrelevant = (1.0 - prop.asy/100) * mean(.data$found_at_exit_irrelevant),
      
      prop_symp_at_entry_relevant = (1.0 - prop.asy/100) * mean(
        (.data$missed_at_exit_relevant & .data$found_at_entry_relevant) |
          (.data$found_at_entry_only_relevant)
      ),
      prop_symp_at_entry_irrelevant = (1.0 - prop.asy/100) * mean(
        (.data$missed_at_exit_irrelevant & .data$found_at_entry_irrelevant) |
          (.data$found_at_entry_only_irrelevant)
      )
    ) %>%
    dplyr::mutate(prop_undetected_relevant = 1.0 - (.data$prop_symp_at_exit_relevant +
                                                      .data$prop_symp_at_entry_relevant))
  
  
  
  
  # return dataframe converted to list object
  return(
    as.list(infection_histories_summary)
  )
}

#' Make confidence interval labels
#'
#' @param x A vector of the central estimate, and the lower and upper confidence
#' limits.
#'
#' @keywords internal
#' @return A string vector of the format `"I (I, I)"``.
make_ci_label <- function(x) {
  x <- round(x)
  return(sprintf("%i (%i, %i)", x[1], x[2], x[3]))
}

#' Generate travellers to have screening applied
#'
#' @param input Input from the Shiny app, giving the duration of the flight,
#' the mean incubation period, the variance of the incubation period, the mean
#' time between infection and symptom onset, the variance of the time to symptom
#' onset, the sensitivity of testing upon departure, the sensitivity of testing
#' upon arrival, and the proportion of asymptomatic infections.
#' @param i The number of simulation runs.
#'
#' @keywords internal
#' @return A data.frame giving the proportion of travellers who are symptomatic
#' upon arrival and departure, given the pathogen parameters and flight duration
#' and the proportions that have severe infections upon arrival, and also the
#' proportion which is infected but undetected upon arrival.

generate_travellers <- function(input, i) {
  as.data.frame(
    do.call(
      calc_probs,
      list(
        input$dur.flight, input$mu_inc, input$sigma_inc,
        input$mu_inf, input$sigma_inf, input$sens.exit,
        input$sens.entry, input$prop.asy, input$prop_fever, input$prop_relevant,input$n_travellers,
        sims = i
      )
    )
  )
}

#' Work out the detection probabilities of travellers
#'
#' @param travellers Output from the [generate_travellers()] function.
#'
#' @importFrom rlang .data
#' @keywords internal
#' @return A data.frame giving the probabilities of travellers who are infected
#' being detected as such at different stages of airline travel.

generate_probabilities <- function(travellers) {
  travellers %>%
    tidyr::pivot_longer(
      cols = c(
        .data$prop_symp_at_exit_relevant,
        .data$prop_symp_at_exit_irrelevant,
        
        .data$prop_symp_at_entry_relevant,
        
        .data$prop_symp_at_entry_irrelevant,
        .data$prop_undetected_relevant
      ),
      names_to = "screening",
      values_to = "prob"
    ) %>%
    dplyr::group_by(.data$screening) %>%
    dplyr::summarise(
      mean_prob = mean(.data$prob * 100),
      lb_prob = stats::quantile(probs = 0.025, x = .data$prob * 100),
      ub_prob = stats::quantile(probs = 0.975, x = .data$prob * 100)
    ) %>%
    tidyr::pivot_longer(cols = c(
      .data$mean_prob,
      .data$lb_prob, .data$ub_prob
    )) %>%
    tidyr::pivot_wider(names_from = .data$screening, values_from = .data$value)
}

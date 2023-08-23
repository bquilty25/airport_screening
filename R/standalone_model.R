pacman::p_load(purrr,furrr,emojifont,gridExtra,knitr,kableExtra,tidyverse,dtplyr,tidyfast,data.table)

source("R/utils.R")
source("data-raw/pathogen_parameters.R")

detect_fun <- function(df){
  #browser()
  travellers <- generate_travellers(df, i = rep(100, df$n_rep))
  probs <- generate_probabilities(travellers)
  
  est_df <- data.frame(CI = apply(X = probs[, -1], 
                                  MARGIN = 2, 
                                  FUN = make_ci_label)) %>%
    rownames_to_column(var = "name") %>%
    mutate(name = factor(name,
                         levels = c("prop_symp_at_exit",
                                    "prop_sev_at_entry",
                                    "prop_symp_at_entry",
                                    "prop_undetected"),
                         labels = c("Detected at exit",
                                    "Severe on flight",
                                    "Detected on entry",
                                    "Not detected"),
                         ordered = TRUE)) %>%
    arrange(name) %>%
    rename(`Detection outcome` = name,
           `Estimate (95% CI)`  = CI)
  
  est_df
  
  waffle_labels <- data.frame(
    desc = factor(c(rep("detected at exit screening",
                        round(probs$prop_symp_at_exit)[1]),
                    rep("detected as severe on flight",
                        round(probs$prop_sev_at_entry)[1]),
                    rep("detected at entry screening",
                        round(probs$prop_symp_at_entry)[1]),
                    rep("not detected",
                        # this is here to make sure we always sum to 100
                        #round(100*(1 - (probs$prob_det_exit + probs$prob_det_entry)),2))),
                        100 - (round(probs$prop_symp_at_exit[1]) + 
                                 round(probs$prop_sev_at_entry)[1] +
                                 round(probs$prop_symp_at_entry[1])))),
                  ordered = T)) %>%
    mutate(
      desc = factor(desc,
                    levels = c("detected at exit screening",
                               "detected as severe on flight",
                               "detected at entry screening",
                               "not detected")))
  
  
  waffle_counts <- count(waffle_labels, desc, .drop = FALSE) %>%
    mutate(desc_comb = paste(n, desc))
  
  waffle_df <- expand.grid(y = -c(1:5), x = 1:20) %>% 
    as_tibble() %>% 
    bind_cols(., waffle_labels) %>%
    mutate(desc_comb = factor(desc, 
                              levels = waffle_counts$desc,
                              labels = waffle_counts$desc_comb))
  
  
  waffle_data <- waffle_df %>% 
    #mutate(group=as.factor(group),
    mutate(label=fontawesome(case_when(
      desc == "detected at exit screening"   ~ "fa-user",
      desc == "detected at entry screening"  ~ "fa-user-circle",
      desc == "detected as severe on flight" ~ "fa-user-circle",
      desc == "not detected"                 ~ "fa-user-circle")))
  
  
  waffle_colors <- RColorBrewer::brewer.pal(4, name = "Set2")[c(1,4,3,2)]
  
  waffle_counts_df <- count(waffle_data, desc, name = "n") 
  waffle_not_detected_on_exit <- waffle_counts_df %>%
    filter(desc != "detected at exit screening") %>%
    summarise(N = sum(n)) %>% pull(N)
  
  waffle_detected_on_entry_or_flight <-
    waffle_counts_df %>%
    filter(desc %in% c("detected as severe on flight",
                       "detected at entry screening")) %>%
    summarise(N = sum(n)) %>% pull(N)
  
  waffle_counts_vec        <- waffle_counts_df$n
  names(waffle_counts_vec) <- waffle_counts_df$desc
  
  # waffle_counts$n <- sum(waffle_counts) - 
  #                      waffle_counts["detected at exit screening"]
  
  waffle_subtitle <- sprintf(
    "%i cases not detected at exit screening.\n%0.2f%% of these %i cases were then detected either during flight or at entry screening",
    waffle_not_detected_on_exit,
    100*waffle_detected_on_entry_or_flight/waffle_not_detected_on_exit,
    waffle_not_detected_on_exit)
  
  names(waffle_colors) <- levels(waffle_data$desc_comb)
  
  
  waffle_plot <- ggplot(waffle_data,
                        aes(x = x, y = y, colour = desc_comb)) + 
    geom_raster(aes(fill = desc_comb),
                alpha = 0.2) +
    geom_text(aes(label=label),
              family="fontawesome-webfont",
              size=6,
              key_glyph="rect") +
    #scale_x_continuous(expand = c(0, 0)) +
    #scale_y_continuous(expand = c(0, 0), trans = 'reverse') +
    scale_colour_manual(values = waffle_colors, drop = FALSE) +
    scale_fill_manual(values = waffle_colors, drop = FALSE) +
    labs(title = "Out of 100 infected travellers:",
         subtitle = waffle_subtitle)+ 
    coord_equal()+
    theme_void()+
    theme(axis.text       = element_blank(),
          axis.title      = element_blank(),
          axis.ticks      = element_blank(),
          legend.title    = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0, size=15),
          plot.subtitle = element_text(hjust = 0, lineheight = 1.5),
          #legend.key.size = unit(3, "line"),
          legend.text     = element_text(size=12)) +
    guides(color=guide_legend(ncol=1)) 
  
  waffle_plot
  
  return(list(res=est_df,
              plot=waffle_plot,
              undetected=0.01*probs %>%filter(name=="mean_prob") %>%  select(prop_undetected)))
}

####################################################################################################

#Incubation period model 

# Initialize scenarios based on pathogen parameters
scenarios <- pathogen_parameters %>%
  filter(name == "Custom") %>%                    # Select scenarios with name "Custom"
  select(-mu_inc) %>%                             # Remove mu_inc column
  mutate(sens.exit = 86,                          # Add sens.exit column with value 86
         sens.entry = 86) %>%                        # Add prop.asy column
  crossing(mu_inc = 1:21, dur.flight = 1:12) %>%  # Create combinations of mu_inc and dur.flight
  mutate(scenario = row_number(),                 # Add scenario column with row numbers
         n_rep = 1000)                            # Add n_rep column with value 1000


# Start timing the simulation
tictoc::tic() 

# Run simulations for each scenario
results <- scenarios %>% 
  group_by(scenario) %>%                         # Group scenarios by scenario ID
  group_split() %>%                              # Split groups into separate data frames
  purrr::map(~ detect_fun(df = .x))              # Apply detect_fun to each split group

# End timing the simulation and display elapsed time
tictoc::toc()

# Process and combine the results of the first simulation element
processed_results_1 <- map(results, 1) %>% 
  bind_rows(.id = "scenario") %>%                # Combine results and add scenario ID column
  mutate(scenario = as.integer(scenario)) %>%    # Convert scenario ID to integer
  left_join(scenarios)                           # Add original scenario parameters

# Extract and display the second simulation element (no further processing mentioned)
processed_results_2 <- map(results, 2)

# Process and combine the results of the third simulation element, and create a heatmap plot
incubation_fig<- map(results, 3) %>% 
  bind_rows(.id = "scenario") %>%                # Combine results and add scenario ID column
  mutate(scenario = as.integer(scenario)) %>%    # Convert scenario ID to integer
  left_join(scenarios) %>%                       # Add original scenario parameters
  ggplot(aes(x = mu_inc, y = dur.flight, fill = prop_undetected)) + 
  geom_tile()+  # Create a heatmap plot using ggplot2
  labs(y = "Flight duration (Hours)", x = "Incubation Period (Days)") +
  theme_classic() +
  scale_x_continuous(breaks=seq(1,21,by=1))+
  scale_y_continuous(breaks=seq(1,12.5,by=1)) +
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 20))

incubation_fig

####################################################################################################

#Asymptomatic period model 

# Initialize scenarios based on pathogen parameters
scenarios <- pathogen_parameters %>%
  filter(name == "Custom") %>%
  select(-prop.asy) %>%
  mutate(sens.exit = 86,
         sens.entry = 86) %>%
  crossing(prop.asy = seq(from=0.01,to=0.90,by=0.05), dur.flight = 1:12) %>%
  mutate(scenario = row_number(),
         n_rep = 1000)

# Start timing the simulation
tictoc::tic() 

# Run simulations for each scenario
results <- scenarios %>% 
  group_by(scenario) %>%                         # Group scenarios by scenario ID
  group_split() %>%                              # Split groups into separate data frames
  purrr::map(~ detect_fun(df = .x))              # Apply detect_fun to each split group

# End timing the simulation and display elapsed time
tictoc::toc()

# Process and combine the results of the first simulation element
processed_results_1 <- map(results, 1) %>% 
  bind_rows(.id = "scenario") %>%                # Combine results and add scenario ID column
  mutate(scenario = as.integer(scenario)) %>%    # Convert scenario ID to integer
  left_join(scenarios)                           # Add original scenario parameters

# Extract and display the second simulation element (no further processing mentioned)
processed_results_2 <- map(results, 2)

# Process and combine the results of the third simulation element, and create a heatmap plot
asym_fig<- map(results, 3) %>% 
  bind_rows(.id = "scenario") %>%                # Combine results and add scenario ID column
  mutate(scenario = as.integer(scenario)) %>%    # Convert scenario ID to integer
  left_join(scenarios) %>%                       # Add original scenario parameters
  ggplot(aes(x = prop.asy, y = dur.flight, fill = prop_undetected)) + 
  geom_tile()+  # Create a heatmap plot using ggplot2
  labs(y = "Flight duration (Hours)", x = "Proportion asymptomatic") +
  theme_classic() +
  scale_x_continuous(breaks=seq(0.0,0.8,by=0.1))+
  scale_y_continuous(breaks=seq(1,12.5,by=1)) +
  scale_fill_viridis_c()+
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 20))

asym_fig


####################################################################################################

#Time to hospitalization period model 

# Initialize scenarios based on pathogen parameters
scenarios <- pathogen_parameters %>%
  filter(name == "Custom") %>%
  select(-mu_inf) %>%
  mutate(sens.exit = 86,
         sens.entry = 86) %>%
  crossing(mu_inf = 1:14, dur.flight = 1:12) %>%
  mutate(scenario = row_number(),
         n_rep = 1000)


# Start timing the simulation
tictoc::tic() 

# Run simulations for each scenario
results <- scenarios %>% 
  group_by(scenario) %>%                         # Group scenarios by scenario ID
  group_split() %>%                              # Split groups into separate data frames
  purrr::map(~ detect_fun(df = .x))              # Apply detect_fun to each split group

# End timing the simulation and display elapsed time
tictoc::toc()

# Process and combine the results of the first simulation element
processed_results_1 <- map(results, 1) %>% 
  bind_rows(.id = "scenario") %>%                # Combine results and add scenario ID column
  mutate(scenario = as.integer(scenario)) %>%    # Convert scenario ID to integer
  left_join(scenarios)                           # Add original scenario parameters

# Extract and display the second simulation element (no further processing mentioned)
processed_results_2 <- map(results, 2)

# Process and combine the results of the third simulation element, and create a heatmap plot
hosp_fig<- map(results, 3) %>% 
  bind_rows(.id = "scenario") %>%                # Combine results and add scenario ID column
  mutate(scenario = as.integer(scenario)) %>%    # Convert scenario ID to integer
  left_join(scenarios) %>%                       # Add original scenario parameters
  ggplot(aes(x = mu_inf, y = dur.flight, fill = prop_undetected)) + 
  geom_tile()+  # Create a heatmap plot using ggplot2
  labs(y = "Flight duration (Hours)", x = "Symptom onset to hospitalisation (Days)") +
  theme_classic() +
  scale_x_continuous(breaks=seq(1,15,by=1))+
  scale_y_continuous(breaks=seq(1,12.5,by=1)) +
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 20))

hosp_fig

####################################################################################################

#Sensitivty of screening

# Initialize scenarios based on pathogen parameters

scenarios <- pathogen_parameters %>%
  filter(name == "Custom") %>%
  mutate( dur.flight = 6,
  ) %>%
  crossing(sens.exit = 1:100, sens.entry = 1:100) %>%
  mutate(scenario = row_number(),
         n_rep = 1000)


# Start timing the simulation
tictoc::tic() 

# Run simulations for each scenario
results <- scenarios %>% 
  group_by(scenario) %>%                         # Group scenarios by scenario ID
  group_split() %>%                              # Split groups into separate data frames
  purrr::map(~ detect_fun(df = .x))              # Apply detect_fun to each split group

# End timing the simulation and display elapsed time
tictoc::toc()

# Process and combine the results of the first simulation element
processed_results_1 <- map(results, 1) %>% 
  bind_rows(.id = "scenario") %>%                # Combine results and add scenario ID column
  mutate(scenario = as.integer(scenario)) %>%    # Convert scenario ID to integer
  left_join(scenarios)                           # Add original scenario parameters

# Extract and display the second simulation element (no further processing mentioned)
processed_results_2 <- map(results, 2)

# Process and combine the results of the third simulation element, and create a heatmap plot
sens.fig<- map(results, 3) %>% 
  bind_rows(.id = "scenario") %>%                # Combine results and add scenario ID column
  mutate(scenario = as.integer(scenario)) %>%    # Convert scenario ID to integer
  left_join(scenarios) %>%                       # Add original scenario parameters
  ggplot(aes(x = sens.entry, y = sens.exit, fill = prop_undetected)) + 
  geom_tile()+  # Create a heatmap plot using ggplot2
  labs(y = "Sensitivity of exit screening", x = "Sensitivity of entry screening") +
  theme_classic() +
  #scale_x_continuous(breaks=seq(1,15,by=1))+
  #scale_y_continuous(breaks=seq(1,12.5,by=1)) +
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 20))

sens.fig

####################################################################################################

#Sensitivity of entry screening vs. flight duration 

#COVID-19

# Initialize scenarios based on pathogen parameters
scenarios <- pathogen_parameters %>%
  filter(name == "nCoV-2019") %>%
  mutate(sens.exit = 0,
  ) %>%
  crossing(sens.entry = 1:100, dur.flight = 1:12) %>%
  mutate(scenario = row_number(),
         n_rep = 1000)



# Start timing the simulation
tictoc::tic() 

# Run simulations for each scenario
results <- scenarios %>% 
  group_by(scenario) %>%                         # Group scenarios by scenario ID
  group_split() %>%                              # Split groups into separate data frames
  purrr::map(~ detect_fun(df = .x))              # Apply detect_fun to each split group

# End timing the simulation and display elapsed time
tictoc::toc()


# Process and combine the results of the third simulation element, and create a heatmap plot
COV_sens.fig<- map(results, 3) %>% 
  bind_rows(.id = "scenario") %>%                # Combine results and add scenario ID column
  mutate(scenario = as.integer(scenario)) %>%    # Convert scenario ID to integer
  left_join(scenarios) %>%                       # Add original scenario parameters
  ggplot(aes(x = sens.entry, y =dur.flight, fill = prop_undetected)) + 
  geom_tile()+  # Create a heatmap plot using ggplot2
  labs(y = "Flight duration (Hours)", x = "Sensitivity of entry screening") +
  theme_classic() +
  #scale_x_continuous(breaks=seq(1,15,by=1))+
  scale_y_continuous(breaks=seq(1,12.5,by=1)) +
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 20))

COV_sens.fig


#H1N1-Like

# Initialize scenarios based on pathogen parameters
scenarios <- pathogen_parameters %>%
  filter(name == "Flu A/H1N1-like (2009)") %>%
  mutate(sens.exit = 0,
  ) %>%
  crossing(sens.entry = 1:100, dur.flight = 1:12) %>%
  mutate(scenario = row_number(),
         n_rep = 1000)



# Start timing the simulation
tictoc::tic() 

# Run simulations for each scenario
results <- scenarios %>% 
  group_by(scenario) %>%                         # Group scenarios by scenario ID
  group_split() %>%                              # Split groups into separate data frames
  purrr::map(~ detect_fun(df = .x))              # Apply detect_fun to each split group

# End timing the simulation and display elapsed time
tictoc::toc()


# Process and combine the results of the third simulation element, and create a heatmap plot
H1N1_sens.fig<- map(results, 3) %>% 
  bind_rows(.id = "scenario") %>%                # Combine results and add scenario ID column
  mutate(scenario = as.integer(scenario)) %>%    # Convert scenario ID to integer
  left_join(scenarios) %>%                       # Add original scenario parameters
  ggplot(aes(x = sens.entry, y =dur.flight, fill = prop_undetected)) + 
  geom_tile()+  # Create a heatmap plot using ggplot2
  labs(y = "Flight duration (Hours)", x = "Sensitivity of entry screening") +
  theme_classic() +
  #scale_x_continuous(breaks=seq(1,15,by=1))+
  scale_y_continuous(breaks=seq(1,12.5,by=1)) +
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 20))

H1N1_sens.fig
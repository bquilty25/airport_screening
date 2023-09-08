pacman::p_load(purrr,furrr,emojifont,gridExtra,knitr,kableExtra,tidyverse,dtplyr,tidyfast,data.table)

source("R/utils.R")

pathogen_parameters <- do.call(
  rbind,
  list(
    data.frame(
      name = "nCoV-2019",
      # (Li et al. (2020) NEJM)
      mu_inc = 5.2,
      sigma_inc = 4.1,
      mu_inf = 9.1,
      sigma_inf = 14.7,
      prop.asy = 0.17
    ),
    data.frame(
      name = "SARS-like (2002)",
      mu_inc = 6.4,
      sigma_inc = 16.7,
      mu_inf = 3.8,
      sigma_inf = 6.0,
      prop.asy = 0.0
    ),
    data.frame(
      name = "Flu A/H1N1-like (2009)",
      mu_inc = 4.3,
      sigma_inc = 1.05,
      mu_inf = 9.3,
      sigma_inf = 0.7,
      prop.asy = 0.16 # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4586318/
    ),
    data.frame(
      name = "MERS-like (2012)",
      mu_inc = 5.5,
      sigma_inc = 6.25,
      # nolint begin
      # https://www.sciencedirect.com/science/article/
      # pii/S1473309913703049?via%3Dihub#sec1
      # nolint end
      mu_inf = 5.0, # https://www.nejm.org/doi/10.1056/NEJMoa1306742
      sigma_inf = 7.5,
      prop.asy = 0.21
      # nolint begin
      # https://doi.org/10.1016/j.tmaid.2018.12.003
      # citing https://www.who.int/csr/disease/coronavirus_infections/
      # risk-assessment-august-2018.pdf?ua=1
      # nolint end
    ),
  data.frame(
    name = "irrelevant",
    mu_inc = 5.0,
    sigma_inc = 5.0,
    mu_inf = 5.0,
    sigma_inf = 5.0,
    prop.asy = 0.5
  ),
    data.frame(
      name = "Custom",
      mu_inc = 5.0,
      sigma_inc = 5.0,
      mu_inf = 5.0,
      sigma_inf = 5.0,
      prop.asy = 0.5
    )
  )
)

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
              prop_undetected=probs %>% slice(1) %>% select(prop_undetected)))
}
 
colnames(scenarios)

# Create Data
scenarios <- pathogen_parameters %>%
  mutate(
    sens.exit = 86,
    sens.entry = 86,
    prop.asy = 17,
    prop_fever = 0.2,         # Example value for prop_fever
    prop_relevant = 0.3,     # Example value for prop_relevant
    n_travellers = 1000      # Example value for n_travellers
  ) %>%
  crossing(., dur.flight = 1:12) %>%
  mutate(scenario = row_number(), n_rep = 1000)

# Run model
tictoc::tic() 
results <- scenarios %>% 
  group_by(scenario) %>% 
  group_split() %>%
  purrr::map(~detect_fun(df=.x)) 
tictoc::toc()

# View results (table and plot) 
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



map(results,1) %>% 
  bind_rows(.id="scenario") %>% 
  mutate(scenario = as.integer(scenario)) %>% 
  left_join(scenarios)

map(results,2)


#########Asymptomatic############
scenarios <- pathogen_parameters %>%
  filter(name == "Custom") %>%
  select(-prop.asy,-mu_inc) %>%
  mutate(sens.exit = 86,
         sens.entry = 86) %>%
  crossing(prop.asy = seq(from=1,to=90,by=5), dur.flight = 6, mu_inc=1:14) %>%
  mutate(scenario = row_number(),
         n_rep = 1000)

tictoc::tic() 

# Run simulations for each scenario
results <- scenarios %>% 
  group_by(scenario) %>%                         # Group scenarios by scenario ID
  group_split() %>%                              # Split groups into separate data frames
  purrr::map(~ detect_fun(df = .x))              # Apply detect_fun to each split group

# End timing the simulation and display elapsed time
tictoc::toc()

asym_fig<- map(results, 3) %>% 
  bind_rows(.id = "scenario") %>%                # Combine results and add scenario ID column
  mutate(scenario = as.integer(scenario)) %>%    # Convert scenario ID to integer
  left_join(scenarios) %>%                       # Add original scenario parameters
  ggplot(aes(x = prop.asy, y = mu_inc, fill = prop_undetected)) + 
  geom_tile()+  # Create a heatmap plot using ggplot2
  labs(x = "Proportion asymptomatic") +
  theme_classic() +
  scale_x_continuous(breaks=seq(0.0,0.8,by=0.1))+
  scale_y_continuous(breaks=seq(1,12.5,by=1)) +
  scale_fill_viridis_c()+
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 20))

asym_fig
 
#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
  library(tidyverse)
  library(dplyr)
  library(furrr)
  #install.packages("emojifont")
  library(emojifont)
  #library(cowplot)
  library(gridExtra)
  library(knitr)
  library(kableExtra)
  library(dtplyr)
  library(tidyfast)
  
})

source("utils.R")


ui <- list(
  tags$style(type="text/css", "
  body { padding-top: 20px; padding-left: 20px; padding-right: 20px}
  .inline label.control-label, .inline .selectize-control.single { 
    display: table-cell; 
    text-align: left; 
    vertical-align: middle; 
  } 
  .inline .form-group { 
    display: table-row;
  }
  .inline .selectize-control.single div.item {
    padding-right: 15px;
  }
  .seroselect .selectize-control.single { width:10em; }
  .demoselect .selectize-control.single { width:15em; }
  .seroselect.inline label.control-label { width: 70px; }
  .demoselect.inline label.control-label { width: 70px; }
  .costs.inline label.control-label { width: 90px; }
  hr { margin:5px; border-top: 1px solid grey; }
"),
  fluidPage(
    titlePanel("Effectiveness of airport screening at detecting infected travellers"),
    shiny::includeMarkdown("date_stamp.md"),
    sidebarLayout(
      sidebarPanel(
        #tags$label(class="h3",)
        sliderInput('dur.flight', 
                    label = "Travel duration (hours)",
                    min = 1, max = 20, step = 1, value = 12),
        sliderInput(inputId = 'sens.exit',
                    label = 'Sensitivity of exit screening',
                    value = 86, min = 0, max = 100, step = 1, post  = " %"),
        sliderInput(inputId = 'sens.entry',
                    label = 'Sensitivity of entry screening',
                    value = 86, min = 0, max = 100, step = 1, post  = " %"),
        sliderInput(inputId = 'prop.asy',
                    label = 'Proportion of cases that are asymptomatic',
                    value = 17, min = 0, max = 100, step = 1),
        hr(),
        div(class="header",
            selectInput('pathogen',label='Pathogen',
                        choices = unique(pathogen$name),
                        selected = pathogen$name[1]),
            numericInput("mu_inc",
                         'Days from infection to symptom onset (mean)',
                         value = 5.2, min = 0.1, max = 20, step = 0.1),
            numericInput("sigma_inc",
                         "Days from infection to symptom onset (variance)",
                         value = 4.1, min = 0.1, max = 20, step = 0.1),
            numericInput("mu_inf",
                         'Days from symptom onset to severe symptoms e.g hospitalisation (mean)',
                         value = 9.1, min = 0.1, max = 20, step = 0.1),
            numericInput("sigma_inf",
                         "Days from symptom onset to severe symptoms e.g hospitalisation (variance)",
                         value = 14.7, min = 0.1, max = 20, step = 0.1)),
        checkboxInput("uncert",
                      label="Show uncertainty (takes longer)",
                      value = FALSE)
        
      ),
      mainPanel(
        tabsetPanel(
          type = "tabs",
          tabPanel(
            title = "Plot",
            (shiny::includeMarkdown("waffle_description.md")),
            fluidRow(
              
              uiOutput("waffle_plot"),
              tableOutput("detailed_estimates"), align = "center"),
            shiny::includeMarkdown("density_description.md"),
            fluidRow(uiOutput("density_plot"))),
          tabPanel(title = "Model",
                   fluidRow(shiny::includeMarkdown("assumptions.md"))),
          tabPanel(title = "References",
                   shiny::includeMarkdown("references.md"))
        )
      )
    )
  )
)


server <- function(input, output, session){
  
  observe({
    default=input$pathogen
    
    updateNumericInput(session,"prop.asy",
                       value = pathogen %>% 
                         filter(name==default) %>% pull(prop.asy))
    updateNumericInput(session,"mu_inc",
                       value = pathogen %>% 
                         filter(name==default) %>% pull(mu_inc))
    updateNumericInput(session,"sigma_inc",
                       value = pathogen %>% 
                         filter(name==default) %>% pull(sigma_inc))
    updateNumericInput(session,"mu_inf",
                       value = pathogen %>% 
                         filter(name==default) %>% pull(mu_inf))
    updateNumericInput(session,"sigma_inf",
                       value = pathogen %>% 
                         filter(name==default) %>% pull(sigma_inf))
  }
  
  )
  
  
  waffle_df <- reactive({
    
    travellers <- generate_travellers(input, i = rep(10000, 1))
    
    
    probs      <- generate_probabilities(travellers)
    
    
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
    return(waffle_df)
  })
  
  nat_hist_periods <- reactive({
    #convert to gamma parameters
    periods <- data.frame(
      inc_period = time_to_event(1e4, input$mu_inc, input$sigma_inc),
      inf_period = time_to_event(1e4, input$mu_inf, input$sigma_inf))
    return(periods)
  })
  
  output$waffleplot <- renderPlot(expr = {
    waffle <- waffle_df()
    waffle_data <- waffle %>% 
      #mutate(group=as.factor(group),
      mutate(label=fontawesome(case_when(
        desc == "detected at exit screening"   ~ "fa-user",
        desc == "detected at entry screening"  ~ "fa-user-circle",
        desc == "detected as severe on flight" ~ "fa-user-circle",
        desc == "not detected"                 ~ "fa-user-circle")))
    
    
    waffle_colors <- RColorBrewer::brewer.pal(4, name = "Set2")[c(1,4,3,2)]
    
    waffle_counts_df <- count(waffle_data, desc, name = "n") 
    waffle_not_detected_on_exit <- waffle_counts_df %>% 
      lazy_dt() %>% 
      filter(desc != "detected at exit screening") %>%
      summarise(N = sum(n)) %>% 
      as.data.frame() %>%  
      pull(N) 
    
    waffle_detected_on_entry_or_flight <-
      waffle_counts_df %>%
      lazy_dt() %>% 
      filter(desc %in% c("detected as severe on flight",
                         "detected at entry screening")) %>%
      summarise(N = sum(n)) %>% 
      as.data.frame() %>% 
      pull(N)
    
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
    
  })
  
  output$densityplot <- renderPlot(expr = {
    
    period_plot_data <- nat_hist_periods()
    
    period_plot_data <- mutate(period_plot_data,
                               severe_period = 
                                 inf_period + inc_period) %>% 
      #lazy_dt() %>% 
      dt_pivot_longer(names_to = "Period",
                          #cols = everything(),
                          values_to = "value") %>%
      mutate(Period = factor(Period,
                             levels = c("inc_period",
                                        "inf_period",
                                        "severe_period"),
                             labels = c("Infection to onset",
                                        "Onset to severe",
                                        "Infection to severe"))) %>% 
      as.data.frame()
    
    
    period_plot <- ggplot(data = period_plot_data,
                          aes(x = value)) +
      geom_density(fill  = "lightskyblue",
                   color = "lightskyblue",
                   alpha = 0.5) +
      labs(x = "Time (days)") +
      theme_minimal() +
      ylab("Density")  +
      facet_grid(. ~ Period) +
      geom_vline(
        data = period_plot_data %>%
          group_by(Period) %>%
          summarise(mean = mean(value)),
        aes(xintercept = mean),
        lty = 2
      ) +
      labs(title = "Vertical lines represent mean time to event") +
      theme(plot.title = element_text(hjust = 0.5))
    
    period_plot
    
  }, execOnResize = FALSE)
  
  
  height = function() {
    session$clientData$output_waffleplot_width
  }
  
  
  output$waffle_plot <- renderUI({
    plotOutput("waffleplot")
  })
  
  output$density_plot <- renderUI({
    plotOutput("densityplot", width = "100%", height = "2in")
  })
  
  output$detailed_estimates <- renderTable(
    
    if (input$uncert==TRUE){
      travellers <- generate_travellers(input, i = rep(100, 200))
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
    }
    
    else{
      NULL
    },
    align = "lr"
  )
  
}

shinyApp(ui=ui,server=server)

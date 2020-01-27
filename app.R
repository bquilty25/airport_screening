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
  library(cowplot)

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
    titlePanel(""),
    sidebarLayout(
      sidebarPanel(
        #tags$label(class="h3",)
        sliderInput('dur.flight', 
                    label = "Travel duration (hours)",
                    min = 0.1, max = 20, step = 1, value = 12),
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
                         value = 1, min = 0.1, max = 20, step = 0.1),
            numericInput("sigma_inc",
                         "Days from infection to symptom onset (variance)",
                         value = 1, min = 0.1, max = 20, step = 0.1),
            numericInput("mu_inf",
                         'Days from symptom onset to severe symptoms e.g hospitalisation (mean)',
                         value = 1, min = 0.1, max = 20, step = 0.1),
            numericInput("sigma_inf",
                         "Days from symptom onset to severe symptoms e.g hospitalisation (variance)",
                         value = 1, min = 0.1, max = 20, step = 0.1)
        )
        
      ),
      mainPanel(
        h1("Effectiveness of airport screening at detecting infected travellers"),
        shiny::includeMarkdown("date_stamp.md"),
        tabsetPanel(type = "tabs",
                    tabPanel(title = "Plot",
                             fluidRow(shiny::includeMarkdown("waffle_description.md")),
                             fluidRow(uiOutput("waffle_plot"))),
                             
                    tabPanel(title = "Model",
                             fluidRow(shiny::includeMarkdown("assumptions.md"))),
                    tabPanel(title = "References",
                             shiny::includeMarkdown("references.md"))
          )
        
      )
    )
  )
  #fluidPage(title = "Key assumptions",
  #          fluidRow(shiny::includeMarkdown("assumptions.md"))
  #),
  #fluidPage(title = "References",
  #          fluidRow(shiny::includeMarkdown("references.md"))
  )


server <- function(input, output, session){
  
  observe({
    default=input$pathogen
    
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
    
    # round(probs$prob_det_exit,2)
    # round(probs$prob_det_entry,2)
    # 100 - (round(probs$prob_det_exit,2) + round(probs$prob_det_entry_only,2))
    # 
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
      mutate(label=fontawesome("fa-user"))
    
    waffle_colors <- RColorBrewer::brewer.pal(4, name = "Set2")[c(1,4,3,2)]
    
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
      labs(title = "Out of 100 infected travellers:")+ 
      coord_equal()+
      theme_void()+
      theme(axis.text       = element_blank(),
            axis.title      = element_blank(),
            axis.ticks      = element_blank(),
            legend.title    = element_blank(),
            legend.position = "bottom",
            plot.title = element_text(hjust = 0.5, size=15),
            #legend.key.size = unit(3, "line"),
            legend.text     = element_text(size=12)) +
      guides(color=guide_legend(ncol=1)) 
    
    period_plot_data <- nat_hist_periods()
    
    inc_period_plot  <- ggplot(period_plot_data,aes(x=inc_period)) + 
      geom_density(fill  = "lightskyblue",
                   color = "lightskyblue",
                   alpha = 0.5) +
      labs(x="Infection to onset (days)") + theme_minimal() +
      ylab("Density") 
    inf_period_plot  <- ggplot(period_plot_data, aes(x=inf_period)) + 
      geom_density(fill  = "lightskyblue",
                   color = "lightskyblue",
                   alpha = 0.5) +
      labs(x="Onset to severe (days)") + theme_minimal() +
      ylab("Density")
    
    
    plots <- align_plots(waffle_plot, inc_period_plot,
                         align = 'v', axis = 'l')
    # then build the bottom row
    bottom_row <- plot_grid(plots[[2]], inf_period_plot)
    
    # then combine with the top row for final plot
    plot_grid(plots[[1]], bottom_row, ncol = 1,
              rel_heights = c(3,2))
    
  },
  

  height = function() {
    session$clientData$output_waffleplot_width
  }
  )
  
  output$waffle_plot <- renderUI({
    plotOutput("waffleplot",height="auto")
  })
}

shinyApp(ui=ui,server=server)

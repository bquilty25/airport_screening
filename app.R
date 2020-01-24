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
  library(dplyr)
  #devtools::install_github("liamgilbey/ggwaffle")
  library(ggwaffle)
  #install.packages("emojifont")
  library(emojifont)
  #install.packages("patchwork")
  library(cowplot)
  
})

# fileInput("file1", "Choose CSV File", multiple = FALSE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
# tags$hr(),
# checkboxInput("header", "Header", TRUE),
# radioButtons("sep", "Separator", choices = c(Comma = ",", Semicolon = ";", Tab = "\t"), selected = ",")

calc_probs <- function(dur.flight,mu_inc,sigma_inc,mu_inf,sigma_inf,sens.exit,sens.entry,prop.asy,sims){
  
  #convert to gamma parameters
  alpha_inc=(mu_inc^2)/sigma_inc
  beta_inc=sigma_inc/mu_inc
  
  #convert to gamma parameters
  alpha_inf=(mu_inf^2)/sigma_inf
  beta_inf=sigma_inf/mu_inf
  
  #convert flight time to days
  dur_flight=dur.flight/24
  
  df = data.frame(incu = rgamma(sims,shape=alpha_inc,scale=beta_inc),
                  inf = rgamma(sims,shape=alpha_inf,scale=beta_inf),
                  flight.at = runif(sims,0,(mu_inc+mu_inf)*2)) %>% 
    mutate(arrive.inf = flight.at < (mu_inc+mu_inf-dur_flight)) %>%
    filter(arrive.inf == T) %>%
    mutate(symp.at.exit = incu<flight.at) %>%
    mutate(symp.at.entry = incu<(flight.at+dur_flight)) %>%
    mutate(found.at.exit = symp.at.exit * (runif(n(),0,1)<sens.exit)) %>%
    mutate(found.at.entry = symp.at.entry * (runif(n(),0,1)<sens.entry)) %>%
    mutate(found.only.at.entry = found.at.entry-found.at.exit)
  
  #results
  prob_det_exit <- ((sum(df$found.at.exit) / dim(df)[1])*(1-prop.asy)) %>% round(2) #prop found at exit
  
  if(sens.entry>=sens.exit){
    prob_det_entry <- ((sum(df$found.only.at.entry) / dim(df)[1])*(1-prop.asy)) %>% round(2)  #additional prop found atentry
    prob_undet <- 1-(prob_det_exit+prob_det_entry)
    return(list("prob_det_exit"=prob_det_exit,"prob_det_entry"=prob_det_entry,"prob_undet"=prob_undet))
  }
  else{
    prob_det_entry <-0
    prob_undet <- 1-(prob_det_exit+prob_det_entry)
    return(list("prob_det_exit"=prob_det_exit,"prob_det_entry"=prob_det_entry,"prob_undet"=prob_undet))
  }
}

ui <- list(
  tags$style(type="text/css", "
  body { padding-top: 70px; }
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
    tags$label(class="h3","Effectiveness of airport screening at detecting infected travellers"),
    sidebarLayout(
      sidebarPanel(
        #tags$label(class="h3",)
        sliderInput('dur.flight', 
                    label = "Flight duration (hours)",
                    min = 0.1, max = 20, step = 1, value = 10),
        sliderInput(inputId = 'sens.exit',
                    label = 'Sensitivity of exit screening',
                    value = 0.9, min = 0, max = 1, step = 0.1),
        sliderInput(inputId = 'sens.entry',
                    label = 'Sensitivity of entry screening',
                    value = 0.9, min = 0, max = 1, step = 0.1),
        hr(),
        div(class="header",
            numericInput("mu_inc",
                         'Days between infection and symptom onset (mean)',
                         value = 6.7, min = 0.1, max = 20, step = 0.1),
            numericInput("sigma_inc",
                         "Days between infection and symptom onset (variance)",
                         value = 16.7, min = 0.1, max = 20, step = 0.1),
            numericInput("mu_inf",
                         'Days between symptom onset and recovery (mean)',
                         value = 9.3, min = 0.1, max = 20, step = 0.1),
            numericInput("sigma_inf",
                         "Days between symptom onset and recovery (variance)",
                         value = 5.2, min = 0.1, max = 20, step = 0.1)
        )
        
      ),
      mainPanel(
        
        fluidRow(shiny::includeMarkdown("waffle_description.md")),
        uiOutput("waffle_plot")
        #tableOutput("testROI")
      )
    )
  ),
  fluidPage(title = "References",
            fluidRow(shiny::includeMarkdown("references.md"))
  ))


server <- function(input, output, session){
  
  waffle_df <- reactive({
    probs <- calc_probs(dur.flight = input$dur.flight,
                        mu_inc=input$mu_inc,
                        sigma_inc=input$sigma_inc,
                        mu_inf=input$mu_inf,
                        sigma_inf=input$sigma_inf,
                        sens.exit=input$sens.exit,
                        sens.entry=input$sens.entry,
                        prop.asy=0,
                        sims=10000)
    
    waffle_df <- expand.grid(y = 1:10, x = 1:10) %>% 
      as_tibble() %>% 
      mutate(desc=forcats::as_factor(c(rep("Case found at exit screening",
                                           round(probs$prob_det_exit*100,2)),
                                       rep("Additional case found at entry screening",
                                           round(probs$prob_det_entry*100,2)),
                                       rep("Case not detected by screening",
                                           round(probs$prob_undet*100,2))))) %>%
      mutate(desc = factor(desc, levels = c("Case found at exit screening",
                                            "Additional case found at entry screening",
                                            "Case not detected by screening"),
                           ordered = TRUE))
    return(waffle_df)
  })
  
  nat_hist_periods <- reactive({
    #convert to gamma parameters
    alpha_inc=(input$mu_inc^2)/input$sigma_inc
    beta_inc=input$sigma_inc/input$mu_inc
    
    #convert to gamma parameters
    alpha_inf=(input$mu_inf^2)/input$sigma_inf
    beta_inf=input$sigma_inf/input$mu_inf
    
    periods <- data.frame(inc_period=rgamma(10000,shape=alpha_inc,scale=beta_inc),
                          inf_period=rgamma(10000,shape=alpha_inf,scale=beta_inf))
    return(periods)
  })
  
  output$waffleplot <- renderPlot(expr = {
    waffle <- waffle_df()
    waffle_data <- waffle %>% 
      #mutate(group=as.factor(group),
      mutate(label=fontawesome("fa-user"))
    
    waffle_colors <- RColorBrewer::brewer.pal(3, name = "Set2")[c(1,3,2)]
    names(waffle_colors) <- c("Case found at exit screening",
                              "Additional case found at entry screening",
                              "Case not detected by screening")
    
    waffle_plot <- ggplot(waffle_data,aes(x = x, y = y, colour = desc)) + 
      geom_raster(aes(fill = desc),
                alpha = 0.2) +
      geom_text(aes(label=label),
                family="fontawesome-webfont",
                size=6, key_glyph="rect") +
      #scale_x_continuous(expand = c(0, 0)) +
      #scale_y_continuous(expand = c(0, 0), trans = 'reverse') +
      scale_colour_manual(values = waffle_colors, drop = FALSE) +
      scale_fill_manual(values = waffle_colors, drop = FALSE) +
      labs(title = paste0("Each person represents 1 out of 100 infected travellers \nwho could infect others at the destination"))+ 
      coord_equal()+
      theme_void()+
      theme(axis.text       = element_blank(),
            axis.title      = element_blank(),
            axis.ticks      = element_blank(),
            legend.title    = element_blank(),
            legend.position = "bottom",
            legend.text     = element_text(size=8)) +
      guides(color=guide_legend(ncol=1)) 
    
    period_plot_data <- nat_hist_periods()
    
    inc_period_plot  <- ggplot(period_plot_data,aes(x=inc_period)) + 
      geom_density(fill = "lightskyblue",
                   color = "lightskyblue",
                   alpha = 0.5) +
      labs(x="Infection to onset (days)") + theme_minimal() +
      ylab("Probability density")
    inf_period_plot  <- ggplot(period_plot_data, aes(x=inf_period)) + 
      geom_density(fill = "lightskyblue",
                   color = "lightskyblue",
                   alpha = 0.5) +
      labs(x="Onset to recovery (days)") + theme_minimal() +
      ylab("Probability density")
    
    
    plots <- align_plots(waffle_plot, inc_period_plot,
                         align = 'v', axis = 'l')
    # then build the bottom row
    bottom_row <- plot_grid(plots[[2]], inf_period_plot)
    
    # then combine with the top row for final plot
    plot_grid(plots[[1]], NULL, bottom_row, ncol = 1,
              rel_heights = c(3,0.25,1))
    
    

  },
  
  
  height = function() {
    session$clientData$output_waffleplot_width
  })
  
  output$waffle_plot <- renderUI({
    plotOutput("waffleplot", height = input$plotSize)
  })
}

shinyApp(ui=ui,server=server)


#' CMMID Shiny for the efficacy of airport screening for infectious diseases
#'
#' @description Runs a local instance of a Shiny app developed by the Centre for
#' the Mathematical Modelling of Infectious Diseases at the London School of
#' Hygiene and tropical Medicine, to model the efficacy of implementing airport
#' screening for infectious diseases, as a part of the CMMID response to the
#' Covid-19 pandemic, in early 2020.
#' @return Runs the `{shiny}` app locally.
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' # Choose options from the app GUI
#' \dontrun{
#' run_app()
#' }
run_app <- function() {
  waffle_description <- system.file(
    "info", "waffle_description.md",
    package = "airportscreening"
  )
  density_description <- system.file(
    "info", "density_description.md",
    package = "airportscreening"
  )
  assumptions <- system.file(
    "info", "assumptions.md",
    package = "airportscreening"
  )
  references <- system.file(
    "info", "references.md",
    package = "airportscreening"
  )

  pathogen_parameters <- airportscreening::pathogen_parameters

  ui <- list(
    shiny::tags$style(type = "text/css", "
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
    shiny::fluidPage(
      shiny::titlePanel(
        paste0(
          "Effectiveness of airport screening at detecting",
          "infected travellers"
        )
      ),
      shiny::markdown(get_date_stamp()),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::sliderInput("dur.flight",
            label = "Travel duration (hours)",
            min = 1, max = 20, step = 1, value = 12
          ),
          shiny::sliderInput(
            inputId = "sens.exit",
            label = "Sensitivity of exit screening",
            value = 86, min = 0, max = 100, step = 1, post = " %"
          ),
          shiny::sliderInput(
            inputId = "sens.entry",
            label = "Sensitivity of entry screening",
            value = 86, min = 0, max = 100, step = 1, post = " %"
          ),
          shiny::sliderInput(
            inputId = "prop.asy",
            label = "Proportion of cases that are asymptomatic",
            value = 17, min = 0, max = 100, step = 1
          ),
          shiny::hr(),
          shiny::div(
            class = "header",
            shiny::selectInput("pathogen",
              label = "Pathogen",
              choices = unique(pathogen_parameters$name),
              selected = pathogen_parameters$name[1]
            ),
            shiny::numericInput("mu_inc",
              "Days from infection to symptom onset (mean)",
              value = 5.2, min = 0.1, max = 20, step = 0.1
            ),
            shiny::numericInput("sigma_inc",
              "Days from infection to symptom onset (variance)",
              value = 4.1, min = 0.1, max = 20, step = 0.1
            ),
            shiny::numericInput("mu_inf",
              paste0(
                "Days from symptom onset to severe symptoms e.g",
                " hospitalisation (mean)"
              ),
              value = 9.1, min = 0.1, max = 20, step = 0.1
            ),
            shiny::numericInput("sigma_inf",
              paste0(
                "Days from symptom onset to severe symptoms e.g",
                " hospitalisation (variance)"
              ),
              value = 14.7, min = 0.1, max = 20, step = 0.1
            )
          ),
          shiny::checkboxInput("uncert",
            label = "Show uncertainty (takes longer)",
            value = FALSE
          )
        ),
        shiny::mainPanel(
          shiny::tabsetPanel(
            type = "tabs",
            shiny::tabPanel(
              title = "Plot",
              (shiny::includeMarkdown(waffle_description)),
              shiny::fluidRow(
                shiny::uiOutput("waffle_plot"),
                shiny::tableOutput("detailed_estimates"),
                align = "center"
              ),
              shiny::includeMarkdown(density_description),
              shiny::fluidRow(shiny::uiOutput("density_plot"))
            ),
            shiny::tabPanel(
              title = "Model",
              shiny::fluidRow(
                shiny::includeMarkdown(assumptions)
              )
            ),
            shiny::tabPanel(
              title = "References",
              shiny::includeMarkdown(references)
            )
          )
        )
      )
    )
  )
  server <- function(input, output, session) {
    shiny::observe({
      pathogen_input <- input$pathogen

      shiny::updateNumericInput(session, "prop.asy",
        value = pathogen_parameters[pathogen_parameters$name ==
          pathogen_input, ]$prop.asy
      )
      shiny::updateNumericInput(session, "mu_inc",
        value = pathogen_parameters[pathogen_parameters$name ==
          pathogen_input, ]$mu_inc
      )
      shiny::updateNumericInput(session, "sigma_inc",
        value = pathogen_parameters[pathogen_parameters$name ==
          pathogen_input, ]$sigma_inc
      )
      shiny::updateNumericInput(session, "mu_inf",
        value = pathogen_parameters[pathogen_parameters$name ==
          pathogen_input, ]$mu_inf
      )
      shiny::updateNumericInput(session, "sigma_inf",
        value = pathogen_parameters[pathogen_parameters$name ==
          pathogen_input, ]$sigma_inf
      )
    })

    waffle_df <- shiny::reactive({
      travellers <- generate_travellers(input, i = rep(10000, 1))

      probs <- generate_probabilities(travellers)

      waffle_labels <- data.frame(
        desc = factor(
          c(
            rep(
              "detected at exit screening",
              round(probs$prop_symp_at_exit)[1]
            ),
            rep(
              "detected as severe on flight",
              round(probs$prop_sev_at_entry)[1]
            ),
            rep(
              "detected at entry screening",
              round(probs$prop_symp_at_entry)[1]
            ),
            rep(
              "not detected",
              100 - (round(probs$prop_symp_at_exit[1]) +
                round(probs$prop_sev_at_entry)[1] +
                round(probs$prop_symp_at_entry[1]))
            )
          ),
          ordered = TRUE
        )
      ) %>%
        dplyr::mutate(
          desc = factor(
            .data$desc,
            levels = c(
              "detected at exit screening",
              "detected as severe on flight",
              "detected at entry screening",
              "not detected"
            )
          )
        )


      waffle_counts <- dplyr::count(waffle_labels, .data$desc) %>%
        dplyr::mutate(desc_comb = paste(.data$n, .data$desc))

      waffle_df <- expand.grid(y = -seq_len(5), x = seq_len(20)) %>%
        tibble::as_tibble() %>%
        dplyr::bind_cols(waffle_labels) %>%
        dplyr::mutate(desc_comb = factor(
          .data$desc,
          levels = waffle_counts$desc,
          labels = waffle_counts$desc_comb
        ))
      return(waffle_df)
    })

    nat_hist_periods <- shiny::reactive({
      # convert to gamma parameters
      periods <- data.frame(
        inc_period = time_to_event(1e4, input$mu_inc, input$sigma_inc),
        inf_period = time_to_event(1e4, input$mu_inf, input$sigma_inf)
      )
      return(periods)
    })

    output$waffleplot <- shiny::renderPlot(expr = {
      waffle <- waffle_df()
      waffle_data <- waffle %>%
        # dplyr::mutate(group=as.factor(group),
        dplyr::mutate(label = emojifont::fontawesome(
          dplyr::case_when(
            .data$desc == "detected at exit screening" ~ "fa-user",
            .data$desc == "detected at entry screening" ~ "fa-user-circle",
            .data$desc == "detected as severe on flight" ~ "fa-user-circle",
            .data$desc == "not detected" ~ "fa-user-circle"
          )
        ))


      waffle_colors <- RColorBrewer::brewer.pal(4, name = "Set2")[c(1, 4, 3, 2)]

      waffle_counts_df <- dplyr::count(waffle_data, .data$desc, name = "n")
      waffle_not_detected_on_exit <- waffle_counts_df %>%
        dplyr::filter(.data$desc != "detected at exit screening") %>%
        dplyr::summarise(N = sum(.data$n)) %>%
        dplyr::pull(.data$N)

      waffle_detected_on_entry_or_flight <-
        waffle_counts_df %>%
        dplyr::filter(
          .data$desc %in% c(
            "detected as severe on flight",
            "detected at entry screening"
          )
        ) %>%
        dplyr::summarise(N = sum(.data$n)) %>%
        as.data.frame() %>%
        dplyr::pull(.data$N)

      waffle_counts_vec <- waffle_counts_df$n
      names(waffle_counts_vec) <- waffle_counts_df$desc

      waffle_subtitle <- sprintf(
        paste0(
          "%i cases not detected at exit screening.\n%0.2f%% of these %i cases",
          " were then detected either during flight or at entry screening"
        ),
        waffle_not_detected_on_exit,
        100 * waffle_detected_on_entry_or_flight / waffle_not_detected_on_exit,
        waffle_not_detected_on_exit
      )

      names(waffle_colors) <- levels(waffle_data$desc_comb)

      waffle_plot <- ggplot2::ggplot(
        data = waffle_data,
        ggplot2::aes(x = .data$x, y = .data$y, colour = .data$desc_comb)
      ) +
        ggplot2::geom_raster(
          ggplot2::aes(fill = .data$desc_comb),
          alpha = 0.2
        ) +
        ggplot2::geom_text(
          ggplot2::aes(label = .data$label),
          family = "fontawesome-webfont",
          size = 6,
          key_glyph = "rect"
        ) +
        ggplot2::scale_colour_manual(values = waffle_colors, drop = FALSE) +
        ggplot2::scale_fill_manual(values = waffle_colors, drop = FALSE) +
        ggplot2::labs(
          title = "Out of 100 infected travellers:",
          subtitle = waffle_subtitle
        ) +
        ggplot2::coord_equal() +
        ggplot2::theme_void() +
        ggplot2::theme(
          axis.text = ggplot2::element_blank(),
          axis.title = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          legend.title = ggplot2::element_blank(),
          legend.position = "bottom",
          plot.title = ggplot2::element_text(hjust = 0, size = 15),
          plot.subtitle = ggplot2::element_text(hjust = 0, lineheight = 1.5),
          # legend.key.size = unit(3, "line"),
          legend.text = ggplot2::element_text(size = 12)
        ) +
        ggplot2::guides(color = ggplot2::guide_legend(ncol = 1))

      waffle_plot
    })

    output$densityplot <- shiny::renderPlot(expr = {
      period_plot_data <- nat_hist_periods()

      period_plot_data <- dplyr::mutate(
        period_plot_data,
        severe_period =
          .data$inf_period + .data$inc_period
      ) %>%
        tidyr::pivot_longer(
          names_to = "Period",
          cols = dplyr::everything(),
          values_to = "value"
        ) %>%
        dplyr::mutate(
          Period = factor(
            .data$Period,
            levels = c(
              "inc_period",
              "inf_period",
              "severe_period"
            ),
            labels = c(
              "Infection to onset",
              "Onset to severe",
              "Infection to severe"
            )
          )
        )

      period_plot <- ggplot2::ggplot(
        data = period_plot_data,
        ggplot2::aes(x = .data$value)
      ) +
        ggplot2::geom_density(
          fill = "lightskyblue",
          color = "lightskyblue",
          alpha = 0.5
        ) +
        ggplot2::labs(x = "Time (days)") +
        ggplot2::theme_minimal() +
        ggplot2::ylab("Density") +
        ggplot2::facet_grid(cols = ggplot2::vars(.data$Period)) +
        ggplot2::geom_vline(
          data = period_plot_data %>%
            dplyr::group_by(.data$Period) %>%
            dplyr::summarise(mean = mean(.data$value)),
          ggplot2::aes(xintercept = .data$mean),
          lty = 2
        ) +
        ggplot2::labs(title = "Vertical lines represent mean time to event") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

      period_plot
    }, execOnResize = FALSE)

    output$waffle_plot <- shiny::renderUI({
      shiny::plotOutput("waffleplot")
    })

    output$density_plot <- shiny::renderUI({
      shiny::plotOutput("densityplot", width = "100%", height = "2in")
    })

    output$detailed_estimates <- shiny::renderTable(

      if (input$uncert == TRUE) {
        travellers <- generate_travellers(input, i = rep(100, 200))
        probs <- generate_probabilities(travellers)

        est_df <- data.frame(
          CI = apply(
            X = probs[, -1],
            MARGIN = 2,
            FUN = make_ci_label
          )
        ) %>%
          tibble::rownames_to_column(var = "name") %>%
          dplyr::mutate(
            name = factor(
              .data$name,
              levels = c(
                "prop_symp_at_exit",
                "prop_sev_at_entry",
                "prop_symp_at_entry",
                "prop_undetected"
              ),
              labels = c(
                "Detected at exit",
                "Severe on flight",
                "Detected on entry",
                "Not detected"
              ),
              ordered = TRUE
            )
          ) %>%
          dplyr::arrange(.data$name) %>%
          dplyr::rename(
            "Detection outcome" = .data$name,
            "Estimate (95% CI)" = .data$CI
          )
        est_df
      } else {
        NULL
      },
      align = "lr"
    )
  }
  shiny::shinyApp(ui = ui, server = server)
}

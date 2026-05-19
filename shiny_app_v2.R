# inst/shiny/app.R
# -----------------
# Shiny interactive app for the lcaStats pedigree module.
# Equivalent of the Python Streamlit app, built for R users.
#
# Run from R console:
#   shiny::runApp("inst/shiny/app.R")
#
# Or if the package is installed:
#   lcaStats::run_pedigree_app()

library(shiny)
library(ggplot2)

# Source the core module (works whether run standalone or via package)
source("pedigree_R.R")
source("pedigree_ilcd.R")
# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

ui <- fluidPage(

  titlePanel(
    windowTitle = "Pedigree Uncertainty Explorer",
    title = div(
      h2("🌿 Pedigree Uncertainty Explorer"),
      p(
        style = "color: #666; font-size: 0.9em;",
        "Convert ecoinvent Data Quality Indicator (DQI) scores into ",
        "lognormal uncertainty parameters for Life Cycle Assessment. ",
        "Based on Weidema et al. (2013)."
      )
    )
  ),

  sidebarLayout(

    # --- Sidebar ---
    sidebarPanel(
      width = 3,

      h4("DQI System"),
      radioButtons("system", NULL,
                   choices = c("ecoinvent pedigree" = "ecoinvent",
                               "ILCD / EU EF"       = "ilcd"),
                   selected = "ecoinvent", inline = TRUE),

      hr(),
      h4("Data Quality Indicators"),
      p(style = "font-size:0.85em; color:#555;",
        "1 = best quality · 5 = worst quality"),

      # ecoinvent-only indicators (hidden for ILCD)
      conditionalPanel(
        condition = "input.system == 'ecoinvent'",
        sliderInput("reliability", "Reliability",
                    min = 1, max = 5, value = 2, step = 1, ticks = TRUE),
        sliderInput("completeness", "Completeness",
                    min = 1, max = 5, value = 2, step = 1, ticks = TRUE)
      ),

      # Shared indicators (shown for both)
      sliderInput("temporal", "Temporal representativeness / correlation",
                  min = 1, max = 5, value = 3, step = 1, ticks = TRUE),
      sliderInput("geographical", "Geographical representativeness / correlation",
                  min = 1, max = 5, value = 2, step = 1, ticks = TRUE),
      sliderInput("technology", "Technology representativeness / correlation",
                  min = 1, max = 5, value = 3, step = 1, ticks = TRUE),

      hr(),
      h4("Basic uncertainty"),
      numericInput("basic_var", "Basic variance (σ² in log-space)",
                   value = 0.0006, min = 0, max = 1, step = 0.0001),
      p(style = "font-size:0.8em; color:#777;",
        "Intrinsic exchange uncertainty. Default 0.0006 is typical for ",
        "many ecoinvent exchanges."),

      hr(),
      h4("Exchange value (optional)"),
      numericInput("exchange_mean", "Inventory amount (any unit)",
                   value = 1.0, min = 0.001, step = 0.1),
      numericInput("mc_n", "Monte Carlo sample size",
                   value = 10000, min = 100, max = 100000, step = 1000)
    ),

    # --- Main panel ---
    mainPanel(
      width = 9,

      fluidRow(

        # Key metrics
        column(4,
          h4("Uncertainty parameters"),
          tableOutput("metrics_table"),
          hr(),
          h5("Lognormal parameters"),
          uiOutput("lognormal_params_text"),
          p(style = "font-size:0.8em; color:#777;",
            "Pass (mu_ln, sigma_ln) to rlnorm() for Monte Carlo sampling.")
        ),

        # Variance breakdown
        column(4,
          h4("Variance breakdown"),
          tableOutput("breakdown_table")
        ),

        # Distribution plot
        column(4,
          h4("Sampled distribution"),
          plotOutput("dist_plot", height = "280px"),
          uiOutput("dist_caption")
        )
      ),

      hr(),

      # Pedigree lookup table
      h4("📋 Pedigree lookup table (Weidema et al., 2013)"),
      tableOutput("pedigree_table"),
      p(style = "font-size:0.8em; color:#777;",
        "Values are additional variance terms (σ²_ln) per indicator and score."),

      hr(),

      # Batch upload
      h4("📂 Batch processing — upload a CSV"),
      p(style = "font-size:0.85em;",
        "Upload a CSV with columns: ",
        code("reliability"), ", ", code("completeness"), ", ",
        code("temporal_correlation"), ", ", code("geographical_correlation"),
        ", ", code("technology_correlation"), ", ", code("basic_var"), ". ",
        "Extra columns (e.g. exchange_name, amount) are preserved."
      ),
      fileInput("csv_upload", "Choose CSV file", accept = ".csv"),
      uiOutput("batch_results"),

      hr(),
      p(style = "font-size:0.8em; color:#888;",
        strong("Reference: "),
        "Weidema, B.P., et al. (2013). Overview and methodology: Data quality ",
        "guideline for the ecoinvent database version 3. Ecoinvent Report 1(v3). ",
        "The ecoinvent Centre, St. Gallen. | ",
        strong("Code: "),
        "MIT Licence · Harper Food Innovation: Digital · ",
        a("github.com/food-innovation-HA-team",
          href = "https://github.com/food-innovation-HA-team", target = "_blank")
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

server <- function(input, output, session) {

  # Reactive: compute uncertainty for whichever system is selected
  result <- reactive({
    if (input$system == "ecoinvent") {
      scores_to_gsd2(
        reliability              = input$reliability,
        completeness             = input$completeness,
        temporal_correlation     = input$temporal,
        geographical_correlation = input$geographical,
        technology_correlation   = input$technology,
        basic_var                = input$basic_var
      )
    } else {
      scores_to_gsd2_ilcd(
        technological_representativeness = input$technology,
        geographical_representativeness  = input$geographical,
        temporal_representativeness      = input$temporal,
        basic_var                        = input$basic_var
      )
    }
  })

  # Reactive: lognormal params
  ln_params <- reactive({
    lognormal_params(mean = input$exchange_mean, sigma_ln = result()$sigma_ln)
  })

  # Reactive: MC samples
  samples <- reactive({
    set.seed(42)
    rlnorm(input$mc_n, meanlog = ln_params()["mu_ln"], sdlog = ln_params()["sigma_ln"])
  })

  # --- Key metrics table ---
  output$metrics_table <- renderTable({
    r <- result()
    rows <- data.frame(
      Parameter = c("GSD² (ecoinvent / adapted)", "GSD", "σ_ln", "Combined variance (σ²_ln)"),
      Value     = c(round(r$gsd2, 4), round(r$gsd, 4),
                    round(r$sigma_ln, 6), round(r$combined_var, 6))
    )
    # Add DQR row for ILCD
    if (input$system == "ilcd") {
      rows <- rbind(rows, data.frame(
        Parameter = sprintf("ILCD DQR (%s)", r$dqr_level),
        Value     = round(r$dqr, 2)
      ))
    }
    rows
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  # --- Lognormal params text ---
  output$lognormal_params_text <- renderUI({
    p <- ln_params()
    HTML(sprintf(
      "For an exchange with mean <strong>%.3f</strong>:<br>
       &mu;<sub>ln</sub> = <strong>%.4f</strong><br>
       &sigma;<sub>ln</sub> = <strong>%.6f</strong>",
      input$exchange_mean, p["mu_ln"], p["sigma_ln"]
    ))
  })

  # --- Variance breakdown table ---
  output$breakdown_table <- renderTable({
    r        <- result()
    label_map <- if (input$system == "ecoinvent") INDICATOR_LABELS else ILCD_LABELS
    sources  <- c("Basic uncertainty", as.character(label_map))
    variances <- c(r$basic_var, r$indicator_vars)
    pct      <- if (r$combined_var > 0) variances / r$combined_var * 100 else rep(0, length(variances))
    data.frame(
      Source            = sources,
      `σ² contribution` = round(variances, 6),
      `% of total`      = paste0(round(pct, 1), "%"),
      check.names       = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  # --- Distribution plot ---
  output$dist_plot <- renderPlot({
    s    <- samples()
    p025 <- quantile(s, 0.025)
    p50  <- quantile(s, 0.500)
    p975 <- quantile(s, 0.975)

    ggplot(data.frame(x = s), aes(x = x)) +
      geom_histogram(aes(y = after_stat(density)),
                     bins = 60, fill = "#2e7d52", alpha = 0.75, colour = NA) +
      geom_vline(xintercept = p50,  linetype = "dashed", colour = "#111", linewidth = 0.8) +
      geom_vline(xintercept = p025, linetype = "dotted", colour = "#b55", linewidth = 0.8) +
      geom_vline(xintercept = p975, linetype = "dotted", colour = "#b55", linewidth = 0.8) +
      labs(x = "Sampled value", y = "Density") +
      theme_minimal(base_size = 11) +
      theme(panel.grid.minor = element_blank())
  })

  output$dist_caption <- renderUI({
    s    <- samples()
    p025 <- round(quantile(s, 0.025), 3)
    p975 <- round(quantile(s, 0.975), 3)
    p(style = "font-size:0.8em; color:#666;",
      sprintf("n = %s draws. 95%% CI: [%s, %s]",
              format(input$mc_n, big.mark = ","), p025, p975))
  })

  # --- Pedigree lookup table ---
  output$pedigree_table <- renderTable({
    if (input$system == "ecoinvent") {
      rows <- lapply(1:5, function(score) {
        row <- list(Score = score)
        for (ind in names(PEDIGREE_TABLE)) {
          row[[INDICATOR_LABELS[ind]]] <- PEDIGREE_TABLE[[ind]][as.character(score)]
        }
        as.data.frame(row)
      })
    } else {
      rows <- lapply(1:5, function(score) {
        row <- list(Score = score)
        for (ind in names(ILCD_TABLE)) {
          row[[ILCD_LABELS[ind]]] <- ILCD_TABLE[[ind]][as.character(score)]
        }
        as.data.frame(row)
      })
    }
    do.call(rbind, rows)
  }, striped = TRUE, hover = TRUE, bordered = TRUE, digits = 4)

  # --- Batch processing ---
  output$batch_results <- renderUI({
    req(input$csv_upload)
    tryCatch({
      df  <- read.csv(input$csv_upload$datapath)
      out <- batch_from_dataframe(df)

      output$batch_table <- renderTable(out, striped = TRUE, hover = TRUE,
                                        bordered = TRUE, digits = 6)
      output$download_batch <- downloadHandler(
        filename = "pedigree_uncertainty_results.csv",
        content  = function(file) write.csv(out, file, row.names = FALSE)
      )

      tagList(
        p(style = "color: #2e7d52; font-weight: bold;",
          sprintf("✓ Processed %d exchanges.", nrow(out))),
        tableOutput("batch_table"),
        downloadButton("download_batch", "⬇️ Download results as CSV")
      )
    }, error = function(e) {
      p(style = "color: red;", sprintf("Could not process file: %s", e$message))
    })
  })
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

shinyApp(ui = ui, server = server)

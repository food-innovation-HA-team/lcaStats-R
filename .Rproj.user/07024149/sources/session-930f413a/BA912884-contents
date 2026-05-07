# monte_carlo.R
# -------------
# Monte Carlo uncertainty propagation for food system LCA.
#
# This module extends the pedigree uncertainty quantification in pedigree.R
# and pedigree_ilcd.R to full Monte Carlo simulation across food items and
# supply chain stages.
#
# Ported and substantially extended from a Python notebook developed within
# the Harper Food Innovation: Digital group.
#
# DESIGN
# ------
# The module is designed to connect directly to batch_from_dataframe() output
# from pedigree.R, so no separate variance file is needed. The full pipeline
# from DQI scores to simulation results runs entirely in R:
#
#   inventory |>
#     batch_from_dataframe() |>      # pedigree.R — generates sigma_ln
#     run_monte_carlo(n = 10000)     # this module
#
# Alternatively, if GSD² values are already available (e.g. from ecoinvent
# or SimaPro exports), they can be supplied directly via run_monte_carlo_gsd2().
#
# REFERENCES
# ----------
# Weidema, B.P., et al. (2013). Overview and methodology: Data quality
# guideline for the ecoinvent database version 3. Ecoinvent Report 1(v3).
#
# Ciroth, A., et al. (2016). Empirically based uncertainty factors for the
# pedigree matrix in ecoinvent. International Journal of Life Cycle
# Assessment, 21(9), 1338-1348. https://doi.org/10.1007/s11367-013-0670-5
#
# Igos, E., et al. (2019). How to treat uncertainties in life cycle
# assessment studies? International Journal of Life Cycle Assessment,
# 24(7), 1261-1277. https://doi.org/10.1007/s11367-018-1477-1


EPSILON <- 1e-10  # Small value to replace zeros or negative means


# ---------------------------------------------------------------------------
# Core simulation function — connects to batch_from_dataframe() output
# ---------------------------------------------------------------------------

#' Run Monte Carlo simulation from pedigree batch output
#'
#' Takes the output of \code{\link{batch_from_dataframe}} (or
#' \code{\link{batch_from_dataframe_ilcd}}) and runs a lognormal Monte Carlo
#' simulation across food items and supply chain stages.
#'
#' Each exchange is treated as an independent lognormal random variable.
#' Stage impacts are summed across exchanges within each stage. Total impacts
#' are summed across all stages.
#'
#' @param df A data frame — output of \code{batch_from_dataframe()}, with
#'   columns: \code{food_item}, \code{stage}, \code{amount}, \code{sigma_ln}.
#'   Any additional columns are ignored.
#' @param n Integer. Number of Monte Carlo iterations. Default 10,000.
#' @param seed Integer or NULL. Random seed for reproducibility. Default 42.
#' @param ci Numeric. Confidence interval width (0-1). Default 0.95
#'   gives a 95% CI (2.5th and 97.5th percentiles).
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{total}{Data frame of summary statistics for total impacts
#'       per food item.}
#'     \item{stages}{Named list of data frames — one per stage — with
#'       summary statistics per food item.}
#'     \item{samples}{Named list of raw sample matrices (food items x n)
#'       for total and each stage. Use for custom analysis or plotting.}
#'   }
#'
#' @examples
#' \dontrun{
#' source("pedigree.R")
#'
#' inventory <- data.frame(
#'   food_item = c("Wheat bread", "Wheat bread", "Beef burger", "Beef burger"),
#'   stage     = c("Agriculture", "Processing", "Agriculture", "Processing"),
#'   amount    = c(0.5, 0.2, 15.0, 2.0),
#'   reliability = c(2L, 2L, 3L, 2L),
#'   completeness = c(2L, 2L, 2L, 2L),
#'   temporal_correlation = c(3L, 2L, 3L, 2L),
#'   geographical_correlation = c(1L, 1L, 2L, 1L),
#'   technology_correlation = c(2L, 2L, 3L, 2L),
#'   basic_var = rep(0.0006, 4)
#' )
#'
#' df_with_uncertainty <- batch_from_dataframe(inventory)
#' results <- run_monte_carlo(df_with_uncertainty, n = 10000)
#' print(results$total)
#' }
#'
#' @export
run_monte_carlo <- function(df, n = 10000, seed = 42, ci = 0.95) {

  # --- Input validation ---
  required <- c("food_item", "stage", "amount", "sigma_ln")
  missing  <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(sprintf(
      "Missing required columns: %s. Did you run batch_from_dataframe() first?",
      paste(missing, collapse = ", ")
    
    ))
  }
  if (n < 100) stop("n must be at least 100 for meaningful results.")
  if (ci <= 0 | ci >= 1) stop("ci must be between 0 and 1 (e.g. 0.95).")

  if (!is.null(seed)) set.seed(seed)

  lo <- (1 - ci) / 2
  hi <- 1 - lo

  food_items <- unique(df$food_item)
  stages     <- unique(df$stage)

  # Storage: raw samples
  total_samples <- matrix(0, nrow = n, ncol = length(food_items),
                          dimnames = list(NULL, food_items))

  stage_samples <- lapply(stages, function(s) {
    matrix(0, nrow = n, ncol = length(food_items),
           dimnames = list(NULL, food_items))
  })
  names(stage_samples) <- stages

  # --- Simulation loop ---
  for (food in food_items) {
    food_rows <- df[df$food_item == food, ]

    for (stage in stages) {
      stage_rows <- food_rows[food_rows$stage == stage, ]
      if (nrow(stage_rows) == 0) next

      # Sum contributions from all exchanges within this food x stage
      stage_draw <- rep(0, n)

      for (i in seq_len(nrow(stage_rows))) {
        mean_val  <- max(stage_rows$amount[i], EPSILON)
        sigma_ln  <- stage_rows$sigma_ln[i]
        mu_ln     <- log(mean_val) - 0.5 * sigma_ln^2

        stage_draw <- stage_draw +
          rlnorm(n, meanlog = mu_ln, sdlog = sigma_ln)
      }

      stage_samples[[stage]][, food] <- stage_draw
      total_samples[, food]          <- total_samples[, food] + stage_draw
    }
  }

  # --- Summary statistics ---
  summarise_samples <- function(mat) {
    do.call(rbind, lapply(colnames(mat), function(food) {
      x <- mat[, food]
      data.frame(
        food_item  = food,
        mean       = mean(x),
        median     = median(x),
        sd         = sd(x),
        ci_lower   = quantile(x, lo),
        ci_upper   = quantile(x, hi),
        cv_pct     = sd(x) / mean(x) * 100,
        row.names  = NULL
      )
    }))
  }

  total_summary <- summarise_samples(total_samples)
  names(total_summary)[5:6] <- c(
    sprintf("ci_%s_lower", ci * 100),
    sprintf("ci_%s_upper", ci * 100)
  )

  stage_summary <- lapply(stage_samples, summarise_samples)
  for (s in names(stage_summary)) {
    names(stage_summary[[s]])[5:6] <- c(
      sprintf("ci_%s_lower", ci * 100),
      sprintf("ci_%s_upper", ci * 100)
    )
  }

  list(
    total   = total_summary,
    stages  = stage_summary,
    samples = c(list(total = total_samples), stage_samples),
    meta    = list(n = n, seed = seed, ci = ci,
                   food_items = food_items, stages = stages)
  )
}


# ---------------------------------------------------------------------------
# Alternative entry point — GSD² supplied directly
# ---------------------------------------------------------------------------

#' Run Monte Carlo from pre-existing GSD² values
#'
#' For cases where GSD² values are already available — for example from an
#' ecoinvent or SimaPro export — rather than derived from DQI scores.
#' Mirrors the structure of the original Python notebook.
#'
#' @param mean_df Data frame of mean impact values. First column must be
#'   \code{food_item}. Remaining columns are supply chain stages.
#' @param gsd2_df Data frame of GSD² values. Must have identical structure
#'   to \code{mean_df}.
#' @param n Integer. Number of Monte Carlo iterations. Default 10,000.
#' @param seed Integer or NULL. Random seed. Default 42.
#' @param ci Numeric. Confidence interval width. Default 0.95.
#'
#' @return Same structure as \code{\link{run_monte_carlo}}.
#'
#' @examples
#' mean_df <- data.frame(
#'   food_item  = c("Wheat bread", "Beef burger"),
#'   Agriculture = c(0.5, 15.0),
#'   Processing  = c(0.2, 2.0)
#' )
#' gsd2_df <- data.frame(
#'   food_item  = c("Wheat bread", "Beef burger"),
#'   Agriculture = c(1.08, 1.22),
#'   Processing  = c(1.05, 1.10)
#' )
#' results <- run_monte_carlo_gsd2(mean_df, gsd2_df, n = 10000)
#' results$total
#'
#' @export
run_monte_carlo_gsd2 <- function(mean_df, gsd2_df,
                                  n = 10000, seed = 42, ci = 0.95) {

  # --- Consistency checks (mirrors Python notebook) ---
  if (!identical(dim(mean_df), dim(gsd2_df))) {
    stop("mean_df and gsd2_df must have the same dimensions.")
  }
  if (!identical(names(mean_df), names(gsd2_df))) {
    stop("Column names of mean_df and gsd2_df must match.")
  }
  if (!identical(mean_df[[1]], gsd2_df[[1]])) {
    stop("Food items in mean_df and gsd2_df must match.")
  }

  # Convert wide format to long, compute sigma_ln from GSD²
  stages     <- names(mean_df)[-1]
  food_items <- mean_df[[1]]

  rows <- do.call(rbind, lapply(seq_along(food_items), function(i) {
    do.call(rbind, lapply(stages, function(stage) {
      mean_val <- max(mean_df[i, stage], EPSILON)
      gsd2_val <- max(gsd2_df[i, stage], 1 + EPSILON)  # GSD² must be > 1
      sigma_ln <- sqrt(log(gsd2_val))  # GSD² = exp(sigma²) → sigma = sqrt(log(GSD²))
      data.frame(
        food_item = food_items[i],
        stage     = stage,
        amount    = mean_val,
        sigma_ln  = sigma_ln,
        stringsAsFactors = FALSE
      )
    }))
  }))

  run_monte_carlo(rows, n = n, seed = seed, ci = ci)
}


# ---------------------------------------------------------------------------
# Export helpers
# ---------------------------------------------------------------------------

#' Export Monte Carlo results to Excel
#'
#' Writes simulation results to a multi-sheet Excel file, matching the
#' output structure of the original Python notebook.
#'
#' Requires the \code{writexl} package.
#'
#' @param results List. Output from \code{\link{run_monte_carlo}} or
#'   \code{\link{run_monte_carlo_gsd2}}.
#' @param file Character. Output file path. Default
#'   \code{"mc_results.xlsx"}.
#' @param include_samples Logical. Whether to include raw sample sheets.
#'   Can produce large files for high n. Default \code{FALSE}.
#'
#' @return Invisibly returns the file path.
#'
#' @examples
#' \dontrun{
#' export_mc_results(results, file = "wheat_beef_mc.xlsx")
#' }
#'
#' @export
export_mc_results <- function(results,
                               file = "mc_results.xlsx",
                               include_samples = FALSE) {

  if (!requireNamespace("writexl", quietly = TRUE)) {
    stop("Package 'writexl' is required. Install with: install.packages('writexl')")
  }

  sheets <- list("Total Impact Stats" = results$total)

  for (stage in names(results$stages)) {
    sheets[[paste(stage, "Stats")]] <- results$stages[[stage]]
  }

  if (include_samples) {
    for (nm in names(results$samples)) {
      sheets[[paste(nm, "Samples")]] <- as.data.frame(results$samples[[nm]])
    }
  }

  writexl::write_xlsx(sheets, path = file)
  message(sprintf("Results saved to '%s'", file))
  invisible(file)
}


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

#' Plot total impact uncertainty across food items
#'
#' Produces a point-range plot showing median and confidence interval
#' for total environmental impact per food item.
#'
#' Requires \code{ggplot2}.
#'
#' @param results List. Output from \code{\link{run_monte_carlo}}.
#' @param title Character. Plot title.
#' @param x_label Character. X axis label (impact category + unit).
#' @param colour Character. Bar colour. Default Harper green.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' \dontrun{
#' plot_mc_total(results,
#'   title = "Monte Carlo results — wheat and beef",
#'   x_label = "GWP100 (kg CO2-eq per kg food)")
#' }
#'
#' @export
plot_mc_total <- function(results,
                           title    = "Monte Carlo uncertainty — total impact",
                           x_label  = "Environmental impact",
                           colour   = "#2e7d52") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  df  <- results$total
  ci_cols <- grep("^ci_", names(df), value = TRUE)
  lo_col  <- ci_cols[grep("lower", ci_cols)]
  hi_col  <- ci_cols[grep("upper", ci_cols)]
  ci_pct  <- results$meta$ci * 100

  ggplot2::ggplot(df,
    ggplot2::aes(x = median,
                 y = reorder(food_item, median))) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data[[lo_col]], xmax = .data[[hi_col]]),
      width = 0.3, colour = colour, linewidth = 0.8,
      orientation = "y"
    ) +
    ggplot2::geom_point(size = 3, colour = colour) +
    ggplot2::labs(
      title    = title,
      subtitle = sprintf("%s%% CI from %s Monte Carlo iterations",
                         ci_pct, format(results$meta$n, big.mark = ",")),
      x = x_label,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor  = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank())
}


#' Plot stage-level uncertainty for a single food item
#'
#' Stacked bar chart showing median stage contributions with uncertainty
#' ranges for a single food item.
#'
#' @param results List. Output from \code{\link{run_monte_carlo}}.
#' @param food Character. Name of the food item to plot.
#' @param title Character. Plot title.
#' @param y_label Character. Y axis label.
#'
#' @return A \code{ggplot2} object.
#'
#' @export
plot_mc_stages <- function(results,
                            food,
                            title   = NULL,
                            y_label = "Environmental impact") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  if (!food %in% results$meta$food_items) {
    stop(sprintf("'%s' not found. Available: %s",
                 food, paste(results$meta$food_items, collapse = ", ")))
  }

  stage_data <- do.call(rbind, lapply(names(results$stages), function(stage) {
    row <- results$stages[[stage]][results$stages[[stage]]$food_item == food, ]
    if (nrow(row) == 0) return(NULL)
    ci_cols <- grep("^ci_", names(row), value = TRUE)
    lo_col  <- ci_cols[grep("lower", ci_cols)]
    hi_col  <- ci_cols[grep("upper", ci_cols)]
    data.frame(
      stage    = stage,
      median   = row$median,
      ci_lower = row[[lo_col]],
      ci_upper = row[[hi_col]],
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(title)) title <- sprintf("Stage contributions — %s", food)
  ci_pct <- results$meta$ci * 100

  ggplot2::ggplot(stage_data,
    ggplot2::aes(x = median, y = reorder(stage, median))) +
    ggplot2::geom_col(fill = "#1565c0", alpha = 0.75, width = 0.6) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
      width = 0.25, colour = "#333", linewidth = 0.7,
      orientation = "y"
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = sprintf("%s%% CI · n = %s",
                         ci_pct, format(results$meta$n, big.mark = ",")),
      x = y_label, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor   = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank())
}


# Two food items × two stages
mean_df <- data.frame(
  food_item   = c("Wheat bread", "Beef burger"),
  Agriculture = c(0.8, 22.0),
  Processing  = c(0.3,  1.5)
)

gsd2_df <- data.frame(
  food_item   = c("Wheat bread", "Beef burger"),
  Agriculture = c(1.05,  1.22),
  Processing  = c(1.08,  1.10)
)

results <- run_monte_carlo_gsd2(mean_df, gsd2_df, n = 10000)

# Inspect outputs
results$total
results$stages$Agriculture
results$stages$Processing

library(ggplot2)

plot_mc_total(results,
              title   = "Monte Carlo — wheat bread vs beef burger",
              x_label = "GWP100 (kg CO2-eq per kg food)")

plot_mc_stages(results, food = "Beef burger",
               y_label = "GWP100 (kg CO2-eq per kg food)")

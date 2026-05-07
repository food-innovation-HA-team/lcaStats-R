# pedigree_R.R
# ----------
# Convert ecoinvent Data Quality Indicator (DQI) scores into lognormal
# uncertainty parameters for Life Cycle Assessment.
#
# Based on the pedigree matrix approach described in:
#
#   Weidema, B.P., et al. (2013). Overview and methodology: Data quality
#   guideline for the ecoinvent database version 3. Ecoinvent Report 1(v3).
#   St. Gallen: The ecoinvent Centre.
#
#   Frischknecht, R., et al. (2005). The ecoinvent database: Overview and
#   methodological framework. International Journal of Life Cycle Assessment,
#   10(1), 3-9. https://doi.org/10.1065/lca2004.10.181.1
#
# The five Data Quality Indicators (DQI):
#   1 = best quality / most representative
#   5 = worst quality / least representative
#
# NOTE ON VALUES
# The sigma2_ln values in PEDIGREE_TABLE are taken from Weidema et al. (2013).
# Users should verify these against Table 5 of the primary source before use
# in peer-reviewed work.


# ---------------------------------------------------------------------------
# Pedigree lookup table
# Values represent additional variance (sigma^2_ln) per indicator and score.
# Source: Weidema et al. (2013), ecoinvent methodology v3, Table 5.
# ---------------------------------------------------------------------------

#' @keywords internal
PEDIGREE_TABLE <- list(
  reliability = c(
    "1" = 0.0000,
    "2" = 0.0001,
    "3" = 0.0006,
    "4" = 0.0006,
    "5" = 0.0208
  ),
  completeness = c(
    "1" = 0.0000,
    "2" = 0.0001,
    "3" = 0.0002,
    "4" = 0.0006,
    "5" = 0.0130
  ),
  temporal_correlation = c(
    "1" = 0.0000,
    "2" = 0.0001,
    "3" = 0.0002,
    "4" = 0.0008,
    "5" = 0.0080
  ),
  geographical_correlation = c(
    "1" = 0.0000,
    "2" = 0.0001,
    "3" = 0.0002,
    "4" = 0.0002,
    "5" = 0.0040
  ),
  technology_correlation = c(
    "1" = 0.0000,
    "2" = 0.0001,
    "3" = 0.0008,
    "4" = 0.0080,
    "5" = 0.0120
  )
)

#' @keywords internal
INDICATOR_LABELS <- c(
  reliability              = "Reliability",
  completeness             = "Completeness",
  temporal_correlation     = "Temporal correlation",
  geographical_correlation = "Geographical correlation",
  technology_correlation   = "Further technology correlation"
)


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

#' Look up additional variance for a single DQI score
#'
#' Returns the additional log-space variance (sigma^2_ln) associated with
#' a given Data Quality Indicator and score, as defined in the ecoinvent
#' pedigree matrix (Weidema et al., 2013).
#'
#' @param indicator Character. One of: \code{"reliability"},
#'   \code{"completeness"}, \code{"temporal_correlation"},
#'   \code{"geographical_correlation"}, \code{"technology_correlation"}.
#' @param score Integer. Data quality score (1 = best, 5 = worst).
#'
#' @return Numeric. Additional variance in log-space (sigma^2_ln).
#'
#' @references
#' Weidema, B.P., et al. (2013). \emph{Overview and methodology: Data quality
#' guideline for the ecoinvent database version 3.} Ecoinvent Report 1(v3).
#' The ecoinvent Centre, St. Gallen.
#'
#' @examples
#' dqi_to_variance("reliability", 3)   # returns 0.0006
#' dqi_to_variance("completeness", 1)  # returns 0.0000
#'
#' @export
dqi_to_variance <- function(indicator, score) {
  if (!indicator %in% names(PEDIGREE_TABLE)) {
    stop(
      sprintf(
        "Unknown indicator '%s'. Valid options: %s",
        indicator,
        paste(names(PEDIGREE_TABLE), collapse = ", ")
      )
    )
  }
  if (!score %in% 1:5) {
    stop(sprintf("Score must be an integer between 1 and 5, got %s.", score))
  }
  PEDIGREE_TABLE[[indicator]][as.character(score)]
}


#' Combine basic uncertainty with all five pedigree indicator variances
#'
#' Computes the total variance in log-space by summing the basic uncertainty
#' variance and the additional variance contributed by each of the five
#' pedigree indicators. The combination assumes independence between
#' uncertainty sources, consistent with the ecoinvent methodology.
#'
#' @param reliability Integer. DQI score for reliability (1-5).
#' @param completeness Integer. DQI score for completeness (1-5).
#' @param temporal_correlation Integer. DQI score for temporal correlation (1-5).
#' @param geographical_correlation Integer. DQI score for geographical
#'   correlation (1-5).
#' @param technology_correlation Integer. DQI score for further technology
#'   correlation (1-5).
#' @param basic_var Numeric. Basic uncertainty expressed as variance in
#'   log-space (sigma^2_ln). A typical default for ecoinvent exchanges is
#'   0.0006. Check the exchange-specific value in your database.
#'
#' @return Numeric. Combined variance in log-space (sigma^2_ln_total).
#'
#' @examples
#' scores_to_combined_variance(2, 2, 3, 1, 2, basic_var = 0.0006)
#'
#' @export
scores_to_combined_variance <- function(reliability,
                                        completeness,
                                        temporal_correlation,
                                        geographical_correlation,
                                        technology_correlation,
                                        basic_var) {
  if (basic_var < 0) {
    stop(sprintf("basic_var must be >= 0, got %s.", basic_var))
  }

  pedigree_var <- sum(
    dqi_to_variance("reliability",              reliability),
    dqi_to_variance("completeness",             completeness),
    dqi_to_variance("temporal_correlation",     temporal_correlation),
    dqi_to_variance("geographical_correlation", geographical_correlation),
    dqi_to_variance("technology_correlation",   technology_correlation)
  )

  basic_var + pedigree_var
}


#' Full pipeline: DQI scores to lognormal uncertainty parameters
#'
#' Converts five DQI scores and a basic uncertainty variance into a complete
#' uncertainty summary including GSD^2 (the ecoinvent uncertainty metric)
#' and the lognormal parameters needed for Monte Carlo simulation.
#'
#' @inheritParams scores_to_combined_variance
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{combined_var}{Total sigma^2_ln (basic + all pedigree contributions).}
#'     \item{gsd2}{GSD^2 = exp(sigma^2_ln), the ecoinvent uncertainty metric.}
#'     \item{gsd}{Geometric standard deviation = sqrt(GSD^2).}
#'     \item{sigma_ln}{Std dev of ln(X); pass directly to \code{rlnorm()}.}
#'     \item{basic_var}{The basic variance as supplied.}
#'     \item{indicator_vars}{Named numeric vector of per-indicator variances.}
#'   }
#'
#' @references
#' Weidema, B.P., et al. (2013). \emph{Overview and methodology: Data quality
#' guideline for the ecoinvent database version 3.} Ecoinvent Report 1(v3).
#' The ecoinvent Centre, St. Gallen.
#'
#' @examples
#' result <- scores_to_gsd2(
#'   reliability = 2, completeness = 2, temporal_correlation = 3,
#'   geographical_correlation = 1, technology_correlation = 2,
#'   basic_var = 0.0006
#' )
#' result$gsd2
#' result$sigma_ln
#'
#' @export
scores_to_gsd2 <- function(reliability,
                            completeness,
                            temporal_correlation,
                            geographical_correlation,
                            technology_correlation,
                            basic_var) {

  scores <- list(
    reliability              = reliability,
    completeness             = completeness,
    temporal_correlation     = temporal_correlation,
    geographical_correlation = geographical_correlation,
    technology_correlation   = technology_correlation
  )

  indicator_vars <- sapply(
    names(scores),
    function(ind) dqi_to_variance(ind, scores[[ind]])
  )

  combined_var <- basic_var + sum(indicator_vars)
  gsd2         <- exp(combined_var)

  list(
    combined_var   = combined_var,
    gsd2           = gsd2,
    gsd            = sqrt(gsd2),
    sigma_ln       = sqrt(combined_var),
    basic_var      = basic_var,
    indicator_vars = indicator_vars
  )
}


#' Compute lognormal distribution parameters for Monte Carlo sampling
#'
#' Returns the mu_ln and sigma_ln parameters of a lognormal distribution
#' given the desired arithmetic mean and log-space standard deviation.
#' These can be passed directly to \code{\link[stats]{rlnorm}}.
#'
#' @param mean Numeric. Arithmetic mean of the exchange (e.g. the inventory
#'   value). Must be > 0.
#' @param sigma_ln Numeric. Standard deviation of ln(X), from
#'   \code{scores_to_gsd2()$sigma_ln}.
#'
#' @return Named numeric vector with elements \code{mu_ln} and
#'   \code{sigma_ln}.
#'
#' @details
#' For a lognormal distribution:
#'   \deqn{\mu_{ln} = \ln(\text{mean}) - 0.5 \cdot \sigma_{ln}^2}
#'
#' @examples
#' result  <- scores_to_gsd2(2, 2, 3, 1, 2, basic_var = 0.0006)
#' params  <- lognormal_params(mean = 1.5, sigma_ln = result$sigma_ln)
#' samples <- rlnorm(10000, meanlog = params["mu_ln"], sdlog = params["sigma_ln"])
#'
#' @export
lognormal_params <- function(mean, sigma_ln) {
  if (mean <= 0) {
    stop(sprintf(
      "mean must be > 0 for a lognormal distribution, got %s.", mean
    ))
  }
  mu_ln <- log(mean) - 0.5 * sigma_ln^2
  c(mu_ln = mu_ln, sigma_ln = sigma_ln)
}


# ---------------------------------------------------------------------------
# Batch processing
# ---------------------------------------------------------------------------

#' Process a data frame of exchanges with DQI scores
#'
#' Applies \code{\link{scores_to_gsd2}} row-wise to a data frame containing
#' DQI scores for multiple exchanges. Any additional columns in the input
#' (e.g. \code{exchange_name}, \code{amount}, \code{unit}) are preserved.
#'
#' @param df A data frame with columns: \code{reliability},
#'   \code{completeness}, \code{temporal_correlation},
#'   \code{geographical_correlation}, \code{technology_correlation}, and a
#'   basic uncertainty column (see \code{basic_var_col}).
#' @param basic_var_col Character. Name of the column containing the basic
#'   uncertainty variance. Default is \code{"basic_var"}.
#'
#' @return The original data frame with four additional columns:
#'   \code{combined_var}, \code{gsd2}, \code{gsd}, \code{sigma_ln}.
#'
#' @examples
#' df <- data.frame(
#'   exchange              = c("N2O field emissions", "Diesel combustion"),
#'   reliability           = c(3L, 2L),
#'   completeness          = c(2L, 2L),
#'   temporal_correlation  = c(3L, 2L),
#'   geographical_correlation = c(2L, 1L),
#'   technology_correlation   = c(3L, 2L),
#'   basic_var             = c(0.0006, 0.0006)
#' )
#' batch_from_dataframe(df)
#'



# ---------------------------------------------------------------------------
# Display helper
# ---------------------------------------------------------------------------

#' Print a readable summary of scores_to_gsd2() output
#'
#' @param result A list returned by \code{\link{scores_to_gsd2}}.
#' @param exchange_name Character. Label shown in the header. Default
#'   \code{"Exchange"}.
#'
#' @return Invisibly returns \code{result}. Called for its side effect of
#'   printing to the console.
#'
#' @examples
#' result <- scores_to_gsd2(2, 2, 3, 1, 2, basic_var = 0.0006)
#' pedigree_summary(result, exchange_name = "N2O field emission")
#'
#' @export
pedigree_summary <- function(result, exchange_name = "Exchange") {
  sep <- strrep("\u2500", 52)
  cat("\n", sep, "\n", sep = "")
  cat("  Pedigree uncertainty summary:", exchange_name, "\n")
  cat(sep, "\n")
  cat(sprintf("  %-36s %.6f\n", "Basic variance", result$basic_var))
  cat("\n")
  for (ind in names(result$indicator_vars)) {
    label <- INDICATOR_LABELS[ind]
    cat(sprintf("  %-36s %.6f\n", label, result$indicator_vars[ind]))
  }
  cat(sep, "\n")
  cat(sprintf("  %-36s %.6f\n", "Combined variance (sigma^2_ln)", result$combined_var))
  cat(sprintf("  %-36s %.4f\n",  "GSD^2",                          result$gsd2))
  cat(sprintf("  %-36s %.4f\n",  "GSD",                            result$gsd))
  cat(sprintf("  %-36s %.6f\n", "sigma_ln (for MC sampling)",      result$sigma_ln))
  cat(sep, "\n")
  invisible(result)
}

#' Apply pedigree uncertainty to every row of an inventory dataframe
#'
#' Takes a long-format inventory with one exchange per row and DQI score
#' columns, calls \code{scores_to_gsd2()} on each row, and returns the
#' original dataframe with uncertainty columns appended.
#'
#' @param df Data frame with columns: \code{food_item}, \code{stage},
#'   \code{amount}, \code{reliability}, \code{completeness},
#'   \code{temporal_correlation}, \code{geographical_correlation},
#'   \code{technology_correlation}, \code{basic_var}.
#'
#' @return The input dataframe with four columns appended:
#'   \code{sigma_ln}, \code{gsd2}, \code{gsd}, \code{combined_var}.
#'
#' @references
#' Weidema, B.P., et al. (2013). \emph{Overview and methodology: Data quality
#' guideline for the ecoinvent database version 3.} Ecoinvent Report 1(v3).
#'
#' @examples
#' inventory <- data.frame(
#'   food_item                = c("Wheat bread", "Beef burger"),
#'   stage                    = c("Agriculture", "Agriculture"),
#'   amount                   = c(0.8, 22.0),
#'   reliability              = c(2L, 3L),
#'   completeness             = c(2L, 2L),
#'   temporal_correlation     = c(3L, 3L),
#'   geographical_correlation = c(1L, 2L),
#'   technology_correlation   = c(2L, 3L),
#'   basic_var                = c(0.0006, 0.0006)
#' )
#' batch_from_dataframe(inventory)
#'
#' @export
batch_from_dataframe <- function(df) {
  
  required <- c("food_item", "stage", "amount",
                "reliability", "completeness", "temporal_correlation",
                "geographical_correlation", "technology_correlation",
                "basic_var")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(sprintf("Missing required columns: %s",
                 paste(missing, collapse = ", ")))
  }
  
  results <- lapply(seq_len(nrow(df)), function(i) {
    scores_to_gsd2(
      reliability              = df$reliability[i],
      completeness             = df$completeness[i],
      temporal_correlation     = df$temporal_correlation[i],
      geographical_correlation = df$geographical_correlation[i],
      technology_correlation   = df$technology_correlation[i],
      basic_var                = df$basic_var[i]
    )
  })
  
  df$sigma_ln     <- sapply(results, `[[`, "sigma_ln")
  df$gsd2         <- sapply(results, `[[`, "gsd2")
  df$gsd          <- sapply(results, `[[`, "gsd")
  df$combined_var <- sapply(results, `[[`, "combined_var")
  
  df
}


# ---------------------------------------------------------------------------
# Shiny / flexible batch processing — no food_item/stage/amount required
# ---------------------------------------------------------------------------

#' @export
batch_from_csv <- function(df, basic_var_col = "basic_var") {
  
  required_cols <- c(names(PEDIGREE_TABLE), basic_var_col)
  missing_cols  <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s",
                 paste(missing_cols, collapse = ", ")))
  }
  
  results <- lapply(seq_len(nrow(df)), function(i) {
    row <- df[i, ]
    out <- scores_to_gsd2(
      reliability              = as.integer(row$reliability),
      completeness             = as.integer(row$completeness),
      temporal_correlation     = as.integer(row$temporal_correlation),
      geographical_correlation = as.integer(row$geographical_correlation),
      technology_correlation   = as.integer(row$technology_correlation),
      basic_var                = as.numeric(row[[basic_var_col]])
    )
    data.frame(combined_var=out$combined_var, gsd2=out$gsd2,
               gsd=out$gsd, sigma_ln=out$sigma_ln)
  })
  
  cbind(df, do.call(rbind, results))
}
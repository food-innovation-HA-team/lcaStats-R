# pedigree_ilcd.R
# ----------------
# ILCD data quality indicator extension for the lcaStats pedigree module.
#
# The International Reference Life Cycle Data System (ILCD) Handbook
# (European Commission / JRC, 2010) defines a Data Quality Rating (DQR)
# system based on THREE representativeness indicators:
#
#   1. Technological representativeness (TeR)
#   2. Geographical representativeness  (GR)
#   3. Temporal representativeness      (TiR)
#
# Each indicator is scored 1 (best) to 5 (worst). The composite DQR is:
#
#   DQR = (TeR + GR + TiR) / 3
#
# DQR is used in the EU Environmental Footprint (EF/PEF) method as a
# data quality label. DQR thresholds map to quality levels as follows
# (European Commission, 2013 EF method):
#
#   1.0 - 1.5  =>  Excellent
#   1.5 - 2.0  =>  Good
#   2.0 - 3.0  =>  Fair
#   3.0 - 4.0  =>  Poor
#   4.0 - 5.0  =>  Very poor
#
# IMPORTANT NOTE ON VARIANCE VALUES
# ----------------------------------
# Unlike the ecoinvent pedigree matrix, the ILCD Handbook does NOT define
# sigma^2_ln variance values for Monte Carlo uncertainty propagation. The
# ILCD DQR system was designed as a data quality label, not a quantitative
# uncertainty model.
#
# To enable uncertainty propagation with ILCD scores, this module adapts
# the empirically validated sigma^2_ln values from Ciroth et al. (2016)
# for the three indicators that overlap between the two systems:
#
#   ILCD indicator                  <- ecoinvent indicator used
#   Technological representativeness <- technology_correlation
#   Geographical representativeness  <- geographical_correlation
#   Temporal representativeness      <- temporal_correlation
#
# Users MUST acknowledge this adaptation in any peer-reviewed publication.
# This approach is consistent with Carlesso et al. (2024) who applied a
# hybrid pedigree/ILCD approach to plastic packaging LCA.
#
# REFERENCES
# ----------
# European Commission - JRC (2010). International Reference Life Cycle Data
# System (ILCD) Handbook - General guide for Life Cycle Assessment - Detailed
# guidance. EUR 24708 EN. Publications Office of the European Union,
# Luxembourg.
#
# European Commission (2013). Recommendations for the use of common methods
# to measure and communicate the life cycle environmental performance of
# products and organisations (Environmental Footprint methods). Official
# Journal of the European Union. 2013/179/EU.
#
# Ciroth, A., Muller, S., Weidema, B., Lesage, P. (2016). Empirically based
# uncertainty factors for the pedigree matrix in ecoinvent. International
# Journal of Life Cycle Assessment, 21(9), 1338-1348.
# https://doi.org/10.1007/s11367-013-0670-5
#
# Carlesso, A., et al. (2024). Data quality assessment of aggregated LCI
# datasets: A case study on fossil-based and bio-based plastic food packaging.
# Journal of Industrial Ecology. https://doi.org/10.1111/jiec.13522


# ---------------------------------------------------------------------------
# ILCD lookup table
# sigma^2_ln values adapted from ecoinvent pedigree (Ciroth et al., 2016)
# for the three ILCD representativeness indicators.
# See header note above regarding this adaptation.
# ---------------------------------------------------------------------------

#' @keywords internal
ILCD_TABLE <- list(
  technological_representativeness = c(
    "1" = 0.0000,
    "2" = 0.0001,
    "3" = 0.0008,
    "4" = 0.0080,
    "5" = 0.0120
  ),
  geographical_representativeness = c(
    "1" = 0.0000,
    "2" = 0.0001,
    "3" = 0.0002,
    "4" = 0.0002,
    "5" = 0.0040
  ),
  temporal_representativeness = c(
    "1" = 0.0000,
    "2" = 0.0001,
    "3" = 0.0002,
    "4" = 0.0008,
    "5" = 0.0080
  )
)

#' @keywords internal
ILCD_LABELS <- c(
  technological_representativeness = "Technological representativeness (TeR)",
  geographical_representativeness  = "Geographical representativeness (GR)",
  temporal_representativeness      = "Temporal representativeness (TiR)"
)

#' @keywords internal
DQR_LEVELS <- data.frame(
  from  = c(1.0, 1.5, 2.0, 3.0, 4.0),
  to    = c(1.5, 2.0, 3.0, 4.0, 5.0),
  level = c("Excellent", "Good", "Fair", "Poor", "Very poor"),
  stringsAsFactors = FALSE
)


# ---------------------------------------------------------------------------
# Core ILCD functions
# ---------------------------------------------------------------------------

#' Compute the ILCD Data Quality Rating (DQR)
#'
#' Calculates the composite Data Quality Rating as defined in the EU
#' Environmental Footprint method: the arithmetic mean of the three
#' representativeness indicator scores.
#'
#' @param technological_representativeness Integer. TeR score (1-5).
#' @param geographical_representativeness Integer. GR score (1-5).
#' @param temporal_representativeness Integer. TiR score (1-5).
#'
#' @return A named list with:
#'   \describe{
#'     \item{dqr}{Numeric. The composite DQR score.}
#'     \item{level}{Character. Quality level label (Excellent/Good/Fair/Poor/Very poor).}
#'   }
#'
#' @references
#' European Commission - JRC (2010). \emph{ILCD Handbook - Detailed guidance.}
#' EUR 24708 EN. Publications Office of the European Union.
#'
#' European Commission (2013). Environmental Footprint methods.
#' Official Journal 2013/179/EU.
#'
#' @examples
#' dqr_score(
#'   technological_representativeness = 2,
#'   geographical_representativeness  = 1,
#'   temporal_representativeness      = 3
#' )
#'
#' @export
dqr_score <- function(technological_representativeness,
                      geographical_representativeness,
                      temporal_representativeness) {

  scores <- c(technological_representativeness,
              geographical_representativeness,
              temporal_representativeness)

  for (s in scores) {
    if (!s %in% 1:5) {
      stop(sprintf("All scores must be integers between 1 and 5, got %s.", s))
    }
  }

  dqr <- mean(scores)

  level <- DQR_LEVELS$level[dqr >= DQR_LEVELS$from & dqr < DQR_LEVELS$to]
  if (length(level) == 0) level <- "Very poor"  # catches dqr == 5.0

  list(dqr = dqr, level = level)
}


#' Full ILCD pipeline: scores to uncertainty parameters
#'
#' Converts three ILCD representativeness scores into lognormal uncertainty
#' parameters, following the same approach as \code{\link{scores_to_gsd2}}
#' but using the ILCD indicator set.
#'
#' The sigma^2_ln variance values are adapted from the ecoinvent pedigree
#' table (Ciroth et al., 2016) for the three overlapping indicators. This
#' adaptation is required because the ILCD Handbook does not define
#' sigma^2_ln values for Monte Carlo uncertainty propagation.
#' See the module header for full methodological justification.
#'
#' @param technological_representativeness Integer. TeR score (1-5).
#' @param geographical_representativeness Integer. GR score (1-5).
#' @param temporal_representativeness Integer. TiR score (1-5).
#' @param basic_var Numeric. Basic uncertainty variance in log-space
#'   (sigma^2_ln). Default 0.0006 is typical for ecoinvent exchanges.
#'
#' @return A named list with:
#'   \describe{
#'     \item{combined_var}{Total sigma^2_ln.}
#'     \item{gsd2}{GSD^2 = exp(sigma^2_ln).}
#'     \item{gsd}{Geometric standard deviation.}
#'     \item{sigma_ln}{Std dev of ln(X); pass to \code{rlnorm()} as sdlog.}
#'     \item{basic_var}{The basic variance as supplied.}
#'     \item{indicator_vars}{Named vector of per-indicator variance contributions.}
#'     \item{dqr}{Composite ILCD Data Quality Rating.}
#'     \item{dqr_level}{Quality level label.}
#'   }
#'
#' @references
#' Ciroth, A., et al. (2016). Empirically based uncertainty factors for the
#' pedigree matrix in ecoinvent. \emph{International Journal of Life Cycle
#' Assessment}, 21(9), 1338-1348. https://doi.org/10.1007/s11367-013-0670-5
#'
#' Carlesso, A., et al. (2024). Data quality assessment of aggregated LCI
#' datasets. \emph{Journal of Industrial Ecology}.
#' https://doi.org/10.1111/jiec.13522
#'
#' @examples
#' result <- scores_to_gsd2_ilcd(
#'   technological_representativeness = 2,
#'   geographical_representativeness  = 1,
#'   temporal_representativeness      = 3,
#'   basic_var = 0.0006
#' )
#' result$gsd2
#' result$dqr
#' result$dqr_level
#'
#' @export
scores_to_gsd2_ilcd <- function(technological_representativeness,
                                 geographical_representativeness,
                                 temporal_representativeness,
                                 basic_var = 0.0006) {

  if (basic_var < 0) stop(sprintf("basic_var must be >= 0, got %s.", basic_var))

  scores <- list(
    technological_representativeness = technological_representativeness,
    geographical_representativeness  = geographical_representativeness,
    temporal_representativeness      = temporal_representativeness
  )

  for (ind in names(scores)) {
    if (!scores[[ind]] %in% 1:5) {
      stop(sprintf("Score for '%s' must be 1-5, got %s.", ind, scores[[ind]]))
    }
  }

  indicator_vars <- sapply(
    names(scores),
    function(ind) ILCD_TABLE[[ind]][as.character(scores[[ind]])]
  )

  combined_var <- basic_var + sum(indicator_vars)
  gsd2         <- exp(combined_var)
  dqr_result   <- dqr_score(
    technological_representativeness,
    geographical_representativeness,
    temporal_representativeness
  )

  list(
    combined_var   = combined_var,
    gsd2           = gsd2,
    gsd            = sqrt(gsd2),
    sigma_ln       = sqrt(combined_var),
    basic_var      = basic_var,
    indicator_vars = indicator_vars,
    dqr            = dqr_result$dqr,
    dqr_level      = dqr_result$level
  )
}


#' Batch process a data frame using ILCD scores
#'
#' Applies \code{\link{scores_to_gsd2_ilcd}} row-wise to a data frame
#' containing ILCD representativeness scores. Any additional columns
#' in the input are preserved.
#'
#' @param df A data frame with columns:
#'   \code{technological_representativeness},
#'   \code{geographical_representativeness},
#'   \code{temporal_representativeness}, and a basic uncertainty column.
#' @param basic_var_col Character. Name of the basic uncertainty column.
#'   Default \code{"basic_var"}.
#'
#' @return The original data frame with additional columns:
#'   \code{combined_var}, \code{gsd2}, \code{gsd}, \code{sigma_ln},
#'   \code{dqr}, \code{dqr_level}.
#'
#' @examples
#' df <- data.frame(
#'   exchange = c("N2O field", "Diesel"),
#'   technological_representativeness = c(2L, 1L),
#'   geographical_representativeness  = c(1L, 1L),
#'   temporal_representativeness      = c(3L, 2L),
#'   basic_var = c(0.0006, 0.0006)
#' )
#' batch_from_dataframe_ilcd(df)
#'
#' @export
batch_from_dataframe_ilcd <- function(df, basic_var_col = "basic_var") {

  required_cols <- c(names(ILCD_TABLE), basic_var_col)
  missing_cols  <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s",
                 paste(missing_cols, collapse = ", ")))
  }

  results <- lapply(seq_len(nrow(df)), function(i) {
    row <- df[i, ]
    out <- scores_to_gsd2_ilcd(
      technological_representativeness = as.integer(row$technological_representativeness),
      geographical_representativeness  = as.integer(row$geographical_representativeness),
      temporal_representativeness      = as.integer(row$temporal_representativeness),
      basic_var                        = as.numeric(row[[basic_var_col]])
    )
    data.frame(
      combined_var = out$combined_var,
      gsd2         = out$gsd2,
      gsd          = out$gsd,
      sigma_ln     = out$sigma_ln,
      dqr          = out$dqr,
      dqr_level    = out$dqr_level,
      stringsAsFactors = FALSE
    )
  })

  cbind(df, do.call(rbind, results))
}


#' Print a readable ILCD uncertainty summary
#'
#' @param result A list returned by \code{\link{scores_to_gsd2_ilcd}}.
#' @param exchange_name Character. Label shown in the header.
#'
#' @return Invisibly returns \code{result}.
#'
#' @examples
#' result <- scores_to_gsd2_ilcd(2, 1, 3, basic_var = 0.0006)
#' ilcd_summary(result, exchange_name = "N2O field emission")
#'
#' @export
ilcd_summary <- function(result, exchange_name = "Exchange") {
  sep <- strrep("\u2500", 52)
  cat("\n", sep, "\n", sep = "")
  cat("  ILCD uncertainty summary:", exchange_name, "\n")
  cat(sep, "\n")
  cat(sprintf("  %-38s %.6f\n", "Basic variance", result$basic_var))
  cat("\n")
  for (ind in names(result$indicator_vars)) {
    label <- ILCD_LABELS[ind]
    cat(sprintf("  %-38s %.6f\n", label, result$indicator_vars[ind]))
  }
  cat(sep, "\n")
  cat(sprintf("  %-38s %.6f\n", "Combined variance (sigma^2_ln)", result$combined_var))
  cat(sprintf("  %-38s %.4f\n",  "GSD^2",                          result$gsd2))
  cat(sprintf("  %-38s %.4f\n",  "GSD",                            result$gsd))
  cat(sprintf("  %-38s %.6f\n", "sigma_ln (for MC sampling)",      result$sigma_ln))
  cat(sep, "\n")
  cat(sprintf("  %-38s %.2f (%s)\n", "ILCD DQR", result$dqr, result$dqr_level))
  cat(sep, "\n")
  cat("  NOTE: sigma^2_ln values adapted from ecoinvent pedigree\n")
  cat("  (Ciroth et al., 2016). See module header for details.\n")
  cat(sep, "\n")
  invisible(result)
}


# ---------------------------------------------------------------------------
# Comparison utility
# ---------------------------------------------------------------------------

#' Compare ecoinvent pedigree and ILCD results side by side
#'
#' A convenience function for studies that need to report both approaches.
#' Takes scores for the shared indicators and shows both GSD^2 values and
#' the ILCD DQR in a single summary.
#'
#' Shared indicators (scored identically in both systems):
#'   - temporal_correlation / temporal_representativeness
#'   - geographical_correlation / geographical_representativeness
#'   - technology_correlation / technological_representativeness
#'
#' @param reliability Integer. Ecoinvent only: reliability score (1-5).
#' @param completeness Integer. Ecoinvent only: completeness score (1-5).
#' @param temporal Integer. Shared temporal score (1-5).
#' @param geographical Integer. Shared geographical score (1-5).
#' @param technology Integer. Shared technology score (1-5).
#' @param basic_var Numeric. Basic uncertainty variance in log-space.
#'
#' @return Invisibly returns a list with both result objects.
#'
#' @examples
#' compare_pedigree_ilcd(
#'   reliability  = 2, completeness = 2,
#'   temporal     = 3, geographical = 1, technology = 2,
#'   basic_var    = 0.0006
#' )
#'
#' @export
compare_pedigree_ilcd <- function(reliability,
                                   completeness,
                                   temporal,
                                   geographical,
                                   technology,
                                   basic_var = 0.0006) {

  # Source pedigree functions if not already loaded
  if (!exists("scores_to_gsd2")) {
    stop("Please source pedigree.R before using compare_pedigree_ilcd().")
  }

  pedigree_result <- scores_to_gsd2(
    reliability              = reliability,
    completeness             = completeness,
    temporal_correlation     = temporal,
    geographical_correlation = geographical,
    technology_correlation   = technology,
    basic_var                = basic_var
  )

  ilcd_result <- scores_to_gsd2_ilcd(
    technological_representativeness = technology,
    geographical_representativeness  = geographical,
    temporal_representativeness      = temporal,
    basic_var                        = basic_var
  )

  sep <- strrep("\u2500", 58)
  cat("\n", sep, "\n", sep = "")
  cat("  Ecoinvent pedigree vs ILCD — side-by-side comparison\n")
  cat(sep, "\n")
  cat(sprintf("  %-30s %12s %12s\n", "Parameter", "Ecoinvent", "ILCD"))
  cat(strrep("\u2500", 58), "\n")
  cat(sprintf("  %-30s %12.6f %12.6f\n",
              "Combined variance (sigma^2_ln)",
              pedigree_result$combined_var, ilcd_result$combined_var))
  cat(sprintf("  %-30s %12.4f %12.4f\n",
              "GSD^2",
              pedigree_result$gsd2, ilcd_result$gsd2))
  cat(sprintf("  %-30s %12.4f %12.4f\n",
              "GSD",
              pedigree_result$gsd, ilcd_result$gsd))
  cat(sprintf("  %-30s %12.6f %12.6f\n",
              "sigma_ln",
              pedigree_result$sigma_ln, ilcd_result$sigma_ln))
  cat(sprintf("  %-30s %12s %12s\n",
              "ILCD DQR", "—",
              sprintf("%.2f (%s)", ilcd_result$dqr, ilcd_result$dqr_level)))
  cat(sep, "\n")
  cat("  Ecoinvent: 5 indicators (reliability, completeness,\n")
  cat("             temporal, geographical, technology)\n")
  cat("  ILCD:      3 indicators (temporal, geographical, technology)\n")
  cat("             sigma^2_ln values adapted from Ciroth et al. (2016)\n")
  cat(sep, "\n")

  invisible(list(pedigree = pedigree_result, ilcd = ilcd_result))
}

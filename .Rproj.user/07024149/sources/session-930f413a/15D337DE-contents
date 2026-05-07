# run_analysis.R
# Full pipeline test — pedigree uncertainty to Monte Carlo
# --------------------------------------------------------

#source("pedigree_R.R")
#source("monte_carlo.R")

# --- Inventory ---
inventory <- data.frame(
  food_item                = c("Wheat bread", "Wheat bread",
                               "Beef burger", "Beef burger"),
  stage                    = c("Agriculture", "Processing",
                               "Agriculture", "Processing"),
  amount                   = c(0.8, 0.3, 22.0, 1.5),
  reliability              = c(2L, 2L, 3L, 2L),
  completeness             = c(2L, 2L, 2L, 2L),
  temporal_correlation     = c(3L, 2L, 3L, 2L),
  geographical_correlation = c(1L, 1L, 2L, 1L),
  technology_correlation   = c(2L, 2L, 3L, 2L),
  basic_var                = rep(0.0006, 4)
)

# --- Full pipeline ---
results <- inventory |>
  batch_from_dataframe() |>
  run_monte_carlo(n = 10000)

# --- Inspect ---
results$total


# --- Runing a conservative scenario ----

inventory_conservative <- data.frame(
  food_item                = c("Wheat bread", "Wheat bread",
                               "Beef burger", "Beef burger"),
  stage                    = c("Agriculture", "Processing",
                               "Agriculture", "Processing"),
  amount                   = c(0.8, 0.3, 22.0, 1.5),
  reliability              = c(2L, 2L, 4L, 2L),   # beef agriculture → 4
  completeness             = c(2L, 2L, 3L, 2L),   # beef agriculture → 3
  temporal_correlation     = c(3L, 2L, 4L, 2L),   # beef agriculture → 4
  geographical_correlation = c(1L, 1L, 3L, 1L),   # beef agriculture → 3
  technology_correlation   = c(2L, 2L, 4L, 2L),   # beef agriculture → 4
  basic_var                = c(0.0006, 0.0006, 0.006, 0.0006)  # higher basic_var for beef
)

results_conservative <- inventory_conservative |>
  batch_from_dataframe() |>
  run_monte_carlo(n = 10000)

results_conservative$total

# --- figures comparison ----

library(ggplot2)

# Add a numeric y position with manual per-scenario offset
comparison$y_pos <- as.numeric(factor(comparison$food_item)) +
  ifelse(comparison$scenario == "Optimistic DQI",       0.20,
         ifelse(comparison$scenario == "Conservative DQI",     0.00,
                -0.20))

food_levels <- levels(factor(comparison$food_item))

ggplot(comparison,
       aes(x = median, y = y_pos, colour = scenario, shape = scenario)) +
  geom_errorbarh(
    aes(xmin = .data[[lo_col]], xmax = .data[[hi_col]]),
    height = 0.12, linewidth = 0.8
  ) +
  geom_point(size = 3) +
  scale_y_continuous(
    breaks = seq_along(food_levels),
    labels = food_levels
  ) +
  scale_colour_manual(values = c(
    "Optimistic DQI"     = "#2e7d52",
    "Conservative DQI"   = "#e65100",
    "GSD² direct (1.22)" = "#1565c0"
  )) +
  labs(
    title    = "Sensitivity to data quality assumptions",
    subtitle = "95% CI from 10,000 Monte Carlo iterations",
    x        = "GWP100 (kg CO2-eq per kg food)",
    y        = NULL,
    colour   = "Scenario",
    shape    = "Scenario"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position    = "bottom"
  ) + facet_wrap(~ food_item, scales = "free", ncol = 1)

# ---------------------------------------------------------------------------
# ILCD walkthrough examples
# ---------------------------------------------------------------------------
source("pedigree_ilcd.R")

# DQR score examples
dqr_score(2, 1, 3)   # DQR = 2.0, "Fair"
dqr_score(1, 1, 1)   # DQR = 1.0, "Excellent"
dqr_score(5, 5, 5)   # DQR = 5.0, "Very poor"

# Full ILCD pipeline
result_ilcd <- scores_to_gsd2_ilcd(
  technological_representativeness = 2,
  geographical_representativeness  = 1,
  temporal_representativeness      = 3,
  basic_var = 0.0006
)
result_ilcd$gsd2
result_ilcd$dqr
result_ilcd$dqr_level

# Formatted summary
ilcd_summary(result_ilcd, exchange_name = "N2O field emission")

# Side-by-side comparison
compare_pedigree_ilcd(
  reliability  = 2, completeness = 2,
  temporal     = 3, geographical = 1, technology = 2,
  basic_var    = 0.0006
)


# Verify this directly
PEDIGREE_TABLE[["reliability"]]["2"]    # 0.0001
PEDIGREE_TABLE[["completeness"]]["2"]   # 0.0001


compare_pedigree_ilcd(
  reliability  = 5, completeness = 3,
  temporal     = 3, geographical = 2, technology = 3,
  basic_var    = 0.006
)


install.packages("writexl")
library(writexl)



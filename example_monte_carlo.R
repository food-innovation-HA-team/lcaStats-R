# example_monte_carlo.R
# ----------------------
# Full pipeline example: DQI scores → pedigree uncertainty → Monte Carlo
#
# This script demonstrates the complete lcaStats-R workflow using a
# simplified food system inventory (wheat bread vs beef burger).
#
# Run line by line in RStudio, or source() the whole file.
#
# Required packages:
#   install.packages(c("ggplot2", "writexl"))

# ---------------------------------------------------------------------------
# 0. Load modules
# ---------------------------------------------------------------------------

source("pedigree_R.R")
source("monte_carlo.R")
library(ggplot2)

# ---------------------------------------------------------------------------
# 1. Define inventory
#    food_item × stage, with DQI scores and mean impact values
#    Impact category: GWP100 (kg CO2-eq per kg food, product stage only)
# ---------------------------------------------------------------------------

inventory <- data.frame(
  food_item = c(
    "Wheat bread", "Wheat bread", "Wheat bread",
    "Beef burger",  "Beef burger",  "Beef burger"
  ),
  stage = c(
    "Agriculture", "Processing", "Packaging",
    "Agriculture", "Processing", "Packaging"
  ),
  # Mean GWP100 values (kg CO2-eq per kg food)
  amount = c(
    0.45, 0.18, 0.05,   # Wheat bread
    14.5, 1.20, 0.15    # Beef burger
  ),
  # DQI scores (1 = best, 5 = worst)
  reliability              = c(2L, 2L, 3L, 3L, 2L, 3L),
  completeness             = c(2L, 2L, 2L, 2L, 2L, 2L),
  temporal_correlation     = c(2L, 2L, 3L, 3L, 2L, 3L),
  geographical_correlation = c(1L, 1L, 2L, 2L, 1L, 2L),
  technology_correlation   = c(2L, 2L, 2L, 3L, 2L, 3L),
  basic_var                = rep(0.0006, 6)
)

# ---------------------------------------------------------------------------
# 2. Step 1 — Pedigree uncertainty
#    Converts DQI scores to sigma_ln for each exchange
# ---------------------------------------------------------------------------

inventory_with_uncertainty <- batch_from_dataframe(inventory)

# Inspect
print(inventory_with_uncertainty[, c("food_item", "stage", "amount",
                                      "gsd2", "sigma_ln")])

# ---------------------------------------------------------------------------
# 3. Step 2 — Monte Carlo simulation
#    Runs 10,000 lognormal draws per exchange, sums within stages and totals
# ---------------------------------------------------------------------------

results <- run_monte_carlo(
  df   = inventory_with_uncertainty,
  n    = 10000,
  seed = 42,
  ci   = 0.95
)

# Total impact summary
cat("\n=== Total GWP100 impacts (kg CO2-eq per kg food) ===\n")
print(results$total)

# Stage-level summary for beef burger
cat("\n=== Stage contributions — Beef burger ===\n")
print(results$stages[["Agriculture"]][
  results$stages[["Agriculture"]]$food_item == "Beef burger", ])

# ---------------------------------------------------------------------------
# 4. Visualisation
# ---------------------------------------------------------------------------

# Total uncertainty comparison
p1 <- plot_mc_total(
  results,
  title   = "Monte Carlo results — product stage GWP100",
  x_label = "GWP100 (kg CO2-eq per kg food)"
)
print(p1)

# Stage contributions for beef burger
p2 <- plot_mc_stages(
  results,
  food    = "Beef burger",
  y_label = "GWP100 (kg CO2-eq per kg food)"
)
print(p2)

# ---------------------------------------------------------------------------
# 5. Export to Excel
#    Same structure as the original Python notebook output
# ---------------------------------------------------------------------------

# install.packages("writexl")  # uncomment if needed
export_mc_results(results, file = "mc_results_food_example.xlsx")

# ---------------------------------------------------------------------------
# 6. Alternative: supply GSD² values directly
#    Use this if you already have GSD² from ecoinvent or SimaPro
# ---------------------------------------------------------------------------

mean_df <- data.frame(
  food_item   = c("Wheat bread", "Beef burger"),
  Agriculture = c(0.45, 14.5),
  Processing  = c(0.18, 1.20),
  Packaging   = c(0.05, 0.15)
)

gsd2_df <- data.frame(
  food_item   = c("Wheat bread", "Beef burger"),
  Agriculture = c(1.08, 1.22),
  Processing  = c(1.05, 1.10),
  Packaging   = c(1.06, 1.08)
)

results_from_gsd2 <- run_monte_carlo_gsd2(mean_df, gsd2_df, n = 10000)
print(results_from_gsd2$total)

# ---------------------------------------------------------------------------
# 7. Real-world usage — loading inventory from CSV
# ---------------------------------------------------------------------------

# In a real study, load your inventory from a CSV file:
# inventory <- read.csv("my_inventory.csv")
#
# Expected columns:
#   food_item, stage, amount,
#   reliability, completeness, temporal_correlation,
#   geographical_correlation, technology_correlation, basic_var
#
# Then run the full pipeline:
# results <- inventory |>
#   batch_from_dataframe() |>
#   run_monte_carlo(n = 10000)

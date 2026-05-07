# lcaStats (R)

**Reproducible statistics for Life Cycle Assessment — R implementation.**

Part of the [lcaStats](https://github.com/food-innovation-HA-team) toolkit.  
A Python version with an equivalent Streamlit app is also available in the same organisation.

---

## What this does

Converts ecoinvent and ILCD **Data Quality Indicator (DQI) scores** into lognormal 
uncertainty parameters using the pedigree matrix approach (Weidema et al., 2013),
and propagates that uncertainty through a full Monte Carlo simulation pipeline.

### The five ecoinvent DQIs

| Indicator | Question |
|---|---|
| Reliability | How well does the source measure what we need? |
| Completeness | How fully does the data cover the relevant population? |
| Temporal correlation | How recent is the data? |
| Geographical correlation | How geographically representative is the data? |
| Technology correlation | How technologically representative is the data? |

Scores run from **1** (best) to **5** (worst). The combined variance is:

```
σ²_total = σ²_basic + Σ σ²_indicator_i
```

converted to GSD² (the ecoinvent uncertainty metric):

```
GSD² = exp(σ²_total)
```

The ILCD / EU Environmental Footprint system uses three indicators
(technological, geographical, and temporal representativeness) and additionally
computes a Data Quality Rating (DQR) for reporting purposes.

---
## Installation

Clone the repo:

    git clone https://github.com/ESique82/lcaStats-R.git

Open `HFI_Digital.Rproj` in RStudio. Then install the required packages if needed:

    install.packages(c("shiny", "ggplot2", "writexl"))

## Running the Shiny app

    shiny::runApp("shiny_app_v2.R")

## Running the pipeline directly

    source("pedigree_R.R")
    source("monte_carlo.R")
    source("example_monte_carlo.R")
---    
## Quickstart

```r
source("pedigree_R.R")

# Single exchange
result <- scores_to_gsd2(
  reliability              = 2,
  completeness             = 2,
  temporal_correlation     = 3,
  geographical_correlation = 1,
  technology_correlation   = 2,
  basic_var                = 0.0006
)

pedigree_summary(result, exchange_name = "N2O field emission")

# Monte Carlo sampling
params  <- lognormal_params(mean = 0.014, sigma_ln = result$sigma_ln)
samples <- rlnorm(10000, meanlog = params["mu_ln"], sdlog = params["sigma_ln"])
hist(samples, breaks = 60, main = "Lognormal uncertainty distribution")
```

### Full Monte Carlo pipeline

```r
source("pedigree_R.R")
source("monte_carlo.R")

inventory <- read.csv("my_inventory.csv")

results <- inventory |>
  batch_from_dataframe() |>
  run_monte_carlo(n = 10000)

results$total
plot_mc_total(results, x_label = "GWP100 (kg CO2-eq per kg food)")
```

**Expected inventory CSV columns:**  
`food_item`, `stage`, `amount`, `reliability`, `completeness`,  
`temporal_correlation`, `geographical_correlation`, `technology_correlation`, `basic_var`

### ILCD / EU EF system

```r
source("pedigree_ilcd.R")

result <- scores_to_gsd2_ilcd(
  technological_representativeness = 2,
  geographical_representativeness  = 1,
  temporal_representativeness      = 3,
  basic_var = 0.0006
)

ilcd_summary(result, exchange_name = "N2O field emission")

# Side-by-side comparison with ecoinvent
compare_pedigree_ilcd(
  reliability = 2, completeness = 2,
  temporal = 3, geographical = 1, technology = 2,
  basic_var = 0.0006
)
```

### Interactive Shiny app

```r
shiny::runApp("shiny_app_v2.R")
```

---

## File structure

```
lcaStats-R/
├── pedigree_R.R          # Core ecoinvent pedigree functions
├── pedigree_ilcd.R       # ILCD / EU EF three-indicator system
├── monte_carlo.R         # Monte Carlo simulation pipeline
├── example_monte_carlo.R # Worked example — wheat bread vs beef burger
├── shiny_app_v2.R        # Interactive Shiny uncertainty explorer
├── run_analysis.R        # Full pipeline + DQI sensitivity analysis
└── README.md
```

---

## Output reference

### `scores_to_gsd2()` / `scores_to_gsd2_ilcd()`

| Field | Description |
|---|---|
| `combined_var` | Total σ²_ln (basic + all pedigree contributions) |
| `gsd2` | GSD² = exp(σ²_ln), the ecoinvent uncertainty metric |
| `gsd` | Geometric standard deviation = √GSD² |
| `sigma_ln` | Std dev of ln(X) — pass to `rlnorm()` as `sdlog` |
| `indicator_vars` | Named vector of per-indicator variance contributions |
| `dqr` | ILCD Data Quality Rating (ILCD system only) |
| `dqr_level` | DQR quality label: Excellent / Good / Fair / Poor / Very poor |

### `run_monte_carlo()`

| Field | Description |
|---|---|
| `results$total` | Summary statistics per food item (mean, median, sd, 95% CI, CV%) |
| `results$stages` | Named list of per-stage summaries |
| `results$samples` | Raw sample matrices for custom analysis |
| `results$meta` | Run metadata (n, seed, ci, food items, stages) |

---

## Notes

- Pedigree variance values are from Weidema et al. (2013), Table 5.
  Verify against the primary source for peer-reviewed work.
- ILCD σ²_ln values are adapted from Ciroth et al. (2016) for the three
  overlapping indicators. This adaptation must be declared in peer-reviewed work.
- Additive variance combination assumes independence between uncertainty
  sources — a known simplification of the ecoinvent approach.
- This module covers **parameter uncertainty** only. Model and scenario
  uncertainty are out of scope.

---

## References

Weidema, B.P., Bauer, C., Hischier, R., Mutel, C., Nemecek, T., Reinhard, J.,
Vadenbo, C.O., & Wernet, G. (2013). *Overview and methodology: Data quality
guideline for the ecoinvent database version 3.* Ecoinvent Report 1(v3).
The ecoinvent Centre, St. Gallen.

Ciroth, A., Muller, S., Weidema, B., & Lesage, P. (2016). Empirically based
uncertainty factors for the pedigree matrix in ecoinvent. *International
Journal of Life Cycle Assessment*, 21(9), 1338–1348.
https://doi.org/10.1007/s11367-013-0670-5

Igos, E., et al. (2019). How to treat uncertainties in life cycle assessment
studies? *International Journal of Life Cycle Assessment*, 24(7), 1261–1277.
https://doi.org/10.1007/s11367-018-1477-1

---

## Licence

MIT Licence · Harper Food Innovation: Digital · Harper Adams University

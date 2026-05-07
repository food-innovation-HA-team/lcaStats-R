# lcaStats (R)

**Reproducible statistics for Life Cycle Assessment — R implementation.**

Part of the [lcaStats](https://github.com/food-innovation-HA-team) toolkit.  
A Python version with an equivalent Streamlit app is also available in the same organisation.

---

## What this does

Converts ecoinvent **Data Quality Indicator (DQI) scores** into lognormal 
uncertainty parameters using the pedigree matrix approach (Weidema et al., 2013).

### The five DQIs

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

---

## Installation

Until the package is on CRAN, source the core file directly:

```r
source("R/pedigree.R")
```

Or install from GitHub using `remotes`:

```r
# install.packages("remotes")
remotes::install_github("food-innovation-HA-team/lcaStats-R")
```

---

## Quickstart

```r
source("R/pedigree.R")

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

### Batch processing

```r
df <- read.csv("my_inventory.csv")
results <- batch_from_dataframe(df)
write.csv(results, "results_with_uncertainty.csv", row.names = FALSE)
```

**Expected CSV columns:**  
`reliability`, `completeness`, `temporal_correlation`,  
`geographical_correlation`, `technology_correlation`, `basic_var`

Extra columns (e.g. `exchange_name`, `amount`, `unit`) are preserved.

### Interactive Shiny app

```r
shiny::runApp("inst/shiny/app.R")
```

---

## File structure

```
lcaStats-R/
├── R/
│   └── pedigree.R          # Core functions — source or library()
├── inst/shiny/
│   └── app.R               # Shiny interactive app
├── vignettes/
│   └── pedigree.Rmd        # Worked example with real data and plots
├── DESCRIPTION             # Package metadata
└── README.md               # This file
```

---

## Output reference

| Field | Description |
|---|---|
| `combined_var` | Total σ²_ln (basic + all pedigree contributions) |
| `gsd2` | GSD² = exp(σ²_ln), the ecoinvent uncertainty metric |
| `gsd` | Geometric standard deviation = √GSD² |
| `sigma_ln` | Std dev of ln(X) — pass to `rlnorm()` as `sdlog` |
| `indicator_vars` | Named vector of per-indicator variance contributions |

---

## Notes

- Values in the pedigree lookup table are from Weidema et al. (2013), Table 5. 
  Verify against the primary source for peer-reviewed work.
- Additive variance combination assumes independence between uncertainty 
  sources — a known simplification of the ecoinvent approach.
- This module covers **parameter uncertainty** only. Model and scenario 
  uncertainty are out of scope.

---

## Reference

Weidema, B.P., Bauer, C., Hischier, R., Mutel, C., Nemecek, T., Reinhard, J.,
Vadenbo, C.O., & Wernet, G. (2013). *Overview and methodology: Data quality
guideline for the ecoinvent database version 3.* Ecoinvent Report 1(v3).
The ecoinvent Centre, St. Gallen.

---

## Licence

MIT Licence · Harper Food Innovation: Digital · Harper Adams University

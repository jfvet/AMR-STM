# Statistical analyses — antibiotic-resistant bacteria carriage in urban wild small mammals

## Study context

This repository contains the R code used for statistical analyses in an M2 internship study investigating the carriage of antibiotic-resistant bacteria in wild small terrestrial mammals captured in the Lyon metropolitan area (France). The study is part of a broader One Health framework that considers humans, animals and the environment as interconnected compartments for antimicrobial resistance (AMR) circulation.

A total of 99 small mammals were sampled across 19 sites between October 2024 and May 2025, including brown rats (*Rattus norvegicus*), wood mice (*Apodemus sylvaticus*), house mice (*Mus musculus*) and field voles (*Microtus arvalis*). Three bacterial groups were targeted by culture on selective chromogenic media: ESBL-producing Enterobacterales, *Staphylococcus aureus* (including MRSA), and *Acinetobacter baumannii*.

The degree of urbanisation at each capture site was assessed using the Human Terrestrial Footprint index (HFP-100, 2020 version), a composite satellite-derived index integrating land use, population density, night-time light and road infrastructure at 100 m resolution.

The associated data and manuscript are described in the companion Zenodo record.

## File structure

```
analyses_AMR_small_mammals.R   main analysis script (this file)
data_individual_records.csv    individual-level dataset (not included — see Data Availability)
README.md                      this file
```

## What the script does

The script is organized into seven sequential sections:

1. **Package loading** — all required R packages are loaded at the top
2. **Data import and recoding** — the raw CSV is imported, variables are renamed for clarity, continuous predictors (weight, HFP) are centred and scaled, and habitat types are grouped into two categories (dense urban vs. green spaces) for univariate tests
3. **Descriptive statistics** — sample sizes by species, site and season; crude prevalence with 95% confidence intervals (Wald method)
4. **Bivariate tests** — Fisher's exact tests (sex, season, habitat type) and Kruskal-Wallis tests (body mass) for *R. norvegicus* carriage of BLSE-Enterobacterales and macrolide-resistant *S. aureus*; Fisher's exact test and Wilcoxon rank-sum test for *A. sylvaticus* and *A. baumannii*
5. **Site-level GLM with Firth correction** — because of quasi-complete separation (ESBL carriage observed exclusively at dense urban sites), a Firth-penalised binomial GLM (`brglm2` package) was used to estimate the effect of HFP-100 on site-level carriage probability; odds ratios and 95% confidence intervals are reported
6. **Individual-level GLMM (appendix)** — generalized linear mixed models with binomial family and site as random intercept (`lme4`); model selection by AICc using `dredge()` (`MuMIn`); these models are included as a methodological appendix documenting convergence failures and identifiability issues
7. **Diagnostics** — collinearity check between HFP and the site random effect (R² of HFP ~ site), events-per-variable calculation for EPV rule, and cluster size diagnostics explaining why a GEE approach was not retained

## How to use the code

### Requirements

R version 4.4 or later is recommended. The following packages are required:

```r
install.packages(c("tidyverse", "janitor", "lme4", "MuMIn",
                   "performance", "knitr", "brglm2"))
```

### Data

The individual-level dataset (`data_individual_records.csv`) must be placed in the working directory. The file uses semicolons as column separators. The minimum required columns are:

| Column | Description |
|--------|-------------|
| `loc_number` | site identifier |
| `loc_type` | habitat type (City, Sewage, Park, Zoo) |
| `hfp_100` | Human Footprint index value at the site |
| `species` | species name |
| `sex` | M / F |
| `season` | Autumn / Winter / Spring |
| `weight` | body mass in grams |
| `entero_blse` | ESBL carriage (0/1) |
| `staph_macro` | macrolide-resistant *S. aureus* carriage (0/1) |
| `acba` | *A. baumannii* carriage (0/1) |

The dataset is not included in this repository because individual capture coordinates and site identifiers could allow localisation of sampling sites. Aggregated data sufficient to reproduce the main results are available upon reasonable request to the corresponding authors.

### Running the script

Open the script in RStudio and run it section by section, or source it entirely:

```r
source("analyses_AMR_small_mammals.R")
```

The script produces printed tables and a bubble plot of ESBL prevalence as a function of HFP-100. No output files are written automatically; tables are formatted with `knitr::kable()` and are intended for use within an R Markdown document or interactive session.

## Key analytical choices and limitations

**Why Firth GLM at the site level rather than GLMM at the individual level for the main analysis?** The HFP index is constant within sites (one raster value per GPS coordinate), which makes it structurally collinear with the site random effect in a GLMM. Including both in the same model leads to non-identifiable fixed effects (R² of HFP ~ site = 1.00, ICC > 0.95). The site-level Firth GLM sidesteps this issue by aggregating data to the site level, at the cost of statistical power (n = 15 sites for the main analysis).

**Why Firth correction?** Complete or quasi-complete separation occurs when ESBL cases are observed exclusively in one habitat category. Standard logistic regression produces infinite or highly unstable maximum likelihood estimates in this situation. The Firth correction penalises the log-likelihood to obtain finite, less biased estimates.

**EPV and model parsimony.** With only 11 ESBL-positive individuals among 71 *R. norvegicus*, the events-per-variable ratio (EPV) falls below the recommended threshold of 10 for models with more than one predictor. As a consequence, only one predictor at a time could be included in individual-level models.

## Citation

If you use this code, please cite the associated manuscript (see Zenodo record DOI).

## License

MIT License — free to reuse and adapt with attribution.

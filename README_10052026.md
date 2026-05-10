# AMR carriage in urban wild small mammals

R analysis code accompanying a veterinary thesis on antimicrobial resistance (AMR)
carriage in wild small mammals in an urban environment.

---

## Repository contents

```
├── analysis.Rmd          # Main analysis script
├── README.md
└── data/                 # Input data (may be available on demand)
    ├── Data_indiv.csv
    └── Data_antibiograms.csv
```

> **Raw data files are not included** in this repository (unpublished data).
> They may be made available upon reasonable request after publication.

---

## Input data

### `Data_indiv.csv` 

| Column | Description |
|--------|-------------|
| `Species` | Host species (*Rattus norvegicus*, *Apodemus sylvaticus*) |
| `Loc_number` | Site identifier |
| `Loc_type` | Site category (`City`, `Sewage`, `Park`, `Zoo`) |
| `Loc_name` | Site name |
| `HFP.100` | Human Footprint Index (100 m radius) |
| `Entero.BLSE` | Rectal ESBL Enterobacteriaceae carriage (0/1) |
| `Staph.macro` | Oro-pharyngeal macrolide-resistant *S. aureus* carriage (0/1) |
| `Sex` | Animal sex (`M` / `F`) |
| `Weight` | Body weight (g) |

### `Data_antibiograms.csv` 

| Column | Description |
|--------|-------------|
| `Specimen_id` | Isolate identifier |
| `Bacterial_species` | Identified bacterial species |
| `Culture_medium` | Selective medium (`BLSE`, `SAID`, `SARM`, `ACBA`) |
| `Host_species` | Host species |
| `*_SIR` | Disk diffusion results (S / I / R) per antibiotic |

---

## Requirements

R 4.4.x. Install packages with:

```r
install.packages(c(
  "tidyverse", "janitor",
  "rstatix", "broom", "geepack", "MuMIn", "brglm2",
  "gtsummary", "flextable",
  "ggplot2", "ggrepel", "ggpubr",
  "pheatmap", "grid", "gtable", "gridExtra", "svglite",
  "AMR"
))
```

---

## Statistical notes

- **GEE instead of GLMM**: a mixed model with a random site effect suffered from
  collinearity with the fixed HFP-100 predictor (site-level variable). GEE with
  exchangeable correlation structure was used to account for within-site
  non-independence while estimating the marginal effect of HFP-100.

- **Firth estimator**: used for site-level logistic regression (small n, near
  complete separation) via the `brglm2` package.

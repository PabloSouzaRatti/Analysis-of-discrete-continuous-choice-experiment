# P2 analysis code and public dataset

This repository contains the public respondent-task dataset and separate scripts for preparing the estimation input, estimating the perceived-risk latent-class MNL, running the acreage intensity and supplementary analyses, and creating figures.

## Files

- `data/Dataset_public.csv` — released dataset, 2,682 task records from 447 respondents.
- `Dataset_Codebook.xlsx` — variable definitions and coding.
- `README_Public_dataset.docx` — dataset documentation.
- `scripts/01_prepare_estimation_data.R` — validates the public CSV and writes `data/estimation_ready.csv`.
- `scripts/02_fit_lcmn_risk.R` — estimates the two-class Apollo model and posterior class probabilities; writes `results/analysis_bundle.rds`.
- `scripts/03_intensity_analysis.R` — runs the conditional acreage models and writes the intensity tables and a model bundle for plotting.
- `scripts/04_extra_analyses.R` — estimates perceived-risk regressions and creates class-profile tables.
- `scripts/05_plot_figures.R` — creates the model figures. It expects the main model bundle; Figure 2 also uses the intensity bundle.

## Run order

From the repository root, run:

```sh
Rscript scripts/01_prepare_estimation_data.R
Rscript scripts/02_fit_lcmn_risk.R
Rscript scripts/03_intensity_analysis.R
Rscript scripts/04_extra_analyses.R
Rscript scripts/05_plot_figures.R
```

The scripts determine the repository path from their own location, so the commands can also be launched from another working directory by supplying the script path.

## R packages

The analysis uses `Apollo`, `dplyr`, `tidyr`, `stringr`, `tibble`, `sandwich`, `lmtest`, `broom`, `gt`, `purrr`, `ggplot2`, `MASS`, `scales`, and `modelsummary`. Install any missing packages before running the scripts. `Apollo` is used for latent-class estimation; the supplied script uses a single core and 30 randomized starting values (seed base 123), then re-estimates from the best log-likelihood. The original analysis reports Apollo 0.3.7; package versions in the current environment may differ.

## Data handling and interpretation notes

`Dataset_public.csv` is already in long respondent-task form and includes the model coding and standardized perceived-risk covariate used by the estimation script. The preparation script validates the released rows and coding, then writes the exact model input without changing or filtering records. It stops if required columns are absent, task rows are duplicated, or respondents do not have six tasks.

The public dataset does not contain `House_LandOwn` or `House_AgriExp`. These fields are omitted from the supplementary perceived-risk regression and class-profile table rather than being inferred from other variables. Consequently, those public-data supplementary specifications are reduced relative to analyses that used the non-public source file. The two scripts do not reconstruct these unavailable variables.

The released minimum coordinated area level of 75 ha is retained as coded in the dataset. The design-generation syntax itself is not needed to rerun the analyses because the final choice-task attributes are included in the public CSV.

The intensity script retains the supplied land-entry correction: when reported enrolled land exceeds owned land, the value is repeatedly divided by ten until it no longer exceeds available owned land. This is a consequential data-cleaning rule from the analysis code and should be reviewed against the questionnaire/export conventions before publication. The scripts print diagnostics for remaining out-of-range records.

## Outputs

Generated files are written beneath `results/` (Apollo estimation outputs, serialized bundles, intensity results, supplementary tables, and figures). The CSV used for estimation is written beneath `data/`. The repository `.gitignore` excludes generated results and that derived CSV; the original public dataset remains included.

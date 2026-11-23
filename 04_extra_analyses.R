# 04_extra_analyses.R
# Supplementary perceived-risk models and respondent class profile tables.
# House_LandOwn and House_AgriExp are omitted because they are not in Dataset_public.csv.

###########################################################################
#----------- Behavioural drivers & Ex-post class profiling------
#
######
# Behavioural drivers regressions

# =========================================================
# SUPPLEMENTARY ANALYSIS:
# Determinants of perceived economic risk
# ---------------------------------------------------------
# This script:
# 1. Builds a clean analysis dataset
# 2. Applies consistent coding across all models
# 3. Rescales land variables for interpretability
# 4. Estimates three nested OLS models with robust SE
# 5. Produces a clean supplementary regression table
# 6. Exports the table to HTML, Word, and LaTeX
# =========================================================

library(dplyr)
library(tidyr)
library(lmtest)
library(sandwich)
library(broom)
library(gt)

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value=TRUE)
if (!length(file_arg)) stop("Run with Rscript scripts/04_extra_analyses.R")
project_root <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash="/")), ".."), winslash="/", mustWork=TRUE)
setwd(project_root)
bundle <- readRDS(file.path(project_root, "results", "analysis_bundle.rds"))
database <- bundle$database
post <- bundle$post
post2 <- bundle$post2
out_dir <- file.path(project_root, "results", "extra")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
setwd(out_dir)

# ---------------------------------------------------------
# 1. Build clean analysis dataset
# ---------------------------------------------------------
analysis_df <- database %>%
  select(
    resp_id,
    econ_risk_std,
    House_Age,
    House_AgriEdu,
    House_AesExp,

    Stru_OwnLand,
    Stru_RentLand,
    Stru_HireEmp,
    Stru_FarType,
    IncSat_IncSource,
    IncSat_Overall,
    Soc_PeerAdoImp
  ) %>%
  distinct()

# ---------------------------------------------------------
# 2. Consistent coding across all models
# ---------------------------------------------------------
analysis_df <- analysis_df %>%
  mutate(
    # rescale land to 100 ha for interpretability
    Stru_OwnLand_100  = Stru_OwnLand / 100,
    Stru_RentLand_100 = Stru_RentLand / 100,
    
    # consistent factor coding
    House_AgriEdu = factor(House_AgriEdu),
    House_AesExp = factor(House_AesExp),
    Stru_FarType = factor(Stru_FarType),
    IncSat_IncSource = factor(IncSat_IncSource),
    IncSat_Overall = factor(IncSat_Overall),
    Soc_PeerAdoImp = factor(Soc_PeerAdoImp)
  )

# Optional: inspect levels
# lapply(analysis_df %>% select(where(is.factor)), levels)

# ---------------------------------------------------------
# 3. Estimate models
# ---------------------------------------------------------
# Model 1: parsimonious income-related drivers
model_1 <- lm(
  econ_risk_std ~
    IncSat_IncSource +
    IncSat_Overall,
  data = analysis_df
)

# Model 2: add structural and socio-demographic controls
model_2 <- lm(
  econ_risk_std ~
    IncSat_IncSource +
    IncSat_Overall +
    Stru_OwnLand_100 +
    Stru_HireEmp +
    Stru_FarType +
    House_Age,
  data = analysis_df
)

# Model 3: extended specification
model_3 <- lm(
  econ_risk_std ~
    IncSat_IncSource +
    IncSat_Overall +
    Stru_OwnLand_100 +
    Stru_RentLand_100 +
    Stru_HireEmp +
    Stru_FarType +
    House_Age +
    House_AgriEdu +
    House_AesExp +
    Soc_PeerAdoImp,
  data = analysis_df
)

# ---------------------------------------------------------
# 4. Robust summary function
# ---------------------------------------------------------
tidy_robust <- function(model) {
  rob <- coeftest(model, vcov = vcovHC(model, type = "HC1"))
  out <- broom::tidy(rob)
  
  out %>%
    mutate(
      stars = case_when(
        p.value < 0.01 ~ "***",
        p.value < 0.05 ~ "**",
        p.value < 0.10 ~ "*",
        TRUE ~ ""
      ),
      estimate_fmt = sprintf("%.3f%s", estimate, stars),
      se_fmt = sprintf("(%.3f)", std.error)
    ) %>%
    select(term, estimate_fmt, se_fmt, p.value)
}

# ---------------------------------------------------------
# 5. Nice variable labels
# ---------------------------------------------------------
label_terms <- function(x) {
  x <- gsub("\\(Intercept\\)", "(Intercept)", x)
  x <- gsub("^IncSat_IncSource", "Income source: ", x)
  x <- gsub("^IncSat_Overall", "Income satisfaction: ", x)
  x <- gsub("^Stru_FarType", "Farm type: ", x)
  x <- gsub("^Soc_PeerAdoImp", "Importance of peer adoption: ", x)
  
  x <- gsub("House_Age", "Age", x)
  x <- gsub("House_AgriEduNej", "Agricultural education: No", x)
  x <- gsub("House_AesExpNej", "Previous AES experience: No", x)
  
  x <- gsub("Stru_OwnLand_100", "Owned land (100 ha)", x)
  x <- gsub("Stru_RentLand_100", "Rented land (100 ha)", x)
  x <- gsub("Stru_HireEmp", "Hired employees", x)
  
  x <- gsub("Income source: Full_time", "Income source: Full-time", x)
  x <- gsub("Income source: Part_full", "Income source: Part-time", x)
  x <- gsub("Income source: Part_hobby", "Income source: Part-hobby", x)
  x <- gsub("Income source: Full_hobby", "Income source: Full-hobby (ref.)", x)
  
  x <- gsub("Income satisfaction: Meget utilfreds", "Income satisfaction: Very unsatisfied", x)
  x <- gsub("Income satisfaction: Utilfreds", "Income satisfaction: Unsatisfied", x)
  x <- gsub("Income satisfaction: Tilfreds", "Income satisfaction: Satisfied", x)
  x <- gsub("Income satisfaction: Meget tilfreds", "Income satisfaction: Very satisfied (ref.)", x)
  
  x <- gsub("Farm type: Agerbrug i øvrigt", "Farm type: Other crops", x)
  x <- gsub("Farm type: Andet husdyrbrug", "Farm type: Other livestock", x)
  x <- gsub("Farm type: Blandet husdyrbrug med planteavl", "Farm type: Mixed livestock and crops", x)
  x <- gsub("Farm type: Frugttræer", "Farm type: Fruit trees", x)
  x <- gsub("Farm type: Gartneri", "Farm type: Horticulture", x)
  x <- gsub("Farm type: Korn og industrifrø", "Farm type: Grain and industrial seed", x)
  x <- gsub("Farm type: Kvæg \\? opdræt/fedning", "Farm type: Beef/fattening", x)
  x <- gsub("Farm type: Malkekvæg", "Farm type: Dairy", x)
  x <- gsub("Farm type: Svin", "Farm type: Pork", x)
  
  x <- gsub("Importance of peer adoption: 2", "Importance of peer adoption: 2", x)
  x <- gsub("Importance of peer adoption: 3", "Importance of peer adoption: 3", x)
  x <- gsub("Importance of peer adoption: 4", "Importance of peer adoption: 4", x)
  x <- gsub("Importance of peer adoption: 5 \\(Meget vigtig\\)", "Importance of peer adoption: 5 (Very important)", x)
  
  x
}

# ---------------------------------------------------------
# 6. Extract robust tables
# ---------------------------------------------------------
t1 <- tidy_robust(model_1) %>%
  rename(m1_est = estimate_fmt, m1_se = se_fmt)

t2 <- tidy_robust(model_2) %>%
  rename(m2_est = estimate_fmt, m2_se = se_fmt)

t3 <- tidy_robust(model_3) %>%
  rename(m3_est = estimate_fmt, m3_se = se_fmt)

# ---------------------------------------------------------
# 7. Merge regression results
# ---------------------------------------------------------
reg_table <- full_join(t1, t2, by = "term") %>%
  full_join(t3, by = "term") %>%
  select(term, m1_est, m1_se, m2_est, m2_se, m3_est, m3_se) %>%
  mutate(
    term = label_terms(term)
  )

# ---------------------------------------------------------
# 8. Add model statistics
# ---------------------------------------------------------
get_nobs <- function(model) nobs(model)
get_r2 <- function(model) summary(model)$r.squared
get_adj_r2 <- function(model) summary(model)$adj.r.squared
get_f_p <- function(model) {
  fs <- summary(model)$fstatistic
  if (is.null(fs)) return(NA_real_)
  pf(fs[1], fs[2], fs[3], lower.tail = FALSE)
}

model_stats <- tibble(
  term = c("Observations", "R-squared", "Adjusted R-squared", "Model p-value"),
  m1_est = c(
    as.character(get_nobs(model_1)),
    sprintf("%.3f", get_r2(model_1)),
    sprintf("%.3f", get_adj_r2(model_1)),
    sprintf("%.5f", get_f_p(model_1))
  ),
  m1_se = c("", "", "", ""),
  m2_est = c(
    as.character(get_nobs(model_2)),
    sprintf("%.3f", get_r2(model_2)),
    sprintf("%.3f", get_adj_r2(model_2)),
    sprintf("%.5f", get_f_p(model_2))
  ),
  m2_se = c("", "", "", ""),
  m3_est = c(
    as.character(get_nobs(model_3)),
    sprintf("%.3f", get_r2(model_3)),
    sprintf("%.3f", get_adj_r2(model_3)),
    sprintf("%.5f", get_f_p(model_3))
  ),
  m3_se = c("", "", "", "")
)

# ---------------------------------------------------------
# 9. Order rows
# ---------------------------------------------------------
preferred_order <- c(
  "(Intercept)",
  "Income source: Full-time",
  "Income source: Part-time",
  "Income source: Part-hobby",
  "Income satisfaction: Very unsatisfied",
  "Income satisfaction: Satisfied",
  "Income satisfaction: Unsatisfied",
  "Owned land (100 ha)",
  "Rented land (100 ha)",
  "Hired employees",
  "Age",
  "Agricultural education: No",
  "Previous AES experience: No",
  "Land ownership: Fully owned",
  "Land ownership: Mixed ownership",
  "Farm type: Other livestock",
  "Farm type: Mixed livestock and crops",
  "Farm type: Fruit trees",
  "Farm type: Horticulture",
  "Farm type: Grain and industrial seed",
  "Farm type: Beef/fattening",
  "Farm type: Dairy",
  "Farm type: Pork",
  "Importance of peer adoption: 2",
  "Importance of peer adoption: 3",
  "Importance of peer adoption: 4",
  "Importance of peer adoption: 5 (Very important)"
)

reg_table <- reg_table %>%
  mutate(order = match(term, preferred_order)) %>%
  mutate(order = ifelse(is.na(order), 999, order)) %>%
  arrange(order, term) %>%
  select(-order)

# ---------------------------------------------------------
# 10. Add statistics at bottom
# ---------------------------------------------------------
final_table <- bind_rows(reg_table, model_stats)

# replace missing with dash
final_table[is.na(final_table)] <- "-"

# ---------------------------------------------------------
# 11. Build GT table
# ---------------------------------------------------------
gt_reg <- final_table %>%
  gt() %>%
  tab_header(
    title = "Determinants of perceived economic risk"
  ) %>%
  tab_spanner(
    label = "Model 1",
    columns = c(m1_est, m1_se)
  ) %>%
  tab_spanner(
    label = "Model 2",
    columns = c(m2_est, m2_se)
  ) %>%
  tab_spanner(
    label = "Model 3",
    columns = c(m3_est, m3_se)
  ) %>%
  cols_label(
    term = "Variable",
    m1_est = "Coef.",
    m1_se = "Std. Error",
    m2_est = "Coef.",
    m2_se = "Std. Error",
    m3_est = "Coef.",
    m3_se = "Std. Error"
  ) %>%
  cols_align(
    align = "left",
    columns = term
  ) %>%
  cols_align(
    align = "center",
    columns = c(m1_est, m1_se, m2_est, m2_se, m3_est, m3_se)
  ) %>%
  tab_source_note(
    source_note = md("Robust standard errors in parentheses. * p < 0.10, ** p < 0.05, *** p < 0.01.")
  )

# show
gt_reg

# ---------------------------------------------------------
# 12. Export
# ---------------------------------------------------------
gtsave(gt_reg, "supplementary_risk_drivers_table.html")
gtsave(gt_reg, "supplementary_risk_drivers_table.tex")
gtsave(gt_reg, "supplementary_risk_drivers_table.docx")

# ---------------------------------------------------------
# 13. Optional: inspect model summaries
# ---------------------------------------------------------
summary(model_1)
summary(model_2)
summary(model_3)




# =========================================================
# TABLE 5. CLASS PROFILE WITH SIGNIFICANCE LEVELS
# Corrected version with category harmonization
# =========================================================

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(gt)

# ---------------------------------------------------------
# 1. BUILD RESPONDENT-LEVEL ANALYSIS DATASET
# ---------------------------------------------------------

resp_df <- post %>%
  distinct(resp_id, most_likely_class, .keep_all = TRUE)

kbau_df <- post2 %>%
  distinct(resp_id, share_KBAU)

covars <- database %>%
  select(
    resp_id,
    econ_risk_std,
    IncSat_Overall,
    IncSat_IncSource,
    Total_Land,
    House_AesExp,
    House_Age,
    House_AgriEdu,


    Stru_FarType,
    Stru_HireEmp,
    Soc_BetEnvPerf,
    Soc_ImpPeer,
    Soc_PeerAdo,
    Soc_PeerAdoImp,
    Soc_PeerAdoInf,
    Soc_WorseInc
  ) %>%
  distinct()

analysis_df <- resp_df %>%
  left_join(kbau_df, by = "resp_id") %>%
  left_join(covars, by = "resp_id")

# ---------------------------------------------------------
# 2. HARMONIZE CATEGORY LABELS / FIX TYPOS
# ---------------------------------------------------------

analysis_df <- analysis_df %>%
  mutate(
    across(where(is.character), str_squish),
    
    # Fix typo in land ownership
    # Harmonize "don't know"
    Soc_PeerAdo = case_when(
      Soc_PeerAdo %in% c("Ved Ikke", "Ved ikke") ~ "Ved ikke",
      TRUE ~ Soc_PeerAdo
    )
  )

# ---------------------------------------------------------
# 3. CONVERT TO FACTORS
# ---------------------------------------------------------

analysis_df <- analysis_df %>%
  mutate(
    most_likely_class = as.integer(most_likely_class),
    
    IncSat_Overall   = as.factor(IncSat_Overall),
    IncSat_IncSource = as.factor(IncSat_IncSource),
    House_AesExp     = as.factor(House_AesExp),
    House_AgriEdu    = as.factor(House_AgriEdu),
    Stru_FarType     = as.factor(Stru_FarType),
    Soc_BetEnvPerf   = as.factor(Soc_BetEnvPerf),
    Soc_ImpPeer      = as.factor(Soc_ImpPeer),
    Soc_PeerAdo      = as.factor(Soc_PeerAdo),
    Soc_PeerAdoImp   = as.factor(Soc_PeerAdoImp),
    Soc_PeerAdoInf   = as.factor(Soc_PeerAdoInf),
    Soc_WorseInc     = as.factor(Soc_WorseInc)
  )

# ---------------------------------------------------------
# 4. OPTIONAL: SET FACTOR ORDER
# ---------------------------------------------------------

if ("IncSat_Overall" %in% names(analysis_df)) {
  lev <- c("Meget tilfreds", "Tilfreds", "Utilfreds", "Meget utilfreds")
  analysis_df$IncSat_Overall <- factor(
    analysis_df$IncSat_Overall,
    levels = lev[lev %in% unique(as.character(analysis_df$IncSat_Overall))]
  )
}

if ("IncSat_IncSource" %in% names(analysis_df)) {
  lev <- c("Full_hobby", "Full_time", "Part_full", "Part_hobby")
  analysis_df$IncSat_IncSource <- factor(
    analysis_df$IncSat_IncSource,
    levels = lev[lev %in% unique(as.character(analysis_df$IncSat_IncSource))]
  )
}

if ("House_AesExp" %in% names(analysis_df)) {
  lev <- c("Ja", "Nej")
  analysis_df$House_AesExp <- factor(
    analysis_df$House_AesExp,
    levels = lev[lev %in% unique(as.character(analysis_df$House_AesExp))]
  )
}

if ("House_AgriEdu" %in% names(analysis_df)) {
  lev <- c("Ja", "Nej")
  analysis_df$House_AgriEdu <- factor(
    analysis_df$House_AgriEdu,
    levels = lev[lev %in% unique(as.character(analysis_df$House_AgriEdu))]
  )
}

if ("Soc_PeerAdo" %in% names(analysis_df)) {
  lev <- c("Ja", "Nej", "Ved ikke")
  analysis_df$Soc_PeerAdo <- factor(
    analysis_df$Soc_PeerAdo,
    levels = lev[lev %in% unique(as.character(analysis_df$Soc_PeerAdo))]
  )
}

# ---------------------------------------------------------
# 5. VARIABLES TO PROFILE
# ---------------------------------------------------------

continuous_vars <- c(
  "econ_risk_std",
  "share_KBAU",
  "Total_Land",
  "House_Age",
  "Stru_HireEmp"
)

categorical_vars <- c(
  "IncSat_Overall",
  "IncSat_IncSource",
  "House_AesExp",
  "House_AgriEdu",
  "Stru_FarType",
  "Soc_BetEnvPerf",
  "Soc_ImpPeer",
  "Soc_PeerAdo",
  "Soc_PeerAdoImp",
  "Soc_PeerAdoInf",
  "Soc_WorseInc"
)

continuous_vars <- continuous_vars[continuous_vars %in% names(analysis_df)]
categorical_vars <- categorical_vars[categorical_vars %in% names(analysis_df)]

# ---------------------------------------------------------
# 6. VARIABLE LABELS
# ---------------------------------------------------------

var_labels <- c(
  econ_risk_std    = "Economic risk (std.)",
  share_KBAU       = "Observed KBAU share",
  Total_Land       = "Total land",
  House_AesExp     = "House: AesExp",
  House_Age        = "House: Age",
  House_AgriEdu    = "House: AgriEdu",
  Stru_FarType     = "Structure: FarType",
  Stru_HireEmp     = "Structure: HireEmp",
  IncSat_Overall   = "Income satisfaction: overall",
  IncSat_IncSource = "Income source",
  Soc_BetEnvPerf   = "Social: BetEnvPerf",
  Soc_ImpPeer      = "Social: ImpPeer",
  Soc_PeerAdo      = "Social: PeerAdo",
  Soc_PeerAdoImp   = "Social: PeerAdoImp",
  Soc_PeerAdoInf   = "Social: PeerAdoInf",
  Soc_WorseInc     = "Social: WorseInc"
)

label_or_name <- function(x) {
  if (x %in% names(var_labels)) var_labels[[x]] else x
}

sig_stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.01 ~ "***",
    p < 0.05 ~ "**",
    p < 0.10 ~ "*",
    TRUE ~ ""
  )
}

# ---------------------------------------------------------
# 7. CONTINUOUS VARIABLES
# ---------------------------------------------------------

cont_results <- map_dfr(continuous_vars, function(v) {
  tmp <- analysis_df %>%
    select(most_likely_class, all_of(v)) %>%
    filter(!is.na(.data[[v]]), !is.na(most_likely_class))
  
  m1 <- tmp %>%
    filter(most_likely_class == 1) %>%
    summarise(val = mean(.data[[v]], na.rm = TRUE)) %>%
    pull(val)
  
  m2 <- tmp %>%
    filter(most_likely_class == 2) %>%
    summarise(val = mean(.data[[v]], na.rm = TRUE)) %>%
    pull(val)
  
  pval <- tryCatch(
    t.test(tmp[[v]] ~ tmp$most_likely_class)$p.value,
    error = function(e) NA_real_
  )
  
  tibble(
    Variable = label_or_name(v),
    `Category / statistic` = "Mean",
    `Class 1` = round(m1, 3),
    `Class 2` = round(m2, 3),
    `p-value` = round(pval, 3),
    stars = sig_stars(pval),
    block_order = match(v, c(continuous_vars, categorical_vars)),
    row_in_block = 1
  )
})

# ---------------------------------------------------------
# 8. CATEGORICAL VARIABLES
# ---------------------------------------------------------

cat_results <- map_dfr(categorical_vars, function(v) {
  tmp <- analysis_df %>%
    select(most_likely_class, all_of(v)) %>%
    filter(!is.na(.data[[v]]), !is.na(most_likely_class))
  
  pval <- tryCatch(
    chisq.test(table(tmp[[v]], tmp$most_likely_class))$p.value,
    error = function(e) NA_real_
  )
  
  levs <- levels(droplevels(tmp[[v]]))
  if (is.null(levs)) levs <- sort(unique(as.character(tmp[[v]])))
  
  map_dfr(seq_along(levs), function(i) {
    lev <- levs[i]
    
    s1 <- tmp %>%
      filter(most_likely_class == 1) %>%
      summarise(val = mean(.data[[v]] == lev, na.rm = TRUE)) %>%
      pull(val)
    
    s2 <- tmp %>%
      filter(most_likely_class == 2) %>%
      summarise(val = mean(.data[[v]] == lev, na.rm = TRUE)) %>%
      pull(val)
    
    tibble(
      Variable = ifelse(i == 1, label_or_name(v), ""),
      `Category / statistic` = lev,
      `Class 1` = round(s1, 3),
      `Class 2` = round(s2, 3),
      `p-value` = ifelse(i == 1, round(pval, 3), NA_real_),
      stars = ifelse(i == 1, sig_stars(pval), ""),
      block_order = match(v, c(continuous_vars, categorical_vars)),
      row_in_block = i
    )
  })
})

# ---------------------------------------------------------
# 9. COMBINE AND FORMAT
# ---------------------------------------------------------

table5_df <- bind_rows(cont_results, cat_results) %>%
  arrange(block_order, row_in_block) %>%
  select(-block_order, -row_in_block) %>%
  mutate(
    `Class 1` = ifelse(is.na(`Class 1`), "", sprintf("%.3f", `Class 1`)),
    `Class 2` = ifelse(is.na(`Class 2`), "", sprintf("%.3f", `Class 2`)),
    # keep NA only for non-first rows of categorical blocks
    `p-value` = ifelse(is.na(`p-value`), "NA", sprintf("%.3f", `p-value`))
  )

# ---------------------------------------------------------
# 10. BUILD GT TABLE
# ---------------------------------------------------------

gt_table5 <- table5_df %>%
  gt() %>%
  tab_header(
    title = "Table 5. Class profile with significance levels"
  ) %>%
  cols_label(
    Variable = "Variable",
    `Category / statistic` = "Category / statistic",
    `Class 1` = "Class 1",
    `Class 2` = "Class 2",
    `p-value` = "p-value",
    stars = ""
  ) %>%
  cols_align("left", columns = c(Variable, `Category / statistic`)) %>%
  cols_align("center", columns = c(`Class 1`, `Class 2`, `p-value`, stars)) %>%
  tab_source_note(
    source_note = md("Continuous variables are reported as class means and tested using Welch two-sample t-tests. Categorical variables are reported as within-class category shares and tested using Pearson chi-square tests. For categorical variables, the p-value is shown on the first row of each variable block. Remaining rows in the same block are marked as NA by construction. Significance levels: * p < 0.10, ** p < 0.05, *** p < 0.01.")
  )

gt_table5

# ---------------------------------------------------------
# 11. EXPORT
# ---------------------------------------------------------

gtsave(gt_table5, "table5_class_profile_corrected.html")
gtsave(gt_table5, "table5_class_profile_corrected.docx")
gtsave(gt_table5, "table5_class_profile_corrected.tex")

# ---------------------------------------------------------
# 12. OPTIONAL CHECKS
# ---------------------------------------------------------

# Check harmonized categories
table(analysis_df$Soc_PeerAdo, useNA = "ifany")

# Check class means
analysis_df %>%
  group_by(most_likely_class) %>%
  summarise(
    risk_mean = mean(econ_risk_std, na.rm = TRUE),
    kbau_mean = mean(share_KBAU, na.rm = TRUE),
    land_mean = mean(Total_Land, na.rm = TRUE),
    aes_yes = mean(House_AesExp %in% c("Ja", "Yes"), na.rm = TRUE),
    age_mean = mean(House_Age, na.rm = TRUE),
    hireemp_mean = mean(Stru_HireEmp, na.rm = TRUE),
    n = n()
  )
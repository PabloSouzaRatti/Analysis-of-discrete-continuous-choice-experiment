# 03_intensity_analysis.R
# Conditional acreage/intensity analysis based on the saved LC model bundle.

##############################################################################
# SECOND-STAGE ANALYSIS FOR model_lcrisk_final
# ---------------------------------------------------------------------------
# Assumes these objects already exist in memory:
#   model_lcrisk_final
#   apollo_inputs
#   post
#
# First-stage LC model features:
#   - 2 classes
#   - class allocation depends on econ_risk_std
#   - utility in WTA space
#   - KBAU utility scaled as b_csa_k * asc_KBAU_k
#
# Second-stage available land:
#   Stru_OwnLand
##############################################################################

library(dplyr)
library(sandwich)
library(lmtest)

##############################################################################
# 0. START FROM THE SAME DATABASE USED IN THE LC MODEL
##############################################################################
args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value=TRUE)
if (!length(file_arg)) stop("Run with Rscript scripts/03_intensity_analysis.R")
project_root <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash="/")), ".."), winslash="/", mustWork=TRUE)
setwd(project_root)
bundle <- readRDS(file.path(project_root, "results", "analysis_bundle.rds"))
model_lcrisk_final <- bundle$model_lcrisk_final
db <- bundle$database
post <- bundle$post
post2 <- bundle$post2
b  <- model_lcrisk_final$estimate
out_dir <- file.path(project_root, "results", "intensity")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
setwd(out_dir)

getb <- function(par_name){
  if(!par_name %in% names(b)) stop("Parameter not found: ", par_name)
  b[[par_name]]
}

##############################################################################
# 0A. CHECK REQUIRED VARIABLES
##############################################################################
required_vars <- c(
  "resp_id", "task", "choice", "econ_risk_std", "Stru_OwnLand",
  "comp_1_s", "comp_2_s",
  "bonus_1_s", "bonus_2_s",
  "min_req_1_s", "min_req_2_s",
  "rec_cons_1", "rec_cons_2",
  "rec_farm_1", "rec_farm_2",
  "land_enrolled_1", "land_enrolled_2", "land_enrolled_3",
  "land_enrolled_4", "land_enrolled_5", "land_enrolled_6"
)

missing_vars <- setdiff(required_vars, names(db))
if(length(missing_vars) > 0){
  stop("Missing variables in db: ", paste(missing_vars, collapse = ", "))
}

required_pars <- c(
  "b_csa_1", "b_csa_2",
  "delta_2", "gamma_risk_2",
  "wta_TB_1", "wta_TB_2",
  "wta_MR_1", "wta_MR_2",
  "wta_recC_1", "wta_recC_2",
  "wta_recF_1", "wta_recF_2",
  "asc_KBAU_1", "asc_KBAU_2"
)

missing_pars <- setdiff(required_pars, names(b))
if(length(missing_pars) > 0){
  stop("Missing parameters in model_lcrisk_final$estimate: ",
       paste(missing_pars, collapse = ", "))
}

##############################################################################
# 1. CREATE TASK-SPECIFIC LAND ENROLLED VARIABLE
##############################################################################
db_dc <- db %>%
  mutate(
    land_enrolled = case_when(
      task == 1 ~ land_enrolled_1,
      task == 2 ~ land_enrolled_2,
      task == 3 ~ land_enrolled_3,
      task == 4 ~ land_enrolled_4,
      task == 5 ~ land_enrolled_5,
      task == 6 ~ land_enrolled_6,
      TRUE ~ NA_real_
    ),
    eligible_land = Stru_OwnLand
  )

##############################################################################
# 2. FIX IMPOSSIBLE LAND VALUES
##############################################################################
# Repeatedly divide by 10 until land_enrolled <= eligible_land
# This handles values like 1000, 10000, 300000, etc.
##############################################################################

fix_land <- function(x, max_land){
  if(is.na(x) || is.na(max_land) || max_land <= 0) return(x)
  
  while(x > max_land){
    x <- x / 10
  }
  
  x
}

db_dc <- db_dc %>%
  rowwise() %>%
  mutate(
    land_enrolled = fix_land(land_enrolled, eligible_land)
  ) %>%
  ungroup()

cat("\nRemaining impossible land values:\n")
print(sum(db_dc$land_enrolled > db_dc$eligible_land, na.rm = TRUE))

##############################################################################
# 3. CREATE SECOND-STAGE DEPENDENT VARIABLE
##############################################################################

db_dc <- db_dc %>%
  mutate(
    share_enrolled = ifelse(
      !is.na(eligible_land) & eligible_land > 0,
      land_enrolled / eligible_land,
      NA_real_
    )
  )

cat("\nSummary of share_enrolled:\n")
print(summary(db_dc$share_enrolled))

cat("\nRange of share_enrolled:\n")
print(range(db_dc$share_enrolled, na.rm = TRUE))



##############################################################################
# 4. COMPUTE FIRST-STAGE CLASS ALLOCATION PROBABILITIES
##############################################################################
# Based on actual estimated class allocation:
#   V_class1 = 0
#   V_class2 = delta_2 + gamma_risk_2 * econ_risk_std
##############################################################################
db_dc <- db_dc %>%
  mutate(
    V_class1 = 0,
    V_class2 = getb("delta_2") + getb("gamma_risk_2") * econ_risk_std
  )

# stable 2-class softmax
m_class <- pmax(db_dc$V_class1, db_dc$V_class2)
e1 <- exp(db_dc$V_class1 - m_class)
e2 <- exp(db_dc$V_class2 - m_class)

db_dc <- db_dc %>%
  mutate(
    pi_class1 = e1 / (e1 + e2),
    pi_class2 = e2 / (e1 + e2)
  )

##############################################################################
# 5. COMPUTE CLASS-SPECIFIC UTILITIES
##############################################################################
# IMPORTANT:
# KBAU is also multiplied by b_csa_k, matching the actual model.
##############################################################################

# Class 1
db_dc <- db_dc %>%
  mutate(
    V1_alt1 = getb("b_csa_1") * (
      comp_1_s +
        getb("wta_TB_1")   * bonus_1_s +
        getb("wta_MR_1")   * min_req_1_s +
        getb("wta_recC_1") * rec_cons_1 +
        getb("wta_recF_1") * rec_farm_1
    ),
    V1_alt2 = getb("b_csa_1") * (
      comp_2_s +
        getb("wta_TB_1")   * bonus_2_s +
        getb("wta_MR_1")   * min_req_2_s +
        getb("wta_recC_1") * rec_cons_2 +
        getb("wta_recF_1") * rec_farm_2
    ),
    V1_KBAU = getb("b_csa_1") * getb("asc_KBAU_1")
  )

# Class 2
db_dc <- db_dc %>%
  mutate(
    V2_alt1 = getb("b_csa_2") * (
      comp_1_s +
        getb("wta_TB_2")   * bonus_1_s +
        getb("wta_MR_2")   * min_req_1_s +
        getb("wta_recC_2") * rec_cons_1 +
        getb("wta_recF_2") * rec_farm_1
    ),
    V2_alt2 = getb("b_csa_2") * (
      comp_2_s +
        getb("wta_TB_2")   * bonus_2_s +
        getb("wta_MR_2")   * min_req_2_s +
        getb("wta_recC_2") * rec_cons_2 +
        getb("wta_recF_2") * rec_farm_2
    ),
    V2_KBAU = getb("b_csa_2") * getb("asc_KBAU_2")
  )

##############################################################################
# 6. COMPUTE CLASS-SPECIFIC AND UNCONDITIONAL CHOICE PROBABILITIES
##############################################################################
softmax3 <- function(v1, v2, v3){
  m <- pmax(v1, v2, v3)
  e1 <- exp(v1 - m)
  e2 <- exp(v2 - m)
  e3 <- exp(v3 - m)
  d  <- e1 + e2 + e3
  list(p1 = e1 / d, p2 = e2 / d, p3 = e3 / d)
}

sm1 <- softmax3(db_dc$V1_alt1, db_dc$V1_alt2, db_dc$V1_KBAU)
sm2 <- softmax3(db_dc$V2_alt1, db_dc$V2_alt2, db_dc$V2_KBAU)

db_dc <- db_dc %>%
  mutate(
    P1_alt1 = sm1$p1,
    P1_alt2 = sm1$p2,
    P1_KBAU = sm1$p3,
    P2_alt1 = sm2$p1,
    P2_alt2 = sm2$p2,
    P2_KBAU = sm2$p3,
    P_alt1  = pi_class1 * P1_alt1 + pi_class2 * P2_alt1,
    P_alt2  = pi_class1 * P1_alt2 + pi_class2 * P2_alt2,
    P_KBAU  = pi_class1 * P1_KBAU + pi_class2 * P2_KBAU,
    P_adopt = P_alt1 + P_alt2,
    P_sum   = P_alt1 + P_alt2 + P_KBAU
  )

cat("\nRange of P_sum (should be ~1):\n")
print(range(db_dc$P_sum, na.rm = TRUE))

##############################################################################
# 7. BUILD SECOND-STAGE SAMPLE: CONTRACT CHOICES ONLY
##############################################################################
d_contract <- db_dc %>%
  filter(choice %in% c(1, 2)) %>%
  mutate(
    comp_chosen   = ifelse(choice == 1, comp_1_s, comp_2_s),
    bonus_chosen  = ifelse(choice == 1, bonus_1_s, bonus_2_s),
    minreq_chosen = ifelse(choice == 1, min_req_1_s, min_req_2_s),
    recC_chosen   = ifelse(choice == 1, rec_cons_1, rec_cons_2),
    recF_chosen   = ifelse(choice == 1, rec_farm_1, rec_farm_2)
  )

cat("\nRows in contract-choice sample:\n")
print(nrow(d_contract))

##############################################################################
# 8. CREATE DUBIN–MCFADDEN CORRECTION TERMS
##############################################################################
eps <- 1e-8

d_contract <- d_contract %>%
  mutate(
    P_alt1 = pmin(pmax(P_alt1, eps), 1 - eps),
    P_alt2 = pmin(pmax(P_alt2, eps), 1 - eps),
    P_KBAU = pmin(pmax(P_KBAU, eps), 1 - eps),
    
    dm_chosen = ifelse(choice == 1, -log(P_alt1), -log(P_alt2)),
    
    dm_other1 = ifelse(
      choice == 1,
      (P_alt2 * log(P_alt2)) / (1 - P_alt2),
      (P_alt1 * log(P_alt1)) / (1 - P_alt1)
    ),
    
    dm_other2 = (P_KBAU * log(P_KBAU)) / (1 - P_KBAU)
  )

cat("\nSummary of DM terms:\n")
print(summary(d_contract[, c("dm_chosen", "dm_other1", "dm_other2")]))

##############################################################################
# 9. MERGE RESPONDENT-LEVEL POSTERIOR PROBABILITIES
##############################################################################
if(!exists("post")) stop("Object 'post' not found. Run LC post-estimation first.")

d_contract <- d_contract %>%
  left_join(
    post %>% select(resp_id, post_class_1, post_class_2, most_likely_class),
    by = "resp_id"
  )

cat("\nMissing post_class_2 values:\n")
print(sum(is.na(d_contract$post_class_2)))

##############################################################################
# 10. DIAGNOSTIC: ZERO LAND AFTER CONTRACT CHOICE
##############################################################################
zero_cases <- d_contract %>%
  filter(land_enrolled == 0)

cat("\nZero-land contract choices:\n")
print(nrow(zero_cases))
cat("\nShare of contract choices with zero land:\n")
print(nrow(zero_cases) / nrow(d_contract))

##############################################################################
# 11. MODEL 1: BENCHMARK INTENSIVE-MARGIN MODEL
##############################################################################
model1_ols <- lm(
  share_enrolled ~
    comp_chosen +
    bonus_chosen +
    minreq_chosen +
    recC_chosen +
    recF_chosen,
  data = d_contract
)

res_m1 <- coeftest(model1_ols, vcov = vcovCL(model1_ols, cluster = ~ resp_id))

model1_frac <- glm(
  share_enrolled ~
    comp_chosen +
    bonus_chosen +
    minreq_chosen +
    recC_chosen +
    recF_chosen,
  family = quasibinomial(link = "logit"),
  data = d_contract
)

res_m1f <- coeftest(model1_frac, vcov = vcovCL(model1_frac, cluster = ~ resp_id))

cat("\n================ MODEL 1 OLS =================\n")
print(res_m1)

cat("\n================ MODEL 1 FRACTIONAL =================\n")
print(res_m1f)

##############################################################################
# 12. MODEL 2: DISCRETE-CONTINUOUS MODEL WITH SELECTION CORRECTION
##############################################################################
model2_ols <- lm(
  share_enrolled ~
    comp_chosen +
    bonus_chosen +
    minreq_chosen +
    recC_chosen +
    recF_chosen +
    dm_chosen +
    dm_other1 +
    dm_other2,
  data = d_contract
)

res_m2 <- coeftest(model2_ols, vcov = vcovCL(model2_ols, cluster = ~ resp_id))

model2_frac <- glm(
  share_enrolled ~
    comp_chosen +
    bonus_chosen +
    minreq_chosen +
    recC_chosen +
    recF_chosen +
    dm_chosen +
    dm_other1 +
    dm_other2,
  family = quasibinomial(link = "logit"),
  data = d_contract
)

res_m2f <- coeftest(model2_frac, vcov = vcovCL(model2_frac, cluster = ~ resp_id))

cat("\n================ MODEL 2 OLS =================\n")
print(res_m2)

cat("\n================ MODEL 2 FRACTIONAL =================\n")
print(res_m2f)

##############################################################################
# 13. MODEL 3: HETEROGENEITY EXTENSION
##############################################################################
model3_ols <- lm(
  share_enrolled ~
    comp_chosen +
    bonus_chosen +
    minreq_chosen +
    recC_chosen +
    recF_chosen +
    post_class_2 +
    econ_risk_std +
    comp_chosen:post_class_2 +
    bonus_chosen:post_class_2 +
    minreq_chosen:post_class_2 +
    dm_chosen +
    dm_other1 +
    dm_other2,
  data = d_contract
)

res_m3 <- coeftest(model3_ols, vcov = vcovCL(model3_ols, cluster = ~ resp_id))

model3_frac <- glm(
  share_enrolled ~
    comp_chosen +
    bonus_chosen +
    minreq_chosen +
    recC_chosen +
    recF_chosen +
    post_class_2 +
    econ_risk_std +
    comp_chosen:post_class_2 +
    bonus_chosen:post_class_2 +
    minreq_chosen:post_class_2 +
    dm_chosen +
    dm_other1 +
    dm_other2,
  family = quasibinomial(link = "logit"),
  data = d_contract
)

res_m3f <- coeftest(model3_frac, vcov = vcovCL(model3_frac, cluster = ~ resp_id))

cat("\n================ MODEL 3 OLS =================\n")
print(res_m3)

cat("\n================ MODEL 3 FRACTIONAL =================\n")
print(res_m3f)

##############################################################################
# 14. ROBUSTNESS: EXCLUDE ZERO-LAND CONTRACT CHOICES
##############################################################################
d_contract_pos <- d_contract %>%
  filter(land_enrolled > 0)

cat("\nRows in positive-land sample:\n")
print(nrow(d_contract_pos))

##############################################################################
# 15. EXPECTED LAND UNDER OBSERVED DESIGN
##############################################################################
d_contract$pred_share_model2 <- predict(model2_ols, newdata = d_contract)
d_contract <- d_contract %>%
  mutate(pred_share_model2 = pmin(pmax(pred_share_model2, 0), 1))

db_dc <- db_dc %>%
  left_join(
    d_contract %>% select(resp_id, task, pred_share_model2),
    by = c("resp_id", "task")
  ) %>%
  mutate(
    expected_land_model2 = eligible_land * P_adopt * pred_share_model2
  )

cat("\nSummary of expected_land_model2:\n")
print(summary(db_dc$expected_land_model2))

##############################################################################
# 16. QUICK CHECKS
##############################################################################
cat("\nMean observed contract adoption rate:\n")
print(mean(db_dc$choice %in% c(1, 2), na.rm = TRUE))

cat("\nMean predicted adoption probability:\n")
print(mean(db_dc$P_adopt, na.rm = TRUE))

cat("\nMean observed share among adopters:\n")
print(mean(d_contract$share_enrolled, na.rm = TRUE))

cat("\nWorkflow complete.\n")


# Save the second-stage objects required by the standalone figure script.
saveRDS(list(d_contract=d_contract, model2_frac=model2_frac),
        file.path(out_dir, "intensity_plot_bundle.rds"))


################################################################################
# SECOND-STAGE TABLE + FIGURE 2
# Uses existing fitted models and d_contract
################################################################################

library(dplyr)
library(ggplot2)
library(sandwich)
library(lmtest)
library(modelsummary)

################################################################################
# PART 1. REGRESSION TABLE
# - contract attributes first
# - correction terms second
# - behavioural heterogeneity last
# - missing values shown as "-"
################################################################################

################################################################################
# UPDATED REGRESSION TABLE WITH SIGNIFICANCE STARS
# Keeps your existing working models unchanged
################################################################################

library(modelsummary)
library(sandwich)
library(dplyr)
library(broom)
library(tidyr)

# --------------------------------------------------
# Clustered vcov
# --------------------------------------------------
vcov_list <- list(
  vcovCL(model1_ols,  cluster = ~ resp_id),
  vcovCL(model1_frac, cluster = ~ resp_id),
  vcovCL(model2_ols,  cluster = ~ resp_id),
  vcovCL(model2_frac, cluster = ~ resp_id),
  vcovCL(model3_ols,  cluster = ~ resp_id),
  vcovCL(model3_frac, cluster = ~ resp_id)
)

# --------------------------------------------------
# Coefficient labels in desired order
# --------------------------------------------------
coef_map <- c(
  "comp_chosen"               = "Compensation",
  "bonus_chosen"              = "Threshold bonus",
  "minreq_chosen"             = "Minimum coordinated area requirement",
  "recC_chosen"               = "Recommendation from consultant",
  "recF_chosen"               = "Recommendation from farmers",
  
  "dm_chosen"                 = "Correction term: chosen alternative",
  "dm_other1"                 = "Correction term: other contract alternative",
  "dm_other2"                 = "Correction term: KBAU",
  
  "post_class_2"              = "Posterior probability of Class 2",
  "econ_risk_std"             = "Economic risk index",
  "comp_chosen:post_class_2"  = "Compensation × Class 2",
  "bonus_chosen:post_class_2" = "Threshold bonus × Class 2",
  "minreq_chosen:post_class_2"= "Minimum requirement × Class 2"
)

# --------------------------------------------------
# Goodness of fit rows
# --------------------------------------------------
gof_map <- tribble(
  ~raw,            ~clean,             ~fmt,
  "nobs",          "Observations",      0,
  "r.squared",     "R-squared",         3,
  "adj.r.squared", "Adj. R-squared",    3,
  "AIC",           "AIC",               1,
  "BIC",           "BIC",               1
)

################################################################################
# MAIN TABLE (with stars)
################################################################################

modelsummary(
  list(
    "Model 1 OLS"        = model1_ols,
    "Model 1 Fractional" = model1_frac,
    "Model 2 OLS"        = model2_ols,
    "Model 2 Fractional" = model2_frac,
    "Model 3 OLS"        = model3_ols,
    "Model 3 Fractional" = model3_frac
  ),
  
  vcov = vcov_list,
  
  coef_map = coef_map,
  gof_map  = gof_map,
  
  coef_omit = "Intercept",
  
  # stars added here
  stars = c(
    "*"   = .10,
    "**"  = .05,
    "***" = .01
  ),
  
  # estimate + stars on first line, SE in parentheses below
  estimate  = "{estimate}{stars}",
  statistic = "({std.error})",
  
  fmt = 3,
  
  output = "second_stage_regression_table_with_stars.docx",
  
  notes = c(
    "* p < 0.10",
    "** p < 0.05",
    "*** p < 0.01",
    "Cluster-robust standard errors at respondent level."
  )
)

cat("\nSaved table: second_stage_regression_table_with_stars.docx\n")

################################################################################
# OPTIONAL VERSION WITH DASHES FOR MISSING COEFFICIENTS + STARS
################################################################################

models <- list(
  "Model 1 OLS"        = model1_ols,
  "Model 1 Fractional" = model1_frac,
  "Model 2 OLS"        = model2_ols,
  "Model 2 Fractional" = model2_frac,
  "Model 3 OLS"        = model3_ols,
  "Model 3 Fractional" = model3_frac
)

star_fun <- function(p){
  case_when(
    p < 0.01 ~ "***",
    p < 0.05 ~ "**",
    p < 0.10 ~ "*",
    TRUE ~ ""
  )
}

tidy_all <- bind_rows(
  lapply(names(models), function(nm){
    
    td <- broom::tidy(
      models[[nm]],
      conf.int = FALSE
    )
    
    td$model <- nm
    td
  })
) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = recode(term, !!!coef_map),
    stars = star_fun(p.value),
    value = sprintf("%.3f%s (%.3f)", estimate, stars, std.error)
  ) %>%
  select(term, model, value) %>%
  pivot_wider(
    names_from = model,
    values_from = value,
    values_fill = "-"
  )

datasummary_df(
  tidy_all,
  output = "second_stage_regression_table_dashes_stars.docx",
  notes = c(
    "* p < 0.10",
    "** p < 0.05",
    "*** p < 0.01"
  )
)

cat("Saved table: second_stage_regression_table_dashes_stars.docx\n")
################################################################################

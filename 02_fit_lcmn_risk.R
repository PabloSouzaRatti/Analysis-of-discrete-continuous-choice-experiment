# 02_fit_lcmn_risk.R
# Estimate the two-class latent-class MNL model with perceived-risk covariate.
# Run scripts/01_prepare_estimation_data.R first.

suppressPackageStartupMessages({
  library(apollo)
  library(dplyr)
})
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]))), "_project_paths.R"))
setwd(project_root)
apollo_initialise()

prepared_file <- file.path(project_root, "data", "estimation_ready.csv")
if (!file.exists(prepared_file)) stop("Run 01_prepare_estimation_data.R first; estimation_ready.csv is missing.")
database <- read.csv(prepared_file, stringsAsFactors=FALSE, check.names=FALSE)
dir.create(file.path(project_root, "results", "apollo"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(project_root, "results"), recursive=TRUE, showWarnings=FALSE)

##############################################################################
# SETTINGS
##############################################################################
nClasses  <- 2
nStarts   <- 30
seed_base <- 123

##############################################################################
# apollo_control
##############################################################################
apollo_control <- list(
  modelName       = paste0("LC_with_riskCov_WTA_", nClasses, "class_b_csa_byClass"),
  modelDescr      = "Latent Class MNL with risk covariate in class allocation; WTA utility; alt1/alt2/KBAU; class-specific b_csa",
  indivID         = "resp_id",
  panelData       = TRUE,
  nCores          = 1,
  mixing          = FALSE,
  analyticGrad    = FALSE,
  outputDirectory = file.path(project_root, "results", "apollo"),
  noValidation    = FALSE
)

##############################################################################
# IMPORTANT: Make sure your database contains a respondent-level covariate:
#   database$risk_construct
# It can be continuous (recommended: mean-centered / standardized) or categorical.
##############################################################################

##############################################################################
# apollo_beta + apollo_fixed
##############################################################################
S <- nClasses

apollo_beta <- c()

# class-specific b_csa
for(k in 1:S){
  apollo_beta[paste0("b_csa_", k)] <- 0.2
}

# class allocation intercepts (delta_1 fixed to 0 for identification)
apollo_beta[paste0("delta_", 1:S)] <- 0

# class allocation covariate effects: risk_construct
# Identification: fix base class effect to 0 (same logic as delta_1).
apollo_beta[paste0("gamma_risk_", 1:S)] <- 0

# class-specific tastes + KBAU ASC
for(k in 1:S){
  apollo_beta[paste0("wta_TB_",   k)] <- 0
  apollo_beta[paste0("wta_MR_",   k)] <- 0
  apollo_beta[paste0("wta_recC_", k)] <- 0
  apollo_beta[paste0("wta_recF_", k)] <- 0
  apollo_beta[paste0("asc_KBAU_", k)] <- 0
}

# FIX base class allocation terms
apollo_fixed <- c("delta_1", "gamma_risk_1")

##############################################################################
# apollo_lcPars
##############################################################################
apollo_lcPars <- function(apollo_beta, apollo_inputs){
  
  lcPars <- list()
  S <- nClasses
  
  lcPars$b_csa    <- lapply(1:S, function(k) get(paste0("b_csa_", k)))
  lcPars$wta_TB   <- lapply(1:S, function(k) get(paste0("wta_TB_",   k)))
  lcPars$wta_MR   <- lapply(1:S, function(k) get(paste0("wta_MR_",   k)))
  lcPars$wta_recC <- lapply(1:S, function(k) get(paste0("wta_recC_", k)))
  lcPars$wta_recF <- lapply(1:S, function(k) get(paste0("wta_recF_", k)))
  lcPars$asc_KBAU <- lapply(1:S, function(k) get(paste0("asc_KBAU_", k)))
  
  # ---- class allocation utilities with covariate ----
  # V_k = delta_k + gamma_risk_k * risk_construct
  
  V_class <- list()
  for(k in 1:S){
    V_class[[paste0("class_", k)]] <-
      get(paste0("delta_", k)) + get(paste0("gamma_risk_", k)) * econ_risk_std
  }
  
  classAlloc_settings <- list(
    classes       = setNames(1:S, paste0("class_", 1:S)),
    utilities     = V_class,
    componentName = "classAlloc"
  )
  
  lcPars$pi_values <- apollo_classAlloc(classAlloc_settings)
  return(lcPars)
}

##############################################################################
# Validate inputs (AFTER apollo_beta exists and database is final)
##############################################################################
apollo_inputs <- apollo_validateInputs()

##############################################################################
# apollo_probabilities
##############################################################################
apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality="estimate"){
  
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  P <- list()
  lcPars <- apollo_lcPars(apollo_beta, apollo_inputs)
  S <- length(lcPars$pi_values)
  
  for(k in 1:S){
    
    b_csa_k <- lcPars$b_csa[[k]]
    
    V <- list(
      alt1 = b_csa_k * (
        comp_1_s +
          lcPars$wta_TB[[k]]   * bonus_1_s +
          lcPars$wta_MR[[k]]   * min_req_1_s +
          lcPars$wta_recC[[k]] * rec_cons_1 +
          lcPars$wta_recF[[k]] * rec_farm_1
      ),
      alt2 = b_csa_k * (
        comp_2_s +
          lcPars$wta_TB[[k]]   * bonus_2_s +
          lcPars$wta_MR[[k]]   * min_req_2_s +
          lcPars$wta_recC[[k]] * rec_cons_2 +
          lcPars$wta_recF[[k]] * rec_farm_2
      ),
      KBAU = b_csa_k * lcPars$asc_KBAU[[k]]
    )
    
    mnl_settings <- list(
      alternatives  = c(alt1=1, alt2=2, KBAU=3),
      avail         = list(alt1=1, alt2=1, KBAU=1),
      choiceVar     = choice,
      utilities     = V,
      componentName = paste0("Class_", k)
    )
    
    P[[paste0("Class_", k)]] <- apollo_mnl(mnl_settings, functionality)
    P[[paste0("Class_", k)]] <- apollo_panelProd(P[[paste0("Class_", k)]], apollo_inputs, functionality)
  }
  
  lc_settings <- list(
    inClassProb = P,
    classProb   = lcPars$pi_values
  )
  
  P[["model"]] <- apollo_lc(lc_settings, apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

##############################################################################
# Random starts (now ALSO randomise gamma_risk_k for k>=2)
##############################################################################
make_random_start_LC <- function(base_beta, nClasses, seed){
  set.seed(seed)
  b <- base_beta
  
  # identification
  b["delta_1"] <- 0
  b["gamma_risk_1"] <- 0
  
  if(nClasses >= 2){
    b["delta_2"] <- runif(1, -1.5, 1.5)
    b["gamma_risk_2"] <- runif(1, -2, 2)   # widen if needed
  }
  
  # randomise class-specific b_csa
  for(k in 1:nClasses){
    b[paste0("b_csa_", k)] <- runif(1, -2, 2)
  }
  
  for(k in 1:nClasses){
    b[paste0("asc_KBAU_", k)] <- runif(1, 2, 10)
    b[paste0("wta_TB_",   k)] <- runif(1, -0.5, 0.5)
    b[paste0("wta_MR_",   k)] <- runif(1, -0.2, 0.2)
    b[paste0("wta_recC_", k)] <- runif(1, -0.8, 0.8)
    b[paste0("wta_recF_", k)] <- runif(1, -0.8, 0.8)
  }
  
  # ensure separation (optional but often helps)
  if(nClasses >= 2){
    b["b_csa_2"]     <- b["b_csa_1"]     + runif(1, -1, 1)
    b["asc_KBAU_2"]  <- b["asc_KBAU_1"]  + runif(1, -2, 2)
    b["wta_TB_2"]    <- b["wta_TB_1"]    + runif(1, -0.4, 0.4)
    b["wta_MR_2"]    <- b["wta_MR_1"]    + runif(1, -0.15, 0.15)
    b["wta_recC_2"]  <- b["wta_recC_1"]  + runif(1, -0.6, 0.6)
    b["wta_recF_2"]  <- b["wta_recF_1"]  + runif(1, -0.6, 0.6)
    b["gamma_risk_2"]<- b["gamma_risk_2"]+ runif(1, -0.5, 0.5)
  }
  
  return(b)
}

##############################################################################
# Multi-start loop
##############################################################################
best_LL   <- -Inf
best_beta <- NULL
best_mod  <- NULL

base_beta <- apollo_beta

for(i in 1:nStarts){
  
  cat("\n============================\n")
  cat("Random start", i, "of", nStarts, "\n")
  cat("============================\n")
  
  apollo_beta_i <- make_random_start_LC(base_beta, nClasses, seed = seed_base + i)
  
  mod_i <- tryCatch(
    apollo_estimate(
      apollo_beta_i, apollo_fixed,
      apollo_probabilities, apollo_inputs,
      estimate_settings = list(estimationRoutine = "bgw")
    ),
    error = function(e){
      cat("Start", i, "FAILED:", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if(is.null(mod_i)) next
  
  LL_i <- mod_i$maximum
  cat("Start", i, "LL =", LL_i, "\n")
  
  if(is.finite(LL_i) && LL_i > best_LL){
    best_LL   <- LL_i
    best_beta <- mod_i$estimate
    best_mod  <- mod_i
  }
}

cat("\nBest LL from random starts:", best_LL, "\n")
if(is.null(best_beta)) stop("All random starts failed. Increase nStarts or widen start ranges.")

##############################################################################
# Final estimation from best start
##############################################################################
model_lcrisk_final <- apollo_estimate(
  best_beta, apollo_fixed,
  apollo_probabilities, apollo_inputs,
  estimate_settings = list(estimationRoutine = "bgw")
)

apollo_modelOutput(model_lcrisk_final)
#Saving output
apollo_saveOutput(model_lcrisk_final)
saveRDS(model_lcrisk_final, file.path(project_root, "results", "model_lcrisk_final.rds"))
##############################################################################
# POST-ESTIMATION TOOLKIT FOR LC MODEL WITH ASC SCALED IN WTA SPACE
# Assumes:
#   model_lc_final
#   apollo_inputs
#   apollo_probabilities() already defined for the ASC-scaled model
#
# Utility now is:
#   U_alt = b_csa_k * (comp + wta_TB*TB + wta_MR*MR + ...)
#   U_KBAU = b_csa_k * asc_KBAU_k
##############################################################################

library(dplyr)
library(ggplot2)

db <- apollo_inputs$database
id <- apollo_inputs$apollo_control$indivID  # usually "resp_id"

# ---------- helper ----------
getb <- function(b, nm){
  if(!nm %in% names(b)) stop("Parameter not found in beta vector: ", nm)
  b[[nm]]
}

# ---------- helper: coerce Apollo prediction output to matrix ----------
as_share_matrix <- function(Pmodel){
  if(is.matrix(Pmodel) || is.data.frame(Pmodel)) return(as.matrix(Pmodel))
  if(is.numeric(Pmodel) && !is.null(names(Pmodel))){
    M <- matrix(Pmodel, nrow = 1)
    colnames(M) <- names(Pmodel)
    return(M)
  }
  if(is.list(Pmodel)){
    nm <- names(Pmodel)
    if(all(c("alt1","alt2","KBAU") %in% nm)){
      M <- cbind(
        alt1 = as.numeric(Pmodel$alt1),
        alt2 = as.numeric(Pmodel$alt2),
        KBAU = as.numeric(Pmodel$KBAU)
      )
      return(M)
    }
  }
  stop("Don't know how to convert P$model to matrix. Run str(P$model).")
}

##############################################################################
# 1) Posterior class probabilities
##############################################################################
P_all <- apollo_probabilities(model_lcrisk_final$estimate, apollo_inputs, functionality = "prediction")
class_names <- grep("^Class_", names(P_all), value = TRUE)
S <- length(class_names)
stopifnot(S >= 2)

chosen_mat_row <- do.call(cbind, lapply(class_names, function(nm){
  as.numeric(P_all[[nm]]$chosen)
}))
colnames(chosen_mat_row) <- class_names
chosen_mat_row <- pmax(chosen_mat_row, 1e-300)

resp_ids <- unique(db[[id]])
N <- length(resp_ids)

logL_ind <- matrix(NA_real_, nrow = N, ncol = S, dimnames = list(NULL, class_names))
for(k in 1:S){
  tmp <- aggregate(log(chosen_mat_row[, k]) ~ db[[id]], FUN = sum)
  tmp <- tmp[match(resp_ids, tmp[[1]]), ]
  logL_ind[, k] <- tmp[[2]]
}

logL_centered <- logL_ind - apply(logL_ind, 1, max, na.rm = TRUE)
Lmat_ind <- exp(logL_centered)

# respondent-level risk covariate
risk_i <- tapply(db$econ_risk_std, db[[id]], function(x) x[1])
risk_i <- risk_i[as.character(resp_ids)]
risk_i <- as.numeric(risk_i)

# class allocation priors
V_prior <- matrix(0, nrow = N, ncol = S)
for(k in 1:S){
  delta_k <- if(k == 1) 0 else getb(model_lcrisk_final$estimate, paste0("delta_", k))
  gamma_k <- if(k == 1) 0 else getb(model_lcrisk_final$estimate, paste0("gamma_risk_", k))
  V_prior[, k] <- delta_k + gamma_k * risk_i
}
V_centered <- V_prior - apply(V_prior, 1, max)
pi_mat <- exp(V_centered)
pi_mat <- pi_mat / rowSums(pi_mat)
colnames(pi_mat) <- paste0("prior_class_", 1:S)

# posterior
post_num <- pi_mat * Lmat_ind
post_den <- rowSums(post_num)
post_ind <- post_num / post_den

post <- as.data.frame(post_ind)
names(post) <- paste0("post_class_", 1:S)
post$most_likely_class <- apply(post[, paste0("post_class_", 1:S), drop = FALSE], 1, which.max)
post[[id]] <- resp_ids

cat("\nPosterior mean class probs:\n")
print(colMeans(post[, paste0("post_class_", 1:S), drop = FALSE]))

cat("\nMost likely class counts:\n")
print(table(post$most_likely_class))

cat("\nMean posterior max-prob (certainty):\n")
print(mean(apply(post[, paste0("post_class_", 1:S), drop = FALSE], 1, max)))

##############################################################################
# 2) Class profiling: observed KBAU share
##############################################################################
kbau_rate <- aggregate(choice == 3 ~ resp_id, data = db, FUN = mean)
names(kbau_rate) <- c(id, "share_KBAU")

post2 <- merge(post, kbau_rate, by = id, all.x = TRUE)

cat("\nMean share_KBAU by most likely class:\n")
print(aggregate(share_KBAU ~ most_likely_class, data = post2, FUN = mean))

cat("\nShare always-KBAU by most likely class:\n")
print(aggregate(I(share_KBAU == 1) ~ most_likely_class, data = post2, FUN = mean))


# Save model objects used by the separate intensity, supplementary, and plotting scripts.
analysis_bundle <- list(
  model_lcrisk_final = model_lcrisk_final,
  database = db,
  post = post,
  post2 = post2
)
saveRDS(analysis_bundle, file.path(project_root, "results", "analysis_bundle.rds"))
cat("Saved analysis bundle to results/analysis_bundle.rds\n")

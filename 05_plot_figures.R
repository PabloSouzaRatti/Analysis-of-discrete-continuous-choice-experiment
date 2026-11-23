##############################################################################
# 05_plot_figures.R
# Standalone figures generated from results/analysis_bundle.rds.

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value=TRUE)
if (!length(file_arg)) stop("Run with Rscript scripts/05_plot_figures.R")
project_root <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash="/")), ".."), winslash="/", mustWork=TRUE)
setwd(project_root)
bundle <- readRDS(file.path(project_root, "results", "analysis_bundle.rds"))
model_lcrisk_final <- bundle$model_lcrisk_final
db <- bundle$database
output_folder <- file.path(project_root, "results", "figures")
dir.create(output_folder, recursive=TRUE, showWarnings=FALSE)

# ============================================================================
# COMMON SETUP FOR FIGURES 3, 4 AND Sx
# ============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)
library(scales)
library(grid)

# ---------------------------------------------------------------------------
# 1. Retrieve estimates, robust covariance matrix and estimation database
# ---------------------------------------------------------------------------

bhat <- model_lcrisk_final$estimate
Vrob <- model_lcrisk_final$robvarcov
db   <- bundle$database

if (is.null(Vrob)) {
  stop("model_lcrisk_final$robvarcov was not found.")
}

Vrob <- as.matrix(Vrob)

# ---------------------------------------------------------------------------
# 2. Helper: draw parameters from asymptotic MVN distribution
#    using ROBUST covariance matrix
# ---------------------------------------------------------------------------

make_robust_draws <- function(par_names,
                              n_draws = 1000,
                              seed = 12345) {
  
  if (!all(par_names %in% names(bhat))) {
    stop(
      "Missing parameter(s): ",
      paste(setdiff(par_names, names(bhat)), collapse = ", ")
    )
  }
  
  if (is.null(rownames(Vrob)) || is.null(colnames(Vrob))) {
    stop("Robust covariance matrix has no parameter names.")
  }
  
  if (!all(par_names %in% rownames(Vrob))) {
    stop(
      "Parameter(s) missing from robust covariance matrix: ",
      paste(setdiff(par_names, rownames(Vrob)), collapse = ", ")
    )
  }
  
  mu <- bhat[par_names]
  V  <- Vrob[par_names, par_names, drop = FALSE]
  
  if (any(!is.finite(V))) {
    stop("Selected robust covariance matrix contains non-finite values.")
  }
  
  # Force exact symmetry
  V <- (V + t(V)) / 2
  
  # Numerical safeguard against tiny negative eigenvalues
  eig <- eigen(V, symmetric = TRUE)
  
  if (min(eig$values) < -1e-6) {
    warning(
      "Selected covariance matrix has a negative eigenvalue: ",
      round(min(eig$values), 8)
    )
  }
  
  eig$values <- pmax(eig$values, 0)
  
  V_psd <- eig$vectors %*%
    diag(eig$values, nrow = length(eig$values)) %*%
    t(eig$vectors)
  
  dimnames(V_psd) <- list(par_names, par_names)
  
  set.seed(seed)
  
  draws <- MASS::mvrnorm(
    n     = n_draws,
    mu    = mu,
    Sigma = V_psd
  )
  
  draws <- as.matrix(draws)
  colnames(draws) <- par_names
  
  draws
}

# ---------------------------------------------------------------------------
# 3. Publication colours
# ---------------------------------------------------------------------------

class_colours <- c(
  "Status quo-oriented (Class 1)" = "#D55E00",
  "Adoption-oriented (Class 2)"   = "#0072B2"
)


################################################################################
### FIGURE 4. PERCEIVED ECONOMIC RISK AND LATENT-CLASS MEMBERSHIP
###
### 0-1 graphical normalization of the OBSERVED standardized risk index
###
### 0 = lowest observed perceived economic risk
### 1 = highest observed perceived economic risk
###
### The fitted model itself is NOT changed.
################################################################################


# ==============================================================================
# 0. PACKAGES
# ==============================================================================

library(dplyr)
library(ggplot2)
library(scales)
library(stringr)
library(MASS)
library(grid)


# ==============================================================================
# 1. MODEL OBJECTS
# ==============================================================================

bhat <- model_lcrisk_final$estimate
Vrob <- as.matrix(
  model_lcrisk_final$robvarcov
)

db <- bundle$database


# ==============================================================================
# 2. OUTPUT FOLDER
# ==============================================================================

output_folder <- file.path(project_root, "results", "figures")
dir.create(output_folder, recursive=TRUE, showWarnings=FALSE)


# ==============================================================================
# 3. RESPONDENT-LEVEL RISK
# ==============================================================================

risk_person <- db %>%
  group_by(resp_id) %>%
  summarise(
    econ_risk_std = first(econ_risk_std),
    .groups = "drop"
  )


# ==============================================================================
# 4. NORMALIZE OBSERVED RISK TO 0-1
#
# IMPORTANT:
# This is ONLY a graphical rescaling.
#
# 0 = lowest observed value
# 1 = highest observed value
# ==============================================================================

risk_min_std <- min(
  risk_person$econ_risk_std,
  na.rm = TRUE
)

risk_max_std <- max(
  risk_person$econ_risk_std,
  na.rm = TRUE
)


risk_person <- risk_person %>%
  mutate(
    risk_01 =
      (
        econ_risk_std -
          risk_min_std
      ) /
      (
        risk_max_std -
          risk_min_std
      )
  )


cat(
  "\nObserved standardized risk range:\n",
  "Minimum = ",
  risk_min_std,
  "\nMaximum = ",
  risk_max_std,
  "\n",
  sep = ""
)


# ==============================================================================
# 5. PREDICTION GRID
#
# Work visually from 0 to 1...
# ==============================================================================

risk_grid_01 <- seq(
  0,
  1,
  length.out = 300
)


# ...but translate it back to standardized risk for the fitted model.
risk_grid_std <-
  risk_min_std +
  risk_grid_01 *
  (
    risk_max_std -
      risk_min_std
  )


# ==============================================================================
# 6. ROBUST PARAMETER DRAWS
# ==============================================================================

par_names <- c(
  "delta_2",
  "gamma_risk_2"
)


mu <- bhat[
  par_names
]

V <- Vrob[
  par_names,
  par_names,
  drop = FALSE
]


# Ensure perfect numerical symmetry
V <- (
  V +
    t(V)
) / 2


set.seed(
  20260812
)


draws_fig4 <- MASS::mvrnorm(
  n     = 2000,
  mu    = mu,
  Sigma = V
)


draws_fig4 <- as.matrix(
  draws_fig4
)

colnames(
  draws_fig4
) <- par_names


# ==============================================================================
# 7. POINT PREDICTION
#
# Class-allocation equation:
#
# V1 = 0
#
# V2 = delta_2 + gamma_risk_2 * risk
#
# P(Class 2) = exp(V2) / [1 + exp(V2)]
# ==============================================================================

eta_point <-
  bhat["delta_2"] +
  bhat["gamma_risk_2"] *
  risk_grid_std


p_point <- plogis(
  eta_point
)


# ==============================================================================
# 8. 95% CONFIDENCE INTERVAL
# ==============================================================================

eta_draws <-
  
  outer(
    draws_fig4[, "delta_2"],
    rep(
      1,
      length(risk_grid_std)
    )
  ) +
  
  outer(
    draws_fig4[, "gamma_risk_2"],
    risk_grid_std
  )


p_draws <- plogis(
  eta_draws
)


fig4_data <- data.frame(
  
  risk_01 =
    risk_grid_01,
  
  estimate =
    p_point,
  
  lower =
    apply(
      p_draws,
      2,
      quantile,
      probs = 0.025,
      na.rm = TRUE
    ),
  
  upper =
    apply(
      p_draws,
      2,
      quantile,
      probs = 0.975,
      na.rm = TRUE
    )
)


# ==============================================================================
# 9. CAPTION
# ==============================================================================

figure_4_caption <- paste(
  
  "Figure 4. Perceived economic risk and probability of membership in the",
  "adoption-oriented latent class.",
  
  "The line shows the predicted probability of membership in Class 2 as",
  "perceived economic risk increases.",
  
  "For graphical interpretation, the standardized economic-risk index is",
  "rescaled over its observed sample range from 0 (lowest observed perceived",
  "economic risk) to 1 (highest observed perceived economic risk).",
  
  "This rescaling affects only the graphical presentation; predictions are",
  "calculated using the standardized risk covariate employed in model estimation.",
  
  "The shaded area represents the 95% confidence interval derived from 2,000",
  "parameter draws using the robust covariance matrix.",
  
  "Rug marks indicate the observed distribution of respondent-level risk values."
)


figure_4_caption_wrapped <- stringr::str_wrap(
  figure_4_caption,
  width = 118
)


# ==============================================================================
# 10. CLEAN JOURNAL FIGURE
# ==============================================================================

p_fig4 <- ggplot(
  fig4_data,
  aes(
    x = risk_01,
    y = estimate
  )
) +
  
  # 95% confidence interval
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    fill  = "#56B4E9",
    alpha = 0.12
  ) +
  
  # Predicted probability
  geom_line(
    colour    = "#0072B2",
    linewidth = 1.15
  ) +
  
  # Observed respondent distribution
  geom_rug(
    data = risk_person,
    
    aes(
      x = risk_01
    ),
    
    inherit.aes = FALSE,
    
    sides = "b",
    
    colour = "grey25",
    
    alpha = 0.40,
    
    linewidth = 0.40,
    
    length = grid::unit(
      3,
      "mm"
    )
  ) +
  
  # X axis
  scale_x_continuous(
    
    limits = c(
      0,
      1
    ),
    
    breaks = seq(
      0,
      1,
      by = 0.2
    ),
    
    labels = scales::label_number(
      accuracy = 0.1
    ),
    
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  
  # Y axis
  scale_y_continuous(
    
    limits = c(
      0,
      1
    ),
    
    breaks = seq(
      0,
      1,
      by = 0.2
    ),
    
    labels = scales::label_percent(
      accuracy = 1
    ),
    
    expand = expansion(
      mult = c(
        0,
        0.01
      )
    )
  ) +
  
  labs(
    
    x =
      "Perceived risk (0 = lowest observed, 1 = highest observed)",
    
    y =
      "Probability of Class 2 membership"
  ) +
  
  theme_bw(
    base_size = 11.5,
    base_family = "sans"
  ) +
  
  theme(
    
    # AXIS TITLES
    axis.title.x = element_text(
      size = 11.5,
      colour = "black",
      margin = margin(t = 10)
    ),
    
    axis.title.y = element_text(
      size = 11.5,
      colour = "black",
      margin = margin(r = 12)
    ),
    
    # AXIS NUMBERS
    axis.text.x = element_text(
      size = 10.5,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 10.5,
      colour = "black",
      margin = margin(r = 4)
    ),
    
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.35
    ),
    
    # GRID
    panel.grid.major.y = element_line(
      colour = "grey88",
      linewidth = 0.30
    ),
    
    panel.grid.major.x = element_blank(),
    
    panel.grid.minor = element_blank(),
    
    # BORDER
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.40
    ),
    
    # IMPORTANT: extra room on LEFT
    plot.margin = margin(
      t = 12,
      r = 14,
      b = 14,
      l = 26
    )
  )


print(
  p_fig4
)


# ==============================================================================
# 11. DISPLAY VERSION WITH TITLE + CAPTION
# ==============================================================================

p_fig4_display <- p_fig4 +
  
  labs(
    
    title =
      "Perceived economic risk and latent-class membership",
    
    subtitle =
      "Probability of belonging to the adoption-oriented class",
    
    caption =
      figure_4_caption_wrapped
  ) +
  
  theme(
    
    plot.title.position =
      "plot",
    
    plot.title = element_text(
      size   = 15,
      face   = "bold",
      hjust  = 0,
      colour = "black",
      margin = margin(
        b = 5
      )
    ),
    
    plot.subtitle = element_text(
      size   = 11,
      hjust  = 0,
      colour = "black",
      margin = margin(
        b = 10
      )
    ),
    
    plot.caption.position =
      "plot",
    
    plot.caption = element_text(
      size       = 9.2,
      hjust      = 0,
      colour     = "black",
      lineheight = 1.08,
      margin = margin(
        t = 13
      )
    ),
    
    plot.margin = margin(
      t = 16,
      r = 14,
      b = 20,
      l = 14
    )
  )


print(
  p_fig4_display
)


# ==============================================================================
# 12. SAVE CAPTION
# ==============================================================================

writeLines(
  figure_4_caption,
  
  file.path(
    output_folder,
    "Figure_4_caption.txt"
  )
)


# ==============================================================================
# 13. EXPORT FUNCTION
# ==============================================================================

save_publication_figure <- function(
    plot_object,
    filename,
    width_mm = 150,
    height_mm = 110
) {
  
  pdf_device <- if (
    capabilities("cairo")
  ) {
    grDevices::cairo_pdf
  } else {
    "pdf"
  }
  
  
  # PDF
  ggsave(
    filename = file.path(
      output_folder,
      paste0(
        filename,
        ".pdf"
      )
    ),
    
    plot = plot_object,
    device = pdf_device,
    
    width  = width_mm,
    height = height_mm,
    units  = "mm",
    
    bg = "white",
    limitsize = FALSE
  )
  
  
  # TIFF
  ggsave(
    filename = file.path(
      output_folder,
      paste0(
        filename,
        "_600dpi.tiff"
      )
    ),
    
    plot = plot_object,
    device = "tiff",
    
    width  = width_mm,
    height = height_mm,
    units  = "mm",
    
    dpi = 600,
    compression = "lzw",
    
    bg = "white",
    limitsize = FALSE
  )
  
  
  # PNG
  ggsave(
    filename = file.path(
      output_folder,
      paste0(
        filename,
        "_600dpi.png"
      )
    ),
    
    plot = plot_object,
    device = "png",
    
    width  = width_mm,
    height = height_mm,
    units  = "mm",
    
    dpi = 600,
    
    bg = "white",
    limitsize = FALSE
  )
}


# ==============================================================================
# 14. EXPORT CLEAN JOURNAL VERSION
# ==============================================================================

save_publication_figure(
  
  plot_object =
    p_fig4,
  
  filename =
    "Figure_4_risk_class_membership_FINAL",
  
  width_mm =
    135,
  
  height_mm =
    105
)


# ==============================================================================
# 15. EXPORT DISPLAY VERSION
# ==============================================================================

save_publication_figure(
  
  plot_object =
    p_fig4_display,
  
  filename =
    "Figure_4_risk_class_membership_DISPLAY",
  
  width_mm =
    150,
  
  height_mm =
    155
)


# ==============================================================================
# 16. REPORT LOCATION
# ==============================================================================

cat(
  "\n============================================================\n",
  "Figure 4 completed.\n",
  
  "\nSaved in:\n",
  normalizePath(
    output_folder
  ),
  
  "\n\nSUBMISSION FILES:\n",
  "Figure_4_risk_class_membership_FINAL.pdf\n",
  "Figure_4_risk_class_membership_FINAL_600dpi.tiff\n",
  "Figure_4_risk_class_membership_FINAL_600dpi.png\n",
  
  "\nDISPLAY / REVIEW FILES:\n",
  "Figure_4_risk_class_membership_DISPLAY.pdf\n",
  "Figure_4_risk_class_membership_DISPLAY_600dpi.tiff\n",
  "Figure_4_risk_class_membership_DISPLAY_600dpi.png\n",
  
  "\nCAPTION:\n",
  "Figure_4_caption.txt\n",
  
  "============================================================\n",
  sep = ""
)

################################################################################
### FIGURE 3. PREDICTED CONTRACT ADOPTION BY LATENT CLASS
###
### Run AFTER:
###   model_lcrisk_final
###   apollo_inputs
###
### FINAL VERSION
################################################################################


# ==============================================================================
# 0. PACKAGES
# ==============================================================================

library(dplyr)
library(ggplot2)
library(MASS)
library(scales)
library(stringr)
library(grid)


# ==============================================================================
# 1. MODEL OBJECTS
# ==============================================================================

bhat <- model_lcrisk_final$estimate
Vrob <- as.matrix(model_lcrisk_final$robvarcov)
db   <- bundle$database


# ==============================================================================
# 2. OUTPUT FOLDER
# ==============================================================================

output_folder <- file.path(project_root, "results", "figures")
dir.create(output_folder, recursive=TRUE, showWarnings=FALSE)


# ==============================================================================
# 3. COLLAPSE DATABASE TO UNIQUE DESIGN PROFILES
# ==============================================================================

design_db <- db %>%
  count(
    comp_1_s,
    comp_2_s,
    bonus_1_s,
    bonus_2_s,
    min_req_1_s,
    min_req_2_s,
    rec_cons_1,
    rec_cons_2,
    rec_farm_1,
    rec_farm_2,
    name = "weight"
  )


# ==============================================================================
# 4. CLASS-SPECIFIC ADOPTION PROBABILITY
# ==============================================================================

predict_adoption_class <- function(beta, class_k, dat) {
  
  b_csa    <- beta[paste0("b_csa_", class_k)]
  wta_TB   <- beta[paste0("wta_TB_", class_k)]
  wta_MR   <- beta[paste0("wta_MR_", class_k)]
  wta_recC <- beta[paste0("wta_recC_", class_k)]
  wta_recF <- beta[paste0("wta_recF_", class_k)]
  asc_KBAU <- beta[paste0("asc_KBAU_", class_k)]
  
  
  V1 <- b_csa * (
    dat$comp_1_s +
      wta_TB   * dat$bonus_1_s +
      wta_MR   * dat$min_req_1_s +
      wta_recC * dat$rec_cons_1 +
      wta_recF * dat$rec_farm_1
  )
  
  
  V2 <- b_csa * (
    dat$comp_2_s +
      wta_TB   * dat$bonus_2_s +
      wta_MR   * dat$min_req_2_s +
      wta_recC * dat$rec_cons_2 +
      wta_recF * dat$rec_farm_2
  )
  
  
  VKBAU <- rep(
    b_csa * asc_KBAU,
    nrow(dat)
  )
  
  
  # Stable softmax
  vmax <- pmax(
    V1,
    V2,
    VKBAU
  )
  
  
  exp1 <- exp(V1 - vmax)
  exp2 <- exp(V2 - vmax)
  exp3 <- exp(VKBAU - vmax)
  
  
  p_adopt <- (
    exp1 + exp2
  ) / (
    exp1 + exp2 + exp3
  )
  
  
  weighted.mean(
    p_adopt,
    w = dat$weight
  )
}


# ==============================================================================
# 5. PARAMETERS USED FOR CONFIDENCE INTERVALS
# ==============================================================================

pars_fig3 <- c(
  "b_csa_1",
  "wta_TB_1",
  "wta_MR_1",
  "wta_recC_1",
  "wta_recF_1",
  "asc_KBAU_1",
  
  "b_csa_2",
  "wta_TB_2",
  "wta_MR_2",
  "wta_recC_2",
  "wta_recF_2",
  "asc_KBAU_2"
)


# Check parameters
stopifnot(
  all(pars_fig3 %in% names(bhat)),
  all(pars_fig3 %in% rownames(Vrob)),
  all(pars_fig3 %in% colnames(Vrob))
)


mu_fig3 <- bhat[
  pars_fig3
]


V_fig3 <- Vrob[
  pars_fig3,
  pars_fig3,
  drop = FALSE
]


# Numerical symmetry
V_fig3 <- (
  V_fig3 +
    t(V_fig3)
) / 2


# Small safeguard if numerical rounding creates tiny negative eigenvalues
eig <- eigen(
  V_fig3,
  symmetric = TRUE
)

eig$values <- pmax(
  eig$values,
  0
)


V_fig3_psd <- eig$vectors %*%
  diag(
    eig$values,
    nrow = length(eig$values)
  ) %*%
  t(eig$vectors)


dimnames(V_fig3_psd) <- list(
  pars_fig3,
  pars_fig3
)


# ==============================================================================
# 6. PARAMETER DRAWS
# ==============================================================================

set.seed(20260811)


draws_fig3 <- MASS::mvrnorm(
  n = 1000,
  mu = mu_fig3,
  Sigma = V_fig3_psd
)


draws_fig3 <- as.matrix(
  draws_fig3
)


colnames(draws_fig3) <- pars_fig3


# ==============================================================================
# 7. FUNCTION TO SIMULATE ONE ATTRIBUTE
#
# The focal attribute is SET to each value.
# Other attributes retain their observed choice-design values.
# ==============================================================================

make_attribute_curve <- function(
    grid_values,
    attribute,
    scenario_label
) {
  
  results <- vector(
    "list",
    length = length(grid_values) * 2
  )
  
  counter <- 1
  
  
  for (class_k in 1:2) {
    
    for (x_value in grid_values) {
      
      d <- design_db
      
      
      # ------------------------------------------------------------------------
      # Set focal attribute
      # ------------------------------------------------------------------------
      
      if (attribute == "compensation") {
        
        d$comp_1_s <- x_value / 10000
        d$comp_2_s <- x_value / 10000
        
      }
      
      
      if (attribute == "minimum") {
        
        d$min_req_1_s <- x_value / 10
        d$min_req_2_s <- x_value / 10
        
      }
      
      
      if (attribute == "bonus") {
        
        d$bonus_1_s <- x_value / 1000
        d$bonus_2_s <- x_value / 1000
        
      }
      
      
      # ------------------------------------------------------------------------
      # Point estimate
      # ------------------------------------------------------------------------
      
      point_estimate <- predict_adoption_class(
        beta = bhat,
        class_k = class_k,
        dat = d
      )
      
      
      # ------------------------------------------------------------------------
      # Predictions from parameter draws
      # ------------------------------------------------------------------------
      
      draw_predictions <- numeric(
        nrow(draws_fig3)
      )
      
      
      for (j in seq_len(nrow(draws_fig3))) {
        
        beta_j <- bhat
        
        beta_j[
          pars_fig3
        ] <- draws_fig3[
          j,
          pars_fig3
        ]
        
        
        draw_predictions[j] <- predict_adoption_class(
          beta = beta_j,
          class_k = class_k,
          dat = d
        )
      }
      
      
      # ------------------------------------------------------------------------
      # Store
      # ------------------------------------------------------------------------
      
      results[[counter]] <- data.frame(
        scenario = scenario_label,
        x = x_value,
        class_k = class_k,
        estimate = point_estimate,
        lower = unname(
          quantile(
            draw_predictions,
            probs = 0.025,
            na.rm = TRUE
          )
        ),
        upper = unname(
          quantile(
            draw_predictions,
            probs = 0.975,
            na.rm = TRUE
          )
        )
      )
      
      
      counter <- counter + 1
    }
  }
  
  
  bind_rows(
    results
  )
}


# ==============================================================================
# 8. ATTRIBUTE RANGES
# ==============================================================================

# Compensation:
# 43,000 to 120,000 DKK/ha
comp_grid <- seq(
  from = 43000,
  to = 120000,
  length.out = 80
)


# Minimum coordinated-area requirement:
# 0 to 100 ha
# 0 = no minimum requirement
mr_grid <- seq(
  from = 0,
  to = 100,
  length.out = 80
)


# Threshold bonus:
# 0 to 25,000 DKK/ha
tb_grid <- seq(
  from = 0,
  to = 25000,
  length.out = 80
)


# ==============================================================================
# 9. RUN SIMULATIONS
# ==============================================================================

curve_comp <- make_attribute_curve(
  grid_values = comp_grid,
  attribute = "compensation",
  scenario_label = "Compensation"
)


curve_mr <- make_attribute_curve(
  grid_values = mr_grid,
  attribute = "minimum",
  scenario_label = "Minimum requirement"
)


curve_tb <- make_attribute_curve(
  grid_values = tb_grid,
  attribute = "bonus",
  scenario_label = "Threshold bonus"
)


# ==============================================================================
# 10. COMBINE RESULTS
#
# IMPORTANT:
#
# Monetary values are divided by 1,000 ONLY FOR DISPLAY:
#
# Compensation:
#   40 = 40,000 DKK/ha
#
# Threshold bonus:
#   5 = 5,000 DKK/ha
#
# Minimum requirement remains in hectares.
# ==============================================================================

fig3_data <- bind_rows(
  curve_comp,
  curve_mr,
  curve_tb
) %>%
  
  mutate(
    
    class = case_when(
      class_k == 1 ~
        "Status quo-oriented (Class 1)",
      
      class_k == 2 ~
        "Adoption-oriented (Class 2)"
    ),
    
    
    class = factor(
      class,
      levels = c(
        "Status quo-oriented (Class 1)",
        "Adoption-oriented (Class 2)"
      )
    ),
    
    
    scenario = factor(
      scenario,
      levels = c(
        "Compensation",
        "Minimum requirement",
        "Threshold bonus"
      )
    ),
    
    
    x_plot = case_when(
      
      scenario == "Compensation" ~
        x / 1000,
      
      scenario == "Threshold bonus" ~
        x / 1000,
      
      TRUE ~
        x
    )
  )


# ==============================================================================
# 11. COLOURS AND LINE TYPES
# ==============================================================================

class_colours <- c(
  "Status quo-oriented (Class 1)" = "#D55E00",
  "Adoption-oriented (Class 2)"   = "#0072B2"
)


class_linetypes <- c(
  "Status quo-oriented (Class 1)" = "longdash",
  "Adoption-oriented (Class 2)"   = "solid"
)


# ==============================================================================
# 12. FIGURE CAPTION
# ==============================================================================

figure_3_caption <- paste(
  "Figure 3. Predicted adoption probability across contract attribute levels, by latent class.",
  "Predicted probabilities of selecting either permanent set-aside contract are shown across",
  "variation in compensation, minimum coordinated area requirement, and threshold bonus.",
  "Monetary values on the x-axes are expressed in thousand DKK/ha; for example, 40 corresponds",
  "to 40,000 DKK/ha. The minimum requirement is expressed in hectares, with 0 representing",
  "the absence of a minimum coordinated-area requirement.",
  "The focal attribute is varied while the remaining contract attributes retain their observed",
  "choice-design values, and predictions are averaged across choice tasks.",
  "Lines show point estimates and shaded areas show 95% confidence intervals derived from",
  "1,000 parameter draws using the robust covariance matrix.",
  "The confidence intervals distinguish more precisely estimated behavioural responses from",
  "more uncertain simulated trajectories; in particular, the wider intervals for the",
  "status quo-oriented class reflect the limited precision of its estimated responses to",
  "the minimum requirement and threshold bonus."
)


figure_3_caption_wrapped <- stringr::str_wrap(
  figure_3_caption,
  width = 118
)


# ==============================================================================
# 13. FINAL CLEAN JOURNAL FIGURE
# ==============================================================================

p_fig3 <- ggplot(
  fig3_data,
  aes(
    x = x_plot,
    y = estimate,
    colour = class,
    fill = class,
    linetype = class
  )
) +
  
  # 95% confidence intervals
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.07,
    colour = NA,
    linetype = 0
  ) +
  
  # Point estimates
  geom_line(
    linewidth = 1.15
  ) +
  
  # Panels
  facet_wrap(
    ~ scenario,
    nrow = 1,
    scales = "free_x"
  ) +
  
  # Colours
  scale_colour_manual(
    values = class_colours,
    drop = FALSE
  ) +
  
  scale_fill_manual(
    values = class_colours,
    drop = FALSE
  ) +
  
  scale_linetype_manual(
    values = class_linetypes,
    drop = FALSE
  ) +
  
  # Y axis
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(
      0,
      1,
      by = 0.20
    ),
    labels = scales::label_percent(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(
        0.005,
        0.015
      )
    )
  ) +
  
  # X axes
  #
  # Because facets use free_x, ggplot generates appropriate
  # breaks separately for each panel:
  #
  # Compensation:
  #   approximately 40, 60, 80, 100, 120
  #
  # Minimum:
  #   approximately 0, 20, 40, 60, 80, 100
  #
  # Threshold bonus:
  #   approximately 0, 5, 10, 15, 20, 25
  #
  scale_x_continuous(
    breaks = scales::breaks_pretty(
      n = 5
    ),
    labels = scales::label_number(
      accuracy = 1,
      big.mark = ""
    )
  ) +
  
  labs(
    x = NULL,
    y = "Predicted adoption probability",
    colour = NULL,
    fill = NULL,
    linetype = NULL
  ) +
  
  theme_bw(
    base_size = 11.5,
    base_family = "sans"
  ) +
  
  theme(
    
    # Axes
    axis.title.y = element_text(
      size = 11.5,
      colour = "black",
      margin = margin(
        r = 8
      )
    ),
    
    axis.text.x = element_text(
      size = 10.5,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 10.5,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.35
    ),
    
    # Grid
    panel.grid.major.x = element_line(
      colour = "grey90",
      linewidth = 0.30
    ),
    
    panel.grid.major.y = element_line(
      colour = "grey88",
      linewidth = 0.30
    ),
    
    panel.grid.minor = element_blank(),
    
    # Borders
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.40
    ),
    
    panel.spacing.x = grid::unit(
      5,
      "mm"
    ),
    
    # Facet labels
    strip.background = element_rect(
      fill = "grey94",
      colour = "black",
      linewidth = 0.40
    ),
    
    strip.text = element_text(
      size = 11,
      face = "bold",
      colour = "black",
      margin = margin(
        5,
        5,
        5,
        5
      )
    ),
    
    # Legend
    legend.position = "bottom",
    legend.direction = "horizontal",
    
    legend.text = element_text(
      size = 10
    ),
    
    legend.key.width = grid::unit(
      15,
      "mm"
    ),
    
    legend.margin = margin(
      t = 5,
      b = 2
    ),
    
    # Margins
    plot.margin = margin(
      t = 10,
      r = 12,
      b = 10,
      l = 10
    )
  )


# ==============================================================================
# 14. DISPLAY CLEAN FIGURE
# ==============================================================================

print(
  p_fig3
)


# ==============================================================================
# 15. DISPLAY VERSION WITH TITLE + FULL CAPTION
# ==============================================================================

p_fig3_display <- p_fig3 +
  
  labs(
    
    title =
      "Predicted adoption responses to permanent set-aside contract attributes",
    
    subtitle =
      "Predicted probabilities by latent class",
    
    caption =
      figure_3_caption_wrapped
  ) +
  
  theme(
    
    plot.title.position =
      "plot",
    
    plot.title = element_text(
      size = 15,
      face = "bold",
      hjust = 0,
      colour = "black",
      margin = margin(
        b = 5
      )
    ),
    
    plot.subtitle = element_text(
      size = 11,
      hjust = 0,
      colour = "black",
      margin = margin(
        b = 10
      )
    ),
    
    plot.caption.position =
      "plot",
    
    plot.caption = element_text(
      size = 9.2,
      hjust = 0,
      colour = "black",
      lineheight = 1.08,
      margin = margin(
        t = 13
      )
    ),
    
    plot.margin = margin(
      t = 16,
      r = 14,
      b = 20,
      l = 14
    )
  )


print(
  p_fig3_display
)


# ==============================================================================
# 16. SAVE CAPTION SEPARATELY
# ==============================================================================

writeLines(
  figure_3_caption,
  file.path(
    output_folder,
    "Figure_3_caption.txt"
  )
)


# ==============================================================================
# 17. EXPORT FUNCTION
# ==============================================================================

save_publication_figure <- function(
    plot_object,
    filename,
    width_mm,
    height_mm
) {
  
  pdf_device <- if (
    capabilities("cairo")
  ) {
    grDevices::cairo_pdf
  } else {
    "pdf"
  }
  
  
  # PDF
  ggsave(
    filename = file.path(
      output_folder,
      paste0(
        filename,
        ".pdf"
      )
    ),
    plot = plot_object,
    device = pdf_device,
    width = width_mm,
    height = height_mm,
    units = "mm",
    bg = "white",
    limitsize = FALSE
  )
  
  
  # TIFF 600 dpi
  ggsave(
    filename = file.path(
      output_folder,
      paste0(
        filename,
        "_600dpi.tiff"
      )
    ),
    plot = plot_object,
    device = "tiff",
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 600,
    compression = "lzw",
    bg = "white",
    limitsize = FALSE
  )
  
  
  # PNG 600 dpi
  ggsave(
    filename = file.path(
      output_folder,
      paste0(
        filename,
        "_600dpi.png"
      )
    ),
    plot = plot_object,
    device = "png",
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
}


# ==============================================================================
# 18. EXPORT CLEAN SUBMISSION VERSION
# ==============================================================================

save_publication_figure(
  plot_object = p_fig3,
  filename = "Figure_3_predicted_adoption_FINAL",
  width_mm = 190,
  height_mm = 105
)


# ==============================================================================
# 19. EXPORT DISPLAY / REVIEW VERSION
# ==============================================================================

save_publication_figure(
  plot_object = p_fig3_display,
  filename = "Figure_3_predicted_adoption_DISPLAY",
  width_mm = 190,
  height_mm = 155
)


# ==============================================================================
# 20. REPORT SAVED FILES
# ==============================================================================

cat(
  "\n============================================================\n",
  "FIGURE 3 COMPLETED\n",
  "============================================================\n",
  "\nSaved in:\n",
  normalizePath(output_folder),
  "\n\nSubmission files:\n",
  "Figure_3_predicted_adoption_FINAL.pdf\n",
  "Figure_3_predicted_adoption_FINAL_600dpi.tiff\n",
  "Figure_3_predicted_adoption_FINAL_600dpi.png\n",
  "\nDisplay/review files:\n",
  "Figure_3_predicted_adoption_DISPLAY.pdf\n",
  "Figure_3_predicted_adoption_DISPLAY_600dpi.tiff\n",
  "Figure_3_predicted_adoption_DISPLAY_600dpi.png\n",
  "\nCaption:\n",
  "Figure_3_caption.txt\n",
  "============================================================\n",
  sep = ""
)

################################################################################
# FIGURE 2. PREDICTED ENROLLED ACREAGE FROM THE INTENSITY MODEL
################################################################################
intensity_bundle_file <- file.path(project_root, "results", "intensity", "intensity_plot_bundle.rds")
if (file.exists(intensity_bundle_file)) {
  intensity_bundle <- readRDS(intensity_bundle_file)
  d_contract <- intensity_bundle$d_contract
  model2_frac <- intensity_bundle$model2_frac
  library(sandwich)
  library(lmtest)
################################################################################
# PART 2. FIGURE 2
# Predicted enrolled acreage conditional on adoption across five attributes
# with shaded 95% confidence intervals
################################################################################

make_second_stage_plot_data_real <- function(attr_name, attr_label){
  
  grid <- data.frame(
    x_scaled = sort(unique(d_contract[[attr_name]]))
  )
  
  grid <- grid %>%
    mutate(
      comp_chosen   = mean(d_contract$comp_chosen, na.rm = TRUE),
      bonus_chosen  = mean(d_contract$bonus_chosen, na.rm = TRUE),
      minreq_chosen = mean(d_contract$minreq_chosen, na.rm = TRUE),
      recC_chosen   = mean(d_contract$recC_chosen, na.rm = TRUE),
      recF_chosen   = mean(d_contract$recF_chosen, na.rm = TRUE),
      dm_chosen     = mean(d_contract$dm_chosen, na.rm = TRUE),
      dm_other1     = mean(d_contract$dm_other1, na.rm = TRUE),
      dm_other2     = mean(d_contract$dm_other2, na.rm = TRUE)
    )
  
  # overwrite focal attribute
  grid[[attr_name]] <- grid$x_scaled
  
  # convert scaled variables to real units for plotting
  grid <- grid %>%
    mutate(
      x_real = case_when(
        attr_name == "comp_chosen"   ~ x_scaled * 10000, # DKK/ha
        attr_name == "bonus_chosen"  ~ x_scaled * 1000,  # DKK/ha
        attr_name == "minreq_chosen" ~ x_scaled * 10,    # ha
        TRUE                         ~ x_scaled
      ),
      attribute = attr_label
    )
  
  # model matrix from your existing model2_frac
  X <- model.matrix(delete.response(terms(model2_frac)), data = grid)
  
  # coefficients + clustered vcov
  beta <- coef(model2_frac)
  V    <- vcovCL(model2_frac, cluster = ~ resp_id)
  
  eta    <- as.numeric(X %*% beta)
  se_eta <- sqrt(diag(X %*% V %*% t(X)))
  
  invlogit <- function(z) exp(z) / (1 + exp(z))
  
  grid <- grid %>%
    mutate(
      pred  = invlogit(eta),
      lower = invlogit(eta - 1.96 * se_eta),
      upper = invlogit(eta + 1.96 * se_eta)
    )
  
  grid
}

# Build data for all five attributes
plot_comp <- make_second_stage_plot_data_real(
  attr_name  = "comp_chosen",
  attr_label = "Compensation (DKK/ha, lump-sum)"
)

plot_bonus <- make_second_stage_plot_data_real(
  attr_name  = "bonus_chosen",
  attr_label = "Threshold bonus (DKK/ha)"
)

plot_minreq <- make_second_stage_plot_data_real(
  attr_name  = "minreq_chosen",
  attr_label = "Minimum coordinated area requirement (ha)"
)

plot_recC <- make_second_stage_plot_data_real(
  attr_name  = "recC_chosen",
  attr_label = "Recommendation from consultant"
)

plot_recF <- make_second_stage_plot_data_real(
  attr_name  = "recF_chosen",
  attr_label = "Recommendation from farmers"
)

plot_all <- bind_rows(
  plot_comp,
  plot_bonus,
  plot_minreq,
  plot_recC,
  plot_recF
)

# Order panels
plot_all$attribute <- factor(
  plot_all$attribute,
  levels = c(
    "Compensation (DKK/ha, lump-sum)",
    "Threshold bonus (DKK/ha)",
    "Minimum coordinated area requirement (ha)",
    "Recommendation from consultant",
    "Recommendation from farmers"
  )
)

# Figure styled like your attached example
figure2 <- ggplot(plot_all, aes(x = x_real, y = pred)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.20, fill = "grey70") +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.6) +
  facet_wrap(~ attribute, scales = "free_x", ncol = 2) +
  labs(
    title = "Figure 2. Predicted enrolled acreage conditional on contract adoption across five contract attributes, with shaded 95% confidence intervals",
    x = NULL,
    y = "Predicted enrolled acreage share | adopt"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(figure2)

ggsave(
  filename = "Figure2_predicted_enrolled_acreage_all_attributes.png",
  plot = figure2,
  width = 12,
  height = 9,
  dpi = 300
)

cat("Saved figure to: Figure2_predicted_enrolled_acreage_all_attributes.png\n")
################################################################################
# END
################################################################################
} else {
  message("Skipping Figure 2: run scripts/03_intensity_analysis.R first to create the intensity model bundle.")
}

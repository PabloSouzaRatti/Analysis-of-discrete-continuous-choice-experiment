# 01_prepare_estimation_data.R
# Validate the public respondent-task data and write the fixed input for estimation.
# No respondent rows are added, removed, or reweighted in this step.

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value=TRUE)
if (!length(file_arg)) stop("Run with Rscript scripts/01_prepare_estimation_data.R")
project_root <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash="/")), ".."), winslash="/", mustWork=TRUE)
input_file <- file.path(project_root, "data", "Dataset_public.csv")
output_file <- file.path(project_root, "data", "estimation_ready.csv")
if (!file.exists(input_file)) stop("Dataset not found: ", input_file)
d <- read.csv(input_file, stringsAsFactors=FALSE, check.names=FALSE)
required <- c("resp_id", "task", "choice", "econ_risk_std", "version",
              "comp_1_s", "comp_2_s", "bonus_1_s", "bonus_2_s",
              "min_req_1_s", "min_req_2_s", "rec_cons_1", "rec_cons_2",
              "rec_farm_1", "rec_farm_2")
missing <- setdiff(required, names(d))
if (length(missing)) stop("Required columns missing: ", paste(missing, collapse=", "))
if (anyNA(d[, required])) stop("Missing values found in required model columns.")
if (!all(d$choice %in% 1:3)) stop("choice must use 1=contract 1, 2=contract 2, 3=KBAU.")
if (!all(d$task %in% 1:6)) stop("task must be numbered 1 to 6.")
counts <- table(d$resp_id)
if (any(counts != 6L)) stop("Each respondent must have exactly six task rows.")
if (anyDuplicated(d[c("resp_id", "task")])) stop("Duplicate respondent-task records found.")
if (length(unique(d$resp_id)) != 447L) warning("Expected 447 respondents based on current release; found ", length(unique(d$resp_id)), ".")
if (any(!is.finite(d$econ_risk_std))) stop("econ_risk_std contains non-finite values.")
# Keep the released columns and coding intact; the source CSV is already long and model-ready.
write.csv(d, output_file, row.names=FALSE, na="")
cat("Validated ", nrow(d), " task rows for ", length(unique(d$resp_id)), " respondents.\n", sep="")
cat("Estimation input written to data/estimation_ready.csv\n")

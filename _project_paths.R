# Shared path helper. Scripts are designed to run from any working directory via Rscript.
args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value=TRUE)
if (length(file_arg)) {
  this_script <- normalizePath(sub("^--file=", "", file_arg[1]), winslash="/", mustWork=TRUE)
  project_root <- normalizePath(file.path(dirname(this_script), ".."), winslash="/", mustWork=TRUE)
} else if (exists("ofile", envir=parent.frame(), inherits=FALSE)) {
  project_root <- normalizePath(file.path(dirname(get("ofile", envir=parent.frame())), ".."), winslash="/", mustWork=TRUE)
} else {
  stop("Could not determine project root. Run this file through Rscript from the scripts folder.")
}

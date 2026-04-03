# Auto-install packages if not present
if (!require("remotes", quietly = TRUE)) install.packages("remotes")

pkgs <- c("httpgd", "languageserver")
missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]

if (length(missing) > 0) {
  cat("Installing missing packages:", paste(missing, collapse = ", "), "\n")
  remotes::install_github("nx10/httpgd", upgrade = "never", quiet = TRUE)
  install.packages("languageserver", upgrade = "never", quiet = TRUE)
}

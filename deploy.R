# =============================================================================
# deploy.R — ADC Peptide Mapper v0.8
# =============================================================================
# Deploy to shinyapps.io (or Posit Connect Cloud).
#
# FIRST-TIME SETUP (run once):
#   install.packages("rsconnect")
#   rsconnect::setAccountInfo(
#     name   = "YOUR_SHINYAPPS_IO_ACCOUNT_NAME",
#     token  = "YOUR_TOKEN",           # from shinyapps.io → Account → Tokens
#     secret = "YOUR_SECRET"
#   )
#
# Then run this script to deploy:
#   source("deploy.R")
#
# NOTES:
#   • The background .rds files in data/ are large (~50–150 MB each).
#     shinyapps.io free tier has a 1 GB bundle limit.
#     Upload only the species you need by setting INCLUDE_SPECIES below.
#   • The ANTHROPIC_API_KEY should be set as a shinyapps.io environment variable
#     (Settings → Environment Variables) — never commit the key to git.
#   • First-run background build: run build_background_db.R locally first,
#     then deploy. The .rds files are included in the bundle automatically.
# =============================================================================

library(rsconnect)

# ── Config ────────────────────────────────────────────────────────────────────

APP_NAME   <- "ADC_Peptide_Mapper"    # appears in the shinyapps.io URL
APP_DIR    <- "."                     # run from project root (contains app.R)
ACCOUNT    <- NULL                    # set to your shinyapps.io account name,
                                      # or NULL to use the default account

# Which background .rds files to include in the bundle
# (delete species you don't need to save bundle space)
INCLUDE_SPECIES <- c("human", "monkey", "rat")

# ── File list ─────────────────────────────────────────────────────────────────
# rsconnect will auto-detect all files in the app directory, but we can
# explicitly list to exclude large test files or dev artefacts.

EXCLUDE_PATTERNS <- c(
  "^tests/",              # unit tests — not needed at runtime
  "^.claude/",
  "^deploy\\.R$",         # this file itself
  "\\.log$",
  "^build_background_db\\.R$"  # build script — not needed at runtime
)

# ── Validate .rds files ───────────────────────────────────────────────────────
for (sp in INCLUDE_SPECIES) {
  rds <- file.path("data", paste0("bg_", sp, ".rds"))
  if (!file.exists(rds)) {
    cat(sprintf("WARNING: %s not found — run build_background_db.R first.\n", rds))
  } else {
    sz_mb <- round(file.size(rds) / 1024 / 1024, 1)
    cat(sprintf("  Found: %-35s  %.1f MB\n", rds, sz_mb))
  }
}

# ── Dependency check ──────────────────────────────────────────────────────────
required_pkgs <- c("shiny", "bs4Dash", "DT", "data.table", "openxlsx",
                   "httr2", "stringr", "dplyr", "shinyjs")
missing <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0L) {
  cat("Installing missing packages:", paste(missing, collapse=", "), "\n")
  install.packages(missing, repos = "https://cloud.r-project.org")
}

# ── Deploy ────────────────────────────────────────────────────────────────────
cat("\n=============================================================\n")
cat(" Deploying ADC Peptide Mapper v0.8 to shinyapps.io\n")
cat("=============================================================\n")

rsconnect::deployApp(
  appDir     = APP_DIR,
  appName    = APP_NAME,
  account    = ACCOUNT,
  launch.browser = FALSE,
  forceUpdate    = TRUE,
  logLevel       = "verbose"
)

cat("\n=============================================================\n")
cat(sprintf(" Done! App URL: https://%s.shinyapps.io/%s\n",
            rsconnect::accounts()$name[1], APP_NAME))
cat(" Next steps:\n")
cat("   1. Go to shinyapps.io → your app → Settings → Environment Variables\n")
cat("   2. Add ANTHROPIC_API_KEY = sk-ant-...\n")
cat("      (This powers the AI Assistant tab)\n")
cat("   3. Add a link from nishi76.github.io to the app URL\n")
cat("=============================================================\n")

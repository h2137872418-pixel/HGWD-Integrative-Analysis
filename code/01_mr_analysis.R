# ====================================================================
# Mendelian Randomization Analysis Script (Core, without dIVW/ConMix)
# Includes: Basic MR methods, MR-PRESSO, sensitivity analysis, visualization
# ====================================================================

# -------------------- 0. Set file paths (modify to your actual paths) --------------------
exposure_file <- "your_exposure_data.csv"        # Exposure GWAS data
outcome_file  <- "your_outcome_data.csv"         # Outcome GWAS data
bfile_path    <- "path/to/1000G_EUR"             # LD reference panel (without extension)
plink_path    <- "path/to/plink.exe"             # PLINK executable file
output_dir    <- "mr_results"                    # Output directory

# -------------------- 1. Install / load required packages --------------------
packages_needed <- c("TwoSampleMR", "data.table", "dplyr", "MRPRESSO", "ggplot2")
for (pkg in packages_needed) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg == "TwoSampleMR") {
      install.packages("TwoSampleMR", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
    } else if (pkg == "MRPRESSO") {
      install.packages("MRPRESSO")
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# -------------------- 2. Create output directory --------------------
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -------------------- 3. Define helper function: calculate F-statistic --------------------
calc_f_stat <- function(beta, se, eaf, n) {
  R2 <- (beta^2 * 2 * eaf * (1 - eaf)) /
    (beta^2 * 2 * eaf * (1 - eaf) + se^2 * n * 2 * eaf * (1 - eaf))
  F_stat <- (R2 * (n - 2)) / (1 - R2)
  return(F_stat)
}

# -------------------- 4. Read exposure data and manually map column names --------------------
cat("===== Step 1: Read exposure data =====\n")
exp_raw <- fread(exposure_file, data.table = FALSE)
cat("Preview first 6 rows of exposure data, modify the mapping below accordingly:\n")
print(head(exp_raw))

# ========== User-defined mapping area (exposure data) ==========
# Please replace the column names below with your actual column names
exp_std <- exp_raw %>%
  transmute(
    SNP = SNP,                              # Replace with SNP column name
    chr.exposure = chr.exposure,            # Replace with chromosome column name
    pos.exposure = pos.exposure,            # Replace with position column name
    other_allele.exposure = other_allele.exposure,  # Non-effect allele
    effect_allele.exposure = effect_allele.exposure, # Effect allele
    beta.exposure = beta.exposure,          # Effect size beta
    se.exposure = se.exposure,              # Standard error
    pval.exposure = pval.exposure,          # P-value
    eaf.exposure = eaf.exposure,            # Effect allele frequency
    samplesize.exposure = 10000,            # Total sample size (please fill in the actual value)
    exposure = "exposure_name",             # Custom exposure name
    mr_keep.exposure = TRUE,
    id.exposure = "exposure_id"            # Custom exposure ID
  )
# ===================================================

cat(sprintf("Exposure data mapping completed, total %d rows\n", nrow(exp_std)))

# -------------------- 5. Filter exposure data by P-value --------------------
cat("\n===== Step 2: P-value filtering (p < 5e-8) =====\n")
pval_threshold <- 5e-8
exp_filtered <- exp_std %>% filter(!is.na(pval.exposure) & pval.exposure < pval_threshold)
cat(sprintf("After filtering, %d SNPs retained\n", nrow(exp_filtered)))
if (nrow(exp_filtered) == 0) stop("No significant SNPs, please lower the P-value threshold or check data.")

# -------------------- 6. Clumping (remove LD) --------------------
cat("\n===== Step 3: Clumping =====\n")
exp_clumped <- clump_data(
  exp_filtered,
  clump_r2 = 0.001,   # Adjust as needed
  clump_kb = 10000,
  bfile = bfile_path,
  plink_bin = plink_path
)
cat(sprintf("After clumping, %d SNPs retained\n", nrow(exp_clumped)))
if (nrow(exp_clumped) == 0) stop("No SNPs after clumping, please adjust parameters.")

# -------------------- 7. Calculate F-statistic and filter weak instruments --------------------
cat("\n===== Step 4: Calculate F-statistic =====\n")
exp_clumped <- exp_clumped %>%
  mutate(F_stat = calc_f_stat(beta.exposure, se.exposure, eaf.exposure, samplesize.exposure))
exp_final <- exp_clumped %>% filter(F_stat > 10)
cat(sprintf("After F > 10 filter, %d strong instruments retained\n", nrow(exp_final)))
if (nrow(exp_final) == 0) stop("No strong instruments.")

# -------------------- 8. Read outcome data and manually map column names --------------------
cat("\n===== Step 5: Read outcome data =====\n")
outcome_raw <- fread(outcome_file, data.table = FALSE)
cat("Preview first 6 rows of outcome data, modify the mapping below accordingly:\n")
print(head(outcome_raw))

# ========== User-defined mapping area (outcome data) ==========
outcome_std <- outcome_raw %>%
  transmute(
    SNP = SNP,                              # Replace with SNP column name
    chr.outcome = chr,                      # Replace with chromosome column name
    pos.outcome = pos,                      # Replace with position column name
    other_allele.outcome = other_allele,    # Non-effect allele
    effect_allele.outcome = effect_allele,  # Effect allele
    beta.outcome = beta,                    # Effect size
    se.outcome = se,                        # Standard error
    pval.outcome = pval,                    # P-value
    eaf.outcome = eaf,                      # Effect allele frequency
    samplesize.outcome = 20000,             # Total sample size (please fill in the actual value)
    outcome = "outcome_name",               # Custom outcome name
    id.outcome = "outcome_id"               # Custom outcome ID
  )
# ===================================================

cat(sprintf("Outcome data mapping completed, total %d rows\n", nrow(outcome_std)))

# -------------------- 9. Data harmonisation --------------------
cat("\n===== Step 6: Harmonise exposure and outcome data =====\n")
harmonized <- harmonise_data(
  exposure_dat = exp_final,
  outcome_dat = outcome_std,
  action = 2
)
cat(sprintf("After harmonisation, %d SNPs retained\n", nrow(harmonized)))
if (nrow(harmonized) == 0) stop("No SNPs after harmonisation.")

# -------------------- 10. Steiger filtering (causal direction test) --------------------
cat("\n===== Step 7: Steiger direction filtering =====\n")
steiger_res <- steiger_filtering(harmonized)
harmonized_steiger <- harmonized[steiger_res$steiger_test, ]
if (nrow(harmonized_steiger) > 0) {
  harmonized <- harmonized_steiger
  cat(sprintf("After Steiger filtering, %d SNPs retained\n", nrow(harmonized)))
} else {
  warning("No SNPs after Steiger filtering, all harmonised data will be used.")
}

# -------------------- 11. Basic MR analysis (IVW, MR-Egger, weighted median, etc.) --------------------
cat("\n===== Step 8: Basic MR analysis =====\n")
mr_results <- mr(harmonized)
mr_or <- generate_odds_ratios(mr_results)
cat("Basic MR methods completed, obtained ", nrow(mr_results), " rows of results\n")

# -------------------- 12. MR-PRESSO test --------------------
cat("\n===== Step 9: MR-PRESSO pleiotropy test =====\n")
presso_data <- harmonized %>%
  select(beta.exposure, beta.outcome, se.exposure, se.outcome) %>%
  na.omit()
if (nrow(presso_data) > 0) {
  presso_res <- tryCatch(
    mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
              SdOutcome = "se.outcome", SdExposure = "se.exposure",
              data = presso_data, NbDistribution = 1000, SignifThreshold = 0.05),
    error = function(e) NULL
  )
  if (!is.null(presso_res)) {
    cat("MR-PRESSO global test P-value =", presso_res$`MR-PRESSO results`$`Global Test`$Pvalue, "\n")
    saveRDS(presso_res, file.path(output_dir, "mr_presso_results.rds"))
  } else {
    cat("MR-PRESSO failed (possibly too few SNPs).\n")
  }
}

# -------------------- 13. Sensitivity analysis (heterogeneity, pleiotropy) --------------------
cat("\n===== Step 10: Sensitivity analysis =====\n")
mr_heter <- mr_heterogeneity(harmonized)
mr_pleio <- mr_pleiotropy_test(harmonized)

# -------------------- 14. Save result tables --------------------
cat("\n===== Step 11: Save result tables =====\n")
fwrite(harmonized, file.path(output_dir, "harmonized_data.csv"))
fwrite(mr_results, file.path(output_dir, "mr_results.csv"))
fwrite(mr_or, file.path(output_dir, "mr_odds_ratios.csv"))
fwrite(mr_heter, file.path(output_dir, "heterogeneity.csv"))
fwrite(mr_pleio, file.path(output_dir, "pleiotropy.csv"))

# Extract IVW result summary
ivw_res <- mr_or %>% filter(method == "Inverse variance weighted")
if (nrow(ivw_res) > 0) {
  ivw_het <- mr_heter %>% filter(method == "IVW")
  summary_tab <- data.frame(
    Exposure = unique(harmonized$exposure)[1],
    Outcome = unique(harmonized$outcome)[1],
    NSnps = ivw_res$nsnp,
    Beta = ivw_res$b,
    SE = ivw_res$se,
    Pval = ivw_res$pval,
    OR = ivw_res$or,
    OR_L95 = ivw_res$or_lci95,
    OR_U95 = ivw_res$or_uci95,
    Heterogeneity_Q = ifelse(nrow(ivw_het) > 0, ivw_het$Q, NA),
    Heterogeneity_P = ifelse(nrow(ivw_het) > 0, ivw_het$Q_pval, NA),
    Egger_Intercept = ifelse(nrow(mr_pleio) > 0, mr_pleio$egger_intercept[1], NA),
    Egger_Intercept_SE = ifelse(nrow(mr_pleio) > 0, mr_pleio$se[1], NA),
    Egger_Intercept_P = ifelse(nrow(mr_pleio) > 0, mr_pleio$pval[1], NA),
    stringsAsFactors = FALSE
  )
  fwrite(summary_tab, file.path(output_dir, "summary_results.csv"))
}

# -------------------- 15. Generate visualization plots --------------------
cat("\n===== Step 12: Generate plots =====\n")
save_plot_multi <- function(plot_obj, filename, width = 8, height = 6) {
  if (is.null(plot_obj) || inherits(plot_obj, "try-error")) return(invisible(NULL))
  ggsave(file.path(output_dir, paste0(filename, ".png")), plot = plot_obj, width = width, height = height, dpi = 300)
  ggsave(file.path(output_dir, paste0(filename, ".pdf")), plot = plot_obj, width = width, height = height)
  cat("  Saved:", filename, "\n")
}

# Scatter plot
p_scatter <- mr_scatter_plot(mr_results, harmonized)
if (length(p_scatter) > 0) save_plot_multi(p_scatter[[1]], "scatter_plot")

# Forest plot
res_single <- mr_singlesnp(harmonized)
p_forest <- mr_forest_plot(res_single)
if (length(p_forest) > 0) save_plot_multi(p_forest[[1]], "forest_plot", width = 10)

# Funnel plot
p_funnel <- mr_funnel_plot(res_single)
if (length(p_funnel) > 0) save_plot_multi(p_funnel[[1]], "funnel_plot")

# Leave-one-out plot
res_loo <- mr_leaveoneout(harmonized)
p_loo <- mr_leaveoneout_plot(res_loo)
if (length(p_loo) > 0) save_plot_multi(p_loo[[1]], "leaveoneout_plot", width = 10)

cat("\n===== All analyses completed! Results saved in:", output_dir, " =====\n")

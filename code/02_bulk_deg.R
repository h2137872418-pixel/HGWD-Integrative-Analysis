# ====================================================================
# Differential Expression Analysis Pipeline (limma-based)
# Core workflow preserved, enhanced for generalizability and code sharing
# ====================================================================

# Load required packages
library(limma)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(cowplot)
library(ggsci)
library(ggpubr)
library(scales)
library(statmod)
library(dplyr)

# 1. Data loading and preprocessing -----------------------------------------
load_and_preprocess_data <- function(expr_file   = "GSE95849_expression_matrix.csv",
                                     sample_file = "pheno_clean.csv",
                                     output_dir  = "Preprocessed_Data",
                                     min_count   = 10) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  if (!file.exists(expr_file))
    stop("Expression file not found: ", expr_file)
  if (!file.exists(sample_file))
    stop("Sample file not found: ", sample_file)
  
  cat("Reading expression data from:", expr_file, "\n")
  
  # Automatically try multiple delimiters
  read_success <- FALSE
  for (sep in c(",", "\t", ";")) {
    tryCatch({
      expr_data <- read.csv(expr_file,
                            sep = sep,
                            check.names = FALSE,
                            stringsAsFactors = FALSE,
                            na.strings = c("NA", ""))
      if (nrow(expr_data) > 1 && ncol(expr_data) > 1) {
        cat("Successfully read data with separator:",
            dQuote(sep), "\n")
        read_success <- TRUE
        break
      }
    }, error = function(e) {
      cat("Attempt with separator", dQuote(sep), "failed:",
          e$message, "\n")
    })
  }
  if (!read_success) {
    stop("Failed to read expression file with common separators")
  }
  
  # Clean column names (remove leading 'X' added by R)
  colnames(expr_data) <- gsub("^X", "", colnames(expr_data))
  
  # Dimension check
  if (nrow(expr_data) < 10 || ncol(expr_data) < 5) {
    warning(sprintf("Suspicious data dimensions: %d rows x %d columns",
                    nrow(expr_data), ncol(expr_data)))
  }
  
  # Explicitly specify the gene column as "gene_name"
  gene_cols <- "gene_name"
  
  if ("gene_name" %in% colnames(expr_data)) {
    cat("Using gene column: gene_name\n")
    gene_symbols  <- expr_data[["gene_name"]]
    # Remove non-sample columns: ID, Phalanx_id, gene_name
    sample_columns <- setdiff(colnames(expr_data),
                              c("ID", "Phalanx_id", "gene_name"))
  } else {
    stop("Required column 'gene_name' not found in expression data")
  }
  
  # Ensure sample columns exist
  if (length(sample_columns) < 3) {
    stop("Insufficient sample columns detected: ",
         paste(sample_columns, collapse = ", "))
  }
  
  expr_matrix <- as.matrix(expr_data[, sample_columns])
  rownames(expr_matrix) <- make.unique(gene_symbols)
  
  # log2 transformation check
  data_transformed <- FALSE
  if (nrow(expr_matrix) > 0 && ncol(expr_matrix) > 0) {
    q90 <- quantile(expr_matrix, 0.9, na.rm = TRUE)
    if (q90 > 500) {
      cat("Data appears to be in raw intensities (90th percentile =",
          round(q90), "). Applying log2 transformation.\n")
      expr_matrix <- log2(expr_matrix + 1)
      data_transformed <- TRUE
    } else {
      cat("Data likely log-transformed (90th percentile =",
          round(q90), "). No transformation applied.\n")
    }
  } else {
    stop("Expression matrix is empty. Columns used: ",
         paste(sample_columns, collapse = ", "))
  }
  
  # Read sample information
  cat("Reading sample information from:", sample_file, "\n")
  sample_info <- read.csv(sample_file, stringsAsFactors = FALSE)
  colnames(sample_info) <- tolower(colnames(sample_info))
  
  required_cols <- c("title", "group")
  if (!all(required_cols %in% colnames(sample_info))) {
    stop("Sample information missing required columns. Found: ",
         paste(colnames(sample_info), collapse = ", "))
  }
  
  # ========== Fix sample processing ==========
  # Reconstruct sample processing logic - avoid type issues in pipes
  sample_info$Sample_ID <- sapply(sample_info$title, function(x) {
    gsm_id <- stringr::str_extract(x, "GSM\\d+")
    ifelse(is.na(gsm_id), gsub(" .*", "", x), gsm_id)
  })
  
  # Process group information
  sample_info$Group <- sapply(sample_info$group, function(g) {
    if (grepl("^DPN", g, ignore.case = TRUE)) "DPN"
    else if (grepl("^DM", g, ignore.case = TRUE)) "DM"
    else if (grepl("^CN", g, ignore.case = TRUE)) "CN"
    else NA_character_
  })
  
  # Explicitly create new data frame
  sample_info_clean <- data.frame(
    Sample_ID = sample_info$Sample_ID,
    Group = sample_info$Group,
    stringsAsFactors = FALSE
  )
  
  # Remove NA values
  sample_info_clean <- sample_info_clean[
    !is.na(sample_info_clean$Group) & 
      !is.na(sample_info_clean$Sample_ID), ]
  
  # Remove duplicates
  sample_info_clean <- sample_info_clean[
    !duplicated(sample_info_clean$Sample_ID), ]
  # ========== End of fix ==========
  
  # Sample matching and reporting
  expr_samples <- colnames(expr_matrix)
  pheno_samples <- sample_info_clean$Sample_ID
  common_samples <- intersect(expr_samples, pheno_samples)
  unmatched_expr  <- setdiff(expr_samples,  pheno_samples)
  unmatched_pheno <- setdiff(pheno_samples, expr_samples)
  
  cat("Sample matching report:\n")
  cat(" - Expression matrix samples:", length(expr_samples), "\n")
  cat(" - Sample info samples:", length(pheno_samples), "\n")
  cat(" - Common samples:", length(common_samples), "\n")
  cat(" - Expression-only samples:", length(unmatched_expr))
  if (length(unmatched_expr) > 0) {
    cat(":", paste(head(unmatched_expr), collapse = ", "), "\n")
  } else {
    cat("\n")
  }
  cat(" - Pheno-only samples:", length(unmatched_pheno))
  if (length(unmatched_pheno) > 0) {
    cat(":", paste(head(unmatched_pheno), collapse = ", "), "\n")
  } else {
    cat("\n")
  }
  
  if (length(common_samples) < 3) {
    stop("Insufficient common samples (n=", length(common_samples),
         "). Analysis cannot proceed.")
  }
  
  expr_matrix <- expr_matrix[, common_samples, drop = FALSE]
  sample_info_final <- sample_info_clean[match(common_samples,
                                               sample_info_clean$Sample_ID), ]
  
  # Gene filtering: remove low-expression genes (lowest 20% by mean)
  original_genes <- nrow(expr_matrix)
  row_means      <- rowMeans(expr_matrix, na.rm = TRUE)
  keep_threshold <- quantile(row_means, 0.20, na.rm = TRUE)
  keep_genes     <- row_means > keep_threshold
  expr_matrix    <- expr_matrix[keep_genes, ]
  
  cat("Genes filtered:\n",
      " - Original genes:",            original_genes, "\n",
      " - Removed low expression:",    original_genes - sum(keep_genes), "\n",
      " - Final genes:",               nrow(expr_matrix), "\n")
  
  # Convert group to factor with desired order
  sample_info_final$Group <- factor(sample_info_final$Group,
                                    levels = c("CN", "DM", "DPN"))
  
  # Save preprocessed data
  save_path <- file.path(output_dir, "preprocessed_data.rds")
  saveRDS(list(expr_matrix      = expr_matrix,
               sample_info      = sample_info_final,
               data_transformed = data_transformed),
          file = save_path)
  cat("Preprocessed data saved to:", save_path, "\n")
  
  invisible(list(expr_matrix      = expr_matrix,
                 sample_info      = sample_info_final,
                 data_transformed = data_transformed))
}

# 2. Group variance analysis (PVCA-style) -----------------------------------
pvca_analysis_improved <- function(expr_matrix, group_info, output_dir,
                                   plot_width  = 8,
                                   plot_height = 6,
                                   plot_dpi    = 300) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  expr_matrix <- apply(expr_matrix, 2, as.numeric)
  rownames(expr_matrix) <- rownames(expr_matrix)
  
  gene_vars <- apply(expr_matrix, 1, var, na.rm = TRUE)
  zero_var_genes <- gene_vars == 0 | is.na(gene_vars)
  if (any(zero_var_genes)) {
    cat("PVCA: Removing zero-variance genes:", sum(zero_var_genes), "\n")
    expr_matrix <- expr_matrix[!zero_var_genes, , drop = FALSE]
  }
  
  pca <- prcomp(t(expr_matrix), scale. = FALSE)
  n_pcs <- min(30, ncol(pca$x), sum(pca$sdev > 1e-5))
  
  var_group <- numeric(n_pcs)
  for (i in seq_len(n_pcs)) {
    pc <- pca$x[, i]
    model_data <- data.frame(pc = pc, group = factor(group_info))
    model      <- lm(pc ~ group, data = model_data)
    anova_res  <- anova(model)
    
    ss_group <- anova_res["group", "Sum Sq"]
    ss_total <- sum(anova_res[, "Sum Sq"], na.rm = TRUE)
    
    var_group[i] <- if (ss_total > 0) ss_group / ss_total else 0
  }
  
  mean_var_group <- mean(var_group, na.rm = TRUE)
  
  results <- data.frame(PC = seq_len(n_pcs),
                        Variance_Group = var_group)
  write.csv(results,
            file.path(output_dir, "pvca_results.csv"),
            row.names = FALSE)
  
  # SCI-level aesthetic parameters
  base_size <- 12  # fixed font size
  
  # Create color gradient
  color_palette <- colorRampPalette(c("#3498db", "#2c3e50"))(n_pcs)
  
  # Generate bar plot showing proportion of variance explained by group for each PC
  pvca_plot <- ggplot(results,
                      aes(x = factor(PC), y = Variance_Group,
                          fill = factor(PC))) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    geom_text(aes(label = paste0(round(Variance_Group * 100, 1), "%")),
              vjust = -0.5, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = color_palette) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, max(results$Variance_Group) * 1.15),
                       expand = c(0, 0)) +
    labs(title = "Group-Explained Variance per Principal Component",
         subtitle = "Proportion of Variance Explained by Group in Each PC",
         x = "Principal Component (PC)",
         y = "Proportion of Variance Explained by Group") +
    theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "black"),
      plot.title = element_text(size = base_size + 2,
                                face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size,
                                   color = "grey40", hjust = 0.5),
      axis.title = element_text(size = base_size + 1, face = "bold"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      legend.position = "none",
      plot.margin = margin(15, 15, 15, 15),
      panel.background = element_blank()
    )
  
  # Add annotation for total group-explained variance
  pvca_plot <- pvca_plot +
    annotate("text", x = Inf, y = Inf,
             label = paste("Total Group-Explained Variance:", 
                           round(mean_var_group * 100, 1), "%"),
             hjust = 1.1, vjust = 1.5, size = 4, color = "red", fontface = "bold")
  
  # Dual-format output (PNG and SVG)
  ggsave(file.path(output_dir, "group_variance_plot.png"),
         pvca_plot,
         width = plot_width, height = plot_height, dpi = plot_dpi)
  ggsave(file.path(output_dir, "group_variance_plot.svg"),
         pvca_plot,
         width = plot_width, height = plot_height)
  
  saveRDS(list(group_variance = mean_var_group,
               detailed_results = results),
          file = file.path(output_dir, "pvca_results.rds"))
  
  invisible(list(group_variance = mean_var_group,
                 detailed_results = results,
                 plot = pvca_plot))
}

# 3. PCA visualization -----------------------------------------------------
visualize_group_effect <- function(expr_matrix, group_info,
                                   output_dir,
                                   detect_outliers = TRUE,
                                   title = "Group Effect Analysis",
                                   plot_width  = 10,
                                   plot_height = 8,
                                   plot_dpi    = 300) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  expr_matrix <- apply(expr_matrix, 2, as.numeric)
  rownames(expr_matrix) <- rownames(expr_matrix)
  
  gene_vars <- apply(expr_matrix, 1, var, na.rm = TRUE)
  zero_var_genes <- gene_vars == 0 | is.na(gene_vars)
  if (any(zero_var_genes)) {
    cat("PCA Visualization: Removing zero-variance genes:",
        sum(zero_var_genes), "\n")
    expr_matrix <- expr_matrix[!zero_var_genes, , drop = FALSE]
  }
  
  pca <- prcomp(t(expr_matrix), scale. = TRUE)
  pca_data <- as.data.frame(pca$x[, 1:2])
  colnames(pca_data) <- c("PC1", "PC2")
  pca_data$Sample <- rownames(pca_data)
  pca_data$Group  <- group_info
  
  variance <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 1)[1:2]
  
  if (detect_outliers) {
    pca_scores <- pca$x[, 1:2]
    center     <- colMeans(pca_scores)
    cov_mat    <- cov(pca_scores)
    mah_dist   <- mahalanobis(pca_scores, center, cov_mat)
    
    outlier_threshold <- quantile(mah_dist, 0.975, na.rm = TRUE)
    pca_data$Outlier  <- mah_dist > outlier_threshold
    pca_data$Mahalanobis_Distance <- mah_dist
  }
  
  # SCI-level aesthetic parameters
  base_size <- 12  # fixed font size
  
  # Create PCA plot
  p_group <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 4, alpha = 0.8) +
    stat_ellipse(level = 0.95, show.legend = FALSE, linewidth = 0.5) +
    labs(title = "PCA of Group Distribution",
         x = paste0("PC1 (", variance[1], "%)"),
         y = paste0("PC2 (", variance[2], "%)"),
         color = "Group") +
    scale_color_manual(values = c("CN" = "#2ecc71",
                                  "DM" = "#3498db",
                                  "DPN" = "#e74c3c")) +
    theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "black"),
      plot.title = element_text(size = base_size + 2,
                                face = "bold", hjust = 0.5),
      axis.title = element_text(size = base_size + 1, face = "bold"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      legend.position = "right",
      legend.background = element_blank(),
      panel.background = element_blank()
    )
  
  # Mark outliers if detected
  if (detect_outliers && any(pca_data$Outlier)) {
    label_size <- 3
    p_group <- p_group +
      geom_point(data = pca_data[pca_data$Outlier, ],
                 shape = 4, size = 5,
                 color = "black", stroke = 1) +
      geom_text_repel(data = pca_data[pca_data$Outlier, ],
                      aes(label = Sample),
                      color = "black", size = label_size,
                      box.padding = 0.8, max.overlaps = 50)
  }
  
  # Dual-format output
  ggsave(file.path(output_dir, "group_effect_pca.png"),
         p_group,
         width = plot_width, height = plot_height, dpi = plot_dpi)
  ggsave(file.path(output_dir, "group_effect_pca.svg"),
         p_group,
         width = plot_width, height = plot_height)
  
  saveRDS(pca_data, file.path(output_dir, "pca_data.rds"))
  
  pca_full_results <- list(coordinates = pca_data,
                           variance_prop = variance)
  saveRDS(pca_full_results,
          file.path(output_dir, "pca_full_results.rds"))
  
  invisible(p_group)
}

# 4. Violin plot of median expression by group ------------------------------
plot_expression_violin <- function(expr_matrix, sample_info, output_dir,
                                   plot_width  = 10,
                                   plot_height = 8,
                                   plot_dpi    = 300) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # SCI-level aesthetic parameters
  base_size <- 12  # fixed font size
  
  median_expr <- apply(expr_matrix, 2, median, na.rm = TRUE)
  plot_data <- data.frame(
    Sample           = names(median_expr),
    MedianExpression = median_expr,
    Group            = sample_info$Group
  )
  
  # Create violin plot (statistical annotations removed for simplicity)
  p <- ggplot(plot_data, aes(x = Group, y = MedianExpression, fill = Group)) +
    geom_violin(alpha = 0.7, trim = FALSE, linewidth = 0.3) +
    geom_boxplot(width = 0.15, fill = "white", alpha = 0.8,
                 outlier.shape = NA, linewidth = 0.3) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.6) +
    labs(title = "Median Expression Distribution by Group",
         x = "Group",
         y = "Median Expression Value") +
    scale_fill_manual(values = c("CN" = "#2ecc71",
                                 "DM" = "#3498db",
                                 "DPN" = "#e74c3c")) +
    theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "black"),
      plot.title = element_text(size = base_size + 2,
                                face = "bold", hjust = 0.5),
      axis.title = element_text(size = base_size + 1, face = "bold"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      legend.position = "none",
      panel.grid.major.y = element_line(color = "grey90",
                                        linewidth = 0.2),
      panel.grid.minor = element_blank()
    )
  
  # Dual-format output
  ggsave(file.path(output_dir, "expression_violin_plot.png"),
         p,
         width = plot_width, height = plot_height, dpi = plot_dpi)
  ggsave(file.path(output_dir, "expression_violin_plot.svg"),
         p,
         width = plot_width, height = plot_height)
  
  saveRDS(plot_data,
          file.path(output_dir, "violin_plot_data.rds"))
  invisible(p)
}

# 5. limma differential expression analysis --------------------------------
run_limma_DE_analysis <- function(expr_matrix,
                                  sample_info,
                                  output_prefix = "DE_Results",
                                  plot_width  = 12,
                                  plot_height = 10,
                                  plot_dpi    = 600,
                                  adj_pval_cutoff = 0.05,
                                  logFC_cutoff    = 1) {
  
  safe_output_prefix <- gsub("[^[:alnum:]._-]", "_", output_prefix)
  if (!dir.exists(safe_output_prefix)) {
    dir.create(safe_output_prefix, recursive = TRUE,
               showWarnings = FALSE)
  }
  
  # Design matrix
  design <- model.matrix(~ 0 + Group, data = sample_info)
  colnames(design) <- gsub("Group", "", colnames(design))
  
  # Define contrasts
  contrasts <- c(
    "DPN-CN",   # DPN vs CN
    "DM-CN",    # DM vs CN
    "DPN-DM"    # DPN vs DM
  )
  
  cm <- makeContrasts(contrasts = contrasts, levels = design)
  
  fit   <- lmFit(expr_matrix, design)
  fit2  <- contrasts.fit(fit, cm)
  fit2  <- eBayes(fit2, trend = TRUE, robust = TRUE)
  
  # Save results for each comparison
  all_results <- list()
  for (i in seq_along(contrasts)) {
    contrast_name <- contrasts[i]
    results <- topTable(fit2, coef = i, number = Inf,
                        adjust.method = "BH", sort.by = "P")
    results$Gene <- rownames(results)
    colnames(results) <- c("logFC", "AveExpr", "t", "P.Value",
                           "adj.P.Val", "B", "Gene")
    
    # Save all results
    all_results[[contrast_name]] <- results
    
    # Filter significant genes
    sig_genes <- results %>%
      filter(!is.na(adj.P.Val),
             adj.P.Val < adj_pval_cutoff,
             abs(logFC) > logFC_cutoff)
    
    up_genes   <- sig_genes %>% filter(logFC >  0) %>%
      arrange(desc(logFC))
    down_genes <- sig_genes %>% filter(logFC <  0) %>%
      arrange(logFC)
    
    # Create contrast-specific directory
    contrast_dir <- file.path(safe_output_prefix, contrast_name)
    if (!dir.exists(contrast_dir))
      dir.create(contrast_dir, recursive = TRUE)
    
    write.csv(results,
              file.path(contrast_dir, "DE_results_all.csv"),
              row.names = FALSE)
    
    write.csv(up_genes,
              file.path(contrast_dir, "DE_upregulated.csv"),
              row.names = FALSE)
    write.csv(down_genes,
              file.path(contrast_dir, "DE_downregulated.csv"),
              row.names = FALSE)
    
    # Generate volcano plot
    volcano_plot <- create_enhanced_volcano(
      res_df    = results,
      contrast_name = contrast_name,
      group_counts = c(
        DPN = sum(sample_info$Group == "DPN"),
        DM  = sum(sample_info$Group == "DM"),
        CN  = sum(sample_info$Group == "CN")
      ),
      output_file = file.path(contrast_dir, "volcano_plot"),
      plot_width  = plot_width,
      plot_height = plot_height,
      plot_dpi    = plot_dpi
    )
    
    cat("\n--- ", contrast_name, " ---\n")
    cat("Significant genes (adj.P.Val<", adj_pval_cutoff,
        " & |logFC|>", logFC_cutoff, "):", nrow(sig_genes), "\n")
    cat("Upregulated:", nrow(up_genes), "\n")
    cat("Downregulated:", nrow(down_genes), "\n")
  }
  
  # Summary report
  cat("\n===== Differential Expression Analysis Summary (limma) =====\n")
  cat("CN  samples:", sum(sample_info$Group == "CN"), "\n")
  cat("DM  samples:", sum(sample_info$Group == "DM"), "\n")
  cat("DPN samples:", sum(sample_info$Group == "DPN"), "\n")
  cat("Total genes analyzed:", nrow(expr_matrix), "\n")
  
  # Save all results to RDS
  saveRDS(list(limma_fit   = fit2,
               de_results  = all_results),
          file = file.path(safe_output_prefix,
                           "de_analysis_results.rds"))
  
  invisible(list(limma_fit   = fit2,
                 de_results  = all_results))
}

# 6. Enhanced volcano plot --------------------------------------------------
create_enhanced_volcano <- function(res_df,
                                    contrast_name,
                                    group_counts,
                                    output_file = "Enhanced_Volcano",
                                    plot_width  = 12,
                                    plot_height = 10,
                                    plot_dpi    = 600) {
  
  logFC_filter   <- 1
  padj_threshold <- 0.05
  top_n          <- 15
  
  volcano_data <- res_df %>%
    mutate(
      group = case_when(
        !is.na(adj.P.Val) & adj.P.Val < padj_threshold &
          logFC >  logFC_filter ~ "Upregulated",
        !is.na(adj.P.Val) & adj.P.Val < padj_threshold &
          logFC < -logFC_filter ~ "Downregulated",
        TRUE ~ "Not significant"
      ),
      neg_log10_padj = -log10(adj.P.Val)
    ) %>%
    mutate(neg_log10_padj = ifelse(is.infinite(neg_log10_padj),
                                   max(neg_log10_padj[
                                     !is.infinite(neg_log10_padj)],
                                     na.rm = TRUE) * 1.1,
                                   neg_log10_padj))
  
  max_logFC          <- max(abs(volcano_data$logFC), na.rm = TRUE) * 1.1
  max_neg_log10_padj <- max(volcano_data$neg_log10_padj, na.rm = TRUE) * 1.1
  
  up_genes   <- volcano_data %>%
    filter(group == "Upregulated") %>%
    arrange(adj.P.Val, desc(abs(logFC))) %>%
    head(top_n)
  
  down_genes <- volcano_data %>%
    filter(group == "Downregulated") %>%
    arrange(adj.P.Val, desc(abs(logFC))) %>%
    head(top_n)
  
  color_gradient <- c("#3288bd", "#66c2a5", "#ffffbf", "#f46d43", "#9e0142")
  base_size      <- 12
  label_size     <- 3.5
  
  # Parse contrast groups
  groups <- strsplit(contrast_name, "-")[[1]]
  group1 <- groups[1]
  group2 <- groups[2]
  
  p <- ggplot(volcano_data, aes(x = logFC, y = neg_log10_padj)) +
    geom_point(aes(color = logFC, size = neg_log10_padj), alpha = 0.6) +
    geom_point(
      data = filter(volcano_data, group != "Not significant"),
      aes(fill = logFC, size = neg_log10_padj),
      shape = 21, color = "black", show.legend = FALSE
    ) +
    geom_text_repel(
      data = up_genes,
      aes(label = Gene),
      box.padding = 0.5,
      nudge_x = max_logFC * 0.05,
      nudge_y = max_neg_log10_padj * 0.05,
      segment.curvature = -0.1,
      segment.ncp = 3,
      segment.size = 0.3,
      direction = "y",
      hjust = 0,
      min.segment.length = 0,
      size = label_size,
      max.overlaps = 100
    ) +
    geom_text_repel(
      data = down_genes,
      aes(label = Gene),
      box.padding = 0.5,
      nudge_x = -max_logFC * 0.05,
      nudge_y = max_neg_log10_padj * 0.05,
      segment.curvature = -0.1,
      segment.ncp = 3,
      segment.size = 0.3,
      direction = "y",
      hjust = 1,
      min.segment.length = 0,
      size = label_size,
      max.overlaps = 100
    ) +
    geom_vline(xintercept = c(-logFC_filter, logFC_filter),
               linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(padj_threshold),
               linetype = "dashed", color = "grey50") +
    scale_color_gradientn(colours = color_gradient,
                          values  = scales::rescale(seq(-max_logFC,
                                                        max_logFC,
                                                        length.out = 5)),
                          limits  = c(-max_logFC, max_logFC)) +
    scale_fill_gradientn(colours = color_gradient,
                         values = scales::rescale(seq(-max_logFC,
                                                      max_logFC,
                                                      length.out = 5)),
                         limits = c(-max_logFC, max_logFC)) +
    scale_size_continuous(range = c(1, 6)) +
    coord_cartesian(xlim = c(-max_logFC, max_logFC),
                    ylim = c(0, max_neg_log10_padj)) +
    labs(
      title = paste("Differentially Expressed Genes:",
                    group1, "vs", group2, "(limma)"),
      subtitle = sprintf(
        "%s: %d samples | %s: %d samples\nSignificant genes: %d (Up: %d, Down: %d)",
        group1, group_counts[[group1]],
        group2, group_counts[[group2]],
        nrow(up_genes) + nrow(down_genes),
        nrow(up_genes), nrow(down_genes)
      ),
      x = paste("Log2 Fold Change (", group1, "vs", group2, ")"),
      y = "-Log10(Adjusted P-value)",
      color = "Log2 Fold Change",
      size  = "-Log10(Adj.P)",
      caption = sprintf("Dashed lines: |logFC| > %.1f & adj.P.Val < %.2f",
                        logFC_filter, padj_threshold)
    ) +
    theme_bw() +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold",
                                   size = base_size * 1.5),
      plot.subtitle = element_text(hjust = 0.5, size = base_size * 1,
                                   margin = margin(b = 15)),
      plot.caption  = element_text(hjust = 0.5,
                                   size = base_size * 0.8,
                                   color = "grey40"),
      axis.title    = element_text(size = base_size * 1.2,
                                   face = "bold"),
      axis.text     = element_text(size = base_size * 0.9),
      legend.title  = element_text(size = base_size * 1.0),
      legend.text   = element_text(size = base_size * 0.8),
      panel.grid    = element_blank(),
      legend.position = "right",
      legend.box      = "vertical",
      legend.background = element_rect(color = "grey80",
                                       fill = "white"),
      panel.border    = element_rect(color = "black", size = 1)
    )
  
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Dual-format output
  ggsave(paste0(output_file, ".png"), p,
         width = plot_width, height = plot_height, dpi = plot_dpi)
  ggsave(paste0(output_file, ".svg"), p,
         width = plot_width, height = plot_height)
  
  saveRDS(volcano_data,
          file = file.path(output_dir, "volcano_plot_data.rds"))
  invisible(p)
}

# 7. Generate downstream analysis matrices ----------------------------------
generate_downstream_matrices <- function(expr_matrix, sample_info, output_dir) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  write.csv(expr_matrix,
            file.path(output_dir, "Expression_Matrix.csv"))
  
  write.csv(expr_matrix,
            file.path(output_dir, "WGCNA_Matrix.csv"))
  
  write.csv(expr_matrix,
            file.path(output_dir, "GSEA_Matrix.csv"))
  
  ml_matrix <- data.frame(t(expr_matrix))
  ml_matrix$Group <- sample_info$Group
  write.csv(ml_matrix,
            file.path(output_dir, "ML_Matrix.csv"),
            row.names = TRUE)
  
  annotation_data <- data.frame(
    Sample = sample_info$Sample_ID,
    Group  = sample_info$Group
  )
  write.csv(annotation_data,
            file.path(output_dir, "sample_annotation.csv"),
            row.names = FALSE)
  
  saveRDS(list(expr_matrix       = expr_matrix,
               wgcna_matrix      = expr_matrix,
               gsea_matrix       = expr_matrix,
               ml_matrix         = ml_matrix,
               sample_annotation = annotation_data),
          file = file.path(output_dir, "downstream_matrices.rds"))
  
  cat("\n===== Downstream Analysis Matrices Generated =====\n")
  cat("Expression Matrix: Expression_Matrix.csv\n")
  cat("WGCNA Matrix:      WGCNA_Matrix.csv\n")
  cat("GSEA Matrix:       GSEA_Matrix.csv\n")
  cat("ML Matrix:         ML_Matrix.csv\n")
  cat("Sample Annotation: sample_annotation.csv\n")
  
  invisible(list(wgcna_matrix      = expr_matrix,
                 gsea_matrix       = expr_matrix,
                 ml_matrix         = ml_matrix,
                 sample_annotation = annotation_data))
}

# 8. Main execution function ------------------------------------------------
run_full_analysis_pipeline <- function(expr_file   = "GSE95849_expression_matrix.csv",
                                       sample_file = "pheno_clean.csv",
                                       output_root = ".") {
  tryCatch({
    data <- load_and_preprocess_data(
      expr_file   = expr_file,
      sample_file = sample_file,
      output_dir  = file.path(output_root, "Preprocessed_Data")
    )
    
    cat("\n===== Group Effect Analysis =====\n")
    pvca_results <- pvca_analysis_improved(
      expr_matrix = data$expr_matrix,
      group_info  = data$sample_info$Group,
      output_dir  = file.path(output_root, "Group_Effect_Analysis"),
      plot_width  = 10,
      plot_height = 8,
      plot_dpi    = 250
    )
    
    pca_plot <- visualize_group_effect(
      expr_matrix   = data$expr_matrix,
      group_info    = data$sample_info$Group,
      output_dir    = file.path(output_root, "Group_Effect_Analysis"),
      detect_outliers = FALSE,
      plot_width    = 10,
      plot_height   = 8,
      plot_dpi      = 250
    )
    
    violin_plot <- plot_expression_violin(
      expr_matrix = data$expr_matrix,
      sample_info = data$sample_info,
      output_dir  = file.path(output_root, "Group_Effect_Analysis"),
      plot_width  = 8,
      plot_height = 6,
      plot_dpi    = 250
    )
    
    cat("\n===== Differential Expression Analysis (limma) =====\n")
    de_results <- run_limma_DE_analysis(
      expr_matrix   = data$expr_matrix,
      sample_info   = data$sample_info,
      output_prefix = file.path(output_root, "DE_Results"),
      plot_width    = 14,
      plot_height   = 12,
      plot_dpi      = 300
    )
    
    cat("\n===== Generating Downstream Matrices =====\n")
    matrices <- generate_downstream_matrices(
      expr_matrix = data$expr_matrix,
      sample_info = data$sample_info,
      output_dir  = file.path(output_root, "Downstream_Matrices")
    )
    
    saveRDS(list(group_effect_analysis = list(pvca = pvca_results,
                                              pca  = pca_plot,
                                              violin = violin_plot),
                 de_analysis          = de_results,
                 downstream_matrices  = matrices),
            file = file.path(output_root, "full_analysis_results.rds"))
    cat("Full analysis results saved to: full_analysis_results.rds\n")
    
    invisible(list(group_effect_analysis = list(pvca = pvca_results,
                                                pca  = pca_plot,
                                                violin = violin_plot),
                   de_analysis          = de_results,
                   downstream_matrices  = matrices))
    
  }, error = function(e) {
    message("\n!!! ANALYSIS FAILED !!!")
    message("Error: ", e$message)
    message("Traceback:")
    traceback(2)
    NULL
  })
}

# Execute the analysis pipeline
analysis_results <- run_full_analysis_pipeline()

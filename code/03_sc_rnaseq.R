# ====================================================================
# Single-cell RNA-seq Integrated Analysis Pipeline
# Includes: Quality Control, Normalization, Clustering
# Core workflow preserved, enhanced for generalizability and code sharing
# ====================================================================

# Load required packages
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(patchwork)
  library(Matrix)
  library(scDblFinder)
  library(BiocParallel)
  library(clustree)
  library(RColorBrewer)
  library(viridis)
  library(scales)
  library(harmony)
  library(future)
})

# ====================================================================
# Part 1: Quality Control Functions
# ====================================================================

# 1.1 Calculate mitochondrial percentage ------------------------------------
calculate_mt_percentage <- function(seurat_obj) {
  # Try mitochondrial gene prefixes for different species
  mt_identifiers <- c(
    "^MT-", "^mt-", "^MT\\.", "^mt\\.", "^MT_", "^mt_", 
    "^MTRNR", "^MTRNR2L", "^MTATP", "^MTCO", "^MTND", 
    "^MTCYB", "^MTTF", "^MTTH", "^MTTL", "^MTTM", "^MTTN", "^MTTS"
  )
  
  # Identify mitochondrial genes
  mt_genes <- unique(unlist(lapply(mt_identifiers, function(prefix) {
    grep(pattern = prefix, 
         x = rownames(seurat_obj), 
         value = TRUE,
         ignore.case = TRUE)
  })))
  
  # Handle case when no mitochondrial genes are found
  if (length(mt_genes) == 0) {
    warning("No mitochondrial genes detected! Using zero values as placeholder.")
    seurat_obj[["percent.mt"]] <- 0
    return(seurat_obj)
  }
  
  # Calculate mitochondrial percentage
  counts <- Seurat::GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  mt_counts <- Matrix::colSums(counts[mt_genes, , drop = FALSE])
  total_counts <- Matrix::colSums(counts)
  percent.mt <- (mt_counts / total_counts) * 100
  
  # Add to metadata
  seurat_obj[["percent.mt"]] <- percent.mt
  
  return(seurat_obj)
}

# 1.2 Calculate doublet risk -------------------------------------------------
calculate_doublet_risk <- function(seurat_obj, sample_id = NULL) {
  # If sample ID is provided, only evaluate that sample
  if (!is.null(sample_id)) {
    if (!sample_id %in% seurat_obj$orig.ident) {
      warning("Sample '", sample_id, "' not found in object. Calculating overall risk.")
    } else {
      sample_cells <- seurat_obj$orig.ident == sample_id
      metadata <- seurat_obj@meta.data[sample_cells, ]
    }
  }
  
  # Fallback to entire object if no sample ID or sample not found
  if (is.null(sample_id) || !exists("metadata")) {
    metadata <- seurat_obj@meta.data
  }
  
  # Calculate median of key metrics
  median_genes <- median(metadata$nFeature_RNA, na.rm = TRUE)
  median_umi <- median(metadata$nCount_RNA, na.rm = TRUE)
  median_mt <- median(metadata$percent.mt, na.rm = TRUE)
  
  # Dynamically adjust thresholds based on sample properties
  gene_threshold <- ifelse(median_genes > 2500, 1.5, 
                           ifelse(median_genes > 1500, 1.8, 2.0))
  umi_threshold <- ifelse(median_umi > 8000, 1.5, 
                          ifelse(median_umi > 4000, 1.8, 2.0))
  mt_threshold <- ifelse(median_mt > 10, 0.6, 0.5)
  
  # Evaluate three key criteria
  gene_outlier <- metadata$nFeature_RNA > gene_threshold * median_genes
  umi_outlier <- metadata$nCount_RNA > umi_threshold * median_umi
  mt_low <- metadata$percent.mt < mt_threshold * median_mt
  
  # Composite risk score
  doublet_risk_score <- (gene_outlier + umi_outlier + mt_low) / 3
  
  # High risk: at least 2 criteria met
  high_risk_cells <- sum(doublet_risk_score >= 2/3, na.rm = TRUE)
  
  # Percentage of high-risk cells
  risk_percent <- high_risk_cells / nrow(metadata) * 100
  
  return(round(risk_percent, 1))
}

# 1.3 Generate initial QC report ---------------------------------------------
generate_qc_report <- function(seurat_obj, 
                               min_genes = 500, 
                               max_genes = 6000,
                               min_umi = 1000,
                               max_umi = 30000,
                               max_mt = 20) {
  
  meta_data <- seurat_obj@meta.data
  
  # Ensure mitochondrial percentage exists
  if (!"percent.mt" %in% colnames(meta_data)) {
    warning("percent.mt column missing. Calculating now...")
    seurat_obj <- calculate_mt_percentage(seurat_obj)
    meta_data <- seurat_obj@meta.data
  }
  
  # Compute QC stats per sample
  qc_stats <- lapply(unique(meta_data$orig.ident), function(sample_id) {
    sample_data <- meta_data[meta_data$orig.ident == sample_id, ]
    
    low_quality <- sum(
      sample_data$nFeature_RNA < min_genes | 
        sample_data$nFeature_RNA > max_genes | 
        sample_data$nCount_RNA < min_umi |
        sample_data$nCount_RNA > max_umi |
        sample_data$percent.mt > max_mt,
      na.rm = TRUE
    )
    
    doublet_risk <- calculate_doublet_risk(seurat_obj, sample_id)
    
    data.frame(
      orig.ident = sample_id,
      Cells = nrow(sample_data),
      Median_Genes = median(sample_data$nFeature_RNA, na.rm = TRUE),
      Median_UMI = median(sample_data$nCount_RNA, na.rm = TRUE),
      Median_MT = median(sample_data$percent.mt, na.rm = TRUE),
      Low_Quality_Cells = low_quality,
      Low_Quality_Pct = round(low_quality / nrow(sample_data) * 100, 1),
      Doublet_Risk = doublet_risk,
      Used_Min_Genes = min_genes,
      Used_Max_Genes = max_genes,
      Used_Max_MT = max_mt
    )
  }) %>% bind_rows()
  
  # Overall statistics
  low_quality_total <- sum(
    meta_data$nFeature_RNA < min_genes | 
      meta_data$nFeature_RNA > max_genes | 
      meta_data$nCount_RNA < min_umi |
      meta_data$nCount_RNA > max_umi |
      meta_data$percent.mt > max_mt,
    na.rm = TRUE
  )
  
  overall_doublet_risk <- calculate_doublet_risk(seurat_obj)
  
  total_stats <- data.frame(
    orig.ident = "Overall",
    Cells = nrow(meta_data),
    Median_Genes = median(meta_data$nFeature_RNA, na.rm = TRUE),
    Median_UMI = median(meta_data$nCount_RNA, na.rm = TRUE),
    Median_MT = median(meta_data$percent.mt, na.rm = TRUE),
    Low_Quality_Cells = low_quality_total,
    Low_Quality_Pct = round(low_quality_total / nrow(meta_data) * 100, 1),
    Doublet_Risk = overall_doublet_risk,
    Used_Min_Genes = min_genes,
    Used_Max_Genes = max_genes,
    Used_Max_MT = max_mt
  )
  
  qc_report <- bind_rows(qc_stats, total_stats)
  return(qc_report)
}

# 1.4 Generate post-filtering QC report --------------------------------------
generate_qc_report_with_thresholds <- function(seurat_obj, sample_thresholds) {
  metadata_df <- seurat_obj@meta.data
  
  if (!"percent.mt" %in% colnames(metadata_df)) {
    warning("percent.mt column missing. Calculating now...")
    seurat_obj <- calculate_mt_percentage(seurat_obj)
    metadata_df <- seurat_obj@meta.data
  }
  
  qc_stats <- lapply(unique(metadata_df$orig.ident), function(sample_id) {
    sample_cells <- metadata_df$orig.ident == sample_id
    sample_data <- metadata_df[sample_cells, ]
    thresholds <- sample_thresholds[[sample_id]]
    
    low_quality <- sum(
      sample_data$nFeature_RNA < thresholds$min_genes |
        sample_data$nFeature_RNA > thresholds$max_genes |
        sample_data$nCount_RNA < thresholds$min_umi |
        sample_data$nCount_RNA > thresholds$max_umi |
        sample_data$percent.mt > thresholds$max_mt,
      na.rm = TRUE
    )
    
    doublet_risk <- calculate_doublet_risk(seurat_obj, sample_id)
    
    data.frame(
      orig.ident = sample_id,
      Cells = nrow(sample_data),
      Median_Genes = median(sample_data$nFeature_RNA, na.rm = TRUE),
      Median_UMI = median(sample_data$nCount_RNA, na.rm = TRUE),
      Median_MT = median(sample_data$percent.mt, na.rm = TRUE),
      Low_Quality_Cells = low_quality,
      Low_Quality_Pct = round(low_quality / nrow(sample_data) * 100, 1),
      Doublet_Risk = doublet_risk,
      Used_Min_Genes = thresholds$min_genes,
      Used_Max_Genes = thresholds$max_genes,
      Used_Max_MT = thresholds$max_mt,
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  
  # Overall metrics using sample-specific thresholds
  total_low_quality <- 0
  for (sample_id in unique(metadata_df$orig.ident)) {
    sample_data <- metadata_df[metadata_df$orig.ident == sample_id, ]
    thresholds <- sample_thresholds[[sample_id]]
    total_low_quality <- total_low_quality + sum(
      sample_data$nFeature_RNA < thresholds$min_genes |
        sample_data$nFeature_RNA > thresholds$max_genes |
        sample_data$nCount_RNA < thresholds$min_umi |
        sample_data$nCount_RNA > thresholds$max_umi |
        sample_data$percent.mt > thresholds$max_mt,
      na.rm = TRUE
    )
  }
  
  overall_doublet_risk <- calculate_doublet_risk(seurat_obj)
  
  total_stats <- data.frame(
    orig.ident = "Overall",
    Cells = nrow(metadata_df),
    Median_Genes = median(metadata_df$nFeature_RNA, na.rm = TRUE),
    Median_UMI = median(metadata_df$nCount_RNA, na.rm = TRUE),
    Median_MT = median(metadata_df$percent.mt, na.rm = TRUE),
    Low_Quality_Cells = total_low_quality,
    Low_Quality_Pct = round(total_low_quality / nrow(metadata_df) * 100, 1),
    Doublet_Risk = overall_doublet_risk,
    Used_Min_Genes = NA_real_,
    Used_Max_Genes = NA_real_,
    Used_Max_MT = NA_real_,
    stringsAsFactors = FALSE
  )
  
  qc_report <- bind_rows(qc_stats, total_stats)
  return(qc_report)
}

# 1.5 Sample-specific filtering ---------------------------------------------
apply_sample_specific_filters <- function(seurat_obj, output_dir) {
  filtered_cells <- c()
  filter_log <- data.frame()
  sample_ids <- unique(seurat_obj$orig.ident)
  sample_thresholds <- list()
  
  for (sample_id in sample_ids) {
    sample_cells <- colnames(seurat_obj)[seurat_obj$orig.ident == sample_id]
    sample_data <- seurat_obj@meta.data[sample_cells, ]
    
    # Calculate sample-specific distribution metrics
    Q_genes <- quantile(sample_data$nFeature_RNA, 
                        probs = c(0.05, 0.25, 0.75, 0.95), 
                        na.rm = TRUE)
    Q_mt <- quantile(sample_data$percent.mt, 
                     probs = c(0.75, 0.95), 
                     na.rm = TRUE)
    IQR_genes <- Q_genes[3] - Q_genes[2]
    
    # Gene count thresholds
    sample_min_genes <- max(200, floor(Q_genes[1]))
    sample_max_genes <- ifelse(IQR_genes > 2000,
                               min(10000, ceiling(Q_genes[4] + IQR_genes * 0.5)),
                               min(8000, ceiling(Q_genes[3] + IQR_genes * 1.5)))
    
    # Mitochondrial threshold
    sample_max_mt <- min(30, ceiling(Q_mt[2] + 2))
    
    # Store sample-specific thresholds
    sample_thresholds[[sample_id]] <- list(
      min_genes = sample_min_genes,
      max_genes = sample_max_genes,
      min_umi = 500,
      max_umi = 50000,
      max_mt = sample_max_mt
    )
    
    # Apply basic filtering
    valid_cells <- rownames(sample_data)[
      sample_data$nFeature_RNA > sample_min_genes &
        sample_data$nFeature_RNA < sample_max_genes &
        sample_data$percent.mt < sample_max_mt
    ]
    
    # Preserve high-gene, low-MT cells (likely doublets but may be biologically relevant)
    high_gene_cells <- rownames(sample_data)[
      sample_data$nFeature_RNA > 7000 & 
        sample_data$percent.mt < 10
    ]
    valid_cells <- unique(c(valid_cells, high_gene_cells))
    
    # Doublet risk estimation
    median_genes <- median(sample_data$nFeature_RNA, na.rm = TRUE)
    doublet_risk <- sum(sample_data$nFeature_RNA > 2 * median_genes) / nrow(sample_data) * 100
    doublet_risk <- round(doublet_risk, 1)
    cat("\nSample", sample_id, "doublet risk:", doublet_risk, "%")
    
    removed_doublets <- 0
    
    # Run doublet detection only for high-risk samples
    if (doublet_risk > 15) {
      temp_obj <- subset(seurat_obj, cells = valid_cells)
      temp_sce <- Seurat::as.SingleCellExperiment(temp_obj)
      sample_labels <- temp_sce$orig.ident
      
      dbr_value <- min(0.3, max(0.05, doublet_risk/100 * 1.5))
      cat(" | Running scDblFinder with dbr =", dbr_value)
      
      # OS-adaptive parallel processing
      if (.Platform$OS.type == "windows") {
        bpparam <- BiocParallel::SnowParam(workers = min(4, parallel::detectCores() - 1))
      } else {
        bpparam <- BiocParallel::MulticoreParam(workers = min(4, parallel::detectCores() - 1))
      }
      
      temp_sce <- scDblFinder::scDblFinder(
        temp_sce,
        samples = sample_labels,
        dbr = dbr_value,
        BPPARAM = bpparam
      )
      
      singlet_cells <- colnames(temp_sce)[temp_sce$scDblFinder.class == "singlet"]
      removed_doublets <- length(valid_cells) - length(singlet_cells)
      
      if (removed_doublets > 0) {
        cat(" | Removed", removed_doublets, "doublets")
        valid_cells <- singlet_cells
      }
    }
    
    # Log detailed statistics
    sample_log <- data.frame(
      Sample = sample_id,
      Initial_Cells = length(sample_cells),
      Filtered_Cells = length(valid_cells),
      Retention_Rate = round(length(valid_cells)/length(sample_cells)*100, 1),
      Min_Genes = sample_min_genes,
      Max_Genes = sample_max_genes,
      Max_MT = sample_max_mt,
      Doublet_Risk = doublet_risk,
      Doublet_Removed = removed_doublets,
      Q1_Genes = Q_genes[2],
      Q3_Genes = Q_genes[3],
      p95_MT = Q_mt[2],
      IQR_Genes = IQR_genes
    )
    filter_log <- rbind(filter_log, sample_log)
    filtered_cells <- c(filtered_cells, valid_cells)
    
    cat("\n - IQR genes:", round(IQR_genes))
    cat("\n - Min genes:", sample_min_genes)
    cat("\n - Max genes:", sample_max_genes)
    cat("\n - Max MT%:", sample_max_mt)
    cat("\n - Cells retained:", length(valid_cells), "/", length(sample_cells),
        "(", round(length(valid_cells)/length(sample_cells)*100), "%)")
  }
  
  write.csv(filter_log, file.path(output_dir, "sample_filter_log.csv"), row.names = FALSE)
  
  return(list(
    filtered_obj = subset(seurat_obj, cells = filtered_cells),
    thresholds = sample_thresholds
  ))
}

# 1.6 QC visualizations ------------------------------------------------------
create_qc_violin_plot <- function(seurat_obj, output_dir, prefix) {
  qc_data <- Seurat::FetchData(
    seurat_obj,
    vars = c("nFeature_RNA", "nCount_RNA", "percent.mt", "orig.ident")
  )
  
  qc_long <- tidyr::pivot_longer(
    qc_data,
    cols = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    names_to = "Metric",
    values_to = "Value"
  )
  
  metric_labels <- c(
    "nFeature_RNA" = "Genes per Cell",
    "nCount_RNA" = "UMIs per Cell",
    "percent.mt" = "Mitochondrial %"
  )
  
  sample_count <- length(unique(seurat_obj$orig.ident))
  color_palette <- scales::hue_pal()(sample_count)
  
  p <- ggplot(qc_long, aes(x = orig.ident, y = Value, fill = orig.ident)) +
    geom_violin(scale = "width", adjust = 1.2, trim = TRUE, alpha = 0.7) +
    geom_boxplot(
      width = 0.15,
      outlier.shape = NA,
      alpha = 0.8,
      position = position_dodge(0.9)
    ) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9),
      size = 0.4,
      alpha = 0.15,
      color = "gray30"
    ) +
    scale_fill_manual(values = color_palette, name = "Sample") +
    scale_y_continuous(labels = scales::comma) +
    facet_wrap(
      ~Metric,
      scales = "free_y",
      nrow = 1,
      labeller = labeller(Metric = metric_labels)
    ) +
    labs(
      title = "Single-cell RNA-seq Quality Metrics",
      x = "Sample",
      y = "Value"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)),
      axis.title = element_text(size = 13, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, vjust = 1),
      axis.text.y = element_text(size = 10),
      legend.position = "none",
      strip.text = element_text(size = 12, face = "bold", margin = margin(b = 8)),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.margin = margin(15, 15, 15, 15),
      strip.background = element_rect(fill = "grey95", color = NA)
    )
  
  ggsave(file.path(output_dir, paste0(prefix, "_qc_violin_plot.png")), 
         plot = p, width = 14, height = 6, dpi = 300)
  ggsave(file.path(output_dir, paste0(prefix, "_qc_violin_plot.svg")), 
         plot = p, width = 14, height = 6)
  
  return(p)
}

create_qc_scatter_plots <- function(seurat_obj, output_dir, prefix) {
  feature_scatter <- Seurat::FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    group.by = "orig.ident",
    pt.size = 0.7
  ) +
    geom_smooth(method = "lm", se = FALSE, color = "darkred", linewidth = 0.8) +
    labs(
      title = "Gene vs UMI Count Correlation",
      x = "UMI Count per Cell",
      y = "Genes Detected per Cell"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  mito_scatter <- Seurat::FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "percent.mt",
    group.by = "orig.ident",
    pt.size = 0.7
  ) +
    labs(
      title = "Mitochondrial vs UMI Count",
      x = "UMI Count per Cell",
      y = "Mitochondrial %"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  combined_scatters <- (feature_scatter + mito_scatter) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom", 
          legend.text = element_text(size = 10))
  
  ggsave(file.path(output_dir, paste0(prefix, "_qc_scatter_plot.png")), 
         plot = combined_scatters, width = 14, height = 7, dpi = 300)
  ggsave(file.path(output_dir, paste0(prefix, "_qc_scatter_plot.svg")), 
         plot = combined_scatters, width = 14, height = 7)
  
  return(combined_scatters)
}

visualize_qc_report <- function(qc_report, output_dir, prefix) {
  plot_data <- qc_report %>%
    tidyr::pivot_longer(
      cols = c(Median_Genes, Median_UMI, Median_MT, Low_Quality_Pct, Doublet_Risk),
      names_to = "Metric",
      values_to = "Value"
    )
  
  metric_labels <- c(
    "Median_Genes" = "Median Genes/Cell",
    "Median_UMI" = "Median UMIs/Cell",
    "Median_MT" = "Median MT %",
    "Low_Quality_Pct" = "Low-Quality Cells (%)",
    "Doublet_Risk" = "Doublet Risk (%)"
  )
  
  p <- ggplot(plot_data, aes(x = orig.ident, y = Metric, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.7, alpha = 0.9) +
    geom_text(
      aes(label = ifelse(
        Metric == "Median_MT" | Metric == "Low_Quality_Pct" | Metric == "Doublet_Risk",
        sprintf("%.1f%%", Value),
        format(round(Value), big.mark = ",")
      )),
      color = "black",
      size = 4.5,
      fontface = "bold"
    ) +
    scale_fill_gradientn(
      colours = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"),
      values = scales::rescale(c(0, 10, 20, 30, 100)),
      limits = c(0, 100),
      na.value = "grey70",
      name = "Value"
    ) +
    scale_y_discrete(labels = metric_labels) +
    labs(
      title = "Single-cell RNA-seq Quality Assessment",
      x = "Sample",
      y = "Quality Metric"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)),
      axis.title = element_text(size = 13, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11, vjust = 1),
      axis.text.y = element_text(size = 11, face = "bold", hjust = 1),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      panel.grid = element_blank(),
      plot.margin = margin(15, 25, 15, 15)
    )
  
  ggsave(file.path(output_dir, paste0(prefix, "_qc_heatmap_report.png")), 
         plot = p, width = 10, height = 6, dpi = 300)
  ggsave(file.path(output_dir, paste0(prefix, "_qc_heatmap_report.svg")), 
         plot = p, width = 10, height = 6)
  
  return(p)
}

visualize_qc_comparison <- function(raw_report, filtered_report, output_dir, prefix) {
  convert_to_char <- function(df) {
    df %>%
      mutate(
        Used_Min_Genes = as.character(Used_Min_Genes),
        Used_Max_Genes = as.character(Used_Max_Genes),
        Used_Max_MT = as.character(Used_Max_MT)
      )
  }
  
  raw_report <- convert_to_char(raw_report)
  filtered_report <- convert_to_char(filtered_report)
  
  comparison_data <- bind_rows(
    mutate(raw_report, Stage = "Raw"),
    mutate(filtered_report, Stage = "Filtered")
  )
  
  p <- ggplot(comparison_data, aes(x = Stage, y = Low_Quality_Pct, fill = Stage)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_text(aes(label = paste0(Low_Quality_Pct, "%")), 
              vjust = -0.8, size = 4, fontface = "bold", color = "darkred") +
    facet_wrap(~orig.ident, scales = "free_y") +
    scale_fill_manual(values = c("Raw" = "#4E79A7", "Filtered" = "#F28E2B")) +
    labs(
      title = "Quality Improvement After Filtering",
      x = "Processing Stage",
      y = "Low-Quality Cells (%)",
      fill = "Stage"
    ) +
    ylim(0, max(comparison_data$Low_Quality_Pct) * 1.2) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank()
    )
  
  ggsave(file.path(output_dir, paste0(prefix, "_qc_comparison.png")), 
         plot = p, width = 10, height = 6, dpi = 300)
  ggsave(file.path(output_dir, paste0(prefix, "_qc_comparison.svg")), 
         plot = p, width = 10, height = 6)
  
  return(p)
}

# 1.7 Integrated QC main function --------------------------------------------
run_integrated_qc <- function(matrix_path,
                              output_dir = "qc_results",
                              project_name = "Integrated_SC",
                              min_genes = 500,
                              max_genes = 6000,
                              min_umi = 1000,
                              max_umi = 30000,
                              max_mt = 20) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  cat("\n===== Starting QC for Integrated Matrix =====")
  cat("\nProject:", project_name)
  cat("\nOutput directory:", output_dir, "\n")
  
  # Step 1: Load integrated matrix
  cat("\n=== Step 1/8: Loading Integrated Matrix ===")
  if (!file.exists(matrix_path)) {
    stop("Matrix file not found: ", matrix_path)
  }
  sc_matrix <- readRDS(matrix_path)
  cat("\nMatrix dimensions:", nrow(sc_matrix), "genes ×", ncol(sc_matrix), "cells")
  
  # Step 2: Create Seurat object and parse metadata
  cat("\n\n=== Step 2/8: Creating Seurat Object with Metadata ===")
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = sc_matrix,
    project = project_name,
    min.cells = 3,
    min.features = 200
  )
  
  cell_ids <- colnames(seurat_obj)
  sample_ids <- sapply(strsplit(cell_ids, "_"), `[`, 1)
  seurat_obj$orig.ident <- sample_ids
  seurat_obj$barcode <- sapply(strsplit(cell_ids, "_"), function(x) paste(x[-1], collapse = "_"))
  
  # Disease group mapping (user-configurable)
  disease_mapping <- c(
    "GSM8239638" = "DPN", "GSM8239639" = "DPN", "GSM8239640" = "DPN",
    "GSM8239642" = "TLA", "GSM8239643" = "TLA", "GSM8239644" = "TLA"
  )
  disease_vector <- disease_mapping[sample_ids]
  names(disease_vector) <- NULL
  seurat_obj$disease <- disease_vector
  
  cat("\nIdentified samples:", toString(unique(sample_ids)))
  cat("\nDisease groups:", toString(unique(disease_vector)))
  
  # Step 3: Calculate mitochondrial percentage
  cat("\n\n=== Step 3/8: Calculating Mitochondrial Percentage ===")
  seurat_obj <- calculate_mt_percentage(seurat_obj)
  
  # Step 4: Initial QC report
  cat("\n\n=== Step 4/8: Generating Initial QC Report ===")
  initial_qc <- generate_qc_report(seurat_obj, min_genes, max_genes, min_umi, max_umi, max_mt)
  write.csv(initial_qc, file.path(output_dir, "initial_qc_report.csv"), row.names = FALSE)
  
  # Step 5: QC visualizations
  cat("\n\n=== Step 5/8: Creating QC Visualizations ===")
  vis_dir <- file.path(output_dir, "visualizations")
  dir.create(vis_dir, showWarnings = FALSE)
  Seurat::Idents(seurat_obj) <- "orig.ident"
  create_qc_violin_plot(seurat_obj, vis_dir, "initial")
  create_qc_scatter_plots(seurat_obj, vis_dir, "initial")
  visualize_qc_report(initial_qc, vis_dir, "initial")
  
  # Step 6: Apply sample-specific filtering
  cat("\n\n=== Step 6/8: Applying Sample-Specific Filters ===")
  filtered_data <- apply_sample_specific_filters(seurat_obj, output_dir)
  seurat_filtered <- filtered_data$filtered_obj
  sample_thresholds <- filtered_data$thresholds
  
  # Step 7: Post-filtering QC report
  cat("\n\n=== Step 7/8: Generating Filtered QC Report ===")
  filtered_qc <- generate_qc_report_with_thresholds(seurat_filtered, sample_thresholds)
  write.csv(filtered_qc, file.path(output_dir, "filtered_qc_report.csv"), row.names = FALSE)
  
  Seurat::Idents(seurat_filtered) <- "orig.ident"
  create_qc_violin_plot(seurat_filtered, vis_dir, "filtered")
  create_qc_scatter_plots(seurat_filtered, vis_dir, "filtered")
  visualize_qc_comparison(initial_qc, filtered_qc, vis_dir, "comparison")
  
  # Step 8: Save results
  cat("\n\n=== Step 8/8: Saving Final Results ===")
  saveRDS(seurat_filtered, file.path(output_dir, "final_filtered_seurat.rds"))
  
  initial_cells <- ncol(seurat_obj)
  filtered_cells <- ncol(seurat_filtered)
  retention_rate <- round(filtered_cells / initial_cells * 100, 1)
  
  summary_data <- data.frame(
    Project = project_name,
    InitialCells = initial_cells,
    FilteredCells = filtered_cells,
    RetentionRate = retention_rate,
    Samples = length(unique(sample_ids))
  )
  write.csv(summary_data, file.path(output_dir, "qc_summary.csv"), row.names = FALSE)
  
  cat("\n\n===== QC Analysis Completed =====")
  cat("\nInitial cells:", initial_cells)
  cat("\nFiltered cells:", filtered_cells)
  cat("\nRetention rate:", retention_rate, "%")
  cat("\nOutput files saved to:", normalizePath(output_dir))
  
  return(list(
    raw_object = seurat_obj,
    filtered_object = seurat_filtered,
    qc_reports = list(initial = initial_qc, filtered = filtered_qc),
    sample_thresholds = sample_thresholds
  ))
}

# ====================================================================
# Part 2: Normalization Functions
# ====================================================================

# 2.1 Mitochondrial percentage (standalone version) --------------------------
calculate_mt_percentage_norm <- function(seurat_obj) {
  mt_identifiers <- c("^MT-", "^mt-", "^MT\\.", "^mt\\.")
  mt_genes <- unique(unlist(lapply(mt_identifiers, function(prefix) {
    grep(prefix, rownames(seurat_obj), value = TRUE, ignore.case = TRUE)
  })))
  
  if (length(mt_genes) == 0) {
    warning("No mitochondrial genes detected! Using zero values.")
    seurat_obj[["percent.mt"]] <- 0
    return(seurat_obj)
  }
  
  counts <- LayerData(seurat_obj, assay = "RNA", layer = "counts")
  mt_counts <- Matrix::colSums(counts[mt_genes, , drop = FALSE])
  total_counts <- Matrix::colSums(counts)
  percent.mt <- (mt_counts / total_counts) * 100
  
  seurat_obj[["percent.mt"]] <- percent.mt
  return(seurat_obj)
}

# 2.2 Identify genes to exclude ---------------------------------------------
identify_excluded_genes <- function(seurat_obj) {
  all_genes <- rownames(seurat_obj)
  
  exclusion_patterns <- c(
    "^MT-",             # mitochondrial genes
    "^IG[HKL]",         # immunoglobulins
    "^HB[^(P)]",        # hemoglobins (except HBP)
    "^MALAT1$|^XIST$|^NEAT1$" # specific lncRNAs
  )
  
  excluded_genes <- unique(unlist(lapply(exclusion_patterns, function(pattern) {
    grep(pattern, all_genes, value = TRUE, ignore.case = TRUE)
  })))
  
  if (file.exists("supplementary_exclude_genes.txt")) {
    supplementary_exclude <- readLines("supplementary_exclude_genes.txt")
    supplementary_exclude <- supplementary_exclude[supplementary_exclude != ""]
    excluded_genes <- unique(c(excluded_genes, supplementary_exclude))
  }
  
  return(excluded_genes[excluded_genes %in% all_genes])
}

# 2.3 Validate normalization quality -----------------------------------------
validate_normalization_quality <- function(seurat_norm, protect_group = NULL) {
  cat("=== Normalization Quality Validation ===\n")
  
  sct_data <- seurat_norm@meta.data
  sct_data$nFeature_SCT <- seurat_norm$nFeature_SCT
  sct_data$nCount_SCT <- seurat_norm$nCount_SCT
  
  umi_correlation <- cor(sct_data$nCount_SCT, sct_data$percent.mt, use = "complete.obs")
  mt_correlation <- cor(sct_data$percent.mt, sct_data$nCount_SCT, use = "complete.obs")
  
  residuals <- tryCatch({
    LayerData(seurat_norm, assay = "SCT", layer = "scale.data")
  }, error = function(e) {
    warning("Cannot access scale.data: ", e$message)
    matrix(0, nrow = nrow(seurat_norm), ncol = ncol(seurat_norm))
  })
  
  residual_mean <- mean(residuals, na.rm = TRUE)
  residual_sd <- sd(residuals, na.rm = TRUE)
  
  var_features <- VariableFeatures(seurat_norm)
  var_feature_count <- length(var_features)
  excluded_in_hvg <- sum(var_features %in% identify_excluded_genes(seurat_norm))
  
  sample_correlation <- NA
  if (length(unique(seurat_norm$orig.ident)) > 1) {
    sample_correlation <- tryCatch({
      sample_expr <- AverageExpression(seurat_norm, assays = "SCT", 
                                       group.by = "orig.ident", layers = "data")$SCT
      cor(sample_expr, method = "pearson")
    }, error = function(e) {
      message("Sample correlation failed: ", e$message)
      return(NA)
    })
  }
  
  harmony_total_var <- NA
  harmony_per_sample_var <- NA
  if ("harmony" %in% names(seurat_norm@reductions)) {
    harmony_dims <- 1:min(20, ncol(Embeddings(seurat_norm, "harmony")))
    harmony_var <- apply(Embeddings(seurat_norm, "harmony")[, harmony_dims], 2, var)
    harmony_total_var <- sum(harmony_var)
    harmony_per_sample_var <- sapply(unique(seurat_norm$orig.ident), function(id) {
      cells <- which(seurat_norm$orig.ident == id)
      sum(apply(Embeddings(seurat_norm, "harmony")[cells, harmony_dims], 2, var))
    })
  }
  
  biological_distance <- NA
  if (!is.null(protect_group) && protect_group %in% colnames(seurat_norm@meta.data)) {
    if ("harmony" %in% names(seurat_norm@reductions)) {
      harmony_embeddings <- Embeddings(seurat_norm, "harmony")[, 1:20]
      group_means <- tryCatch({
        aggregate(harmony_embeddings, 
                  by = list(seurat_norm@meta.data[[protect_group]]), 
                  mean)
      }, error = function(e) {
        message("Biological distance calculation failed: ", e$message)
        return(NULL)
      })
      
      if (!is.null(group_means) && nrow(group_means) > 1) {
        rownames(group_means) <- group_means[,1]
        biological_distance <- dist(group_means[, -1])
      }
    }
  }
  
  report <- list(
    umi_correlation = umi_correlation,
    mt_correlation = mt_correlation,
    residual_mean = residual_mean,
    residual_sd = residual_sd,
    var_feature_count = var_feature_count,
    excluded_in_hvg = excluded_in_hvg,
    sample_correlation = sample_correlation,
    harmony_total_variance = harmony_total_var,
    harmony_per_sample_variance = harmony_per_sample_var,
    biological_distance = biological_distance
  )
  
  cat("UMI-MT Correlation:", round(umi_correlation, 4), "\n")
  cat("MT-UMI Correlation:", round(mt_correlation, 4), "\n")
  cat("Residual Mean:", format(residual_mean, scientific = TRUE, digits = 3), "\n")
  cat("Residual SD:", round(residual_sd, 4), "\n")
  cat("Variable Features:", var_feature_count, "\n")
  cat("Excluded Genes in HVGs:", excluded_in_hvg, "\n")
  
  if (!is.na(harmony_total_var)) {
    cat("Harmony Integration Metrics:\n")
    cat(" - Total Variance in Harmony Space:", round(harmony_total_var, 2), "\n")
  }
  
  if (!is.null(biological_distance)) {
    cat("Biological Group Distances:\n")
    print(biological_distance)
  }
  
  return(report)
}

# 2.4 Generate normalization report ------------------------------------------
generate_normalization_report <- function(validation_report, output_dir, prefix, excluded_genes) {
  report_text <- c(
    "Enhanced Normalization Quality Report",
    "=====================================",
    paste("Date:", Sys.Date()),
    paste("Group:", prefix),
    "",
    "Key Metrics:",
    paste("- UMI-MT Correlation:", round(validation_report$umi_correlation, 4)),
    paste("- MT-UMI Correlation:", round(validation_report$mt_correlation, 4)),
    paste("- Residual Mean:", format(validation_report$residual_mean, scientific = TRUE, digits = 3)),
    paste("- Residual SD:", round(validation_report$residual_sd, 4)),
    paste("- Variable Features:", validation_report$var_feature_count),
    paste("- Excluded Genes in HVGs:", validation_report$excluded_in_hvg),
    paste("- Total Excluded Genes:", length(excluded_genes))
  )
  
  if (!is.na(validation_report$harmony_total_variance)) {
    report_text <- c(
      report_text,
      "",
      "Harmony Integration Metrics:",
      paste("- Total Variance in Harmony Space:", round(validation_report$harmony_total_variance, 2))
    )
  }
  
  if (!is.null(validation_report$biological_distance)) {
    report_text <- c(
      report_text,
      "",
      "Biological Group Distances:",
      capture.output(print(validation_report$biological_distance))
    )
  }
  
  writeLines(report_text, file.path(output_dir, paste0(prefix, "_normalization_report.txt")))
  
  metrics_df <- data.frame(
    Metric = c(
      "UMI_MT_Correlation",
      "MT_UMI_Correlation",
      "Residual_Mean",
      "Residual_SD",
      "Variable_Features",
      "Excluded_Genes_in_HVGs",
      "Total_Excluded_Genes"
    ),
    Value = c(
      validation_report$umi_correlation,
      validation_report$mt_correlation,
      validation_report$residual_mean,
      validation_report$residual_sd,
      validation_report$var_feature_count,
      validation_report$excluded_in_hvg,
      length(excluded_genes)
    )
  )
  
  if (!is.na(validation_report$harmony_total_variance)) {
    metrics_df <- rbind(metrics_df, 
                        data.frame(Metric = "Harmony_Total_Variance", 
                                   Value = validation_report$harmony_total_variance))
  }
  
  write.csv(metrics_df, file.path(output_dir, paste0(prefix, "_normalization_metrics.csv")), row.names = FALSE)
}

# 2.5 Enhanced normalization main function -----------------------------------
perform_enhanced_normalization <- function(filtered_seurat_path,
                                           output_dir = "normalization_results",
                                           group_name = "Unknown",
                                           n_genes = 5000,
                                           harmony_theta = 1.0,
                                           harmony_lambda = 0.1,
                                           conserve.memory = TRUE) {
  
  norm_dir <- file.path(output_dir, group_name)
  dir.create(norm_dir, showWarnings = FALSE, recursive = TRUE)
  cat("\n===== Starting Normalization for:", group_name, "Group =====\n")
  
  # Memory optimization settings
  options(future.globals.maxSize = 50 * 1024^3)
  options(Seurat.memsafe = TRUE)
  
  monitor_memory <- function() {
    if (.Platform$OS.type == "windows") {
      mem_usage <- system("wmic OS get FreePhysicalMemory /Value", intern = TRUE)
      free_mem <- as.numeric(gsub("[^0-9]", "", mem_usage[grep("FreePhysicalMemory", mem_usage)]))
      cat("Free physical memory:", round(free_mem/1e6, 1), "GB\n")
    } else {
      mem_info <- system('free -m | grep Mem', intern = TRUE)
      mem_values <- strsplit(mem_info, "\\s+")[[1]]
      free_mem <- as.numeric(mem_values[4])
      cat("Free memory:", free_mem, "MB\n")
    }
  }
  
  cat("Current memory status: ")
  monitor_memory()
  
  # Step 1: Load filtered Seurat object
  cat("\n=== Step 1/8: Loading Filtered Seurat Object ===\n")
  if (!file.exists(filtered_seurat_path)) {
    stop("Filtered Seurat object not found: ", filtered_seurat_path)
  }
  seurat_obj <- readRDS(filtered_seurat_path)
  
  if (!"disease" %in% colnames(seurat_obj@meta.data)) {
    cat("\nDisease grouping not found. Creating based on sample names...\n")
    disease_mapping <- c(
      "GSM8239638" = "DPN", "GSM8239639" = "DPN", "GSM8239640" = "DPN",
      "GSM8239642" = "TLA", "GSM8239643" = "TLA", "GSM8239644" = "TLA"
    )
    sample_ids <- as.character(seurat_obj$orig.ident)
    disease_group <- disease_mapping[sample_ids]
    disease_group[is.na(disease_group)] <- "Unknown"
    seurat_obj$disease <- disease_group
    cat("Created disease grouping:\n")
    print(table(seurat_obj$disease))
  } else {
    cat("\nExisting disease grouping found:\n")
    print(table(seurat_obj$disease))
  }
  
  cat("Loaded object dimensions:", ncol(seurat_obj), "cells ×", nrow(seurat_obj), "genes\n")
  
  # Step 2: Verify mitochondrial percentage
  cat("\n=== Step 2/8: Verifying Mitochondrial Percentage ===\n")
  if (!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj <- calculate_mt_percentage_norm(seurat_obj)
    cat("Calculated percent.mt for", ncol(seurat_obj), "cells\n")
  }
  
  # Step 3: Create gene whitelist
  cat("\n=== Step 3/8: Creating Gene Whitelist ===\n")
  excluded_genes <- identify_excluded_genes(seurat_obj)
  all_genes <- rownames(seurat_obj)
  whitelist_genes <- setdiff(all_genes, excluded_genes)
  cat("Whitelist genes count:", length(whitelist_genes), "\n")
  
  # Step 4: Stratified SCTransform (gold standard)
  cat("\n=== Step 4/8: Merging and Stratified SCTransform ===\n")
  new_cellnames <- paste0(seurat_obj$orig.ident, "_", colnames(seurat_obj))
  seurat_obj <- RenameCells(seurat_obj, new.names = new_cellnames)
  
  merged_seurat <- SCTransform(
    object = seurat_obj,
    assay = "RNA",
    new.assay.name = "SCT",
    variable.features.n = n_genes,
    vars.to.regress = c("percent.mt"),
    vst.flavor = "v2",
    method = "glmGamPoi",
    conserve.memory = conserve.memory,
    verbose = TRUE,
    residual.features = NULL,
    do.correct.umi = TRUE,
    return.only.var.genes = TRUE,
    batch_var = "orig.ident"
  )
  
  var_features <- VariableFeatures(merged_seurat)
  cat("Selected variable features:", length(var_features), "\n")
  
  rm(seurat_obj)
  gc(full = TRUE)
  cat("Post-normalization memory status: ")
  monitor_memory()
  
  # Step 5: Run PCA
  cat("\n=== Step 5/8: Running PCA ===\n")
  DefaultAssay(merged_seurat) <- "SCT"
  merged_seurat <- RunPCA(
    merged_seurat, 
    npcs = 50, 
    verbose = FALSE,
    reduction.name = "pca",
    reduction.key = "PC_",
    approx = (ncol(merged_seurat) > 10000)
  )
  cat("PCA dimensions:", dim(merged_seurat@reductions$pca), "\n")
  
  # Step 6: Harmony integration (protect biological differences)
  cat("\n=== Step 6/8: Harmony Integration (Protecting Biological Differences) ===\n")
  if ("disease" %in% colnames(merged_seurat@meta.data)) {
    cat("Protecting biological group: disease (DPN vs TLA)\n")
    merged_seurat$disease <- as.factor(merged_seurat$disease)
    merged_seurat <- RunHarmony(
      object = merged_seurat,
      group.by.vars = "orig.ident",
      covariates = "disease",
      reduction = "pca",
      reduction.save = "harmony",
      dims = 1:30,
      theta = harmony_theta,
      lambda = harmony_lambda,
      verbose = TRUE,
      max.iter = 20
    )
  } else {
    cat("Warning: Disease group not found, using standard Harmony integration\n")
    merged_seurat <- RunHarmony(
      object = merged_seurat,
      group.by.vars = "orig.ident",
      reduction = "pca",
      reduction.save = "harmony",
      dims = 1:30,
      theta = harmony_theta,
      lambda = harmony_lambda,
      verbose = TRUE
    )
  }
  cat("Harmony integration complete, dimensions:", dim(merged_seurat@reductions$harmony), "\n")
  
  # Step 7: Validate normalization quality
  cat("\n=== Step 7/8: Validating Normalization Quality ===\n")
  validation_report <- validate_normalization_quality(merged_seurat, "disease")
  print(validation_report)
  
  # Step 8: Save results
  cat("\n=== Step 8/8: Saving Final Results ===\n")
  merged_seurat@commands <- list()
  merged_seurat@graphs <- list()
  
  saveRDS(merged_seurat, file.path(norm_dir, "normalized_seurat.rds"))
  generate_normalization_report(validation_report, norm_dir, group_name, excluded_genes)
  saveRDS(excluded_genes, file.path(norm_dir, "excluded_genes.rds"))
  
  harmony_params <- data.frame(
    theta = harmony_theta,
    lambda = harmony_lambda,
    n_genes = n_genes,
    protected_group = "disease"
  )
  write.csv(harmony_params, file.path(norm_dir, "harmony_parameters.csv"), row.names = FALSE)
  
  cat("\n===== Normalization Summary =====")
  cat("\nGroup:", group_name)
  cat("\nFinal cell count:", ncol(merged_seurat))
  cat("\nVariable features:", length(var_features))
  cat("\nHarmony parameters (theta/lambda):", harmony_theta, "/", harmony_lambda)
  cat("\nProtected biological group: disease (DPN vs TLA)")
  cat("\nOutput directory:", normalizePath(norm_dir))
  cat("\n=================================\n")
  
  return(list(
    normalized_object = merged_seurat,
    validation_report = validation_report,
    excluded_genes = excluded_genes
  ))
}

# ====================================================================
# Part 3: Clustering Functions
# ====================================================================

# 3.1 Generate clustering report ---------------------------------------------
generate_clustering_report <- function(output_dir, 
                                       project_name,
                                       resolution,
                                       cluster_count,
                                       reduction_used,
                                       n_cells,
                                       n_samples) {
  
  report_file <- file.path(output_dir, "clustering_summary_report.txt")
  
  report_content <- paste(
    "===== Clustering Analysis Report =====",
    paste("Project Name:", project_name),
    paste("Date:", Sys.Date()),
    "",
    "===== Key Parameters =====",
    paste("Selected Resolution:", resolution),
    paste("Number of Clusters:", cluster_count),
    paste("Dimensionality Reduction Used:", reduction_used),
    paste("Number of Cells:", n_cells),
    paste("Number of Samples:", n_samples),
    "",
    "===== Quality Assessment =====",
    "This section would include quality metrics if implemented",
    sep = "\n"
  )
  
  writeLines(report_content, report_file)
  cat("Saved clustering report to:", report_file, "\n")
}

# 3.2 Comprehensive clustering main function ---------------------------------
run_comprehensive_clustering <- function(seurat_obj, 
                                         output_dir = "comprehensive_clustering",
                                         project_name = "MultiSample_Analysis",
                                         dims = 1:20,
                                         resolutions = seq(0.1, 1.5, by = 0.1),
                                         selected_resolution = NULL,
                                         algorithm = 2,
                                         clustree_width = 12,
                                         clustree_height = 10,
                                         pt_size = 0.5,
                                         label_size = 4,
                                         # UMAP parameters
                                         umap_n_neighbors = 30,
                                         umap_min_dist = 0.3,
                                         umap_metric = "cosine",
                                         umap_seed = 42,
                                         # t-SNE parameters
                                         tsne_perplexity = 30,
                                         tsne_max_iter = 1000,
                                         tsne_seed = 42,
                                         # Visualization parameters
                                         max_cells_visualization = 15000) {
  
  suppressPackageStartupMessages({
    require(ggplot2)
    require(clustree)
    require(Seurat)
    require(patchwork)
    require(viridis)
    require(RColorBrewer)
    require(scales)
  })
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  cat("\n===== Starting Comprehensive Clustering for", project_name, "=====\n")
  cat("Cells:", ncol(seurat_obj), "| Genes:", nrow(seurat_obj), "\n")
  
  n_samples <- length(unique(seurat_obj$orig.ident))
  cat("Detected", n_samples, "samples\n")
  
  if (!"SCT" %in% names(seurat_obj@assays)) {
    stop("Input Seurat object must be normalized using SCTransform. SCT assay not found.")
  }
  
  required_meta <- c("nCount_SCT", "nFeature_SCT", "percent.mt", "orig.ident")
  missing_meta <- setdiff(required_meta, colnames(seurat_obj@meta.data))
  if (length(missing_meta) > 0) {
    stop("Missing critical metadata columns: ", paste(missing_meta, collapse = ", "))
  }
  
  DefaultAssay(seurat_obj) <- "SCT"
  
  # Select dimensionality reduction method (prefer Harmony)
  cat("\n[1/5] Selecting dimensionality reduction method...\n")
  reduction_use <- "pca"
  reduction_name <- "PCA"
  
  if ("harmony" %in% names(seurat_obj@reductions)) {
    cat("Using Harmony-corrected PCA for clustering\n")
    reduction_use <- "harmony"
    reduction_name <- "Harmony PCA"
    
    available_dims <- ncol(seurat_obj@reductions$harmony)
    if (max(dims) > available_dims) {
      warning("Requested dims (", max(dims), ") exceed available harmony dimensions (", 
              available_dims, "). Using max available.")
      dims <- 1:min(max(dims), available_dims)
    }
  } else if ("pca" %in% names(seurat_obj@reductions)) {
    cat("Using standard PCA for clustering\n")
    reduction_use <- "pca"
    reduction_name <- "PCA"
    
    available_dims <- ncol(seurat_obj@reductions$pca)
    if (max(dims) > available_dims) {
      warning("Requested dims (", max(dims), ") exceed available PCA dimensions (", 
              available_dims, "). Using max available.")
      dims <- 1:min(max(dims), available_dims)
    }
  } else {
    cat("No PCA reduction found. Running PCA...\n")
    seurat_obj <- RunPCA(seurat_obj, npcs = max(dims), verbose = FALSE)
    reduction_use <- "pca"
    reduction_name <- "PCA"
  }
  
  # Visualize dimensionality reduction
  cat("\n[2/5] Visualizing dimensionality reduction...\n")
  reduction_dir <- file.path(output_dir, "reduction_analysis")
  dir.create(reduction_dir, showWarnings = FALSE)
  
  elbow_plot <- ElbowPlot(seurat_obj, ndims = max(dims), reduction = reduction_use) +
    geom_vline(xintercept = max(dims), color = "red", linetype = "dashed", linewidth = 1) +
    labs(title = paste(reduction_name, "Elbow Plot"),
         subtitle = paste("Selected", length(dims), "dimensions (red dashed line)"),
         caption = "Red line indicates maximum dimensions used in downstream analysis") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  ggsave(file.path(reduction_dir, paste0(tolower(reduction_use), "_elbow_plot.png")), 
         elbow_plot, width = 10, height = 8, dpi = 300)
  ggsave(file.path(reduction_dir, paste0(tolower(reduction_use), "_elbow_plot.svg")), 
         elbow_plot, width = 10, height = 8)
  
  # Multi-resolution clustering
  cat("\n[3/5] Performing multi-resolution clustering...\n")
  prefix <- "SCT_snn_res."
  
  if (!"SCT_snn" %in% names(seurat_obj@graphs)) {
    cat("Building nearest neighbor graph...\n")
    seurat_obj <- FindNeighbors(seurat_obj, reduction = reduction_use, dims = dims)
  }
  
  cat("Running clustering at", length(resolutions), "resolutions...\n")
  seurat_obj <- FindClusters(
    seurat_obj,
    resolution = resolutions,
    algorithm = algorithm,
    group.singletons = FALSE,
    verbose = FALSE
  )
  
  # UMAP
  cat("\n[4/5] Running UMAP...\n")
  if (!"umap" %in% names(seurat_obj@reductions)) {
    seurat_obj <- RunUMAP(
      seurat_obj,
      reduction = reduction_use,
      dims = dims,
      n.neighbors = umap_n_neighbors,
      min.dist = umap_min_dist,
      metric = umap_metric,
      seed.use = umap_seed,
      verbose = FALSE
    )
  }
  
  # t-SNE
  cat("\n[5/5] Running t-SNE...\n")
  if (!"tsne" %in% names(seurat_obj@reductions)) {
    seurat_obj <- RunTSNE(
      seurat_obj,
      reduction = reduction_use,
      dims = dims,
      perplexity = tsne_perplexity,
      max_iter = tsne_max_iter,
      seed.use = tsne_seed,
      verbose = FALSE
    )
  }
  
  # Resolution selection
  cat("\n[6/6] Selecting resolution and generating results...\n")
  cluster_cols <- grep(paste0("^", prefix), colnames(seurat_obj@meta.data), value = TRUE)
  cluster_counts <- sapply(cluster_cols, function(col) {
    length(unique(seurat_obj@meta.data[[col]]))
  })
  res_values <- as.numeric(gsub(prefix, "", names(cluster_counts)))
  
  if (!is.null(selected_resolution)) {
    cat("Using user-specified resolution:", selected_resolution, "\n")
    rec_res <- selected_resolution
    if (!rec_res %in% res_values) {
      closest_res <- res_values[which.min(abs(res_values - rec_res))]
      warning("Specified resolution ", rec_res, " not tested. Using closest tested resolution: ", closest_res)
      rec_res <- closest_res
    }
  } else {
    cat("Automatically selecting optimal resolution...\n")
    stable_ranges <- list()
    current_start <- 1
    for (i in 2:length(cluster_counts)) {
      if (cluster_counts[i] != cluster_counts[i - 1]) {
        range_length <- i - 1 - current_start + 1
        if (range_length >= 1) {
          stable_ranges[[length(stable_ranges) + 1]] <- list(
            start_res = res_values[current_start],
            end_res = res_values[i - 1],
            length = range_length,
            clusters = cluster_counts[current_start]
          )
        }
        current_start <- i
      }
    }
    
    if (current_start <= length(cluster_counts)) {
      range_length <- length(cluster_counts) - current_start + 1
      if (range_length >= 1) {
        stable_ranges[[length(stable_ranges) + 1]] <- list(
          start_res = res_values[current_start],
          end_res = res_values[length(res_values)],
          length = range_length,
          clusters = cluster_counts[current_start]
        )
      }
    }
    
    if (length(stable_ranges) > 0) {
      valid_ranges <- stable_ranges[sapply(stable_ranges, function(x) x$length >= 2)]
      if (length(valid_ranges) > 0) {
        longest_idx <- which.max(sapply(valid_ranges, function(x) x$length))
        rec_range <- valid_ranges[[longest_idx]]
        stable_resolutions <- res_values[
          res_values >= rec_range$start_res & 
            res_values <= rec_range$end_res
        ]
        mid_value <- mean(c(rec_range$start_res, rec_range$end_res))
        rec_res <- stable_resolutions[which.min(abs(stable_resolutions - mid_value))]
        cat("Detected stable resolution range: ", 
            rec_range$start_res, "to", rec_range$end_res, "\n")
        cat("Selected resolution: ", rec_res, "\n")
      } else {
        rec_res <- median(resolutions)
        cat("No suitable stable range found. Using median resolution:", rec_res, "\n")
      }
    } else {
      changes <- abs(diff(cluster_counts))
      min_change_idx <- which.min(changes)
      rec_res <- res_values[min_change_idx + 1]
      cat("No stable range found. Using minimal change resolution:", rec_res, "\n")
    }
  }
  
  cat("Final selected resolution:", rec_res, "\n")
  res_col <- paste0(prefix, rec_res)
  Idents(seurat_obj) <- res_col
  seurat_obj$recommended_clusters <- seurat_obj[[res_col]]
  
  # Visualizations
  vis_dir <- file.path(output_dir, "visualizations")
  dir.create(vis_dir, showWarnings = FALSE)
  
  clustree_plot <- clustree(seurat_obj, prefix = prefix) +
    ggtitle("Cluster Stability Across Resolutions") +
    theme_minimal(base_size = 12)
  
  ggsave(file.path(vis_dir, "clustree_plot.png"), 
         clustree_plot, width = clustree_width, height = clustree_height, dpi = 300)
  ggsave(file.path(vis_dir, "clustree_plot.svg"), 
         clustree_plot, width = clustree_width, height = clustree_height)
  
  sample_cluster_plot <- ggplot(seurat_obj@meta.data, 
                                aes(x = recommended_clusters, fill = orig.ident)) +
    geom_bar(position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_viridis_d(option = "viridis") +
    labs(title = "Sample Composition per Cluster",
         x = "Cluster",
         y = "Percentage",
         fill = "Sample") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(vis_dir, "sample_cluster_composition.png"), 
         sample_cluster_plot, width = 12, height = 8, dpi = 300)
  ggsave(file.path(vis_dir, "sample_cluster_composition.svg"), 
         sample_cluster_plot, width = 12, height = 8)
  
  # Generate separate dimension plots per grouping
  create_dim_plots <- function(seurat_obj, reduction, group_name, output_dir, max_cells) {
    if (ncol(seurat_obj) > max_cells) {
      set.seed(123)
      sampled_cells <- sample(colnames(seurat_obj), max_cells)
      seurat_sub <- subset(seurat_obj, cells = sampled_cells)
      message("Subsampled to ", max_cells, " cells for ", reduction, " visualization")
    } else {
      seurat_sub <- seurat_obj
    }
    
    group_vars <- c("orig.ident")
    if ("disease" %in% colnames(seurat_sub@meta.data)) {
      group_vars <- c(group_vars, "disease")
    }
    group_vars <- c(group_vars, "recommended_clusters")
    
    for (group_var in group_vars) {
      if (group_var == "disease") {
        cols <- RColorBrewer::brewer.pal(3, "Set1")
      } else if (group_var == "orig.ident") {
        cols <- RColorBrewer::brewer.pal(8, "Set2")
      } else {
        cols <- viridis_pal(option = "viridis")(length(unique(seurat_sub$recommended_clusters)))
      }
      
      dim_plot <- DimPlot(
        seurat_sub,
        reduction = reduction,
        group.by = group_var,
        pt.size = pt_size,
        label = if (group_var == "recommended_clusters") TRUE else FALSE,
        label.size = if (group_var == "recommended_clusters") label_size else 3,
        repel = if (group_var == "recommended_clusters") TRUE else FALSE,
        cols = cols
      ) +
        labs(title = paste0(toupper(reduction), " - ", group_var),
             subtitle = paste("Colored by", group_var)) +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold", hjust = 0.5),
              legend.position = "bottom")
      
      file_prefix <- paste0(group_name, "_", reduction, "_", group_var)
      ggsave(file.path(output_dir, paste0(file_prefix, ".png")), 
             dim_plot, width = 10, height = 8, dpi = 300)
      ggsave(file.path(output_dir, paste0(file_prefix, ".svg")), 
             dim_plot, width = 10, height = 8)
    }
  }
  
  create_dim_plots(seurat_obj, "umap", project_name, vis_dir, max_cells_visualization)
  create_dim_plots(seurat_obj, "tsne", project_name, vis_dir, max_cells_visualization)
  
  # Harmony integration plots (if used)
  if (reduction_use == "harmony") {
    cat("Generating Harmony integration plots...\n")
    if (ncol(seurat_obj) > max_cells_visualization) {
      set.seed(123)
      sampled_cells <- sample(colnames(seurat_obj), max_cells_visualization)
      seurat_sub <- subset(seurat_obj, cells = sampled_cells)
      message("Subsampled to ", max_cells_visualization, " cells for Harmony visualization")
    } else {
      seurat_sub <- seurat_obj
    }
    
    p1 <- DimPlot(seurat_sub,
                  reduction = "harmony",
                  group.by = "orig.ident",
                  pt.size = pt_size,
                  cols = RColorBrewer::brewer.pal(8, "Set2")) +
      labs(title = "Harmony Integration by Sample",
           subtitle = "Colored by sample origin") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            legend.position = "bottom")
    
    ggsave(file.path(vis_dir, paste0(project_name, "_harmony_samples.png")), 
           p1, width = 10, height = 8, dpi = 300)
    ggsave(file.path(vis_dir, paste0(project_name, "_harmony_samples.svg")), 
           p1, width = 10, height = 8)
    
    if ("disease" %in% colnames(seurat_sub@meta.data)) {
      p2 <- DimPlot(seurat_sub,
                    reduction = "harmony",
                    group.by = "disease",
                    pt.size = pt_size,
                    cols = RColorBrewer::brewer.pal(3, "Set1")) +
        labs(title = "Harmony Integration by Disease",
             subtitle = "Colored by biological group") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold", hjust = 0.5),
              legend.position = "bottom")
      
      ggsave(file.path(vis_dir, paste0(project_name, "_harmony_disease.png")), 
             p2, width = 10, height = 8, dpi = 300)
      ggsave(file.path(vis_dir, paste0(project_name, "_harmony_disease.svg")), 
             p2, width = 10, height = 8)
    }
  }
  
  # Save results
  saveRDS(seurat_obj, file.path(output_dir, "comprehensive_clustered_seurat.rds"))
  
  res_summary <- data.frame(
    Resolution = res_values,
    Clusters = cluster_counts
  )
  write.csv(res_summary, file.path(output_dir, "resolution_summary.csv"), row.names = FALSE)
  write.csv(seurat_obj@meta.data, file.path(output_dir, "cell_metadata.csv"))
  
  generate_clustering_report(
    output_dir = output_dir,
    project_name = project_name,
    resolution = rec_res,
    cluster_count = length(unique(Idents(seurat_obj))),
    reduction_used = reduction_use,
    n_cells = ncol(seurat_obj),
    n_samples = n_samples
  )
  
  cat("\n===== Comprehensive Clustering Completed =====\n")
  cat("Results saved to:", normalizePath(output_dir), "\n")
  cat("Selected resolution:", rec_res, "with", 
      length(unique(Idents(seurat_obj))), "clusters\n")
  cat("Dimensionality reduction used:", reduction_use, "\n")
  cat("UMAP parameters: n.neighbors =", umap_n_neighbors, 
      "| min.dist =", umap_min_dist, "| metric =", umap_metric, "\n")
  cat("t-SNE parameters: perplexity =", tsne_perplexity, 
      "| max_iter =", tsne_max_iter, "\n")
  
  return(seurat_obj)
}

# ====================================================================
# Execute the Pipeline (example usage)
# ====================================================================

# --- Quality Control ---
integrated_qc_results <- run_integrated_qc(
  matrix_path = "integrated_sc_matrix.rds",
  output_dir = "qc_results/Integrated_v3",
  project_name = "GEO_Integrated_SC"
)

# --- Normalization ---
normalized_results <- perform_enhanced_normalization(
  filtered_seurat_path = "qc_results/Integrated_v3/final_filtered_seurat.rds",
  output_dir = "30harmony_Results",
  group_name = "30harmony",
  n_genes = 5000,
  harmony_theta = 1.0,
  harmony_lambda = 0.1,
  conserve.memory = TRUE
)

# --- Clustering ---
clustered_obj <- run_comprehensive_clustering(
  seurat_obj = normalized_results$normalized_object,
  output_dir = "0.9clustering",
  project_name = "Custom_Analysis",
  dims = 1:20,
  resolutions = seq(0.1, 1.6, by = 0.1),
  selected_resolution = 0.9,
  umap_n_neighbors = 40,
  umap_min_dist = 0.3,
  tsne_perplexity = 60,
  max_cells_visualization = 15000
)

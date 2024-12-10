# Set the working directory
setwd("/home/akkey/Akshay/Phd/RNAseq/Genesis_paper")
getwd()

# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(Rtsne)  # For t-SNE
library(umap)
library(pheatmap)

# Load the data
Counts <- read.delim("GSE71318_HS_count.csv", header = TRUE, sep = ",")

# Inspect the data structure
str(Counts)

# Aggregate counts by 'Geneid'
Counts <- aggregate(. ~ Geneid, data = Counts, FUN = sum)

# Set 'Geneid' as row names and remove the column
rownames(Counts) <- Counts$Geneid
Counts <- Counts[, -1]

# Filter genes with non-zero counts in at least 50% of samples
min_samples_with_counts <- ceiling(ncol(Counts) * 0.5)
Counts <- Counts[rowSums(Counts > 0) >= min_samples_with_counts, ]

# Further filter to remove low-expression genes
Counts <- Counts[rowSums(Counts) > 50, ]

# Display first few rows of the filtered data
head(Counts)

# Calculate zero values per sample
zero_counts_per_sample <- colSums(Counts == 0)
print(zero_counts_per_sample)

# Define experimental conditions
condition <- factor(c("8cell", "8cell", "8cell", 
                      "Morula", "Morula", "Morula", 
                      "Blastocyst", "Blastocyst", "Blastocyst"))

# Create coldata for DESeq2
coldata <- data.frame(row.names = colnames(Counts), condition)
print(coldata)

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = Counts, colData = coldata, design = ~ condition)

# Run DESeq pipeline
dds <- DESeq(dds)

# Perform variance stabilizing transformation
vsdata <- vst(dds, blind = FALSE)

# --- Overall Analysis: PCA, Dispersion, and Heatmap ---
overall_output_folder <- "Overall_Analysis"
if (!dir.exists(overall_output_folder)) {
  dir.create(overall_output_folder)
}

# PCA Plot
pca_data <- plotPCA(vsdata, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4, alpha = 0.8) +
  labs(
    title = "PCA Plot: Overall",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  ) +
  theme_minimal() +
  theme(legend.title = element_text(face = "bold"), legend.position = "top")

ggsave(file.path(overall_output_folder, "PCA_Plot_Overall.png"), plot = pca_plot, width = 12, height = 10, dpi = 600)
print(pca_plot)

# Dispersion Estimates
png(file.path(overall_output_folder, "Dispersion_Estimates_Overall.png"), width = 2400, height = 1800, res = 600)
plotDispEsts(dds, main = "Dispersion Estimates: Overall")
dev.off()

# Heatmap
count_data <- assay(vsdata)
stage_labels <- coldata$condition

calculate_pcc <- function(data, stage) {
  stage_data <- data[, stage_labels == stage]
  cor(stage_data, method = "pearson")
}

stages <- unique(stage_labels)
pcc_list <- lapply(stages, function(stage) {
  list(stage = stage, pcc = calculate_pcc(count_data, stage))
})

combined_pcc <- sapply(pcc_list, function(stage_data) colMeans(stage_data$pcc, na.rm = TRUE))
rownames(combined_pcc) <- stages
colnames(combined_pcc) <- stages

png(file.path(overall_output_folder, "Heatmap_Overall.png"), width = 2400, height = 1800, res = 600)
pheatmap(
  combined_pcc,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = FALSE,
  main = "Heatmap: Overall",
  color = colorRampPalette(c("blue", "white", "red"))(50)
)
dev.off()

# --- Compare Multiple Conditions ---
comparisons <- list(
  c("8cell", "Morula"),
  c("Morula", "Blastocyst"),
  c("8cell", "Blastocyst")
)

# Loop through comparisons
for (comparison in comparisons) {
  # Extract conditions for the current comparison
  condition1 <- comparison[1]
  condition2 <- comparison[2]
  comparison_name <- paste(condition1, "vs", condition2, sep = "_")
  
  # Create a folder for the current comparison
  output_folder <- paste0("DEGs_", comparison_name)
  if (!dir.exists(output_folder)) {
    dir.create(output_folder)
  }
  
  # Perform differential expression analysis
  res <- results(dds, contrast = c("condition", condition2, condition1))
  res_clean <- na.omit(res)
  DEGs <- res_clean[res_clean$padj < 0.05, ]
  
  # Save DEGs
  write.csv(DEGs, file.path(output_folder, paste0("DEGs_", comparison_name, ".csv")), row.names = TRUE)
  
  # Identify upregulated and downregulated genes
  upregulated_genes <- DEGs[DEGs$log2FoldChange > 1, ]
  downregulated_genes <- DEGs[DEGs$log2FoldChange < -1, ]
  
  # Save upregulated and downregulated genes
  write.csv(upregulated_genes, file.path(output_folder, paste0("Upregulated_Genes_", comparison_name, ".csv")), row.names = TRUE)
  write.csv(downregulated_genes, file.path(output_folder, paste0("Downregulated_Genes_", comparison_name, ".csv")), row.names = TRUE)
  
  # Summary
  cat("Comparison:", comparison_name, "\n")
  cat("Total DEGs:", nrow(DEGs), "\n")
  cat("Upregulated genes:", nrow(upregulated_genes), "\n")
  cat("Downregulated genes:", nrow(downregulated_genes), "\n")
  
  # --- Volcano Plot ---
  res_clean$Significance <- ifelse(res_clean$padj < 0.05 & abs(res_clean$log2FoldChange) > 1,
                                   ifelse(res_clean$log2FoldChange > 0, "Upregulated", "Downregulated"),
                                   "Not Significant")
  
  volcano_plot <- ggplot(res_clean, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
    theme_minimal() +
    labs(
      title = paste("Volcano Plot:", condition1, "vs", condition2),
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    ) +
    theme(legend.title = element_blank(), legend.position = "top")
  
  ggsave(file.path(output_folder, paste0("Volcano_Plot_", comparison_name, ".png")), plot = volcano_plot, width = 8, height = 6, dpi = 600)
  print(volcano_plot)
}

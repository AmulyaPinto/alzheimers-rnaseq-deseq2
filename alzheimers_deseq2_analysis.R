# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differential expression analysis: Alzheimer's disease (AD) vs Healthy brain
# Dataset: GEO GSE203206 (occipital cortex, Brodmann area 17)
# Methods: DESeq2, volcano plots (ggplot2 + ggrepel), heatmap (pheatmap)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# SECTION 1 - PACKAGES


# Install once (remove the # to run)
#install.packages("BiocManager")
#BiocManager::install(c("DESeq2", "org.Hs.eg.db"))
#install.packages(c("ggplot2", "ggrepel", "pheatmap", "RColorBrewer"))

# Run the installed packages (needs to be done every time you want to run the following analysis code)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(AnnotationDbi)
library(org.Hs.eg.db)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# SECTION 2 - LOADING THE DATA


# Set the working directory to the folder holding the two downloaded files
setwd("E:/R Language/Alzheimers bulk RNA-seq study")

# Quick check
getwd()

# Create folders for saved tables and figures
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Loading the count matrix
counts <- read.table("GSE203206_Subramaniam.ADRC_brain.counts.tsv",
                     header = TRUE,
                     row.names = 1,
                     sep = "\t")

# DESeq2 requires whole-number counts, so convert to integers
counts <- round(counts)

# Quick check
head(counts)
dim(counts)

# Loading the metadata
metadata <- read.table("GSE203206_Subramaniam.ADRC_brain.metadata.tsv",
                       header = TRUE,
                       row.names = 1,
                       sep = "\t")

# Quick check
head(metadata)
table(metadata$PATHDX)  # samples per group: Alzheimer's and Healthy

# Metadata rows must be in the same order as the count matrix columns
metadata <- metadata[colnames(counts), ]
stopifnot(all(colnames(counts) == rownames(metadata)))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# SECTION 3 - RUN DESeq2


dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ PATHDX
)

# Healthy is the reference group, so a positive log2 fold change means higher expression in Alzheimer's disease
dds$PATHDX <- relevel(dds$PATHDX,
                      ref = "Healthy")

# Filter low counts and keep genes with at least 10 reads in total across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Run the DESeq2 model
dds <- DESeq(dds)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# SECTION 4 - EXTRACT AND SUMMARISE RESULTS


res <- as.data.frame(results(dds))
res$gene <- rownames(res)
res <- res[!is.na(res$padj), ]

# Total genes after filtering
nrow(res)

# Genes with adjusted p-value < 0.05
sum(res$padj < 0.05)

# The most significant gene
res[which.min(res$padj), ]

# Quick check
head(res[order(res$padj), ])

# Direction column (up/down): adjusted p-value < 0.05 and |log2 fold change| > 1
res$direction <- ifelse(
  res$padj < 0.05 & res$log2FoldChange > 1, "up",
  ifelse(res$padj < 0.05 & res$log2FoldChange < -1, "down", "ns")
)

# Quick check
table(res$direction)

# Key numbers for the README
summary_table <- data.frame(
  measure = c("Samples: Alzheimers",
              "Samples: Healthy",
              "Genes tested (after filtering)",
              "Significant genes (adj. p < 0.05)",
              "Upregulated in AD (adj. p < 0.05, log2FC > 1)",
              "Downregulated in AD (adj. p < 0.05, log2FC < -1)"),
  value = c(sum(metadata$PATHDX == "Alzheimers"),
            sum(metadata$PATHDX == "Healthy"),
            nrow(res),
            sum(res$padj < 0.05),
            sum(res$direction == "up"),
            sum(res$direction == "down"))
)

summary_table
write.csv(summary_table, "results/summary_counts.csv", row.names = FALSE)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# SECTION 5 - GENE SYMBOLS AND TOP 10 GENES


# Converting Ensembl IDs to gene symbols (unmapped IDs give NA)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys = res$gene,
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")

# Quick check
head(res[, c("gene", "symbol")])

# Ten most significant up/down genes, that will be used for labelling the volcano plot
sig_genes <- res[res$direction != "ns", ]
top10 <- head(sig_genes[order(sig_genes$padj), ], 10)

# Quick check
top10

# Save results table
write.csv(res, "results/deseq2_all_results.csv", row.names =  FALSE)
write.csv(sig_genes[order(sig_genes$padj), ],
          "results/significant_genes.csv", row.names = FALSE)
write.csv(top10, "results/top10_genes.csv", row.names = FALSE)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# SECTION 6 - VOLCANO PLOTS


direction_colours <- c(
  "up" = "forestgreen",
  "down" = "red",
  "ns" = "orange")

volcano_plot <- ggplot(res, aes(
  x = log2FoldChange,
  y = -log10(padj),
  fill = direction)) +
  geom_point(shape = 21, color = "black", alpha = 0.6, size = 2) +
  scale_fill_manual(values = direction_colours, name = "Direction") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             linewidth = 0.5, color = "grey25") +
  geom_hline(yintercept = 1.3, linetype = "dashed",
             linewidth = 0.5, color = "navy") +
  theme_classic() +
  labs(title = "Significance of Changes in Genes",
       subtitle = "Based on amount of regulation",
       x = "Log2 Fold Change",
       y = "Adjusted P-value (-log10)") +
  xlim(-10, 26) +  # axis limits wide enough to show all genes
  ylim(0, 28)


volcano_plot

# Saving plot in the figures folder for README
ggsave("figures/volcano_plot.png", plot = volcano_plot,
       width = 8, height = 6, dpi = 300)

# Same plot with the top 5 genes labelled
volcano_labelled <- volcano_plot +
  geom_text_repel(data = head(top10, 5),
                  aes(label = symbol, color = direction),
                  size = 3,
                  max.overlaps = Inf) +
  scale_color_manual(values = direction_colours, guide = "none")


volcano_labelled

# Saving plot in the figures folder for README
ggsave("figures/volcano_labelled.png", plot = volcano_labelled,
       width = 8, height = 6, dpi = 300)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# SECTION 7 - HEATMAP OF THE 30 MOST VARIABLE GENES


# Variance stabilising transformation (normalises counts for visualisation)
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)

top_var <- head(order(rowVars(mat), decreasing = TRUE), 30)
mat30 <- mat[top_var, ]

# Z-score each gene across samples
mat_scaled <- t(scale(t(mat30)))

# Quick check
dim(mat_scaled)

pheatmap(mat_scaled)

# Adding a diagnosis annotation bar and a clearer colour scale (drop = FALSE keeps the results as a data frame, which pheatmap requires)
col_ann <- as.data.frame(metadata[, "PATHDX", drop = FALSE])

heatmap_annotated <- pheatmap(mat_scaled,
                              annotation_col = col_ann,
                              annotation_colors = list(PATHDX = c("Alzheimers" = "firebrick",
                                                                  "Healthy" = "steelblue")),
                              color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
                              show_colnames = FALSE,  # sample names are too crowded, not showing them for neater heatmap
                              main = "Top 30 Variable Genes\nAD vs Healthy")

# Saving heatmap in the figures folder for README
png("figures/heatmap_top30.png",
    width = 2400, height = 2000, res = 300)
grid::grid.newpage()
grid::grid.draw(heatmap_annotated$gtable)
dev.off()

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# SECTION 8 - CHECK THE FOUR GENES WITH EXTREME FOLD CHANGES


# Normalised counts in each group for the four genes with log2 fold change above 20
png("figures/extreme_genes_counts.png",
    width = 2000, height = 1600, res = 250)
par(mfrow = c(2,2))
for (g in c("ENSG00000204389", "ENSG00000230230",
            "ENSG00000231247", "ENSG00000275246")) {
  plotCounts(dds, gene = g, intgroup = "PATHDX",
             main = paste0(res[g, "symbol"], " (", g, ")"))
}
dev.off()

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Save package versions
writeLines(capture.output(sessionInfo()), "results/session_info.txt")

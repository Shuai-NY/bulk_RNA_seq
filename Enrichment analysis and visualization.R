getwd()
setwd('C:/R_temp/RNA-Seq/shP2RY1')
data<- read.table(file ='count.csv', sep = ',',row.names = 1,header = T)
head(data)
#Option 2: Aggregate duplicates by summing counts
library(dplyr)

# Create a column for Ensembl IDs without version
data$ensembl_id <- sub("\\..*", "", rownames(data))

# Group by Ensembl ID and sum counts
data_combined <- data %>%
  group_by(ensembl_id) %>%
  summarise(across(where(is.numeric), sum))

# Convert back to data frame
data_combined <- as.data.frame(data_combined)
rownames(data_combined) <- data_combined$ensembl_id
data_combined$ensembl_id <- NULL
data <- data_combined

library(biomaRt)

# Connect to Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Get mapping for protein-coding genes only
mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                 filters = "ensembl_gene_id",
                 values = rownames(data),
                 mart = mart)

# Filter for protein-coding
mapping <- mapping[mapping$gene_biotype == "protein_coding", ]

# Merge with your count matrix
count_matrix_filtered <- data[rownames(data) %in% mapping$ensembl_gene_id, ]


# ----  Merge mapping with counts ----
count_matrix_filtered$ensembl_gene_id <- rownames(count_matrix_filtered)
merged <- merge(mapping, count_matrix_filtered, by = "ensembl_gene_id")


# Replace Ensembl IDs with gene symbols and sum duplicates ----
final_matrix <- merged %>%
  select(-ensembl_gene_id, -gene_biotype) %>%
  group_by(hgnc_symbol) %>%
  summarise(across(everything(), sum))


final_matrix <- as.data.frame(final_matrix)  # Convert tibble to data frame
rownames(final_matrix) <- final_matrix$hgnc_symbol
final_matrix$hgnc_symbol <- NULL


# ---- 11. Save final matrix ----
write.csv(final_matrix, "counts_protein_coding.csv")


library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)

#---- 2. Load raw count matrix ----
counts <- read.csv("counts_protein_coding.csv", row.names = 1)
counts<-counts[,-2]
counts<-counts[,-5]
counts['P2RY1',]

# ---- 3. Create sample metadata ----
coldata <- data.frame(
  row.names = colnames(counts),
  condition = c("shP2RY1", "shP2RY1", "shScramble", "shScramble")  # Example
)

# Set factor levels so shScramble is the reference
coldata$condition <- factor(coldata$condition, levels = c("shScramble", "shP2RY1"))

# ---- 4. Create DESeq2 dataset ----
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

# ---- 5. Filter low counts ----
dds <- dds[rowSums(counts(dds)) > 10, ]

# ---- 6. Run DESeq2 ----
dds <- DESeq(dds)

# ---- 7. Results ----
res <- results(dds, contrast = c("condition", "shScramble", "shP2RY1"))
resLFC <- lfcShrink(dds, coef = 2, type = "apeglm")

# Save results
write.csv(as.data.frame(resLFC), "DESeq2_results.csv")

# ---- 8. Normalized counts ----
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, "dseq2normalized_counts.csv")

# ---- 9. Volcano Plot with Highlighted Significant Genes and Labels ----
res_df <- as.data.frame(resLFC)
res_df$padj[is.na(res_df$padj)] <- 1

# Create status column for color coding
res_df$status <- "nosig"
res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange >= 1] <- "up"
res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange <= -1] <- "down"
write.csv(res_df, "dseq2res_df.csv")

# Count up/downregulated genes
up_count <- sum(res_df$status == "up")
down_count <- sum(res_df$status == "down")
cat("Upregulated genes:", up_count, "\n")
cat("Downregulated genes:", down_count, "\n")


# Select top 10 upregulated and top 10 downregulated genes
sig_genes <- res_df[res_df$padj < 0.05, ]
top_up <- head(sig_genes[order(-sig_genes$log2FoldChange), ], 10)
top_down <- head(sig_genes[order(sig_genes$log2FoldChange), ], 10)
label_genes <- c(rownames(top_up), rownames(top_down))


# Add P2RY1 if present
extra_label <- if ("P2RY1" %in% rownames(res_df)) "P2RY1" else NULL
label_genes <- unique(c(rownames(top_up), rownames(top_down), extra_label))

# Volcano plot




ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "nosig" = "grey")) +
  labs(title = paste("Volcano Plot\nUp:", up_count, " | Down:", down_count),
       x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = subset(res_df, rownames(res_df) %in% label_genes),
                  aes(label = rownames(subset(res_df, rownames(res_df) %in% label_genes))),
                  size = 3, max.overlaps = Inf) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  xlim(-5, 5)  # Fix x-axis range





# ---- 10. PCA Plot ----
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  labs(title = "PCA Plot")

# ---- 11. Heatmap with Top Up/Downregulated Genes ----
selected_genes <- label_genes
mat <- assay(vsd)[selected_genes, ]
mat <- mat - rowMeans(mat)

pheatmap(mat,
         annotation_col = coldata,
         show_rownames = TRUE,
         cluster_cols = TRUE,
         main = "Top 10 Upregulated & Downregulated Genes Heatmap")



# Get normalized counts
norm_counts <- counts(dds, normalized = TRUE)

# Identify all significant up/downregulated genes
sig_genes <- rownames(res_df[res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1, ])

# Subset normalized counts
mat <- norm_counts[sig_genes, ]

# Center rows for better visualization
mat <- mat - rowMeans(mat)

# Heatmap
library(pheatmap)
pheatmap(mat,
         annotation_col = coldata,
         show_rownames = FALSE,  # Hide if too many genes
         show_colnames = FALSE,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         scale = "row",  # Standardizes each gene
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         treeheight_row = 0,   # Default ~50
         treeheight_col = 0,
         main = paste("Heatmap of", length(sig_genes), "Up/Downregulated Genes"))

# Use DESeq2 results
res_df <- as.data.frame(resLFC)

# Remove NA values
res_df <- res_df[!is.na(res_df$padj), ]

# Create ranked list: gene names and log2FoldChange
gene_list <- res_df$log2FoldChange
names(gene_list) <- rownames(res_df)

# Sort in decreasing order
gene_list <- sort(gene_list, decreasing = TRUE)



library(clusterProfiler)
library(org.Hs.eg.db)


# Convert gene symbols to Entrez IDs
gene_ids <- bitr(names(gene_list), fromType = "SYMBOL",
                 toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Remove duplicates: keep first mapping
gene_ids <- gene_ids[!duplicated(gene_ids$SYMBOL), ]

# Filter gene_list to match unique SYMBOLs
gene_list <- gene_list[names(gene_list) %in% gene_ids$SYMBOL]

# Reorder gene_ids to match gene_list order
gene_ids <- gene_ids[match(names(gene_list), gene_ids$SYMBOL), ]

# Assign Entrez IDs as names
names(gene_list) <- gene_ids$ENTREZID

# Load MSigDB Hallmark gene sets (downloaded separately)
gmt_file <- "h.all.v2025.1.Hs.entrez.gmt"  # Example path
hallmark_sets <- read.gmt(gmt_file)
go_gmt_file <- "c5.go.v2025.1.Hs.entrez.gmt"  # Example path
go_sets <- read.gmt(go_gmt_file)
kegg_gmt_file <- "c2.cp.kegg_medicus.v2025.1.Hs.entrez.gmt"  # Example path
kegg_sets <- read.gmt(kegg_gmt_file)
# Run GSEA
gsea_res <- GSEA(gene_list, TERM2GENE = hallmark_sets, pvalueCutoff = 0.05)
gsea_res_go <- GSEA(gene_list, TERM2GENE = go_sets, pvalueCutoff = 0.05)
gsea_res_kegg <- GSEA(gene_list, TERM2GENE = kegg_sets, pvalueCutoff = 0.05)

# View results
head(gsea_res)


write.csv(gsea_res, "hallmark_results.csv", row.names = FALSE)
write.csv(gsea_res_go, "go_results.csv", row.names = FALSE)
write.csv(gsea_res_kegg, "kegg_results.csv", row.names = FALSE)



library(fgsea)



dotplot(gsea_res_go, showCategory = 10, x = "NES") +
  ggtitle("Top 10 Enriched GO Terms") +
  theme_minimal()
dotplot(gsea_res, showCategory = 10, x = "NES") +
  ggtitle("Top 10 Enriched Hallmark Terms") +
  theme_minimal()

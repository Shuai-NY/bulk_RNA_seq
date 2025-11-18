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



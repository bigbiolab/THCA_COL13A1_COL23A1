# Differential Gene Expression Analysis of GSE165724

# 01. Load required packages
library(tidyverse)
library(DESeq2)
library(biomaRt)
library(conflicted)
library(readr)
library(dplyr)

# 02. Load count data
count_data <- read_tsv("data/GSE165724/GSE165724_counts.tsv")

glimpse(count_data)


# 03. Load metadata
meta_data <- read.csv("data/GSE165724/GSE165724_metadata.csv")

glimpse(meta_data)

# 04. Create a matrix & add gene ids as row names
count_data <- GSE165724 |> select(2:73)  |> as.matrix()
rownames(count_data) <- GSE165724$GeneID

# 05. Match metadata with count data
meta_data <- meta_data |> 
  filter(Samples %in% colnames(count_data)) |> 
  arrange(match(Samples, colnames(count_data)))

# 07. Prepare Sample information
colData <- data.frame( condition = as.factor(meta_data$Conditions), 
                       row.names = colnames(count_data))

#$ 08. Create DESeq2 data set object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = colData,
                              design = ~ condition)

# 09. filter any counts less than 10
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# 10. Run the below to see the DESeq2 object elements
dds

# 12. Run DESeq2
dds <- DESeq(dds)

# 13. Analyze the data set and compile result
res <- results(dds)
res$Gene_ID <- rownames(res)

# 14. Export results to a CSV file
dir.create("outputs/GSE165724")

write.csv(res, "outputs/GSE165724/GSE165724.csv",
          row.names = FALSE)

# ------------------------------------------------------------------------------
# 16. load results 
GSE165724 <- read.csv("outputs/GSE165724/GSE165724.csv")

# 17. Set up biomaRt for annotation
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 18. Fetch gene annotations for Entrez IDs
annotations <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol", "description"),
  filters = "entrezgene_id",
  values = GSE165724$Gene_ID,
  mart = mart
)

# 19. Merge annotations directly with the combined results
annotated_results <- merge(GSE165724, 
                           annotations, by.x = "Gene_ID", 
                           by.y = "entrezgene_id", 
                           all.x = TRUE)

# 20. Rename columns for clarity
colnames(annotated_results)[colnames(annotated_results) == "hgnc_symbol"] <- "Gene_Symbol"
colnames(annotated_results)[colnames(annotated_results) == "description"] <- "Gene_Description"

# 21. Remove NA's in the gene symbol column
annotated_results <- annotated_results |> drop_na(Gene_Symbol) |>
  distinct()

# 22. Export results to a CSV file
write.csv(annotated_results, "outputs/GSE165724/GSE165724_annotated.csv",
          row.names = FALSE)
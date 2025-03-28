# Differential Gene Expression Analysis of GSE153659

# 01. Load required packages
library(tidyverse)
library(DESeq2)
library(biomaRt)
library(conflicted)
library(readr)
library(dplyr)

# 02. Load count data
GSE153659 <- read_tsv("data/GSE153659/GSE153659_counts.tsv")

glimpse(GSE153659)

# 03. Load metadata
meta_data <- read.csv("data/GSE153659/GSE153659_metadata.csv")

glimpse(meta_data)

# 04. Create a matrix & add gene ids as row names
count_data <- GSE153659 |> select(2:32)  |> as.matrix()
rownames(count_data) <- GSE153659$GeneID

# 05. Match metadata with count data
meta_data <- meta_data |> 
  filter(Samples %in% colnames(count_data)) |> 
  arrange(match(Samples, colnames(count_data)))

# 07. Prepare Sample information
colData <- data.frame( condition = as.factor(meta_data$Conditions), 
                       row.names = colnames(count_data))

# 08. Create DESeq2 data set object
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
dir.create("outputs/GSE153659")

write.csv(res, "outputs/GSE153659/GSE153659.csv",
          row.names = FALSE)

# ------------------------------------------------------------------------------

# 15. Load results
GSE153659 <- read.csv("outputs/GSE153659/GSE153659.csv")

glimpse(GSE153659)

# 16. Set up biomaRt for annotation
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 17. Annotate the results
annotations <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol", "description"),
  filters = "entrezgene_id",
  values = GSE153659$Gene_ID,
  mart = mart
)


# 19. Merge annotations directly with the combined results
annotated_results <- merge(GSE153659, 
                           annotations, by.x = "Gene_ID", 
                           by.y = "entrezgene_id", 
                           all.x = TRUE)

# 20. Rename columns for clarity
colnames(annotated_results)[colnames(annotated_results) == "hgnc_symbol"] <- "Gene_Symbol"
colnames(annotated_results)[colnames(annotated_results) == "description"] <- "Gene_Description"

# 21. Remove NA's in the gene symbol column
annotated_results <- annotated_results |> drop_na(Gene_Symbol) |>
  distinct()

# 22. Export the annotated results
write.csv(annotated_results, "outputs/GSE153659/GSE153659_annotated.csv",
          row.names = FALSE)

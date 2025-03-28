# Volcano plot for visualization of differentially expressed genes

# Load required pacakges
library(EnhancedVolcano)
library(tidyverse)

# Load data 
GSE153659 <- read.csv("outputs/GSE153659/GSE153659_annotated.csv")
GSE165724 <- read.csv("outputs/GSE165724/GSE165724_annotated.csv")
GSE150899 <- read.csv("outputs/GSE150899/GSE150899_annotated.csv")

# Combine data
degenes <- rbind(GSE153659, GSE165724, GSE150899)

# 22. Export the results
dir.create("outputs/Combined Genes")

write.csv(degenes, "outputs/Combined Genes/degenes.csv",
          row.names = FALSE)

# Generate volcano plot
volcano <- EnhancedVolcano(degenes,
                           lab = as.character(degenes$Gene_Symbol),
                           x = 'log2FoldChange',
                           y = 'pvalue',
                           pointSize = 2.0,
                           labSize = 5.0,
                           col=c('#636363', '#9ecae1', '#a1d99b', '#e6550d'),
                           colAlpha = 1,
                           legendLabels=c('NS','Log2FC','p-value',
                                          'p-Value & Log2FC'),
                           legendPosition = 'top',
                           legendLabSize = 16,
                           legendIconSize = 5.0,
                           title = "",  
                           titleLabSize = 16,
                           subtitle = "",
                           subtitleLabSize = 18,
                           pCutoff = 10e-6,
                           FCcutoff = 2,
                           cutoffLineType = "dashed",
                           border = "partial")

# Export volcano plot
ggsave(filename = ("figures/volcano_plot.png"), 
       plot = volcano, 
       width = 8.5, 
       height = 8, 
       dpi = 300)
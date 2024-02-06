# ===== Differential Gene Expression Analysis =====

# Human Airway Smooth Muscle (HASM) Transcriptome Changes in Response to Asthma Medications
# GSE52778

# Identify transcriptional changes due to treatment with dexamethasone


# ----- Load packages
library(DESeq2)
library(airway)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(ggplot2)
library(EnhancedVolcano)
library(clusterProfiler)

# ----- Load data
data(airway)
countsData <- as.data.frame(assay(airway))
head(countsData)
# 63677 genes as rows and 8 samples as columns

sample_metadata <- as.data.frame(colData(airway))
sample_metadata <- sample_metadata[,c(2,3)]
sample_metadata$dex <- gsub('trt', 'treated', sample_metadata$dex)
sample_metadata$dex <- gsub('untrt', 'untreated', sample_metadata$dex)
colnames(sample_metadata) <- c('cellLine', 'dexamethasone')

# Check the row names in sample_metadata match the column names in countsData
all(row.names(sample_metadata) %in% colnames(countsData))
all(row.names(sample_metadata) == colnames(countsData))

# ----- Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = countsData,
                              colData = sample_metadata,
                              design = ~ dexamethasone)
dds

# ----- Prefiltering
# Remove genes with low counts across samples (< 10)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds
# We keep 22369 genes

# Set untreated condition as control (reference)
dds$dexamethasone <- relevel(dds$dexamethasone, ref = 'untreated')
dds <- DESeq(dds)
# ----- Data Exploration

# Variance-Stabilization Transformation
dds_norm <- vst(dds)
norm_counts <- as.data.frame(assay(dds_norm))
plotPCA(dds_norm, intgroup = 'dexamethasone')

plotDispEsts(dds)

# ----- Differential Expression Analysis
res <- results(dds)
res.df <- as.data.frame(res)
head(res.df)
res.df <- res.df[complete.cases(res.df), ]

# ----- Volcano plot

# Add a column with gene symbols
res.df$symbol <- mapIds(org.Hs.eg.db,
                        keys=row.names(res.df),
                        keytype='ENSEMBL',
                        column='SYMBOL')
head(res.df)

EnhancedVolcano(res.df,
                x = 'log2FoldChange',
                y = 'padj',
                lab = res.df$symbol,
                pCutoff = 1e-10,
                FCcutoff = 1,
                pointSize = 3.0)

known_responsive_genes <- c('DUSP1', 'FKBP5', 'KLF15', 'TSC22D3', 'PER1')
EnhancedVolcano(res.df,
                x = 'log2FoldChange',
                y = 'padj',
                lab = res.df$symbol,
                selectLab = known_responsive_genes,
                pCutoff = 1e-10,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')


data_to_plot <- as.data.frame(t(norm_counts[c('ENSG00000120129', 'ENSG00000152583'),]))
data_to_plot$dexamethasone <- sample_metadata$dexamethasone

# DUSP1
ggplot(data_to_plot, aes(x=dexamethasone, y=ENSG00000120129)) + 
  geom_boxplot(fill='#A4A4A4', color="black") +
  theme_classic() + ylab('DUSP1')

# SPARCL1
ggplot(data_to_plot, aes(x=dexamethasone, y=ENSG00000152583)) + 
  geom_boxplot(fill='#A4A4A4', color="black") +
  theme_classic() + ylab('SPARCL1')

# ----- Results Exploration

# Top 10 differentially expressed genes
head(res.df[order(res.df$padj), ], 10)

# Filter up and downregulated genes that are significantly differentially expressed
sigs <- res.df[abs(res.df$log2FoldChange) > 1 & res.df$padj < 1e-10, ]
sigs <- sigs[order(sigs$log2FoldChange, decreasing = TRUE), ]
head(sigs)

mat <- counts(dds, normalized=T)[row.names(sigs), ]

# Compute the Z-score for each gene
mat.scaled <- t(apply(mat, 1, scale))
colnames(mat.scaled) <- colnames(mat)

# Heatmap of all the differentially expressed genes
ha <- HeatmapAnnotation(condition = sample_metadata$dexamethasone,
                        col = list(condition = c('untreated' = 'green', 'treated' = 'purple')))

Heatmap(mat.scaled,
        top_annotation = ha,
        cluster_rows = T,
        show_row_dend = F,
        show_row_names = F,
        name = 'Z-score')

# Detailed heatmap of top 15 up and downregulated genes
num_genes <- 15
keep_rows <- c(1:num_genes, (nrow(sigs)-num_genes+1):nrow(sigs))


lfc_values <- as.matrix(sigs$log2FoldChange[keep_rows])
colnames(lfc_values) <- 'log2FC'
lfc_col <- colorRamp2(c(min(lfc_values),0,max(lfc_values)), c('blue', 'white', 'red'))


dh1 <- Heatmap(mat.scaled[keep_rows, ],
              top_annotation = ha,
              cluster_rows = F,
              name = 'Z-score')
dh2 <- Heatmap(lfc_values,
               row_labels = sigs$symbol[keep_rows],
               cluster_rows = F,
               name = 'log2FC',
               col = lfc_col,
               cell_fun = function(j, i, x, y, w, h, f) {
                 grid.text(round(lfc_values[i, j], 2), x, y)
               })

hm <- dh1 + dh2

hm

# ----- Over Representation Analysis
upregulated_genes <- row.names(sigs)[sigs$log2FoldChange > 1 & sigs$padj < 1e-10]
upreg_enrich_BP <- enrichGO(gene = upregulated_genes,
                         OrgDb = org.Hs.eg.db,
                         ont = 'BP',
                         pvalueCutoff = 0.05,
                         keyType = 'ENSEMBL')

dotplot(upreg_enrich_BP)
cnetplot(upreg_enrich_BP)



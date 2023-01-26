# this script will perform differential gene expression analysis using DESeq2

# set working directory
setwd("~/Library/CloudStorage/OneDrive-NewcastleUniversity/Molecular Microbiology/MMB8052")

# installing and loading packages 
install.packages('BiocManager')
library(BiocManager)
install(c('tximport', 'DESeq2', 'biomaRt', 'pheatmap', 'RColorBrewer'))
install.packages(
  pkgs = "DESeqAnalysis",
  repos = c(
    "https://r.acidgenomics.com",
    BiocManager::repositories()
  ),
  dependencies = TRUE
)
BiocManager::install("geneplotter")
BiocManager::install("vidger")
library(tximport)
library(DESeq2)
library(biomaRt)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(DESeqAnalysis)
library(vidger)
library(geneplotter)

# installation problems on mac - run code on unsuccessful packages shown below
## install('[name of package]', type='binary', force=TRUE)

# step 1: preparing count data
sample_table = read_csv('https://raw.githubusercontent.com/sjcockell/mmb8052/main/practicals/practical_08/data/sample_table.csv')
files = pull(sample_table, Run)
files = paste0('counts/', files, '/quant.sf')
names(files) = pull(sample_table, Run)
gene_map = read_csv('https://github.com/sjcockell/mmb8052/raw/main/practicals/practical_08/extdata/gene_map.csv')
txi = tximport(files, 
                 type='salmon',
                 tx2gene=gene_map,
                 ignoreTxVersion=TRUE)

# step 2: creating a DESeq dataset 
dds = DESeqDataSetFromTximport(txi, colData = sample_table, design = ~ Group)
dds = DESeq(dds)

# plotting dispersion estimates
plotDispEsts(dds)

# colours used for figures
cols = RColorBrewer::brewer.pal(9, "PuBu")
# MA plot - RColorBrewer::brewer.pal(3, "Dark2")
# heatmap -  RColorBrewer::brewer.pal(3,"Paired")
# volcano - RColorBrewer::brewer.pal(3, "Set2")

# data quality control: PCA and heatmap

# principle component analysis 
rld = rlog(dds)
pcaData <- plotPCA(rld, intgroup=c('Group'), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2)) + 
      geom_point(size=3, aes(col=Group, fill=Group)) +
      scale_color_manual(values=c("#A6CEE3", "#1F78B4", "#82DF8A")) +
      geom_text_repel(aes(label=name), direction="both", nudge_y=1.1, 
                      point.padding = 0.6, box.padding=0.55,
                      min.segment.length = unit(0.2, 'lines'), size=2.5) + 
      ggtitle("Principle Component Analysis for GSE116538") + 
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance"))

# heatmap
sample_distance = dist(t(assay(rld)), method='euclidian')
sample_distance_matrix = as.matrix(sample_distance)
heatmap_annotation = data.frame(group=colData(dds)[,c('Group')], 
                                row.names=rownames(colData(dds)))
pheatmap(sample_distance_matrix,
         clustering_distance_rows=sample_distance,
         clustering_distance_cols=sample_distance,
         cutree_rows = 3,
         cutree_cols = 3,
         color = cols,
         annotation_col = heatmap_annotation,
         annotation_colors = list(group=c(Allo24h="#A6CEE3", 
                                          Allo2h="#1F78B4", Naive="#82DF8A")))

# step 3: running DESeq function to perform differentiation analysis

# creating results table for the three comparisons
results_table1 = results(dds, contrast= c('Group', 'Allo24h', 'Naive'))
results_table2 = results(dds, contrast= c('Group', 'Allo2h', 'Naive'))
results_table3 = results(dds, contrast= c('Group', 'Allo24h', 'Allo2h'))


results_tibble1 = as_tibble(results_table1, rownames='ensembl_gene_id')
results_tibble2 = as_tibble(results_table2, rownames='ensembl_gene_id')
results_tibble3 = as_tibble(results_table3, rownames='ensembl_gene_id')

filtered_results1 = filter(results_tibble1, complete.cases(results_tibble1))
filtered_results2 = filter(results_tibble2, complete.cases(results_tibble2))
filtered_results3 = filter(results_tibble3, complete.cases(results_tibble3))

# creating annotated dataframe with all attributes
ensembl108 = useEnsembl(biomart="ensembl", version=108)
ensembl108 = useDataset("mmusculus_gene_ensembl", mart=ensembl108)

annotation = getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 
                                'start_position', 'end_position', 
                                'strand', 'gene_biotype', 'external_gene_name', 
                                'description'), 
                   filters = 'ensembl_gene_id', 
                   values = filtered_results1$ensembl_gene_id, 
                   mart = ensembl108)
annot_results1 = left_join(filtered_results1, annotation)
annot_results1 = arrange(annot_results1, padj)

annotation = getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 
                                'start_position', 'end_position', 
                                'strand', 'gene_biotype', 'external_gene_name', 
                                'description'), 
                   filters = 'ensembl_gene_id', 
                   values = filtered_results2$ensembl_gene_id, 
                   mart = ensembl108)
annot_results2 = left_join(filtered_results2, annotation)
annot_results2 = arrange(annot_results2, padj)

annotation = getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 
                                'start_position', 'end_position', 
                                'strand', 'gene_biotype', 'external_gene_name', 
                                'description'), 
                   filters = 'ensembl_gene_id', 
                   values = filtered_results3$ensembl_gene_id, 
                   mart = ensembl108)
annot_results3 = left_join(filtered_results3, annotation)
annot_results3 = arrange(annot_results3, padj)


# creating column for neither, pval and both based on pval and lfc
annot_results1 = mutate(annot_results1, logPval= -log10(padj))
annot_results1$Significance <- 'Neither'
annot_results1$Significance[annot_results1$padj < 0.05] <- 'Pval only'
annot_results1$Significance[annot_results1$log2FoldChange < -1 & 
                              annot_results1$padj < 0.05] <- 'Both'
annot_results1$Significance[annot_results1$log2FoldChange > 1 & 
                              annot_results1$padj < 0.05] <- 'Both'

annot_results2 = mutate(annot_results2, logPval= -log10(padj))
annot_results2$Significance <- 'Neither'
annot_results2$Significance[annot_results2$padj < 0.05] <- 'Pval only'
annot_results2$Significance[annot_results2$log2FoldChange < -1 & 
                              annot_results2$padj < 0.05] <- 'Both'
annot_results2$Significance[annot_results2$log2FoldChange > 1 & 
                              annot_results2$padj < 0.05] <- 'Both'

annot_results3 = mutate(annot_results3, logPval= -log10(padj))
annot_results3$Significance <- 'Neither'
annot_results3$Significance[annot_results3$padj < 0.05] <- 'Pval only'
annot_results3$Significance[annot_results3$log2FoldChange < -1 & 
                              annot_results3$padj < 0.05] <- 'Both'
annot_results3$Significance[annot_results3$log2FoldChange > 1 & 
                              annot_results3$padj < 0.05] <- 'Both'

# creating column containing the top 10 genes: used for labelling on volcano plot
top10degs1 = head(arrange(annot_results1, abs(log2FoldChange) < 1 & padj), 10)
annot_results1$delabel = ifelse(annot_results1$external_gene_name %in% 
                                  top10degs1$external_gene_name, 
                               annot_results1$external_gene_name, NA)

top10degs2 = head(arrange(annot_results2, abs(log2FoldChange) < 1 & padj), 10)
annot_results2$delabel = ifelse(annot_results2$external_gene_name %in% 
                                  top10degs2$external_gene_name, 
                                annot_results2$external_gene_name, NA)

top10degs3 = head(arrange(annot_results3, abs(log2FoldChange) < 1 & padj), 10)
annot_results3$delabel = ifelse(annot_results3$external_gene_name %in% 
                                  top10degs3$external_gene_name, 
                                annot_results3$external_gene_name, NA)


# theme for volcano plot
theme_set(theme_classic(base_size = 20) + 
            theme(
              axis.title.y = element_text(face = "bold", 
                                          margin = margin(0,20,0,0), 
                                          size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", 
                                          margin = margin(20,0,0,0), 
                                          size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

# visualising data with volcano plot
ggplot(annot_results1, aes(x=log2FoldChange, 
                        y=logPval, col = Significance, label = delabel,
                        name = external_gene_name)) + 
       geom_vline(xintercept = c(-1, 1), 
                  col = 'black', alpha = 0.5, linetype = 'dashed') +
       geom_hline(yintercept =-log10(0.05), 
                  col = 'black', alpha =0.5, linetype = 'dashed') +
       geom_point(size = 1, alpha = 0.3) + 
       scale_color_manual(values = c("#FC8D62", "#66C2A5", "#8DA0CB"), 
                     labels = c("Both", "Neither", "Pval only")) + 
       coord_cartesian(ylim = c(0, 200), xlim = c(-20, 20)) + 
       scale_x_continuous(breaks = seq(-20, 20, 4)) +
       ggtitle('Allo24h vs Naive') +
       geom_text_repel(aes(label=delabel), direction="both", nudge_y=1.1, 
                  point.padding = 0.6, box.padding=0.55,
                  min.segment.length = unit(0.2, 'lines'), size=2.5)

ggplot(annot_results2, aes(x=log2FoldChange, 
                           y=logPval, col = Significance, label = delabel,
                           name = external_gene_name)) + 
  geom_vline(xintercept = c(-1, 1), 
             col = 'black', alpha = 0.5, linetype = 'dashed') +
  geom_hline(yintercept =-log10(0.05), 
             col = 'black', alpha =0.5, linetype = 'dashed') +
  geom_point(size = 1, alpha = 0.3) + 
  scale_color_manual(values = c("#FC8D62", "#66C2A5", "#8DA0CB"), 
                     labels = c("Both", "Neither", "Pval only")) + 
  coord_cartesian(ylim = c(0, 90), xlim = c(-16, 16)) + 
  scale_x_continuous(breaks = seq(-16, 16, 4)) +
  ggtitle('Allo2h vs Naive') +
  geom_text_repel(aes(label=delabel), direction="both", nudge_y=1.1, 
                  point.padding = 0.6, box.padding=0.55,
                  min.segment.length = unit(0.2, 'lines'), size=2.5)

ggplot(annot_results3, aes(x=log2FoldChange, 
                           y=logPval, col = Significance, label = delabel,
                           name = external_gene_name)) + 
  geom_vline(xintercept = c(-1, 1), 
             col = 'black', alpha = 0.5, linetype = 'dashed') +
  geom_hline(yintercept =-log10(0.05), 
             col = 'black', alpha =0.5, linetype = 'dashed') +
  geom_point(size = 1, alpha = 0.3) + 
  scale_color_manual(values = c("#FC8D62", "#66C2A5", "#8DA0CB"), 
                     labels = c("Both", "Neither", "Pval only")) + 
  coord_cartesian(ylim = c(0, 220), xlim = c(-20, 20)) + 
  scale_x_continuous(breaks = seq(-20, 20, 4)) +
  ggtitle('Allo24h vs Allo2h') +
  geom_text_repel(aes(label=delabel), direction="both", nudge_y=1.1, 
                  point.padding = 0.6, box.padding=0.55,
                  min.segment.length = unit(0.2, 'lines'), size=2.5)

# visualising data with MA plot
DESeqAnalysis::plotMA(
  results_table1,
  direction = c("both", "up", "down"),
  alphaThreshold = 0.05,
  lfcThreshold = 1,
  pointColor = c(downregulated = "#FC8D62", upregulated = "#66C2A5", 
                 nonsignificant = "#8DA0CB"),
  pointSize = 2L,
  pointAlpha = 0.3,
  labels = list(title = "Allo24h vs Naive", 
  subtitle = "n = 1510; up (green): 808, down (orange): 702, 
  P-value = 0.05")
)

DESeqAnalysis::plotMA(
  results_table2,
  direction = c("both", "up", "down"),
  alphaThreshold = 0.05,
  lfcThreshold = 1,
  pointColor = c(downregulated = "#FC8D62", upregulated = "#66C2A5", 
                 nonsignificant = "#8DA0CB"),
  pointSize = 2L,
  pointAlpha = 0.3,
  labels = list(title = "Allo2h vs Naive", 
  subtitle = "n = 1268; up (green): 788, down (orange): 480, 
  P-value = 0.05")
)

DESeqAnalysis::plotMA(
  results_table3,
  direction = c("both", "up", "down"),
  alphaThreshold = 0.05,
  lfcThreshold = 1,
  pointColor = c(downregulated = "#FC8D62", upregulated = "#66C2A5", 
                 nonsignificant = "#8DA0CB"),
  pointSize = 2L,
  pointAlpha = 0.3,
  labels = list(title = "Allo24h vs Allo2h", 
  subtitle = "n = 2495; up (green): 1072, down (orange): 1423, 
  P-value = 0.05")
)

# vidger four-way plot 
vsFourWay(x = "Allo2h", y = "Allo24h", control = "Naive", lfc = 1, data = dds, 
          d.factor = 'Group', 
          type = c('deseq'), padj = 0.05)
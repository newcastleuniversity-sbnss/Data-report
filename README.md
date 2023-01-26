# data-report

# author
Lewis Isaac Chan // source code for MMB8052 Bioinformatics data report
180030622 
26/01/23

# process
One script containing the entire process. Outline of sections: 
- Installation of data
- Preparing count data 
- Creating a DESeq dataset
- Data quality control
- DESeq function to perform differential analysis
- Generating volcano plot
- Generating MA plot
- Generating four-way plot

# data
Information about the dataset: [https://doi.org/10.1111/ajt.15751](https://doi.org/10.1111/ajt.15751)

Downloading the RNA-Seq data for GSE116538: [https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116583/soft/GSE116583_family.soft.gz](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116583/soft/GSE116583_family.soft.gz) 

Downloading the count data used in this report: [https://github.com/sjcockell/mmb8052/raw/main/practicals/practical_08/results/counts.zip](https://github.com/sjcockell/mmb8052/raw/main/practicals/practical_08/results/counts.zip)

# software
Software
- R (v4.2.2)
- Bioconductor (v3.16)
- tximport (v1.26.1)
- DESeq2 (v1.38.3)
- biomaRt (v2.54.0)
- pheatmap (v1.0.12)
- RColorBrewer (v1.1-3)
- ggplot2 (v3.4.0)
- ggrepel (v0.9.2)
- DESeqAnalysis (0.6.6)
- geneplotter (v1.76.0)
- Vidger (1.18.0)

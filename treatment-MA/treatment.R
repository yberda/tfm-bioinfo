# load libraries
library("readr")
library("dplyr")

library("tximeta")
library("biomaRt")
library("DESeq2")
library("limma")
library("edgeR")
library("meta")
library("metafor")
library("TFEA.ChIP")
library("clusterProfiler")
library("enrichplot")

library("ggplot2")
library("ggfortify")
library("RColorBrewer")
library("pheatmap")
library("corrplot")
library("ggVennDiagram")
library("egg")
library("gridExtra")
library("grid")
library("ggrepel")
library("plotly")


setwd("./tfm-bioinfo/treatment-MA")

# load salmon output (gse objects)
load("./salmon-output-treatment.RData")

# in case of running the code, load the data above and skip the data loading 
# lines ("load transcriptomes" section) up to the definition of gse objects 
# (excluding the cell lines and the treatment definition)


# ----------------------------- plot functions ---------------------------------

# plot log(CPM) density before and after filtering
plot.filt <- function(y, dge, ylim = c(0, 0.65)) {
    lcpm.plot <- cpm(y, log = TRUE)
    lcpm.filt.plot <- cpm(dge, log = TRUE)
    col <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Set1"))
    lcpm.cutoff <- log2(10/(mean(y$samples$lib.size) * 1e-06) + 
                                         2/(median(y$samples$lib.size) * 1e-06))
    par(mfrow = c(1, 2))
    plot(density(lcpm.plot[, 1]), main = "raw counts", lwd = 2, col = col[1], 
    						   ylim = ylim, xlab = "log(cpm)")
    abline(v = lcpm.cutoff, lty = 2)
    for (i in 2:dim(y)[2]) {
        lines(density(lcpm.plot[, i]), col = col[i], lwd = 2)
    }
    plot(density(lcpm.filt.plot[, 1]), main = "filtered counts", lwd = 2, 
        col = col[1], ylim = ylim, xlab = "log(cpm)")
    abline(v = lcpm.cutoff, lty = 2)
    for (i in 2:dim(y)[2]) {
        lines(density(lcpm.filt.plot[, i]), col = col[i], lwd = 2)
    }
    par(mfrow = c(1, 1))
}


# volcano plots using limma results
volcanoplot.colors.limma <- function(genes, adjpval = 0.05, logfc = 1.5, 
						        title = "Volcano plot") {
    # add a column of NAs
    genes$diffexpressed <- "NO"
    # if log2Foldchange > logfc and pvalue < adjpval, set as 'UP'
    genes$diffexpressed[genes$logFC >= logfc & 
    					       genes$adj.P.Val < adjpval] <- "UP"
    # if log2Foldchange < -logfc and pvalue < adjpval, set as 'DOWN'
    genes$diffexpressed[genes$logFC <= -logfc & 
    					     genes$adj.P.Val < adjpval] <- "DOWN"

    # plot and color the points with 'diffexpressed'
    p <- ggplot(data = genes, aes(x = logFC, y = -log10(adj.P.Val), 
                                                         col = diffexpressed)) +
        geom_point() + theme_minimal()

    # Add lines
    p2 <- p + geom_vline(xintercept = c(-logfc, logfc), col = "red") + 
                           geom_hline(yintercept = -log10(adjpval), col = "red")

    ## Change point color

    # create a named vector: the values are the colors to be used, the names
    # are the categories they will be assigned to:
    mycolors <- c("blue", "red", "black")
    names(mycolors) <- c("DOWN", "UP", "NO")
    p3 <- p2 + scale_colour_manual(values = mycolors)
    p3 + xlab("log(Fold Change)") + ylab("log(Adj p-value)") + 
        guides(color = "none", size = "none") + labs(title = title) + 
        theme(plot.title = element_text(hjust = 0.5, size = 20), 
              axis.title = element_text(size = 17), 
              axis.text = element_text(size = 17))
    return(p3)
}



# This script contains the code used to do the following:

# 1. Data preparation: data are loaded and processed to be used for 
#    metaanalyses.
# 2. Metaanalyses using 6-48h treatment data: this section includes the 
#    metaanalyses, the exploration of the results (including 
#    the clusterProfiler, GSEA, TFEA.ChIP, actomyosin network and cell cycle 
#    analysis) and a subgroup and metaregression analysis.
# 3. Metaanalyses using 24-48h treatment data: this section includes the 
#    metaanalyses and the exploration of the results (including the GSEA, 
#    TFEA.ChIP and actomyosin network analysis)



################################################################################
#                            1. DATA PREPARATION                               #
################################################################################

# ----------------------------- load transcriptomes ----------------------------

## obenauf dataset ----

# load data
files.obenauf.lines <- read.table("../output-salmon/files-treat/obenauf.txt")
files.obenauf.lines <- file.path(unlist(files.obenauf.lines))
all(file.exists(files.obenauf.lines))
tdata.obenauf.lines <- as.data.frame(files.obenauf.lines)
colnames(tdata.obenauf.lines) <- "files"
cell.lines.obenauf.lines <- c(rep("A375", 9), rep("Colo800", 4), 
                              rep("UACC62", 4))
tdata.obenauf.lines$cell <- factor(cell.lines.obenauf.lines)
treat.obenauf.lines <- relevel(factor(c("P", "P", "P", "BRAFi6h", "BRAFi6h", 
        "BRAFi6h", "BRAFi48h", "BRAFi48h", "BRAFi48h", "P", "P", "BRAFi48h", 
        "BRAFi48h", "P", "P", "BRAFi48h", "BRAFi48h"), 
    levels = c("P", "BRAFi6h", "BRAFi48h")), ref = "P")
tdata.obenauf.lines$treatment <- treat.obenauf.lines
tdata.obenauf.lines$names <- paste("Sample", 1:17, sep = "")

# reference transcriptome checksum recognized by tximeta
se.obenauf.lines <- tximeta(tdata.obenauf.lines)  
# summarize quantifications to gene-level
gse.obenauf.lines <- summarizeToGene(se.obenauf.lines, 
                                       countsFromAbundance = "lengthScaledTPM")  
y.obenauf.lines <- makeDGEList(gse.obenauf.lines)

# filter low or null counts
dim(y.obenauf.lines)
group.obenauf.lines <- paste(cell.lines.obenauf.lines, 
                                                 treat.obenauf.lines, sep = ".")
design.obenauf.lines <- model.matrix(~0 + group.obenauf.lines)
colnames(design.obenauf.lines) <- c("A375.BRAFi48h", "A375.BRAFi6h", "A375.P", 
    "Colo800.BRAFi48h", "Colo800.P", "UACC62.BRAFi48h", "UACC62.P")
keep.obenauf.lines <- filterByExpr(y.obenauf.lines, design.obenauf.lines)
dge.obenauf.lines <- y.obenauf.lines[keep.obenauf.lines, , 
                                                         keep.lib.sizes = FALSE]
dim(dge.obenauf.lines)[1]/dim(y.obenauf.lines)[1] * 100  

# plots
plot.filt(y.obenauf.lines, dge.obenauf.lines, c(0, 0.7))
plotMDS(dge.obenauf.lines, col = as.numeric(as.factor(treat.obenauf.lines)))

# differential expression analysis
contrast.obenauf.lines <- makeContrasts(
                             A375_6h_vs_cont= A375.BRAFi48h - A375.P,
                             A375_48h_vs_cont = A375.BRAFi6h - A375.P, 
                             Colo800_48h_vs_cont = Colo800.BRAFi48h - Colo800.P, 
                             UACC62_48h_vs_cont = UACC62.BRAFi48h - UACC62.P, 
                             levels = design.obenauf.lines)

v.obenauf.lines <- voom(dge.obenauf.lines, design.obenauf.lines, plot = TRUE)

fit.obenauf.lines <- lmFit(v.obenauf.lines, design.obenauf.lines)
fit.obenauf.lines <- contrasts.fit(fit.obenauf.lines, contrast.obenauf.lines)
fit.obenauf.lines <- eBayes(fit.obenauf.lines, trend = FALSE)

# A375 6h treatment genes
genes.obenauf.lines.A375.6h <- topTable(fit.obenauf.lines, 
    coef = "A375_6h_vs_cont", adjust.method = "BH", sort.by = "p", 
    number = dim(fit.obenauf.lines)[1])
sum(genes.obenauf.lines.A375.6h$adj.P.Val < 0.05)
sum(genes.obenauf.lines.A375.6h$adj.P.Val < 0.001)

# A375 48h treatment
genes.obenauf.lines.A375.48h <- topTable(fit.obenauf.lines, 
                                         coef = "A375_48h_vs_cont",
                                         adjust.method = "BH", sort.by = "p", 
                                         number = dim(fit.obenauf.lines)[1])
sum(genes.obenauf.lines.A375.48h$adj.P.Val < 0.05)
sum(genes.obenauf.lines.A375.48h$adj.P.Val < 0.001)

# Colo800 48h treatment
genes.obenauf.lines.colo.48h <- topTable(fit.obenauf.lines, 
                                         coef = "Colo800_48h_vs_cont", 
                                         adjust.method = "BH", sort.by = "p", 
                                         number = dim(fit.obenauf.lines)[1])
sum(genes.obenauf.lines.colo.48h$adj.P.Val < 0.05)
sum(genes.obenauf.lines.colo.48h$adj.P.Val < 0.001)

# UACC62 48h treatment
genes.obenauf.lines.uacc.48h <- topTable(fit.obenauf.lines, 
                                         coef = "UACC62_48h_vs_cont",
                                         adjust.method = "BH", sort.by = "p", 
                                         number = dim(fit.obenauf.lines)[1])
sum(genes.obenauf.lines.colo.48h$adj.P.Val < 0.05)
sum(genes.obenauf.lines.colo.48h$adj.P.Val < 0.001)

# volcano plots
volcanoplot.colors.limma(genes = genes.obenauf.lines.A375.6h, adjpval = 0.05,
    logfc = 1.5, title = "Obenauf, A375 treated 6h")
volcanoplot.colors.limma(genes = genes.obenauf.lines.A375.48h, adjpval = 0.05, 
    logfc = 1.5, title = "Obenauf, A375 treated 48h")
volcanoplot.colors.limma(genes = genes.obenauf.lines.colo.48h, adjpval = 0.05, 
   logfc = 1.5,  title = "Obenauf, Colo800 treated 48h")
volcanoplot.colors.limma(genes = genes.obenauf.lines.uacc.48h, adjpval = 0.05, 
    logfc = 1.5, title = "Obenauf, UACC62 treated 48h")


## fallahi dataset ----

# load data
files.fallahi.lines <- read.table("../output-salmon/files-treat/fallahi.txt")
files.fallahi.lines <- file.path(unlist(files.fallahi.lines))
all(file.exists(files.fallahi.lines))
tdata.fallahi.lines <- as.data.frame(files.fallahi.lines)
colnames(tdata.fallahi.lines) <- "files"
cell.lines.fallahi.lines <- c(rep("Colo858", 6), rep("MMACSF", 6))
tdata.fallahi.lines$cell <- factor(cell.lines.fallahi.lines)
treat.fallahi.lines <- relevel(factor(c("P", "P", "BRAFi24h", "BRAFi24h", 
    "BRAFi48h", "BRAFi48h", "P", "P", "BRAFi24h", "BRAFi24h", "BRAFi48h", 
    "BRAFi48h"), levels = c("P", "BRAFi24h", "BRAFi48h")), ref = "P")
tdata.fallahi.lines$treatment <- treat.fallahi.lines
tdata.fallahi.lines$names <- paste("Sample", 1:12, sep = "")

# reference transcriptome checksum recognized by tximeta
se.fallahi.lines <- tximeta(tdata.fallahi.lines)

# summarize quantifications to gene-level
gse.fallahi.lines <- summarizeToGene(se.fallahi.lines, 
                                     countsFromAbundance = "lengthScaledTPM")  
y.fallahi.lines <- makeDGEList(gse.fallahi.lines)

# filter low or null counts
dim(y.fallahi.lines)
group.fallahi.lines <- paste(cell.lines.fallahi.lines, treat.fallahi.lines, 
                             sep = ".")
design.fallahi.lines <- model.matrix(~0 + group.fallahi.lines)
colnames(design.fallahi.lines) <- c("Colo858.BRAFi24h", "Colo858.BRAFi48h", 
    "Colo858.P", "MMACSF.BRAFi24h", "MMACSF.BRAFi48h", "MMACSF.P")
keep.fallahi.lines <- filterByExpr(y.fallahi.lines, design.fallahi.lines)
dge.fallahi.lines <- y.fallahi.lines[keep.fallahi.lines, , 
                                     keep.lib.sizes = FALSE]
dim(dge.fallahi.lines)[1]/dim(y.fallahi.lines)[1] * 100 

# plots
plot.filt(y.fallahi.lines, dge.fallahi.lines, c(0, 1.7))
plotMDS(dge.fallahi.lines, col = as.numeric(as.factor(treat.fallahi.lines)))

# differential expression analysis
contrast.fallahi.lines <- makeContrasts(
                             Colo858_24h_vs_cont = Colo858.BRAFi24h - Colo858.P, 
                             Colo858_48h_vs_cont = Colo858.BRAFi48h - Colo858.P, 
                             MMACSF_24h_vs_cont = MMACSF.BRAFi24h - MMACSF.P, 
                             MMACSF_48h_vs_cont = MMACSF.BRAFi48h - MMACSF.P, 
                             levels = design.fallahi.lines)

v.fallahi.lines <- voom(dge.fallahi.lines, design.fallahi.lines, plot = TRUE)

fit.fallahi.lines <- lmFit(v.fallahi.lines, design.fallahi.lines)
fit.fallahi.lines <- contrasts.fit(fit.fallahi.lines, contrast.fallahi.lines)
fit.fallahi.lines <- eBayes(fit.fallahi.lines, trend = FALSE)

# Colo858 24h treatment genes
genes.fallahi.lines.colo.24h <- topTable(fit.fallahi.lines, 
                                         coef = "Colo858_24h_vs_cont", 
                                         adjust.method = "BH", sort.by = "p", 
                                         number = dim(fit.fallahi.lines)[1])
sum(genes.fallahi.lines.colo.24h$adj.P.Val < 0.05)
sum(genes.fallahi.lines.colo.24h$adj.P.Val < 0.001)

# Colo858 48h treatment genes
genes.fallahi.lines.colo.48h <- topTable(fit.fallahi.lines, 
                                         coef = "Colo858_48h_vs_cont",
                                         adjust.method = "BH", sort.by = "p", 
                                         number = dim(fit.fallahi.lines)[1])
sum(genes.fallahi.lines.colo.48h$adj.P.Val < 0.05)
sum(genes.fallahi.lines.colo.48h$adj.P.Val < 0.001)

# MMACSF 24h treatment genes
genes.fallahi.lines.mmacsf.24h <- topTable(fit.fallahi.lines, 
                                           coef = "MMACSF_24h_vs_cont", 
                                           adjust.method = "BH", sort.by = "p", 
                                           number = dim(fit.fallahi.lines)[1])
sum(genes.fallahi.lines.mmacsf.24h$adj.P.Val < 0.05)
sum(genes.fallahi.lines.mmacsf.24h$adj.P.Val < 0.001)

# MMACSF 48h treatment genes
genes.fallahi.lines.mmacsf.48h <- topTable(fit.fallahi.lines, 
                                           coef = "MMACSF_48h_vs_cont",
                                           adjust.method = "BH", sort.by = "p", 
                                           number = dim(fit.fallahi.lines)[1])
sum(genes.fallahi.lines.mmacsf.48h$adj.P.Val < 0.05)
sum(genes.fallahi.lines.mmacsf.48h$adj.P.Val < 0.001)

# volcano plots
volcanoplot.colors.limma(genes = genes.fallahi.lines.colo.24h, adjpval = 0.05, 
    logfc = 1.5, title = "Fallahi, Colo858 treated 24h")
volcanoplot.colors.limma(genes = genes.fallahi.lines.colo.48h, adjpval = 0.05, 
    logfc = 1.5, title = "Fallahi, Colo858 treated 48h")
volcanoplot.colors.limma(genes = genes.fallahi.lines.mmacsf.24h, adjpval = 0.05,
    logfc = 1.5, title = "Fallahi, MMACSF treated 24h")
volcanoplot.colors.limma(genes = genes.fallahi.lines.mmacsf.48h, adjpval = 0.05,
    logfc = 1.5, title = "Fallahi, MMACSF treated 48h")


## corre dataset ----

# load data
files.corre <- read.table("../output-salmon/files-treat/corre.txt")
files.corre <- file.path(unlist(files.corre))
all(file.exists(files.corre))
tdata.corre <- as.data.frame(files.corre)
colnames(tdata.corre) <- "files"
cell.lines.corre <- c(rep("501Mel", 4))
tdata.corre$cell <- factor(cell.lines.corre)
treat.corre <- relevel(factor(c("P", "P", "BRAFi48h", "BRAFi48h"), 
                              levels = c("P", "BRAFi48h")), ref = "P")
tdata.corre$treatment <- treat.corre
tdata.corre$names <- paste("Sample", 1:4, sep = "")

# reference transcriptome checksum recognized by tximeta
se.corre <- tximeta(tdata.corre)  

# summarize quantifications to gene-level
gse.corre <- summarizeToGene(se.corre, countsFromAbundance = "lengthScaledTPM")  
y.corre <- makeDGEList(gse.corre)

# filter null and low counts
dim(y.corre)
design.corre <- model.matrix(~treat.corre)
keep.corre <- filterByExpr(y.corre, design.corre)
dge.corre <- y.corre[keep.corre, , keep.lib.sizes = FALSE]
dim(dge.corre)[1]/dim(y.corre)[1] * 100  

# plots
plot.filt(y.corre, dge.corre, c(0, 1.4))
plotMDS(dge.corre, col = as.numeric(as.factor(treat.corre)))

# differential expression analysis
v.corre <- voom(dge.corre, design.corre, plot = TRUE)

fit.corre <- lmFit(v.corre, design.corre)
fit.corre <- eBayes(fit.corre, trend = FALSE)

# 501Mel genes
genes.corre <- topTable(fit.corre, coef = 2, adjust.method = "BH", 
                        sort.by = "p", number = dim(fit.corre)[1])
sum(genes.corre$adj.P.Val < 0.05)
sum(genes.corre$adj.P.Val < 0.001)

# volcano plot
volcanoplot.colors.limma(genes = genes.corre, adjpval = 0.05, logfc = 1.5, 
                         title = "Corre, 501Mel")


## gerosa dataset ----

# load data
files.gerosa <- read.table("../output-salmon/files-treat/gerosa.txt")
files.gerosa <- file.path(unlist(files.gerosa))
all(file.exists(files.gerosa))
tdata.gerosa <- as.data.frame(files.gerosa)
colnames(tdata.gerosa) <- "files"
cell.lines.gerosa <- c(rep("A375", 8))
tdata.gerosa$cell <- factor(cell.lines.gerosa)
treat.gerosa <- relevel(factor(c("BRAFi24h", "BRAFi24h", "P", "P", "BRAFi24h",
    "BRAFi24h","P", "P"), levels = c("P", "BRAFi24h")), ref = "P")
tdata.gerosa$treatment <- treat.gerosa
tdata.gerosa$names <- paste("Sample", 1:8, sep = "")

# reference transcriptome checksum recognized by tximeta
se.gerosa <- tximeta(tdata.gerosa)  

# summarize quantifications to gene-level
gse.gerosa <- summarizeToGene(se.gerosa, 
                                countsFromAbundance = "lengthScaledTPM")  
y.gerosa <- makeDGEList(gse.gerosa)

# filter null and low counts
dim(y.gerosa)
design.gerosa <- model.matrix(~treat.gerosa)
keep.gerosa <- filterByExpr(y.gerosa, design.gerosa)
dge.gerosa <- y.gerosa[keep.gerosa, , keep.lib.sizes = FALSE]
dim(dge.gerosa)[1]/dim(y.gerosa)[1] * 100  

# plots
plot.filt(y.gerosa, dge.gerosa, c(0, 1.4))
plotMDS(dge.gerosa, col = as.numeric(as.factor(treat.gerosa)))

# differential expression analysis
v.gerosa <- voom(dge.gerosa, design.gerosa, plot = TRUE)

fit.gerosa <- lmFit(v.gerosa, design.gerosa)
fit.gerosa <- eBayes(fit.gerosa, trend = FALSE)


# A375 genes
genes.gerosa <- topTable(fit.gerosa, coef = 2, adjust.method = "BH", 
                         sort.by = "p", number = dim(fit.gerosa)[1])
sum(genes.gerosa$adj.P.Val < 0.05)
sum(genes.gerosa$adj.P.Val < 0.001)

# volcano plot
volcanoplot.colors.limma(genes = genes.gerosa, adjpval = 0.05, logfc = 1.5, 
                         title = "Gerosa, A375")

## smalley dataset ----

# laod data
files.smalley <- read.table("../output-salmon/files-treat/smalley.txt")
files.smalley <- file.path(unlist(files.smalley))
all(file.exists(files.smalley))
tdata.smalley <- as.data.frame(files.smalley)
colnames(tdata.smalley) <- "files"
cell.lines.smalley <- rep("WM164", 6)
tdata.smalley$cell <- factor(cell.lines.smalley)
treat.smalley <- relevel(factor(c(rep("P", 3), rep("BRAFi8h", 3)), 
                                levels = c("P", "BRAFi8h")), ref = "P")
tdata.smalley$treatment <- treat.smalley
tdata.smalley$names <- paste("Sample", 1:6, sep = "")

# reference transcriptome checksum recognized by tximeta
se.smalley <- tximeta(tdata.smalley)  

# summarize quantifications to gene-level
gse.smalley <- summarizeToGene(se.smalley, 
                               countsFromAbundance = "lengthScaledTPM")  
y.smalley <- makeDGEList(gse.smalley)

# filter low or null counts
dim(y.smalley)
design.smalley <- model.matrix(~treat.smalley)
keep.smalley <- filterByExpr(y.smalley, design.smalley)
dge.smalley <- y.smalley[keep.smalley, , keep.lib.sizes = FALSE]
dim(dge.smalley)[1]/dim(y.smalley)[1] * 100 

# plots
plot.filt(y.smalley, dge.smalley, c(0, 0.8))
plotMDS(dge.smalley, col = as.numeric(as.factor(treat.smalley)))

# differential expression analysis
v.smalley <- voom(dge.smalley, design.smalley, plot = TRUE)

fit.smalley <- lmFit(v.smalley, design.smalley)
fit.smalley <- eBayes(fit.smalley, trend = FALSE)

# WM164 8h treatment genes
genes.smalley <- topTable(fit.smalley, coef = 2, adjust.method = "BH",
                          sort.by = "p", number = dim(fit.smalley)[1])
sum(genes.smalley$adj.P.Val < 0.05)
sum(genes.smalley$adj.P.Val < 0.001)

# volcano plot
volcanoplot.colors.limma(genes = genes.smalley, adjpval = 0.05, logfc = 1.5, 
                         title = "Smalley, WM164")


## reganfendt dataset ----

# load data
files.reganfendt <- read.table("../output-salmon/files-treat/reganfendt.txt")
files.reganfendt <- file.path(unlist(files.reganfendt))
all(file.exists(files.reganfendt))
tdata.reganfendt <- as.data.frame(files.reganfendt)
colnames(tdata.reganfendt) <- "files"
cell.lines.reganfendt <- rep("A375", 6)
tdata.reganfendt$cell <- factor(cell.lines.reganfendt)
treat.reganfendt <- relevel(factor(c(rep("P", 3), rep("BRAFi8h", 3)), 
                                   levels = c("P", "BRAFi8h")), ref = "P")
tdata.reganfendt$treatment <- treat.reganfendt
tdata.reganfendt$names <- paste("Sample", 1:6, sep = "")

# reference transcriptome checksum recognized by tximeta
se.reganfendt <- tximeta(tdata.reganfendt)  

# summarize quantifications to gene-level
gse.reganfendt <- summarizeToGene(se.reganfendt, 
                                  countsFromAbundance = "lengthScaledTPM")  
y.reganfendt <- makeDGEList(gse.reganfendt)

# filter low or null counts
dim(y.reganfendt)
design.reganfendt <- model.matrix(~treat.reganfendt)
keep.reganfendt <- filterByExpr(y.reganfendt, design.reganfendt)
dge.reganfendt <- y.reganfendt[keep.reganfendt, , keep.lib.sizes = FALSE]
dim(dge.reganfendt)[1]/dim(y.reganfendt)[1] * 100  

# plots
plot.filt(y.reganfendt, dge.reganfendt, c(0, 1.4))
plotMDS(dge.reganfendt, col = as.numeric(as.factor(treat.reganfendt)))

# differential expression analysis
v.reganfendt <- voom(dge.reganfendt, design.reganfendt, plot = TRUE)

fit.reganfendt <- lmFit(v.reganfendt, design.reganfendt)
fit.reganfendt <- eBayes(fit.reganfendt, trend = FALSE)

# A375 8h genes
genes.reganfendt <- topTable(fit.reganfendt, coef = 2, adjust.method = "BH", 
                             sort.by = "p", number = dim(fit.reganfendt)[1])
sum(genes.reganfendt$adj.P.Val < 0.05)
sum(genes.reganfendt$adj.P.Val < 0.001)

# volcano plot
volcanoplot.colors.limma(genes = genes.reganfendt, adjpval = 0.05, logfc = 1.5, 
                         title = "Reganfendt, A375")


## song dataset ----

# load data
files.song <- read.table("../output-salmon/files-treat/song_treatment.txt")
files.song <- file.path(unlist(files.song))
all(file.exists(files.song))
tdata.song <- as.data.frame(files.song)
colnames(tdata.song) <- "files"
cell.lines.song <- c(rep("M229", 4), rep("M238", 4))
tdata.song$cell <- factor(cell.lines.song)
treat.song <- relevel(factor(rep(c(rep("P", 2), rep("BRAFi48h", 2)), 2), 
                             levels = c("P", "BRAFi48h")), ref = "P")
tdata.song$treatment <- treat.song
tdata.song$names <- paste("Sample", 1:8, sep = "")

# reference transcriptome checksum recognized by tximeta
se.song <- tximeta(tdata.song)  

# summarize quantifications to gene-level
gse.song <- summarizeToGene(se.song, countsFromAbundance = "lengthScaledTPM")  
y.song <- makeDGEList(gse.song)

# filter low or null counts
dim(y.song)
group.song.lines <- paste(cell.lines.song, treat.song, sep = ".")
design.song.lines <- model.matrix(~0 + group.song.lines)
colnames(design.song.lines) <- c("M229.BRAFi48h", "M229.P", 
                                 "M238.BRAFi48h", "M238.P")
keep.song.lines <- filterByExpr(y.song, design.song.lines)
dge.song.lines <- y.song[keep.song.lines, , keep.lib.sizes = FALSE]
dim(dge.song.lines)[1]/dim(y.song)[1] * 100 

# plots
plot.filt(y.song, dge.song.lines, c(0, 0.6))
plotMDS(y.song, col = as.numeric(as.factor(treat.song)), 
        labels = paste(dge.song.lines$samples$cell,
                       dge.song.lines$samples$treatment, sep = "."), 
        top = dim(dge.song.lines$counts)[1])

# differential expression analysis
contrast.song.lines <- makeContrasts(M229_48h_vs_cont = M229.BRAFi48h - M229.P, 
                                     M238_48h_vs_cont = M238.BRAFi48h - M238.P, 
                                     levels = design.song.lines)

v.song.lines <- voom(dge.song.lines, design.song.lines, plot = TRUE)

fit.song.lines <- lmFit(v.song.lines, design.song.lines)
fit.song.lines <- contrasts.fit(fit.song.lines, contrast.song.lines)
fit.song.lines <- eBayes(fit.song.lines, trend = FALSE)

# M229 genes
genes.song.lines.M229.48h <- topTable(fit.song.lines, coef = "M229_48h_vs_cont",
    adjust.method = "BH", sort.by = "p", number = dim(fit.song.lines)[1])
sum(genes.song.lines.M229.48h$adj.P.Val < 0.05)
sum(genes.song.lines.M229.48h$adj.P.Val < 0.001)

# M238 genes
genes.song.lines.M238.48h <- topTable(fit.song.lines, coef = "M238_48h_vs_cont",
    adjust.method = "BH", sort.by = "p", number = dim(fit.song.lines)[1])
sum(genes.song.lines.M238.48h$adj.P.Val < 0.05)
sum(genes.song.lines.M238.48h$adj.P.Val < 0.001)

# volcano plots
volcanoplot.colors.limma(genes = genes.song.lines.M229.48h, adjpval = 0.05, 
                         logfc = 1.5, title = "Song, M229 48h treatment")
volcanoplot.colors.limma(genes = genes.song.lines.M238.48h, adjpval = 0.05, 
                         logfc = 1.5, title = "Song, M238 48h treatment")


# ----------------------------- logFC and SE -----------------------------------

# obenauf A375 6h
lfcSE.obenauf.lines.A375.6h <- (sqrt(fit.obenauf.lines$s2.post) * 
                           fit.obenauf.lines$stdev.unscaled)[,"A375_6h_vs_cont"]
genes.obenauf.lines.A375.6h <- cbind(genes.obenauf.lines.A375.6h, 
                                     lfcSE = lfcSE.obenauf.lines.A375.6h[
                                       rownames(genes.obenauf.lines.A375.6h)])
genes.obenauf.lines.A375.6h.lfc <- genes.obenauf.lines.A375.6h[, 
                                                            c("logFC", "lfcSE")]

# obenauf A375 48h
lfcSE.obenauf.lines.A375.48h <- (sqrt(fit.obenauf.lines$s2.post) * 
                         fit.obenauf.lines$stdev.unscaled)[, "A375_48h_vs_cont"]
genes.obenauf.lines.A375.48h <- cbind(genes.obenauf.lines.A375.48h, 
                                      lfcSE = lfcSE.obenauf.lines.A375.48h[
                                        rownames(genes.obenauf.lines.A375.48h)])
genes.obenauf.lines.A375.48h.lfc <- genes.obenauf.lines.A375.48h[, 
                                                            c("logFC", "lfcSE")]

# obenauf Colo800 48h
lfcSE.obenauf.lines.colo.48h <- (sqrt(fit.obenauf.lines$s2.post) * 
                      fit.obenauf.lines$stdev.unscaled)[, "Colo800_48h_vs_cont"]
genes.obenauf.lines.colo.48h <- cbind(genes.obenauf.lines.colo.48h, 
                                      lfcSE = lfcSE.obenauf.lines.colo.48h[
                                        rownames(genes.obenauf.lines.colo.48h)])
genes.obenauf.lines.colo.48h.lfc <- genes.obenauf.lines.colo.48h[, 
                                                            c("logFC", "lfcSE")]

# obenauf UACC62 48h
lfcSE.obenauf.lines.uacc.48h <- (sqrt(fit.obenauf.lines$s2.post) * 
                       fit.obenauf.lines$stdev.unscaled)[, "UACC62_48h_vs_cont"]
genes.obenauf.lines.uacc.48h <- cbind(genes.obenauf.lines.uacc.48h, 
                                      lfcSE = lfcSE.obenauf.lines.uacc.48h[
                                        rownames(genes.obenauf.lines.uacc.48h)])
genes.obenauf.lines.uacc.48h.lfc <- genes.obenauf.lines.uacc.48h[, 
                                                            c("logFC", "lfcSE")]

# fallahi Colo848 24h
lfcSE.fallahi.lines.colo.24h <- (sqrt(fit.fallahi.lines$s2.post) * 
                      fit.fallahi.lines$stdev.unscaled)[, "Colo858_24h_vs_cont"]
genes.fallahi.lines.colo.24h <- cbind(genes.fallahi.lines.colo.24h, 
                                      lfcSE = lfcSE.fallahi.lines.colo.24h[
                                        rownames(genes.fallahi.lines.colo.24h)])
genes.fallahi.lines.colo.24h.lfc <- genes.fallahi.lines.colo.24h[, 
                                                            c("logFC", "lfcSE")]

# fallahi Colo848 48h
lfcSE.fallahi.lines.colo.48h <- (sqrt(fit.fallahi.lines$s2.post) * 
                      fit.fallahi.lines$stdev.unscaled)[, "Colo858_48h_vs_cont"]
genes.fallahi.lines.colo.48h <- cbind(genes.fallahi.lines.colo.48h, 
                                      lfcSE = lfcSE.fallahi.lines.colo.48h[
                                        rownames(genes.fallahi.lines.colo.48h)])
genes.fallahi.lines.colo.48h.lfc <- genes.fallahi.lines.colo.48h[, 
                                                            c("logFC", "lfcSE")]

# fallahi MMACSF 24h
lfcSE.fallahi.lines.mmacsf.24h <- (sqrt(fit.fallahi.lines$s2.post) * 
                       fit.fallahi.lines$stdev.unscaled)[, "MMACSF_24h_vs_cont"]
genes.fallahi.lines.mmacsf.24h <- cbind(genes.fallahi.lines.mmacsf.24h, 
                                        lfcSE = lfcSE.fallahi.lines.mmacsf.24h[
                                      rownames(genes.fallahi.lines.mmacsf.24h)])
genes.fallahi.lines.mmacsf.24h.lfc <- genes.fallahi.lines.mmacsf.24h[, 
                                                            c("logFC", "lfcSE")]

# fallahi MMACSF 48h
lfcSE.fallahi.lines.mmacsf.48h <- (sqrt(fit.fallahi.lines$s2.post) * 
                       fit.fallahi.lines$stdev.unscaled)[, "MMACSF_48h_vs_cont"]
genes.fallahi.lines.mmacsf.48h <- cbind(genes.fallahi.lines.mmacsf.48h, 
                                        lfcSE = lfcSE.fallahi.lines.mmacsf.48h[
                                      rownames(genes.fallahi.lines.mmacsf.48h)])
genes.fallahi.lines.mmacsf.48h.lfc <- genes.fallahi.lines.mmacsf.48h[, 
                                                            c("logFC", "lfcSE")]

# corre 501Mel 48h
lfcSE.corre <- (sqrt(fit.corre$s2.post) * 
                  fit.corre$stdev.unscaled)[, "treat.correBRAFi48h"]
genes.corre <- cbind(genes.corre, lfcSE = lfcSE.corre[rownames(genes.corre)])
genes.corre.lfc <- genes.corre[, c("logFC", "lfcSE")]

# gerosa A375 24h
lfcSE.gerosa <- (sqrt(fit.gerosa$s2.post) * 
                   fit.gerosa$stdev.unscaled)[, "treat.gerosaBRAFi24h"]
genes.gerosa <- cbind(genes.gerosa, 
                           lfcSE = lfcSE.gerosa[rownames(genes.gerosa)])
genes.gerosa.lfc <- genes.gerosa[, c("logFC", "lfcSE")]

# smalley WM164 8h
lfcSE.smalley <- (sqrt(fit.smalley$s2.post) * 
                    fit.smalley$stdev.unscaled)[, "treat.smalleyBRAFi8h"]
genes.smalley <- cbind(genes.smalley, 
                       lfcSE = lfcSE.smalley[rownames(genes.smalley)])
genes.smalley.lfc <- genes.smalley[, c("logFC", "lfcSE")]

# reganfendt A375 8h
lfcSE.reganfendt <- (sqrt(fit.reganfendt$s2.post) * 
               fit.reganfendt$stdev.unscaled)[, "treat.reganfendtBRAFi8h"]
genes.reganfendt <- cbind(genes.reganfendt, 
                          lfcSE = lfcSE.reganfendt[rownames(genes.reganfendt)])
genes.reganfendt.lfc <- genes.reganfendt[, c("logFC", "lfcSE")]

# song M229 48h
lfcSE.song.lines.M229.48h <- (sqrt(fit.song.lines$s2.post) * 
                           fit.song.lines$stdev.unscaled)[, "M229_48h_vs_cont"]
genes.song.lines.M229.48h <- cbind(genes.song.lines.M229.48h,
                                   lfcSE = lfcSE.song.lines.M229.48h[
                                           rownames(genes.song.lines.M229.48h)])
genes.song.lines.M229.48h.lfc <- genes.song.lines.M229.48h[,c("logFC", "lfcSE")]

# song M238 48h
lfcSE.song.lines.M238.48h <- (sqrt(fit.song.lines$s2.post) * 
                            fit.song.lines$stdev.unscaled)[, "M238_48h_vs_cont"]
genes.song.lines.M238.48h <- cbind(genes.song.lines.M238.48h, 
                                   lfcSE = lfcSE.song.lines.M238.48h[
                                     rownames(genes.song.lines.M238.48h)])
genes.song.lines.M238.48h.lfc <- genes.song.lines.M238.48h[,c("logFC", "lfcSE")]


################################################################################
#              2. METAANALYSES USING 6-48H TREATMENT DATA                      #
################################################################################

# These meta-analyses combine the effect size (LFC) of 6, 8, 24 or 48 h 
# treated samples (TE-MA)

# save gene names
all.genes <- rownames(y.obenauf.lines)

# function to format the results of the differential expression analysis  
# and prepare the data frame for the metaanalyses
search.gene <- function(gene) {
    mat <- matrix(ncol = 6, nrow = 0)
    df <- data.frame(mat)

    # obenauf A375 6h
    if (nrow(genes.obenauf.lines.A375.6h.lfc[
                   rownames(genes.obenauf.lines.A375.6h.lfc) == gene, ]) == 1) {
        a <- genes.obenauf.lines.A375.6h.lfc[
               rownames(genes.obenauf.lines.A375.6h.lfc) == gene, ][1, 1]
        b <- genes.obenauf.lines.A375.6h.lfc[
               rownames(genes.obenauf.lines.A375.6h.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Obenauf, A375 6h treatment", "Obenauf", "6h", 
            "A375", "short"))
    } else {
        df <- rbind(df, c(NA, NA, "Obenauf, A375 6h treatment", "Obenauf", "6h",
            "A375", "short"))
    }

    # obenauf A375 48h
    if (nrow(genes.obenauf.lines.A375.48h.lfc[
                  rownames(genes.obenauf.lines.A375.48h.lfc) == gene, ]) == 1) {
        a <- genes.obenauf.lines.A375.48h.lfc[
              rownames(genes.obenauf.lines.A375.48h.lfc) == gene, ][1, 1]
        b <- genes.obenauf.lines.A375.48h.lfc[
              rownames(genes.obenauf.lines.A375.48h.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Obenauf, A375 48h treatment", "Obenauf", "48h",
            "A375", "long"))
    } else {
        df <- rbind(df, c(NA, NA, "Obenauf, A375 48h treatment", "Obenauf", 
            "48h", "A375", "long"))
    }

    # obenauf Colo800 48h
    if (nrow(genes.obenauf.lines.colo.48h.lfc[
                  rownames(genes.obenauf.lines.colo.48h.lfc) == gene, ]) == 1) {
        a <- genes.obenauf.lines.colo.48h.lfc[
              rownames(genes.obenauf.lines.colo.48h.lfc) == gene, ][1, 1]
        b <- genes.obenauf.lines.colo.48h.lfc[
              rownames(genes.obenauf.lines.colo.48h.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Obenauf, Colo800 48h treatment", "Obenauf", 
            "48h", "Colo800", "long"))
    } else {
        df <- rbind(df, c(NA, NA, "Obenauf, Colo800 48h treatment", "Obenauf", 
            "48h", "Colo800", "long"))
    }

    # obenauf UACC62 48h
    if (nrow(genes.obenauf.lines.uacc.48h.lfc[
                  rownames(genes.obenauf.lines.uacc.48h.lfc) == gene, ]) == 1) {
        a <- genes.obenauf.lines.uacc.48h.lfc[
              rownames(genes.obenauf.lines.uacc.48h.lfc) == gene, ][1, 1]
        b <- genes.obenauf.lines.uacc.48h.lfc[
              rownames(genes.obenauf.lines.uacc.48h.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Obenauf, UACC62 48h treatment", "Obenauf",
            "48h", "UACC62", "long"))
    } else {
        df <- rbind(df, c(NA, NA, "Obenauf, UACC62 48h treatment", "Obenauf", 
            "48h", "UACC62", "long"))
    }

    # fallahi Colo858 24h
    if (nrow(genes.fallahi.lines.colo.24h.lfc[
                  rownames(genes.fallahi.lines.colo.24h.lfc) == gene, ]) == 1) {
        a <- genes.fallahi.lines.colo.24h.lfc[
              rownames(genes.fallahi.lines.colo.24h.lfc) == gene, ][1, 1]
        b <- genes.fallahi.lines.colo.24h.lfc[
              rownames(genes.fallahi.lines.colo.24h.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Fallahi, Colo858 24h treatment", "Fallahi", 
            "24h", "Colo858", "long"))
    } else {
        df <- rbind(df, c(NA, NA, "Fallahi, Colo858 24h treatment", "Fallahi", 
            "24h", "Colo858", "long"))
    }

    # fallahi Colo858 48h
    if (nrow(genes.fallahi.lines.colo.48h.lfc[
                  rownames(genes.fallahi.lines.colo.48h.lfc) == gene, ]) == 1) {
        a <- genes.fallahi.lines.colo.48h.lfc[
              rownames(genes.fallahi.lines.colo.48h.lfc) == gene, ][1, 1]
        b <- genes.fallahi.lines.colo.48h.lfc[
              rownames(genes.fallahi.lines.colo.48h.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Fallahi, Colo858 48h treatment", "Fallahi", 
            "48h", "Colo858", "long"))
    } else {
        df <- rbind(df, c(NA, NA, "Fallahi, Colo858 48h treatment", "Fallahi", 
            "48h", "Colo858", "long"))
    }

    # fallahi MMACSF 24h
    if (nrow(genes.fallahi.lines.mmacsf.24h.lfc[
                rownames(genes.fallahi.lines.mmacsf.24h.lfc) == gene, ]) == 1) {
        a <- genes.fallahi.lines.mmacsf.24h.lfc[
              rownames(genes.fallahi.lines.mmacsf.24h.lfc) == gene, ][1, 1]
        b <- genes.fallahi.lines.mmacsf.24h.lfc[
              rownames(genes.fallahi.lines.mmacsf.24h.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Fallahi, MMACSF 24h treatment", "Fallahi", 
            "24h", "MMACSF", "long"))
    } else {
        df <- rbind(df, c(NA, NA, "Fallahi, MMACSF 24h treatment", "Fallahi", 
            "24h", "MMACSF", "long"))
    }

    # fallahi MMACSF 48h
    if (nrow(genes.fallahi.lines.mmacsf.48h.lfc[
                rownames(genes.fallahi.lines.mmacsf.48h.lfc) == gene, ]) == 1) {
        a <- genes.fallahi.lines.mmacsf.48h.lfc[
              rownames(genes.fallahi.lines.mmacsf.48h.lfc) == gene, ][1, 1]
        b <- genes.fallahi.lines.mmacsf.48h.lfc[
              rownames(genes.fallahi.lines.mmacsf.48h.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Fallahi, MMACSF 48h treatment", "Fallahi", 
            "48h", "MMACSF", "long"))
    } else {
        df <- rbind(df, c(NA, NA, "Fallahi, MMACSF 48h treatment", "Fallahi", 
            "48h", "MMACSF", "long"))
    }

    # corre 501Mel 48h
    if (nrow(genes.corre.lfc[rownames(genes.corre.lfc) == gene, ]) == 1) {
        a <- genes.corre.lfc[rownames(genes.corre.lfc) == gene, ][1, 1]
        b <- genes.corre.lfc[rownames(genes.corre.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Corre, 501Mel 48h treatment", "Corre", "48h",
            "501Mel", "long"))
    } else {
        df <- rbind(df, c(NA, NA, "Corre, 501Mel 48h treatment", "Corre", "48h",
            "501Mel", "long"))
    }

    # gerosa A375 24h
    if (nrow(genes.gerosa.lfc[rownames(genes.gerosa.lfc) == gene, ]) == 1) {
        a <- genes.gerosa.lfc[rownames(genes.gerosa.lfc) == gene, ][1, 1]
        b <- genes.gerosa.lfc[rownames(genes.gerosa.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Gerosa, A375 24h treatment", "Gerosa", "24h", 
            "A375", "long"))
    } else {
        df <- rbind(df, c(NA, NA, "Gerosa, A375 24h treatment", "Gerosa", "24h",
            "A375", "long"))
    }

    # smalley WM164 8h
    if (nrow(genes.smalley.lfc[rownames(genes.smalley.lfc) == gene, ]) == 1) {
        a <- genes.smalley.lfc[rownames(genes.smalley.lfc) == gene, ][1, 1]
        b <- genes.smalley.lfc[rownames(genes.smalley.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Smalley, WM164 8h treatment", "Smalley", "8h", 
            "WM164", "short"))
    } else {
        df <- rbind(df, c(NA, NA, "Smalley, WM164 8h treatment", "Smalley", 
            "8h", "WM164", "short"))
    }

    # reganfendt A375 8h
    if (nrow(genes.reganfendt.lfc[rownames(genes.reganfendt.lfc) == gene, ]) 
                                                                         == 1) {
        a <- genes.reganfendt.lfc[rownames(genes.reganfendt.lfc) == gene, ][1,1]
        b <- genes.reganfendt.lfc[rownames(genes.reganfendt.lfc) == gene, ][1,2]
        df <- rbind(df, c(a, b, "Reganfendt, A375 8h treatment", "Reganfendt", 
            "8h", "A375", "short"))
    } else {
        df <- rbind(df, c(NA, NA, "Reganfendt, A375 8h treatment", "Reganfendt",
            "8h", "A375", "short"))
    }

    # song M229 48h
    if (nrow(genes.song.lines.M229.48h.lfc[
                     rownames(genes.song.lines.M229.48h.lfc) == gene, ]) == 1) {
        a <- genes.song.lines.M229.48h.lfc[
              rownames(genes.song.lines.M229.48h.lfc) == gene, ][1, 1]
        b <- genes.song.lines.M229.48h.lfc[
              rownames(genes.song.lines.M229.48h.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Song, M229 48h treatment", "Song", "48h", 
            "M229", "long"))
    } else {
        df <- rbind(df, c(NA, NA, "Song, M229 48h treatment", "Song", "48h", 
            "M229", "long"))
    }

    # song M238 48h
    if (nrow(genes.song.lines.M238.48h.lfc[
                     rownames(genes.song.lines.M238.48h.lfc) == gene, ]) == 1) {
        a <- genes.song.lines.M238.48h.lfc[
              rownames(genes.song.lines.M238.48h.lfc) == gene, ][1, 1]
        b <- genes.song.lines.M238.48h.lfc[
              rownames(genes.song.lines.M238.48h.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Song, M238 48h treatment", "Song", "48h", 
            "M238", "long"))
    } else {
        df <- rbind(df, c(NA, NA, "Song, M238 48h treatment", "Song", "48h", 
            "M238", "long"))
    }

    # end of search
    colnames(df) <- c("logFC", "lfcSE", "studylab", "author", "hours", "line", 
                      "treatment")
    return(df)
}

# example of use
search.gene("ENSG00000000003.15")

# loop to iterate over all the genes using the function
data.meta.lines <- lapply(all.genes, FUN = search.gene)
names(data.meta.lines) <- all.genes

length(data.meta.lines)

# remove null results
for (g in all.genes) {
    if (all(is.na(data.meta.lines[[g]][, c(1, 2)]))) {
        data.meta.lines[[g]] <- NULL
    }
}
length(data.meta.lines)

# remove genes with info from only 1 study
for (g in names(data.meta.lines)) {
    if (sum(is.na(data.meta.lines[[g]])) == 26) {
        data.meta.lines[[g]] <- NULL
    }
}
length(data.meta.lines)

# remove digits after dots and save rds
names(data.meta.lines) <- as.vector(mapply(names(data.meta.lines), 
    FUN = function(x) {gsub("\\.\\d+", "", x)}))

# saveRDS(data.meta.lines, file = 'data-meta-lines.rds')


## metaanalyses

# metaanalysis function
metaanalyses <- function(data.metaanalyses, gene) {
    data.m <- data.metaanalyses[[gene]]
    data.m[, c("logFC", "lfcSE")] <- sapply(data.m[c("logFC", "lfcSE")], 
                                                                   as.numeric)
    m.gen <- metagen(TE = logFC, 
                     seTE = lfcSE, 
                     studlab = studylab, 
                     data = data.m,
                     sm = "MD", 
                     fixed = FALSE, 
                     random = TRUE, 
                     method.tau = "REML", 
                     id = author,
                     title = gene, 
                     prediction = TRUE, 
                     control = list(optimizer = "nmk"))
    summary(m.gen)
    return(m.gen)
}

# use a loop to iterate over all the genes using the function
counter <- 0
total <- length(data.meta.lines)
genes.meta.lines <- lapply(names(data.meta.lines), FUN = function(x) {
    counter <<- counter + 1
    print(paste0("iteration number ", counter, "; total ", total))
    metaanalyses(data.meta.lines, x)
})

names(genes.meta.lines) <- names(data.meta.lines)

# remove null results and save rds
length(genes.meta.lines)

for (g in names(data.meta.lines)) {
    if (is.null(genes.meta.lines[[g]])) {
        genes.meta.lines[[g]] <- NULL
    }
}

length(genes.meta.lines)

# saveRDS(genes.meta.lines, file = 'genes-meta-lines.rds') 


# summarize results in a data frame
results.meta <- function(meta, gene) {
    res <- data.frame(
      num.studies = meta[[gene]]$k.study, 
      num.data = meta[[gene]]$k,
      lfc = meta[[gene]]$TE.random, 
      lfcSE = meta[[gene]]$seTE.random, 
      CI95lower = meta[[gene]]$lower.random,
      CI95upper = meta[[gene]]$upper.random, 
      zval = meta[[gene]]$zval.random, 
      stat = meta[[gene]]$statistic.random,
      pval = meta[[gene]]$pval.random, 
      I2 = meta[[gene]]$I2, 
      Q = meta[[gene]]$Q,
      dfQ = meta[[gene]]$df.Q, 
      pvalQ = meta[[gene]]$pval.Q)
    return(res)
}

# example of use
results.meta(genes.meta.lines, "ENSG00000000003")

# loop to iterate over all the genes using the function
res.meta <- lapply(names(genes.meta.lines), FUN = function(x) {
    results.meta(genes.meta.lines, x)
})

# add names and transform list into a data frame
names(res.meta) <- names(genes.meta.lines)
res.meta.df <- do.call(rbind.data.frame, res.meta)
dim(res.meta.df)

# multiple testing correction: FDR/BH
res.meta.df$adj.pval <- p.adjust(res.meta.df$pval, method = "BH", 
                                 n = length(res.meta.df$pval))

# multiple testing correction: Bonferroni
res.meta.df$adj.pval.bonfer <- p.adjust(res.meta.df$pval, method = "bonferroni",
    n = length(res.meta.df$pval))

# saveRDS(res.meta.df, file='res-meta-df.rds')


# ----------------------- metaanalyses results ---------------------------------

# metaanalyses vs. experiments correlation 

## experiments

# save the LFC values of the experiments in a new variable
lfc.genes <- lapply(data.meta.lines, '[[', 'logFC')
lfc.genes <- do.call('rbind.data.frame', lfc.genes)
lfc.genes <- sapply(lfc.genes, as.numeric)
colnames(lfc.genes) <- data.meta.lines[[1]]$studylab
rownames(lfc.genes) <- names(data.meta.lines)
lfc.genes <- lfc.genes[names(genes.meta.lines), ]
dim(lfc.genes)


## metaanalyses

# save the LFC values of the metaanalyses in a new variable
lfc.meta <- lapply(genes.meta.lines, '[[', 'TE.random')
lfc.meta <- do.call('rbind.data.frame', lfc.meta)
lfc.meta <- sapply(lfc.meta, as.numeric)
colnames(lfc.meta) <- 'lfc.meta'
rownames(lfc.meta) <- names(genes.meta.lines)

# merge both
lfc.genes <- merge(lfc.genes, lfc.meta, by = 0)
rownames(lfc.genes) <- lfc.genes$Row.names
lfc.genes <- lfc.genes[ , -1]

# plot correlations using all the genes: metaanalyses vs. each experiment
lfc.genes.cor <- lfc.genes

# change downregulated genes sign first
lfc.genes.cor[lfc.genes[, "lfc.meta"] < 0, "lfc.meta"] <- -1 * 
  lfc.genes.cor[lfc.genes[, "lfc.meta"] < 0, "lfc.meta"]

# and plot
par(mfrow = c(4, 4))
for (i in 1:14) {
  lfc.genes.cor[lfc.genes[, "lfc.meta"] < 0, i] <- -1 * 
    lfc.genes.cor[lfc.genes[, "lfc.meta"] < 0, i]
  plot(lfc.genes.cor[, i], lfc.genes.cor[, "lfc.meta"], 
       main = colnames(lfc.genes.cor)[i], ylab = "meta lfc", 
       xlab = "individual lfc", pch = 21, col = "blue4", bg = "blue1",
       cex = 1, lwd = 1)
}
par(mfrow = c(1, 1))

# compute correlations using all the genes
first.cor <- cor.test(lfc.genes.cor[, 1], lfc.genes.cor[, "lfc.meta"], 
                                                            method = "pearson")
cor.res.all <- data.frame(cor = first.cor$estimate, pval = first.cor$p.value)
for (i in 2:14) {
    cor <- cor.test(lfc.genes.cor[, i], lfc.genes.cor[, "lfc.meta"],
                                                            method = "pearson")
    cor.res.all <- rbind(cor.res.all, c(cor$estimate, cor$p.value))
}
dim(cor.res.all)
cor.res.all

# now take only differentially expressed genes (DEGs) in order to reduce noise 
# in correlations and plots
dim(lfc.genes)
lfc.genes.05 <- lfc.genes[rownames(res.meta.df[
                                     res.meta.df$adj.pval <= 0.05, ]),]
dim(lfc.genes.05)

# plots: down in blue and up in red
par(mfrow = c(7, 4))
for (i in 1:14) {
    plot(lfc.genes.05[lfc.genes.05[, "lfc.meta"] > 0, i], 
         lfc.genes.05[lfc.genes.05[, "lfc.meta"] > 0, "lfc.meta"], 
         main = paste(colnames(lfc.genes.05)[i], "- UP",
        sep = " "), ylab = "meta lfc", xlab = "individual lfc", pch = 21, 
        col = "darkred", bg = "red", cex = 1, lwd = 1)
    plot(-1 * (lfc.genes.05[lfc.genes.05[, "lfc.meta"] < 0, i]), -1 * 
           (lfc.genes.05[lfc.genes.05[, "lfc.meta"] < 0, "lfc.meta"]), 
         main = paste(colnames(lfc.genes.05)[i], "- DOWN",
         sep = " "), ylab = "meta lfc", xlab = "individual lfc", pch = 21, 
         col = "blue4", bg = "blue1", cex = 1, lwd = 1)
}
par(mfrow = c(1, 1))

# DEGs correlations
lfc.genes.05.cor <- lfc.genes.05
for (i in 1:14) {
    lfc.genes.05.cor[lfc.genes.05.cor[, "lfc.meta"] < 0, i] <- -1 * 
      (lfc.genes.05.cor[lfc.genes.05.cor[, "lfc.meta"] < 0, i])
}

# change downregulated DEGs sign first
lfc.genes.05.cor[lfc.genes.05.cor[, "lfc.meta"] < 0, "lfc.meta"] <- -1 * 
  (lfc.genes.05.cor[lfc.genes.05.cor[, "lfc.meta"] < 0, "lfc.meta"])

# compute correlations
first.cor <- cor.test(lfc.genes.05.cor[, 1], lfc.genes.05.cor[, "lfc.meta"], 
                                                            method = "pearson")
cor.res <- data.frame(cor = first.cor$estimate, pval = first.cor$p.value)
for (i in 2:14) {
    cor <- cor.test(lfc.genes.05.cor[, i], lfc.genes.05.cor[, "lfc.meta"], 
                                                            method = "pearson")
    cor.res <- rbind(cor.res, c(cor$estimate, cor$p.value))
}
dim(cor.res)
cor.res

# histogram of correlations
p.hist <- ggplot(data = cor.res, aes(x = cor)) + geom_histogram(color = "black",
    fill = "lightcyan3", bins = 15) + xlim(c(0, 1)) +
  xlab("Pearson correlation coefficient") +
  ylab("Frequency") + theme_light() + 
  theme(plot.title = element_text(hjust = 0.5, size = 20), 
        axis.title = element_text(size = 8), axis.text = element_text(size = 8))

p.hist

# plot DEGs in one plot for each experiment
par(mfrow = c(4, 4))
for (i in 1:14) {
    plot(lfc.genes.05.cor[, i], lfc.genes.05.cor[, "lfc.meta"], 
         main = colnames(lfc.genes.05.cor)[i],
         ylab = "meta lfc", xlab = "individual lfc", 
         pch = 21, col = "darkorchid4",
         bg = "darkorchid2", cex = 1, lwd = 1)
}
par(mfrow = c(1, 1))

# and now plot using ggplot
for (i in 1:14) {
    dd <- as.data.frame(lfc.genes.05.cor[, c(i, 15)])
    colnames(dd) <- c("LFC.individual", "LFC.MA")
    plot.cor <- ggplot(dd, aes(x = LFC.individual, y = LFC.MA)) + 
      geom_point(aes(fill = "A"), colour = "blue4", pch = 21, size = 2.5) + 
      scale_fill_manual(values = "blue1") + theme_bw() + 
      ylab("LFC metaanalyses") + xlab("LFC individual") + 
      ggtitle(sub("treatment", "MAPKi", colnames(lfc.genes.05.cor)[i])) + 
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5))
    assign(paste0("p.cor.", i), plot.cor)
}

figure <- ggarrange(p.cor.1, p.cor.2, p.cor.3, p.cor.4, p.cor.5, p.cor.6, 
    p.cor.7, p.cor.8, p.cor.9, p.cor.10, p.cor.11, p.cor.12, p.cor.13, p.cor.14, 
    ncol = 4, nrow = 4)

figure


## comparison with RankProd and RankSum

# load rankprod and ranksum results
rp.results.tables <- readRDS("data/rp-results-tables.rds")
rs.results.tables <- readRDS("data/rs-results-tables.rds")

RP.genes <- c(rownames(rp.results.tables$Table1), 
              rownames(rp.results.tables$Table2))
RS.genes <- c(rownames(rs.results.tables$Table1), 
              rownames(rs.results.tables$Table2))

# now take DEGs that were obtained with values of expression in all samples 
# in TE-MA
res.meta.df.05.14st <- res.meta.df[res.meta.df$num.data == 14 &
                                     res.meta.df$adj.pval < 0.05, ]

# now compare TE-MA, RankProd and RankSum DEGs
ggVennDiagram(list(TE_MA = rownames(res.meta.df.05.14st), 
                   RankProd_MA = RP.genes,
                   RankSum_MA = RS.genes),
              color = "black", lwd = 0.8, lty = 1, label = "count") +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  scale_color_manual(values = rep("darkgray", 3)) + 
  theme_void() + theme(legend.position = "none")

# now compute correlations selecting only DEGs detected by the 3 MA
degs.14st <- intersect(rownames(res.meta.df.05.14st), 
                       intersect(RP.genes, RS.genes))
lfc.genes.05.cor.14st <- lfc.genes.05.cor[degs.14st, ]
dim(lfc.genes.05.cor.14st)

first.cor <- cor.test(lfc.genes.05.cor.14st[, 1], 
                      lfc.genes.05.cor.14st[, "lfc.meta"], method = "pearson")
cor.res.degs.14st <- data.frame(cor = first.cor$estimate, 
                                pval = first.cor$p.value)

for (i in 2:14) {
    cor <- cor.test(lfc.genes.05.cor.14st[, i], 
                    lfc.genes.05.cor.14st[, "lfc.meta"], method = "pearson")
    cor.res.degs.14st <- rbind(cor.res.degs.14st, c(cor$estimate, cor$p.value))
}
dim(cor.res.degs.14st)
cor.res.degs.14st

# plot
par(mfrow = c(4, 4))
for (i in 1:14) {
    plot(lfc.genes.05.cor.14st[, i], lfc.genes.05.cor.14st[, "lfc.meta"], 
         main = colnames(lfc.genes.05.cor.14st)[i], ylab = "meta lfc", 
         xlab = "experiment lfc", pch = 21, col = "darkorchid4",
         bg = "darkorchid2", cex = 1, lwd = 1)
}
par(mfrow = c(1, 1))


# metaanalysis and experiments values comparison

## metaanalyses

res.meta.df.05 <- res.meta.df[degs.14st, ]
dim(res.meta.df.05)[1]/dim(res.meta.df)[1] * 100
dim(res.meta.df.05)

# upregulated genes
sum(res.meta.df.05$lfc > 0)  # n upregulated genes
median(res.meta.df.05[res.meta.df.05$lfc > 0, ]$lfc)  # median upregulated g

# downregulated genes
sum(res.meta.df.05$lfc < 0)  # n downregulated genes
median(res.meta.df.05[res.meta.df.05$lfc < 0, ]$lfc)  # median downregulated g

## experiments
degs.lfc.indiv <- merge(genes.obenauf.lines.A375.6h[
          genes.obenauf.lines.A375.6h$adj.P.Val < 0.05, "logFC", drop = FALSE], 
                        genes.obenauf.lines.A375.48h[
          genes.obenauf.lines.A375.48h$adj.P.Val < 0.05, "logFC", drop = FALSE], 
          all = TRUE, by = 0) %>%
    merge(genes.obenauf.lines.colo.48h[
      genes.obenauf.lines.colo.48h$adj.P.Val < 0.05, "logFC", drop = FALSE], 
      all = TRUE, by.x = "Row.names", by.y = 0) %>%
    merge(genes.obenauf.lines.uacc.48h[
      genes.obenauf.lines.uacc.48h$adj.P.Val < 0.05, "logFC", drop = FALSE], 
      all = TRUE, by.x = "Row.names", by.y = 0) %>%
    merge(genes.fallahi.lines.colo.24h[
      genes.fallahi.lines.colo.24h$adj.P.Val < 0.05, "logFC", drop = FALSE], 
      all = TRUE, by.x = "Row.names", by.y = 0) %>%
    merge(genes.fallahi.lines.colo.48h[
      genes.fallahi.lines.colo.48h$adj.P.Val < 0.05, "logFC", drop = FALSE], 
      all = TRUE, by.x = "Row.names", by.y = 0) %>%
    merge(genes.fallahi.lines.mmacsf.24h[
      genes.fallahi.lines.mmacsf.24h$adj.P.Val < 0.05, "logFC", drop = FALSE], 
      all = TRUE, by.x = "Row.names", by.y = 0) %>%
    merge(genes.fallahi.lines.mmacsf.48h[
      genes.fallahi.lines.mmacsf.48h$adj.P.Val < 0.05, "logFC", drop = FALSE], 
      all = TRUE, by.x = "Row.names", by.y = 0) %>%
    merge(genes.corre[
      genes.corre$adj.P.Val < 0.05, "logFC", drop = FALSE], 
      all = TRUE, by.x = "Row.names", by.y = 0) %>%
    merge(genes.gerosa[
      genes.gerosa$adj.P.Val < 0.05, "logFC", drop = FALSE], 
      all = TRUE, by.x = "Row.names", by.y = 0) %>%
    merge(genes.smalley[
      genes.smalley$adj.P.Val < 0.05, "logFC", drop = FALSE], 
      all = TRUE, by.x = "Row.names", by.y = 0) %>%
    merge(genes.reganfendt[
      genes.reganfendt$adj.P.Val < 0.05, "logFC", drop = FALSE], 
      all = TRUE, by.x = "Row.names", by.y = 0) %>%
    merge(genes.song.lines.M229.48h[
      genes.song.lines.M229.48h$adj.P.Val < 0.05, "logFC", drop = FALSE], 
      all = TRUE, by.x = "Row.names", by.y = 0) %>%
    merge(genes.song.lines.M238.48h[
      genes.song.lines.M238.48h$adj.P.Val < 0.05, "logFC", drop = FALSE], 
      all = TRUE, by.x = "Row.names", by.y = 0)
rownames(degs.lfc.indiv) <- degs.lfc.indiv$Row.names
degs.lfc.indiv <- degs.lfc.indiv[, -1]
rownames(degs.lfc.indiv) <- as.vector(mapply(rownames(degs.lfc.indiv), 
    FUN = function(x) {gsub("\\.\\d+", "", x)}))
dim(degs.lfc.indiv)

# take the median from each experiment and the compute mean of medians
md.up <- c()
md.down <- c()
for (i in 1:14) {
    md.up <- c(md.up, median(degs.lfc.indiv[degs.lfc.indiv[, i] > 0, i], 
                                                                 na.rm = TRUE))
    md.down <- c(md.down, median(degs.lfc.indiv[degs.lfc.indiv[, i] < 0, i], 
                                                                 na.rm = TRUE))
}
mean(md.up)
mean(md.down)

# plot distributions
p.degs.df <- rbind(data.frame(LFC = res.meta.df.05$lfc, 
                              COND = rep("MA", length(res.meta.df.05$lfc))),
    data.frame(LFC = degs.lfc.indiv[, 1], 
               COND = rep("EXP1", length(degs.lfc.indiv[, 1]))), 
    data.frame(LFC = degs.lfc.indiv[, 2], 
               COND = rep("EXP2", length(degs.lfc.indiv[, 2]))), 
    data.frame(LFC = degs.lfc.indiv[, 3], 
               COND = rep("EXP3", length(degs.lfc.indiv[, 3]))), 
    data.frame(LFC = degs.lfc.indiv[, 4], 
               COND = rep("EXP4", length(degs.lfc.indiv[, 4]))), 
    data.frame(LFC = degs.lfc.indiv[, 5], 
               COND = rep("EXP5", length(degs.lfc.indiv[, 5]))), 
    data.frame(LFC = degs.lfc.indiv[, 6], 
               COND = rep("EXP6", length(degs.lfc.indiv[, 6]))), 
    data.frame(LFC = degs.lfc.indiv[, 7], 
               COND = rep("EXP7", length(degs.lfc.indiv[, 7]))), 
    data.frame(LFC = degs.lfc.indiv[, 8], 
               COND = rep("EXP8", length(degs.lfc.indiv[, 8]))), 
    data.frame(LFC = degs.lfc.indiv[, 9], 
               COND = rep("EXP9", length(degs.lfc.indiv[, 9]))), 
    data.frame(LFC = degs.lfc.indiv[, 10], 
               COND = rep("EXP10", length(degs.lfc.indiv[, 10]))), 
    data.frame(LFC = degs.lfc.indiv[, 11], 
               COND = rep("EXP11", length(degs.lfc.indiv[, 11]))), 
    data.frame(LFC = degs.lfc.indiv[, 12], 
               COND = rep("EXP12", length(degs.lfc.indiv[, 12]))), 
    data.frame(LFC = degs.lfc.indiv[, 13], 
               COND = rep("EXP13", length(degs.lfc.indiv[, 13]))), 
    data.frame(LFC = degs.lfc.indiv[, 14], 
               COND = rep("EXP14", length(degs.lfc.indiv[, 14]))))
p.degs.df$COND <- factor(p.degs.df$COND, levels = c("EXP1", "EXP2", "EXP3", 
    "EXP4", "EXP5", "EXP6", "EXP7", "EXP8", "EXP9", "EXP10", "EXP11", "EXP12", 
    "EXP13", "EXP14", "MA"), ordered = TRUE)
p.degs.df$EXPRESSION <- NA
p.degs.df$EXPRESSION[p.degs.df$LFC > 0] <- 'UP'
p.degs.df$EXPRESSION[p.degs.df$LFC < 0] <- 'DOWN'

p.degs <- ggplot(subset(p.degs.df, !is.na(LFC)), 
                 aes(EXPRESSION, LFC, fill = COND)) +
  geom_boxplot(alpha = 1, outlier.shape = NA) + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8), 
        legend.title = element_blank()) +
  labs(y = "LFC",
       x = "") +
  ylim (c(-4, 4)) +
  scale_x_discrete(labels = c("DOWN", "UP")) +
  scale_fill_manual(values = c(rep('gray68', 14), 'limegreen')) +
  geom_hline(yintercept = c(mean(md.up), mean(md.down)), col="red", lty = 5)

p.degs


## heatmap: short vs long-term treatment
# compare DEGs in long and short treatment data
lfc.genes.05.14st <- lfc.genes[degs.14st,]
colnames(lfc.genes.heatmap) <- sub("treatment", "MAPKi",
                                   colnames(lfc.genes.heatmap))

col.hm <- colorRampPalette(c("green", "black", "red"))(50)
pheatmap(lfc.genes.05.14st[, -15], col = col.hm, show_rownames = FALSE, 
         show_colnames = TRUE, breaks = seq(-4, 4, length = 50), 
         fontsize_row = 7, border_color = NA)


# ------------------------ subgroup analysis -----------------------------------

# select genes expressed in all the studies
genes.meta.lines.14st <- genes.meta.lines[res.meta.df$num.data == 14]

# example of subgroup analysis
update.meta(genes.meta.lines.14st[[1]], tau.common = TRUE, subgroup = treatment)

# now define function and do the subgroup analysis on all the genes
i <- 0
error.subg.genes <- c()
sub.lines.fx <- function(x, meta.list) {
    print(paste0("iteration number ", i))
    i <<- i + 1
    res <- try(update.meta(meta.list[[x]], subgroup = treatment, 
                           tau.common = TRUE))
    if (class(res) == "try-error") {
        res <- try(update.meta(meta.list[[x]], subgroup = as.factor(treatment), 
                               tau.common = TRUE))
        if (class(res) == "try-error") {
            res <- NULL
            error.genes <<- c(error.subg.genes, x)
        }
    }
    return(res)
}

# example of use
sub.lines.fx("ENSG00000000003", genes.meta.lines.14st)

# loop to iterate over all the genes using the function
subg.meta.lines <- lapply(names(genes.meta.lines.8st), FUN = function(x) {
    sub.lines.fx(x, genes.meta.lines.8st)
})

names(subg.meta.treatment) <- names(genes.meta.lines.8st)
# saveRDS(subg.meta.treatment, 'subg-meta-treatment.rds')

# null results
null.subg <- lapply(names(subg.meta.treatment), function(x) {
    is.null(subg.meta.treatment[[x]])
})
null.subg.names <- names(subg.meta.treatment)[
           which(do.call(rbind.data.frame, null.subg) == TRUE)]
subg.meta.treatment[null.subg.names]

# remove those items (2) and and continue analyzing results
length(subg.meta.treatment)
subg.meta.treatment[null.subg.names] <- NULL
length(subg.meta.treatment)

# function to retrieve subgroup analysis significances
diff.subg.meta <- function(meta, gene) {
    res <- data.frame(pval.btw = meta[[gene]]$pval.Q.b.random, 
                      pval.wth = meta[[gene]]$pval.Q.w.random)
    return(res)
}

# example of use
diff.subg.meta(subg.meta.treatment, "ENSG00000000003")

# loop to iterate over all the genes using the function
diff.subg.data <- lapply(names(subg.meta.treatment), 
                         function(x) {diff.subg.meta(subg.meta.treatment, x)})

names(diff.subg.data) <- names(subg.meta.treatment)
diff.subg.data <- do.call(rbind.data.frame, diff.subg.data)

# select genes with the same expression within groups but not between them
diff.genes.btw <- diff.subg.data[diff.subg.data$pval.btw <= 0.05 & 
                                   diff.subg.data$pval.wth > 0.05, ]

# add symbols 
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
btw.sy <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), 
                filters = "ensembl_gene_id",
                values = rownames(diff.genes.btw), 
                mart = mart)

diff.genes.btw <- merge(btw.sy, diff.genes.btw, 
                        by.x = "ensembl_gene_id", by.y = 0, all = TRUE)
rownames(diff.genes.btw) <- diff.genes.btw$ensembl_gene_id
diff.genes.btw <- diff.genes.btw[, -1]
dim(diff.genes.btw)

# heatmap 
lfc.genes.subgdiff <- lfc.genes[rownames(diff.genes.btw), -15]
colnames(lfc.genes.subgdiff) <- sub("treatment", "MAPKi", 
                                    colnames(lfc.genes.subgdiff))
p.subgdiff <- pheatmap(lfc.genes.subgdiff, col = col.hm, 
                       cluster_rows = TRUE, cluster_cols = TRUE,
                       breaks = seq(-2, 2, length = 50), 
                       show_rownames = FALSE, fontsize_row = 7)

p.subgdiff


# ---------------------- metaregression analysis -------------------------------

hours.reg <- as.numeric(sub("h", "", genes.meta.lines[[1]]$data$hours))

# loop to iterate over genes with differences btw groups performing 
# a metaregression analysis
subg.reg <- lapply(rownames(diff.genes.btw), 
                FUN = function(x) {metareg(genes.meta.lines[[x]], ~hours.reg)})
names(subg.reg) <- rownames(diff.genes.btw)

# function to retrieve info from the results
subg.reg.fx <- function(x) {
  res <- data.frame(intrcpt = subg.reg[[x]]$beta[1], 
                    pval.intrcpt = subg.reg[[x]]$pval[1],
                    hours.var = subg.reg[[x]]$beta[2], 
                    pval.hours.var = subg.reg[[x]]$pval[2])
    return(res)
}

# function to iterate over the genes using the function
subg.reg.data <- lapply(rownames(diff.genes.btw), 
                        FUN = function(x) {subg.reg.fx(x)})
subg.reg.data <- do.call(rbind.data.frame, subg.reg.data)
rownames(subg.reg.data) <- rownames(diff.genes.btw)

# significant estimates of the regression weight for 'hours'
dim(subg.reg.data)
sum(subg.reg.data$pval.hours.var < 0.05)
subg.reg.data.05 <- subg.reg.data[subg.reg.data$pval.hours.var < 0.05, ]

# plot one result
bubble(subg.reg[[rownames(subg.reg.data.05)[3]]], studlab = TRUE)


# ----------------- enrichment analysis - clusterProfiler ----------------------

# GO biological process
ans.go <- enrichGO(gene = rownames(res.meta.df.05), 
                   keyType = "ENSEMBL", ont = "BP",
                   OrgDb = "org.Hs.eg.db", 
                   universe = rownames(res.meta.df), 
                   readable = TRUE, 
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05)
tab.go <- as.data.frame(ans.go)
head(tab.go)

# plot
p1 <- barplot(ans.go, showCategory = 25)
p1
p1.2 <- dotplot(ans.go, showCategory = 25)
p1.2
p1.3 <- upsetplot(ans.go, n = 10)
p1.3

# GO cellular component
ans.go2 <- enrichGO(gene = rownames(res.meta.df.05), keyType = "ENSEMBL", 
    ont = "CC", OrgDb = "org.Hs.eg.db", universe = rownames(res.meta.df), 
    readable = TRUE, pAdjustMethod = "BH", pvalueCutoff = 0.05)
tab.go2 <- as.data.frame(ans.go2)
head(tab.go2)

# plot
p2 <- barplot(ans.go2, showCategory = 10)
p2
p2.2 <- dotplot(ans.go2, showCategory = 25)
p2.2
p2.3 <- upsetplot(ans.go2, n = 10)
p2.3


# subgroup analysis results - GO terms
# no results
ans.go.subg <- enrichGO(gene = rownames(diff.genes.btw),
                        keyType = "ENSEMBL", ont = "MF",
                        OrgDb = "org.Hs.eg.db", 
                        universe = rownames(res.meta.df), 
                        readable = TRUE, pAdjustMethod = "BH",
                        pvalueCutoff = 0.05)
tab.go.subg <- as.data.frame(ans.go.subg)
head(tab.go.subg)


# ------------------------------- TFEA.ChIP ------------------------------------

tfea.data <- res.meta.df[, c("lfc", "pval", "adj.pval")]
list.genes.entrez <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id"), 
                           values = rownames(tfea.data), mart = mart)
tfea.data <- merge(list.genes.entrez, tfea.data, by.x = "ensembl_gene_id", 
                   by.y = 0, all.y = TRUE)
tfea.data <- tfea.data[order(tfea.data$lfc, decreasing = TRUE), ]
tfea.data <- tfea.data[, -1]
colnames(tfea.data) <- c("Genes", "log2FoldChange", "pvalue", "pval.adj")

# UP-regulated genes extract vector with names of upregulated genes
genes.upreg <- Select_genes(tfea.data, min_LFC = 1)

# extract vector with names of non-responsive genes
genes.control <- Select_genes(tfea.data, min_pval = 0.5, max_pval = 1, 
                              min_LFC = -0.25, max_LFC = 0.25)

## association analysis

# generates list of contingency tables, one per dataset
CM_list_UP <- contingency_matrix(genes.upreg) 

# generates list of p-values and OR from association test
pval_mat_UP <- getCMstats(CM_list_UP)  
head(pval_mat_UP)

TF_ranking_up <- rankTFs(pval_mat_UP, rankMethod = "gsea", makePlot = TRUE)
TF_ranking_up[["TFranking_plot"]]

head(TF_ranking_up[["TF_ranking"]])
tf.ranking.up <- na.omit(TF_ranking_up[["TF_ranking"]][
                                TF_ranking_up[["TF_ranking"]]$pVal < 0.05, ])

plot_CM(pval_mat_UP)  # plot p-values against ORs

tfs.up <- rownames(tf.ranking.up)
names(tfs.up) <- rownames(tf.ranking.up)
col.up <- sample(viridis::turbo(length(tfs.up)), length(tfs.up))

# plot p-values against ORs highlighting indicated TFs
plot_CM(pval_mat_UP, specialTF = tfs.up, TF_colors = col.up)  

# DOWN-regulated genes extract vector with names of upregulated genes
genes.downreg <- Select_genes(tfea.data, max_LFC = -1)
# extract vector with names of non-responsive genes
genes.control <- Select_genes(tfea.data, min_pval = 0.5, max_pval = 1, 
                              min_LFC = -0.25, max_LFC = 0.25)

## association analysis

# generates list of contingency tables, one per dataset
CM_list_DOWN <- contingency_matrix(genes.downreg) 

# generates list of p-values and OR from association test
pval_mat_DOWN <- getCMstats(CM_list_DOWN) 
head(pval_mat_DOWN)

TF_ranking_down <- rankTFs(pval_mat_DOWN, rankMethod = "gsea", makePlot = TRUE)
TF_ranking_down[["TFranking_plot"]]
tf.ranking.down <- na.omit(TF_ranking_down[["TF_ranking"]][
                                TF_ranking_down[["TF_ranking"]]$pVal < 0.05, ])

head(TF_ranking_down[["TF_ranking"]])

plot_CM(pval_mat_DOWN)  # plot p-values against ORs

tfs.down <- rownames(tf.ranking.down)
names(tfs.down) <- rownames(tf.ranking.down)
col.down <- sample(viridis::turbo(length(tfs.down)), length(tfs.down))

# plot p-values against ORs highlighting indicated TFs
plot_CM(pval_mat_DOWN, specialTF = tfs.down, TF_colors = col.down)  


# -------------------- actomyosin network analysis -----------------------------

## MA volcano plot

res.meta.df$diffexp <- "NO"
res.meta.df$diffexp[res.meta.df$lfc > 1 & res.meta.df$adj.pval < 0.05] <- "UP"
res.meta.df$diffexp[res.meta.df$lfc < -1 & 
                                         res.meta.df$adj.pval < 0.05] <- "DOWN"
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p <- ggplot(data = res.meta.df, 
            aes(x = lfc, y = -log10(adj.pval), col = diffexp)) +
    geom_point() + theme_minimal() + geom_vline(xintercept = c(-1, 1), 
                                                col = "red") +
    geom_hline(yintercept = -log10(0.05), col = "red") + 
  scale_colour_manual(values = mycolors) +
    theme(legend.position = "none", axis.title = element_text(size = 8), 
          axis.text = element_text(size = 8)) +
    ylab("-log10(adj. p-value)") + xlab("LFC")
p


## MA volcano plot selecting the actomyosin genes

# read list containing actomyosin genes
actomyosin <- read.table("data/actomyosin.txt")
actomyosin <- actomyosin$V1

# add symbols
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
list.genes <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), 
                    values = rownames(res.meta.df),
                    mart = mart)
dim(list.genes)
res.meta.df.sy <- merge(list.genes, res.meta.df, all.y = TRUE, 
                        by.x = "ensembl_gene_id", by.y = 0)

# add a column indicating if the gene is down, up or not regulated
res.meta.df.sy$actomy <- "NO"
res.meta.df.sy$actomy[(res.meta.df.sy$hgnc_symbol %in% actomyosin) & 
                        (res.meta.df.sy$lfc > 0)] <- "UP"
res.meta.df.sy$actomy[(res.meta.df.sy$hgnc_symbol %in% actomyosin) & 
                        (res.meta.df.sy$lfc < 0)] <- "DOWN"

# choose colors
mycolors <- c("darkgray", "blue", "red")
names(mycolors) <- c("DOWN", "UP", "NO")
res.meta.df.sy$actomy <- factor(res.meta.df.sy$actomy, 
                                levels = c("NO", "DOWN", "UP"))

# order data
res.meta.df.sy <- res.meta.df.sy[order(res.meta.df.sy$actomy), ]
res.meta.df.sy$delabel <- NA
res.meta.df.sy$delabel[
  res.meta.df.sy$actomy != "NO"] <- res.meta.df.sy$hgnc_symbol[
                                           res.meta.df.sy$actomy != "NO"]

# plot
p2 <- ggplot(data = res.meta.df.sy, aes(x = lfc, y = -log10(adj.pval), 
                                        col = actomy, label = delabel)) + 
  geom_point(shape = 21, size = 2.2, fill = mycolors[res.meta.df.sy$actomy]) +
  theme_light() + 
  scale_colour_manual(values = c("darkgray", "black", "black")) +
    labs(x = "log2(Fold Change)", y = "log10(Adj. p-value)") + 
  theme(legend.position = "none") +
    scale_x_continuous(breaks = seq(-8, 8, by = 1))
p2 + geom_label_repel(data = subset(res.meta.df.sy, adj.pval < 0.05 & 
                                  (lfc > 1.5 | lfc < -1.5)), max.overlaps = 100)

# zoom in removing two points with the smallest pvalue
res.meta.df.sy.volcano <- res.meta.df.sy[
                                  -c(which(res.meta.df.sy$adj.pval < 1e-45)), ]
p3 <- ggplot(data = res.meta.df.sy.volcano, aes(x = lfc, y = -log10(adj.pval), 
                                               col = actomy, label = delabel)) + 
  geom_point(shape = 21, size = 2.2, 
             fill = mycolors[res.meta.df.sy.volcano$actomy]) +
  theme_light() + 
  scale_colour_manual(values = c("darkgray", "black", "black")) +
  labs(x = "LFC", y = "-log10(p-valor ajustado)") + 
  theme(legend.position = "none",
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12)) + 
  scale_x_continuous(breaks = seq(-8, 8, by = 1))
p3 + xlim(c(-6, 6)) + 
  geom_label_repel(data = subset(res.meta.df.sy.volcano, adj.pval < 0.05 & 
                                   (lfc > 1.5 | lfc < -1.5)), 
                   max.overlaps = 100, 
                   position = position_dodge(width = 6), size = 2.5)


## forest plot of some genes

# MYL9
forest.meta(genes.meta.lines[['ENSG00000101335']], 
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2', 
            studlab = sub('treatment', 'MAPKi', 
                          genes.meta.lines[[1]]$data$studylab),
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# MYH7
forest.meta(genes.meta.lines[['ENSG00000092054']], 
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2', 
            studlab = sub('treatment', 'MAPKi', 
                          genes.meta.lines[[1]]$data$studylab),
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# MYO15B
forest.meta(genes.meta.lines[['ENSG00000266714']], 
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2', 
            studlab = sub('treatment', 'MAPKi', 
                          genes.meta.lines[[1]]$data$studylab),
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# MYOT
forest.meta(genes.meta.lines[['ENSG00000120729']], 
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2', 
            studlab = sub('treatment', 'MAPKi', 
                          genes.meta.lines[[1]]$data$studylab),
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# MYOZ1
forest.meta(genes.meta.lines[['ENSG00000177791']], 
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2', 
            studlab = sub('treatment', 'MAPKi', 
                          genes.meta.lines[[1]]$data$studylab),
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# ACTA2
forest.meta(genes.meta.lines[['ENSG00000107796']], 
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2', 
            studlab = sub('treatment', 'MAPKi', 
                          genes.meta.lines[[1]]$data$studylab),
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# ACTBL2
forest.meta(genes.meta.lines[['ENSG00000169067']], 
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2', 
            studlab = sub('treatment', 'MAPKi', 
                          genes.meta.lines[[1]]$data$studylab),
            print.I2 = FALSE,
            print.pval.Q = FALSE)

## actomyosin network heatmap
lfc.genes.sy <- merge(list.genes, lfc.genes[, -15], all.y = TRUE, 
                      by.x = "ensembl_gene_id", by.y = 0)
lfc.genes.sy.actomy <- lfc.genes.sy[which(lfc.genes.sy$hgnc_symbol %in% 
                                                             actomyosin), ]
lfc.genes.sy.actomy <- lfc.genes.sy.actomy[, -1]

dim(lfc.genes.sy.actomy)
rownames(lfc.genes.sy.actomy) <- lfc.genes.sy.actomy$hgnc_symbol
lfc.genes.sy.actomy <- lfc.genes.sy.actomy[, -1]
colnames(lfc.genes.sy.actomy) <- sub("treatment", "MAPKi", 
                                     colnames(lfc.genes.sy.actomy))

pheatmap(na.omit(lfc.genes.sy.actomy), col = col.hm, show_rownames = FALSE, 
         breaks = seq(-2, 2, length = 50), fontsize_row = 7, border_color = NA)


# ----------------------- cell cycle analysis ----------------------------------

# load genes
kegg.cell.cycle <- read.table("data/KEGG_CELL_CYCLE.v7.5.1.grp", )

# merge
kegg.cell.cycle.ensemblid <- merge(kegg.cell.cycle, list.genes, by.x = "V1", 
                                   by.y = "hgnc_symbol")
kegg.cell.cycle.ensemblid <-intersect(kegg.cell.cycle.ensemblid$ensembl_gene_id,
                                                            rownames(lfc.genes))
lfc.cell.cycle <- lfc.genes[kegg.cell.cycle.ensemblid, -15]
colnames(lfc.cell.cycle) <- sub("treatment", "MAPKi", colnames(lfc.cell.cycle))

# plot
col.hm <- colorRampPalette(c("green", "black", "red"))(50)
pheat.cell.cycle <- pheatmap(na.omit(lfc.cell.cycle), col = col.hm, 
                             show_rownames = FALSE, 
                             breaks = seq(-5, 5, length = 50), 
                             fontsize_row = 7, border_color = NA, fontsize = 15)
pheat.cell.cycle

# now select only the second half of the list to get a better heatmap
hcluster <- pheat.cell.cycle$tree_row
lbl <- cutree(hcluster, 2)
which(lbl == 2)
kegg.cell.cycle.ensemblid.clean <- names(lbl)[lbl == 2]

# plot
pheat.cell.cycle.clean <- pheatmap(na.omit(lfc.cell.cycle[
                                           kegg.cell.cycle.ensemblid.clean, ]), 
                                   col = col.hm, show_rownames = FALSE, 
                                   breaks = seq(-5, 5, length = 50),
                                   fontsize_row = 7, 
                                   border_color = NA, fontsize = 15)
pheat.cell.cycle.clean


################################################################################
#              3. METAANALYSES USING 24-48H TREATMENT DATA                     #
################################################################################

# These metaanalyses combine the effect size (LFC) of 24 or 48 h 
# treated samples (TE-MA)
# 6 and 8 h treatment data were excluded

# function to prepare the data to metaanalyze
search.gene.long <- function(gene) {
    mat <- matrix(ncol = 6, nrow = 0)
    df <- data.frame(mat)

    # obenauf A375 48h
    if (nrow(genes.obenauf.lines.A375.48h.lfc[
                  rownames(genes.obenauf.lines.A375.48h.lfc) == gene, ]) == 1) {
        a <- genes.obenauf.lines.A375.48h.lfc[
              rownames(genes.obenauf.lines.A375.48h.lfc) == gene, ][1, 1]
        b <- genes.obenauf.lines.A375.48h.lfc[
              rownames(genes.obenauf.lines.A375.48h.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Obenauf, A375 48h treatment", "Obenauf", "48h",
            "A375", "long"))
    } else {
        df <- rbind(df, c(NA, NA, "Obenauf, A375 48h treatment", "Obenauf", 
            "48h", "A375", "long"))
    }

    # obenauf Colo800 48h
    if (nrow(genes.obenauf.lines.colo.48h.lfc[
                  rownames(genes.obenauf.lines.colo.48h.lfc) == gene, ]) == 1) {
        a <- genes.obenauf.lines.colo.48h.lfc[
              rownames(genes.obenauf.lines.colo.48h.lfc) == gene, ][1, 1]
        b <- genes.obenauf.lines.colo.48h.lfc[
              rownames(genes.obenauf.lines.colo.48h.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Obenauf, Colo800 48h treatment", "Obenauf", 
            "48h", "Colo800", "long"))
    } else {
        df <- rbind(df, c(NA, NA, "Obenauf, Colo800 48h treatment", "Obenauf", 
            "48h", "Colo800", "long"))
    }

    # obenauf UACC62 48h
    if (nrow(genes.obenauf.lines.uacc.48h.lfc[
                  rownames(genes.obenauf.lines.uacc.48h.lfc) == gene, ]) == 1) {
        a <- genes.obenauf.lines.uacc.48h.lfc[
              rownames(genes.obenauf.lines.uacc.48h.lfc) == gene, ][1, 1]
        b <- genes.obenauf.lines.uacc.48h.lfc[
              rownames(genes.obenauf.lines.uacc.48h.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Obenauf, UACC62 48h treatment", "Obenauf", 
            "48h", "UACC62", "long"))
    } else {
        df <- rbind(df, c(NA, NA, "Obenauf, UACC62 48h treatment", "Obenauf", 
            "48h", "UACC62", "long"))
    }

    # fallahi Colo858 24h
    if (nrow(genes.fallahi.lines.colo.24h.lfc[
                  rownames(genes.fallahi.lines.colo.24h.lfc) == gene, ]) == 1) {
        a <- genes.fallahi.lines.colo.24h.lfc[
              rownames(genes.fallahi.lines.colo.24h.lfc) == gene, ][1, 1]
        b <- genes.fallahi.lines.colo.24h.lfc[
              rownames(genes.fallahi.lines.colo.24h.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Fallahi, Colo858 24h treatment", "Fallahi",
            "24h", "Colo858", "long"))
    } else {
        df <- rbind(df, c(NA, NA, "Fallahi, Colo858 24h treatment", "Fallahi", 
            "24h", "Colo858", "long"))
    }

    # fallahi Colo858 48h
    if (nrow(genes.fallahi.lines.colo.48h.lfc[
                 rownames(genes.fallahi.lines.colo.48h.lfc) == gene, ]) == 1) {
        a <- genes.fallahi.lines.colo.48h.lfc[
              rownames(genes.fallahi.lines.colo.48h.lfc) == gene, ][1, 1]
        b <- genes.fallahi.lines.colo.48h.lfc[
              rownames(genes.fallahi.lines.colo.48h.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Fallahi, Colo858 48h treatment", "Fallahi", 
            "48h", "Colo858", "long"))
    } else {
        df <- rbind(df, c(NA, NA, "Fallahi, Colo858 48h treatment", "Fallahi", 
            "48h", "Colo858", "long"))
    }

    # fallahi MMACSF 24h
    if (nrow(genes.fallahi.lines.mmacsf.24h.lfc[
                rownames(genes.fallahi.lines.mmacsf.24h.lfc) == gene, ]) == 1) {
        a <- genes.fallahi.lines.mmacsf.24h.lfc[
              rownames(genes.fallahi.lines.mmacsf.24h.lfc) == gene, ][1, 1]
        b <- genes.fallahi.lines.mmacsf.24h.lfc[
              rownames(genes.fallahi.lines.mmacsf.24h.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Fallahi, MMACSF 24h treatment", "Fallahi", 
                          "24h", "MMACSF", "long"))
    } else {
        df <- rbind(df, c(NA, NA, "Fallahi, MMACSF 24h treatment", "Fallahi", 
                          "24h", "MMACSF", "long"))
    }

    # fallahi MMACSF 48h
    if (nrow(genes.fallahi.lines.mmacsf.48h.lfc[
                rownames(genes.fallahi.lines.mmacsf.48h.lfc) == gene, ]) == 1) {
        a <- genes.fallahi.lines.mmacsf.48h.lfc[
              rownames(genes.fallahi.lines.mmacsf.48h.lfc) == gene, ][1, 1]
        b <- genes.fallahi.lines.mmacsf.48h.lfc[
              rownames(genes.fallahi.lines.mmacsf.48h.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Fallahi, MMACSF 48h treatment", "Fallahi", 
            "48h", "MMACSF", "long"))
    } else {
        df <- rbind(df, c(NA, NA, "Fallahi, MMACSF 48h treatment", "Fallahi", 
                          "48h", "MMACSF", "long"))
    }

    # corre 501Mel 48h
    if (nrow(genes.corre.lfc[rownames(genes.corre.lfc) == gene, ]) == 1) {
        a <- genes.corre.lfc[rownames(genes.corre.lfc) == gene, ][1, 1]
        b <- genes.corre.lfc[rownames(genes.corre.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Corre, 501Mel 48h treatment", "Corre", "48h", 
                          "501Mel", "long"))
    } else {
        df <- rbind(df, c(NA, NA, "Corre, 501Mel 48h treatment", "Corre", "48h",
            "501Mel", "long"))
    }

    # gerosa A375 24h
    if (nrow(genes.gerosa.lfc[rownames(genes.gerosa.lfc) == gene, ]) == 1) {
        a <- genes.gerosa.lfc[rownames(genes.gerosa.lfc) == gene, ][1, 1]
        b <- genes.gerosa.lfc[rownames(genes.gerosa.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Gerosa, A375 24h treatment", "Gerosa", "24h", 
            "A375", "long"))
    } else {
        df <- rbind(df, c(NA, NA, "Gerosa, A375 24h treatment", "Gerosa", "24h",
            "A375", "long"))
    }

    # song M229 48h
    if (nrow(genes.song.lines.M229.48h.lfc[
                    rownames(genes.song.lines.M229.48h.lfc) == gene, ]) == 1) {
        a <- genes.song.lines.M229.48h.lfc[
              rownames(genes.song.lines.M229.48h.lfc) == gene, ][1, 1]
        b <- genes.song.lines.M229.48h.lfc[
              rownames(genes.song.lines.M229.48h.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Song, M229 48h treatment", "Song", "48h", 
            "M229", "long"))
    } else {
        df <- rbind(df, c(NA, NA, "Song, M229 48h treatment", "Song", "48h", 
                          "M229", "long"))
    }

    # song M238 48h
    if (nrow(genes.song.lines.M238.48h.lfc[
                    rownames(genes.song.lines.M238.48h.lfc) == gene, ]) == 1) {
        a <- genes.song.lines.M238.48h.lfc[
              rownames(genes.song.lines.M238.48h.lfc) == gene, ][1, 1]
        b <- genes.song.lines.M238.48h.lfc[
              rownames(genes.song.lines.M238.48h.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Song, M238 48h treatment", "Song", "48h", 
            "M238", "long"))
    } else {
        df <- rbind(df, c(NA, NA, "Song, M238 48h treatment", "Song", "48h", 
            "M238", "long"))
    }

    # end of search
    colnames(df) <- c("logFC", "lfcSE", "studylab", "author", "hours", "line", 
                      "treatment")
    return(df)
}

# example of use
search.gene.long("ENSG00000000003.15")

# loop to iterate over all the genes using the function
data.meta.lines.long <- lapply(all.genes, 
                               FUN = function(x) {search.gene.long(x)})
names(data.meta.lines.long) <- all.genes

length(data.meta.lines.long)

# remove null results
for (g in all.genes) {
    if (all(is.na(data.meta.lines.long[[g]][, c(1, 2)]))) {
        data.meta.lines.long[[g]] <- NULL
    }
}
length(data.meta.lines.long)

# remove genes with info from only 1 study
for (g in names(data.meta.lines.long)) {
    if (sum(is.na(data.meta.lines.long[[g]])) == 20) {
        data.meta.lines.long[[g]] <- NULL
    }
}
length(data.meta.lines.long)

# remove digits after dots and save
names(data.meta.lines.long) <- as.vector(mapply(names(data.meta.lines.long), 
                                   FUN = function(x) {gsub("\\.\\d+", "", x)}))

# saveRDS(data.meta.lines.long, file = 'data-meta-lines-long.rds')


## metaanalyses
counter <- 0
total <- length(data.meta.lines.long)
genes.meta.lines.long <- lapply(names(data.meta.lines.long), FUN = function(x) {
    counter <<- counter + 1
    print(paste0("iteration number ", counter, "; total ", total))
    metaanalyses(data.meta.lines.long, x)
})

names(genes.meta.lines.long) <- names(data.meta.lines.long)

# remove null results and save
length(genes.meta.lines.long)

for (g in names(data.meta.lines.long)) {
    if (is.null(genes.meta.lines.long[[g]])) {
        genes.meta.lines.long[[g]] <- NULL
    }
}
length(genes.meta.lines.long)

# saveRDS(genes.meta.lines.long, file='genes-meta-lines-long.rds')

# summarize results in a data frame
res.meta.long <- lapply(names(genes.meta.lines.long), FUN = function(x) {
    results.meta(genes.meta.lines.long, x)
})
names(res.meta.long) <- names(genes.meta.lines.long)
res.meta.long.df <- do.call(rbind.data.frame, res.meta.long)
dim(res.meta.long.df)

# multiple testing correction: FDR/BH
res.meta.long.df$adj.pval <- p.adjust(res.meta.long.df$pval, method = "BH", 
                                      n = length(res.meta.long.df$pval))

# multiple testing correction: Bonferroni
res.meta.long.df$adj.pval.bonfer <- p.adjust(res.meta.long.df$pval, 
                                             method = "bonferroni",
                                             n = length(res.meta.long.df$pval))

# saveRDS(res.meta.long.df, file='res-meta-long-df.rds') res.meta.long.df <-


# ----------------------- metaanalyses results ---------------------------------

# meta-analysis and experiments correlation

## experiments

# save the LFC values of the experiments in a new variable
lfc.genes.long <- lapply(data.meta.lines.long, '[[', 'logFC')
lfc.genes.long <- do.call('rbind.data.frame', lfc.genes.long)
lfc.genes.long <- sapply(lfc.genes.long, as.numeric)
colnames(lfc.genes.long) <- data.meta.lines.long[[1]]$studylab
rownames(lfc.genes.long) <- names(data.meta.lines.long)
lfc.genes.long <- lfc.genes.long[names(genes.meta.lines.long), ]
dim(lfc.genes.long)


## metaanalyses

# save the LFC values of the metaanalyses in a new variable
lfc.meta.long <- lapply(genes.meta.lines.long, '[[', 'TE.random')
lfc.meta.long <- do.call('rbind.data.frame', lfc.meta.long)
lfc.meta.long <- sapply(lfc.meta.long, as.numeric)
colnames(lfc.meta.long) <- 'lfc.meta.long'
rownames(lfc.meta.long) <- names(genes.meta.lines.long)

# merge both
lfc.genes.long <- merge(lfc.genes.long, lfc.meta.long, by = 0)
rownames(lfc.genes.long) <- lfc.genes.long$Row.names
lfc.genes.long <- lfc.genes.long[ , -1]


# plot correlations
lfc.genes.long.cor <- lfc.genes.long
lfc.genes.long.cor[lfc.genes.long[, "lfc.meta.long"] < 0, "lfc.meta.long"] <- -1*
    lfc.genes.long.cor[lfc.genes.long[, "lfc.meta.long"] < 0, "lfc.meta.long"]

par(mfrow = c(3, 4))
for (i in 1:11) {
    lfc.genes.long.cor[lfc.genes.long[, "lfc.meta.long"] < 0, i] <- -1 * 
      lfc.genes.long.cor[lfc.genes.long[, "lfc.meta.long"] < 0, i]
    plot(lfc.genes.long.cor[, i], lfc.genes.long.cor[, "lfc.meta.long"], 
         main = colnames(lfc.genes.long.cor)[i],
         ylab = "meta lfc", xlab = "individual lfc", pch = 21, 
         col = "blue4", bg = "blue1", cex = 1, lwd = 1)
}
par(mfrow = c(1, 1))

# compute correlations
first.cor <- cor.test(lfc.genes.long.cor[, 1], 
                      lfc.genes.long.cor[, "lfc.meta.long"],
                      method = "pearson")
cor.res.long.all <- data.frame(cor = first.cor$estimate, 
                               pval = first.cor$p.value)

for (i in 2:11) {
    cor <- cor.test(lfc.genes.long.cor[, i], 
                    lfc.genes.long.cor[, "lfc.meta.long"], method = "pearson")
    cor.res.long.all <- rbind(cor.res.long.all, c(cor$estimate, cor$p.value))
}

dim(cor.res.long.all)
cor.res.long.all


# now take only DEGs in order to reduce noise in correlations and plots
dim(lfc.genes.long)
lfc.genes.long.05 <- lfc.genes.long[
             rownames(res.meta.long.df[res.meta.long.df$adj.pval <= 0.05, ]), ]
dim(lfc.genes.long.05)
lfc.genes.long.05.cor <- lfc.genes.long.05

# change downregulated DEGs sign first
for (i in 1:11) {
    lfc.genes.long.05.cor[lfc.genes.long.05.cor[,"lfc.meta.long"] < 0, i] <- -1*
        (lfc.genes.long.05.cor[lfc.genes.long.05.cor[, "lfc.meta.long"] < 0, i])
}

lfc.genes.long.05.cor[
  lfc.genes.long.05.cor[, "lfc.meta.long"] < 0, "lfc.meta.long"] <- -1 * 
  (lfc.genes.long.05.cor[
               lfc.genes.long.05.cor[, "lfc.meta.long"] < 0, "lfc.meta.long"])

# compute correlations
first.cor <- cor.test(lfc.genes.long.05.cor[, 1], 
                      lfc.genes.long.05.cor[, "lfc.meta.long"], 
                      method = "pearson")
cor.res.long <- data.frame(cor = first.cor$estimate, pval = first.cor$p.value)

for (i in 2:11) {
    cor <- cor.test(lfc.genes.long.05.cor[, i], 
                    lfc.genes.long.05.cor[, "lfc.meta.long"],
                    method = "pearson")
    cor.res.long <- rbind(cor.res.long, c(cor$estimate, cor$p.value))
}

dim(cor.res.long)
cor.res.long

# plot DEGs 
par(mfrow = c(3, 4))
for (i in 1:11) {
    plot(lfc.genes.long.05.cor[, i], lfc.genes.long.05.cor[, "lfc.meta.long"], 
         main = colnames(lfc.genes.long.05.cor)[i],
        ylab = "meta lfc", xlab = "individual lfc", pch = 21, 
        col = "darkorchid4", bg = "darkorchid2", cex = 1, lwd = 1)
}
par(mfrow = c(1, 1))

# plot using ggplot2
for (i in 1:11) {
    dd <- as.data.frame(lfc.genes.long.05.cor[, c(i, 12)])
    colnames(dd) <- c("LFC.individual", "LFC.MA")
    plot.cor <- ggplot(dd, aes(x = LFC.individual, y = LFC.MA)) + 
      geom_point(aes(fill = "A"), colour = "blue4", pch = 21, size = 2.5) + 
      scale_fill_manual(values = "blue1") +
      theme_bw() + ylab("LFC en meta-analisis") + 
      xlab("LFC individual") + ggtitle(sub("treatment", "MAPKi", 
                                          colnames(lfc.genes.long.05.cor)[i])) + 
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5))
    assign(paste0("p.cor.long.", i), plot.cor)
}

figure2 <- ggarrange(p.cor.long.1, p.cor.long.2, p.cor.long.3, p.cor.long.4, 
    p.cor.long.5, p.cor.long.6, p.cor.long.7, p.cor.long.8, p.cor.long.9, 
    p.cor.long.10, p.cor.long.11, ncol = 4, nrow = 3)

figure2

# correlation between the data two-by-two 
# ! these correlations were computed taking into account all the genes, 
# not only DEGs
cor.matrix.long.lfc <- cor(lfc.genes.long[, -12], use = "pairwise.complete.obs")
colnames(cor.matrix.long.lfc) <- sub("treatment", "MAPKi", 
                                                 colnames(cor.matrix.long.lfc))
rownames(cor.matrix.long.lfc) <- sub("treatment", "MAPKi", 
                                                 rownames(cor.matrix.long.lfc))

col4 <- colorRampPalette(c("yellow2", "white", "darkorchid4"))(50)
pheatmap(cor.matrix.long.lfc, col = col4, display_numbers = TRUE, 
         breaks = seq(-1, 1, length = 50), fontsize = 12)


## comparison with RankProd and RankSum

# load rankprod and ranksum results
rp.results.long.tables <- readRDS("data/rp-results-long-tables.rds")
rs.results.long.tables <- readRDS("data/rs-results-long-tables.rds")

RP.genes.long <- c(rownames(rp.results.long.tables$Table1), 
                     rownames(rp.results.long.tables$Table2))
RS.genes.long <- c(rownames(rs.results.long.tables$Table1), 
                     rownames(rs.results.long.tables$Table2))

# now take DEGs that were obtained with values of expression in all samples 
# in TE-MA
res.meta.long.df.05.11st <- res.meta.long.df[res.meta.long.df$num.data == 11 & 
                                          res.meta.long.df$adj.pval < 0.05, ]

# compare (venn diagram)
ggVennDiagram(list(ES_MA = rownames(res.meta.long.df.05.11st), 
                   RankProd_MA = RP.genes.long,
                   RankSum_MA = RS.genes.long), 
              color = "black", lwd = 0.8, lty = 1, label = "count",
              label_alpha = 0, 
              category.names = c("TE-MA", "RankProd", "RankSum"), 
              label_size = 7) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  scale_color_manual(values = rep("darkgray", 3)) + 
  theme_void() + theme(legend.position = "none")


# compute correlations using DEGs selected by the 3 MA
degs.11st <- intersect(rownames(res.meta.long.df.05.11st), 
                       intersect(RP.genes.long, RS.genes.long))
lfc.genes.long.05.cor.11st <- lfc.genes.long.05.cor[degs.11st, ]
dim(lfc.genes.long.05.cor.11st)

first.cor <- cor.test(lfc.genes.long.05.cor.11st[, 1], 
                      lfc.genes.long.05.cor.11st[, "lfc.meta.long"], 
                      method = "pearson")
cor.res.degs.11st <- data.frame(cor = first.cor$estimate, 
                                pval = first.cor$p.value)
for (i in 2:11) {
    cor <- cor.test(lfc.genes.long.05.cor.11st[, i], 
                    lfc.genes.long.05.cor.11st[, "lfc.meta.long"], 
                    method = "pearson")
    cor.res.degs.11st <- rbind(cor.res.degs.11st, c(cor$estimate, cor$p.value))
}
dim(cor.res.degs.11st)
cor.res.degs.11st


# plot DEGs: experiment vs. MA LFC
par(mfrow = c(3, 4))
for (i in 1:11) {
    plot(lfc.genes.long.05.cor.11st[, i], 
         lfc.genes.long.05.cor.11st[, "lfc.meta.long"],
         main = colnames(lfc.genes.long.05.cor.11st)[i], 
         ylab = "meta lfc", xlab = "individual lfc",
         pch = 21, col = "darkorchid4", bg = "darkorchid2", cex = 1, lwd = 1)
}
par(mfrow = c(1, 1))


# histogram to compare correlations obtained with 24-48h data vs. with 6-48h
cor.res.compare <- data.frame(cor = cor.res.degs.14st[, 1], 
                              group = rep("6-48 h"))
cor.res.compare <- rbind(cor.res.compare, 
                         data.frame(cor = cor.res.degs.11st[, 1],
                                    group = rep("24-48 h")))

p.hist.long <- ggplot(data = cor.res.compare, aes(x = cor, fill = group)) + 
  geom_histogram(aes(y = 0.15 * ..density..), color = "black", alpha = 0.5, 
                 position = "identity", binwidth = 0.15) +
    scale_fill_brewer(palette = "Set1") + xlim(c(0, 1)) + theme_light() + 
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15), 
        legend.title = element_blank(), legend.position = "top",
    legend.text = element_text(size = 12)) + 
  xlab("Pearson correlation coefficient") +
  ylab("Relative frequency")

p.hist.long


# select DEGs from TE-MA and from the experiments and compare medians

## metaanalyses
res.meta.long.df.05 <- res.meta.long.df[degs.11st, ]
dim(res.meta.long.df.05)[1]/dim(res.meta.long.df)[1] * 100
dim(res.meta.long.df.05)

# upregulated genes
sum(res.meta.long.df.05$lfc > 0)  
median(res.meta.long.df.05[res.meta.long.df.05$lfc > 0, ]$lfc)

# downregulated genes
sum(res.meta.long.df.05$lfc < 0)  
median(res.meta.long.df.05[res.meta.long.df.05$lfc < 0, ]$lfc)  

## experiments
degs.lfc.indiv.long <- merge(genes.obenauf.lines.A375.48h[
            genes.obenauf.lines.A375.48h$adj.P.Val<0.05, 'logFC', drop = FALSE], 
                             genes.obenauf.lines.colo.48h[
            genes.obenauf.lines.colo.48h$adj.P.Val<0.05, 'logFC', drop = FALSE],
             all = TRUE, by = 0) %>%
  merge(genes.obenauf.lines.uacc.48h[
          genes.obenauf.lines.uacc.48h$adj.P.Val < 0.05, 'logFC', drop = FALSE], 
           all = TRUE, by.x = 'Row.names', by.y = 0) %>%
  merge(genes.fallahi.lines.colo.24h[
          genes.fallahi.lines.colo.24h$adj.P.Val < 0.05, 'logFC', drop = FALSE], 
           all = TRUE, by.x = 'Row.names', by.y = 0) %>%
  merge(genes.fallahi.lines.colo.48h[
          genes.fallahi.lines.colo.48h$adj.P.Val < 0.05, 'logFC', drop = FALSE], 
           all = TRUE, by.x = 'Row.names', by.y = 0) %>%
  merge(genes.fallahi.lines.mmacsf.24h[
        genes.fallahi.lines.mmacsf.24h$adj.P.Val < 0.05, 'logFC', drop = FALSE], 
         all = TRUE, by.x = 'Row.names', by.y = 0) %>%
  merge(genes.fallahi.lines.mmacsf.48h[
        genes.fallahi.lines.mmacsf.48h$adj.P.Val < 0.05, 'logFC', drop = FALSE], 
         all = TRUE, by.x = 'Row.names', by.y = 0) %>%
  merge(genes.corre[genes.corre$adj.P.Val < 0.05, 'logFC', drop = FALSE], 
        all = TRUE, by.x = 'Row.names', by.y = 0) %>%
  merge(genes.gerosa[genes.gerosa$adj.P.Val < 0.05, 'logFC', drop = FALSE], 
        all = TRUE, by.x = 'Row.names', by.y = 0) %>%
  merge(genes.song.lines.M229.48h[
         genes.song.lines.M229.48h$adj.P.Val < 0.05, 'logFC', drop = FALSE], 
          all = TRUE, by.x = 'Row.names', by.y = 0) %>%
  merge(genes.song.lines.M238.48h[
         genes.song.lines.M238.48h$adj.P.Val<0.05, 'logFC', drop = FALSE], 
          all = TRUE, by.x = 'Row.names', by.y = 0) 

rownames(degs.lfc.indiv.long) <- degs.lfc.indiv.long$Row.names
degs.lfc.indiv.long <- degs.lfc.indiv.long[ , -1]
rownames(degs.lfc.indiv.long) <- as.vector(mapply(rownames(degs.lfc.indiv.long), 
                                      FUN=function(x){gsub("\\.\\d+", "", x)}))
dim(degs.lfc.indiv.long)

# compute mean of medians
md.up.long <- c()
md.down.long <- c()
for (i in 1:11) {
    md.up.long <- c(md.up.long, 
     median(degs.lfc.indiv.long[degs.lfc.indiv.long[, i] > 0, i], na.rm = TRUE))
    md.down.long <- c(md.down.long, 
     median(degs.lfc.indiv.long[degs.lfc.indiv.long[, i] < 0, i], na.rm = TRUE))
}
mean(md.up.long)
mean(md.down.long)

# plot distributions
p.degs.long.df <- rbind(data.frame(LFC = res.meta.long.df.05$lfc, 
                                   COND = rep('MA', 
                                              length(res.meta.long.df.05$lfc))),
                        data.frame(LFC = degs.lfc.indiv.long[ , 1], 
                                   COND = rep('EXP1', 
                                            length(degs.lfc.indiv.long[ , 1]))),
                        data.frame(LFC = degs.lfc.indiv.long[ , 2], 
                                   COND = rep('EXP2', 
                                            length(degs.lfc.indiv.long[ , 2]))),
                        data.frame(LFC = degs.lfc.indiv.long[ , 3], 
                                   COND = rep('EXP3', 
                                            length(degs.lfc.indiv.long[ , 3]))),
                        data.frame(LFC = degs.lfc.indiv.long[ , 4], 
                                   COND = rep('EXP4', 
                                            length(degs.lfc.indiv.long[ , 4]))),
                        data.frame(LFC = degs.lfc.indiv.long[ , 5], 
                                   COND = rep('EXP5', 
                                            length(degs.lfc.indiv.long[ , 5]))),
                        data.frame(LFC = degs.lfc.indiv.long[ , 6], 
                                   COND = rep('EXP6', 
                                            length(degs.lfc.indiv.long[ , 6]))),
                        data.frame(LFC = degs.lfc.indiv.long[ , 7], 
                                   COND = rep('EXP7', 
                                            length(degs.lfc.indiv.long[ , 7]))),
                        data.frame(LFC = degs.lfc.indiv.long[ , 8], 
                                   COND = rep('EXP8', 
                                            length(degs.lfc.indiv.long[ , 8]))),
                        data.frame(LFC = degs.lfc.indiv.long[ , 9], 
                                   COND = rep('EXP9', 
                                            length(degs.lfc.indiv.long[ , 9]))),
                        data.frame(LFC = degs.lfc.indiv.long[ , 10], 
                                   COND = rep('EXP10',
                                           length(degs.lfc.indiv.long[ , 10]))),
                        data.frame(LFC = degs.lfc.indiv.long[ , 11], 
                                   COND = rep('EXP11', 
                                           length(degs.lfc.indiv.long[ , 11]))))

p.degs.long.df$COND <- factor(p.degs.long.df$COND, 
                              levels = c("EXP1", "EXP2", "EXP3", "EXP4",  
                                         "EXP5",  "EXP6", "EXP7",  "EXP8", 
                                         "EXP9",  "EXP10", "EXP11", "MA"), 
                              ordered = TRUE)

p.degs.long.df$EXPRESSION <- NA
p.degs.long.df$EXPRESSION[p.degs.long.df$LFC > 0] <- "UP"
p.degs.long.df$EXPRESSION[p.degs.long.df$LFC < 0] <- "DOWN"

# plot
p.degs.long <- ggplot(subset(p.degs.long.df, !is.na(LFC)), 
                      aes(EXPRESSION, LFC, fill = COND)) +
  geom_boxplot(alpha = 1, outlier.shape = NA) + 
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12), 
        legend.position = 'none') +
  labs(y = "LFC",
       x = "") +
  ylim (c(-4, 4)) +
  scale_x_discrete(labels = c("DOWN", "UP")) +
  scale_fill_manual(values = c(rep('gray68', 11), 'limegreen')) +
  geom_hline(yintercept = c(mean(md.up.long), mean(md.down.long)),
             col = "red", lty = 5)

p.degs.long


# ----------------- enrichment analysis - clusterProfiler ----------------------

# GO BP with upregulated genes
# no results
go.up.genes <- rownames(res.meta.long.df.05[res.meta.long.df.05$lfc > 0, ])[
  rownames(res.meta.long.df.05[res.meta.long.df.05$lfc > 0, ]) %in% degs.11st]
ans.go3 <- enrichGO(gene = go.up.genes, 
                    keyType = "ENSEMBL", 
                    ont = "BP", 
                    OrgDb = "org.Hs.eg.db",
                    universe = rownames(res.meta.long.df), 
                    readable = TRUE, 
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05)

# GO CC with upregulated genes
ans.go3.2 <- enrichGO(gene = go.up.genes,
                      keyType = 'ENSEMBL',
                      ont = "CC",
                      OrgDb ="org.Hs.eg.db",
                      universe = rownames(res.meta.long.df),
                      readable=TRUE,
                      pAdjustMethod = 'BH',
                      pvalueCutoff = 0.05)

# plots
p1.long.2 <- barplot(ans.go3.2, showCategory = 20)
p1.long.2
p1.2.long.2 <- dotplot(ans.go3.2, showCategory = 20, font.size = 15)
p1.2.long.2

# GO BP with downregulated genes
go.down.genes <- rownames(res.meta.long.df.05[res.meta.long.df.05$lfc < 0, ])[
  rownames(res.meta.long.df.05[res.meta.long.df.05$lfc < 0, ]) %in% degs.11st]
ans.go4 <- enrichGO(gene = go.down.genes,
                    keyType = 'ENSEMBL',
                    ont = "BP",
                    OrgDb ="org.Hs.eg.db",
                    universe = rownames(res.meta.long.df),
                    readable=TRUE,
                    pAdjustMethod = 'BH',
                    pvalueCutoff = 0.05)

# plots
p2.long <- barplot(ans.go4, showCategory = 20)
p2.long
p2.2.long <- dotplot(ans.go4, showCategory = 20, font.size = 15)
p2.2.long

# GO CC with downregulated genes
ans.go4.2 <- enrichGO(gene = go.down.genes,
                      keyType = 'ENSEMBL',
                      ont = "cc",
                      OrgDb ="org.Hs.eg.db",
                      universe = rownames(res.meta.long.df),
                      readable=TRUE,
                      pAdjustMethod = 'BH',
                      pvalueCutoff = 0.05)

# plots
p2.long.2 <- barplot(ans.go4.2, showCategory = 20)
p2.long.2
p2.2.long.2 <- dotplot(ans.go4.2, showCategory = 20, font.size = 15)
p2.2.long.2


# --------------------------- GSEA analysis ------------------------------------

# plot results tft
tft.gsea.pos <- read_tsv("data/gsea_report_pos_tft.tsv")
tft.gsea.neg <- read_tsv("data/gsea_report_neg_tft.tsv")
tft.gsea.pos <- tft.gsea.pos[1:5, c(1, 4, 6, 8)]
tft.gsea.neg <- tft.gsea.neg[1:50, c(1, 4, 6, 8)]

# positive NES
tft.pos <- ggplot(tft.gsea.pos, aes(x = NES, y = `FDR q-val`, 
                                    color = `FDR q-val` < 0.05, label = NAME)) + 
  geom_point(size = tft.gsea.pos$SIZE/100) + theme_light() +
  scale_color_manual(values = c("black", "red")) + 
  geom_hline(yintercept = 0.05, col = "gray10", lty = 2) + ylab("") + 
  xlab("NES") + ylim(c(-0.02, 0.08)) + 
  theme(legend.position = "none", 
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17)) + 
  geom_label_repel(position = position_dodge(width = 6), 
                   size = 4, max.overlaps = 100)

tft.pos

# negative NES
tft.neg <- ggplot(tft.gsea.neg, aes(x = NES, y = `FDR q-val`, 
                                    color = `FDR q-val` < 0.05, label = 1:50)) + 
  geom_point(size = tft.gsea.neg$SIZE/100) + theme_light() + 
  scale_color_manual(values = c("black", "darkblue")) + 
  geom_hline(yintercept = 0.05, col = "gray10", lty = 2) + 
  ylab("FDR") + xlab("NES") + ylim(c(-0.02, 0.08)) + 
  theme(legend.position = "none", axis.title = element_text(size = 17),
        axis.text = element_text(size = 17)) + 
  geom_label_repel(position = position_dodge(width = 6), size = 4,
                   max.overlaps = 100)

tft.neg

# plot
tft.plot <- ggarrange(tft.neg, tft.pos, ncol = 2, nrow = 1)

tft.plot

# save terms that are not shown 
write.table(paste(1:50, tft.gsea.neg$NAME[1:50], sep = ". "), 
            quote = FALSE, file = "gsea-tft-treat-list.txt",
            col.names = FALSE, row.names = FALSE)


# plot results hallmarks
hallm.gsea.pos <- read_tsv("data/gsea_report_pos_hallmarks.tsv")
hallm.gsea.neg <- read_tsv("data/gsea_report_neg_hallmarks.tsv")
hallm.gsea.pos <- hallm.gsea.pos[1:5, c(1, 4, 6, 8)]
hallm.gsea.neg <- hallm.gsea.neg[1:25, c(1, 4, 6, 8)]

# positive NES
hallm.pos <- ggplot(hallm.gsea.pos, aes(x = NES, y = `FDR q-val`, 
                                        color = `FDR q-val` < 0.05, 
                                        label = sub("HALLMARK_", "", NAME))) +
  geom_point(size = hallm.gsea.pos$SIZE/100) +
  theme_light() + scale_color_manual(values = c("black", "red")) + 
  geom_hline(yintercept = 0.05, col = "gray10", lty = 2) + ylab("") + 
  xlab("Puntuaci?n normalizada") + ylim(c(-0.02, 0.125)) + 
  theme(legend.position = "none", axis.title = element_text(size = 17), 
        axis.text = element_text(size = 17)) + 
  geom_label_repel(position = position_dodge(width = 6), size = 4, 
                   max.overlaps = 100)
hallm.pos

# negative NES
hallm.neg <- ggplot(hallm.gsea.neg, aes(x = NES, y = `FDR q-val`,
                                        color = `FDR q-val` < 0.05, 
                                        label = 1:25)) + 
  geom_point(size = hallm.gsea.neg$SIZE/100) + theme_light() +
  scale_color_manual(values = c("black", "darkblue")) + 
  geom_hline(yintercept = 0.05, col = "gray10", lty = 2) + ylab("FDR") + 
  xlab("NES") + ylim(c(-0.02, 0.125)) + 
  theme(legend.position = "none", axis.title = element_text(size = 17),
        axis.text = element_text(size = 17)) + 
  geom_label_repel(position = position_dodge(width = 6),
    size = 4, max.overlaps = 100)

hallm.neg

# plot
hallm.plot <- ggarrange(hallm.neg, hallm.pos, ncol = 2, nrow = 1)

hallm.plot

# save terms that are not shown 
write.table(paste(1:25, sub("HALLMARK_", "", hallm.gsea.neg$NAME[1:25]), 
                                                     sep = ". "), quote = FALSE, 
            file = "gsea-hallm-treat-list.txt", 
            col.names = FALSE, row.names = FALSE)


# plot results kegg
kegg.gsea.pos <- read_tsv("data/gsea_report_pos_kegg.tsv")
kegg.gsea.neg <- read_tsv("data/gsea_report_neg_kegg.tsv")
kegg.gsea.pos <- kegg.gsea.pos[1:5, c(1, 4, 6, 8)]
kegg.gsea.neg <- kegg.gsea.neg[1:20, c(1, 4, 6, 8)]

# positive NES
kegg.pos <- ggplot(kegg.gsea.pos, aes(x = NES, y = `FDR q-val`,
                                      color = `FDR q-val` < 0.05, 
                                      label = sub("KEGG_", "", NAME))) + 
  geom_point(size = kegg.gsea.pos$SIZE/100) + theme_light() + 
  scale_color_manual(values = c("black", "red")) + 
  geom_hline(yintercept = 0.05, col = "gray10", lty = 2) + ylab("") + 
  xlab("NES") + ylim(c(-0.02, 0.16)) + 
  theme(legend.position = "none", axis.title = element_text(size = 17),
        axis.text = element_text(size = 17)) + 
  geom_label_repel(position = position_dodge(width = 6), 
                   size = 4, max.overlaps = 100)
kegg.pos


# negative NES
kegg.neg <- ggplot(kegg.gsea.neg, aes(x = NES, y = `FDR q-val`, 
                                      color = `FDR q-val` < 0.05, 
                                      label = 1:20)) + 
  geom_point(size = kegg.gsea.neg$SIZE/100) + theme_light() +
  scale_color_manual(values = c("black", "darkblue")) + 
  geom_hline(yintercept = 0.05, col = "gray10", lty = 2) + ylab("FDR") + 
  xlab("NES") + ylim(c(-0.02, 0.16)) + 
  theme(legend.position = "none", axis.title = element_text(size = 17), 
        axis.text = element_text(size = 17)) + 
  geom_label_repel(position = position_dodge(width = 6),
    size = 4, max.overlaps = 100)

kegg.neg

# plot
kegg.plot <- ggarrange(kegg.neg, kegg.pos, ncol = 2, nrow = 1)

kegg.plot

# save terms that are not shown 
write.table(paste(1:20, sub("KEGG_", "", kegg.gsea.neg$NAME[1:20]), sep = ". "),
    quote = FALSE, file = "gsea-kegg-treat-list.txt", 
    col.names = FALSE, row.names = FALSE)


# ---------------------------------- TFEA.ChIP ---------------------------------

# select the data to be used
tfea.data.long <- res.meta.long.df[, c("lfc", "pval", "adj.pval")]
list.genes.entrez <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id"), 
                           values = rownames(tfea.data), mart = mart)
tfea.data.long <- merge(list.genes.entrez, tfea.data.long, 
                        by.x = "ensembl_gene_id", by.y = 0, all.y = TRUE)
tfea.data.long <- tfea.data.long[order(tfea.data.long$lfc, decreasing = TRUE), ]
tfea.data.long <- tfea.data.long[, -1]
colnames(tfea.data.long) <- c("Genes", "log2FoldChange", "pvalue", "pval.adj")

# UP-regulated genes extract vector with names of upregulated genes
genes.upreg.long <- Select_genes(tfea.data.long, min_LFC = 1)

# extract vector with names of non-responsive genes
genes.control.long <- Select_genes(tfea.data.long, min_pval = 0.5, max_pval = 1,
    min_LFC = -0.25, max_LFC = 0.25)

## association analysis

# generates list of contingency tables, one per dataset
CM_list_UP_long <- contingency_matrix(genes.upreg.long)  
# generates list of p-values and OR from association test
pval_mat_UP_long <- getCMstats(CM_list_UP_long)  
head(pval_mat_UP_long)

TF_ranking_up_long <- rankTFs(pval_mat_UP_long, rankMethod = "gsea", 
                              makePlot = TRUE)
TF_ranking_up_long[["TFranking_plot"]]

head(TF_ranking_up_long[["TF_ranking"]])
tf.ranking.up.long <- na.omit(TF_ranking_up_long[["TF_ranking"]][
                             TF_ranking_up_long[["TF_ranking"]]$pVal < 0.05, ])

plot_CM(pval_mat_UP_long)  # plot p-values against ORs

tfs.up.long <- rownames(tf.ranking.up.long)
names(tfs.up.long) <- rownames(tf.ranking.up.long)
col.up.long <- sample(viridis::turbo(length(tfs.up.long)), length(tfs.up.long))

# plot p-values against ORs highlighting indicated TFs
plot_CM(pval_mat_UP_long, specialTF = tfs.up.long, TF_colors = col.up.long)  
TF_ranking_up_long[["TFranking_plot"]]


# DOWN-regulated genes extract vector with names of upregulated genes
genes.downreg.long <- Select_genes(tfea.data.long, max_LFC = -1)
# extract vector with names of non-responsive genes
genes.control.long <- Select_genes(tfea.data.long, min_pval = 0.5, max_pval = 1,
    min_LFC = -0.25, max_LFC = 0.25)

## association analysis

# generates list of contingency tables, one per dataset
CM_list_DOWN_long <- contingency_matrix(genes.downreg.long)  
# generates list of p-values and OR from association test
pval_mat_DOWN_long <- getCMstats(CM_list_DOWN_long)  
head(pval_mat_DOWN_long)

TF_ranking_down_long <- rankTFs(pval_mat_DOWN_long, rankMethod = "gsea", 
                                makePlot = TRUE)
TF_ranking_down_long[["TFranking_plot"]]
tf.ranking.down.long <- na.omit(TF_ranking_down_long[["TF_ranking"]][
  TF_ranking_down_long[["TF_ranking"]]$pVal < 0.05, ])

head(TF_ranking_down_long[["TF_ranking"]])

plot_CM(pval_mat_DOWN_long)  # plot p-values against ORs

tfs.down.long <- rownames(tf.ranking.down.long)
names(tfs.down.long) <- rownames(tf.ranking.down.long)
col.down.long <- sample(viridis::turbo(length(tfs.down.long)), 
                        length(tfs.down.long))

# plot p-values against ORs highlighting indicated TFs
plot_CM(pval_mat_DOWN_long, specialTF = tfs.down.long, 
        TF_colors = col.down.long)  
TF_ranking_down_long[["TFranking_plot"]]


# ------------------------- actomyosin network analysis ------------------------

## MA volcano plot
res.meta.long.df$diffexp <- "NO"
res.meta.long.df$diffexp[res.meta.long.df$lfc > 1 &
                                       res.meta.long.df$adj.pval < 0.05] <- "UP"
res.meta.long.df$diffexp[res.meta.long.df$lfc < -1 & 
                                     res.meta.long.df$adj.pval < 0.05] <- "DOWN"
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")

p.long <- ggplot(data = res.meta.long.df, 
                 aes(x = lfc, y = -log10(adj.pval), fill = diffexp)) +
    geom_point(shape = 21, size = 2.2) + theme_minimal() + 
  geom_vline(xintercept = c(-1, 1), col = "darkgray", lty = 2) + 
  geom_hline(yintercept = -log10(0.05), col = "darkgray", lty = 2) + 
scale_fill_manual(values = c("blue", "black", "red")) + 
  theme(legend.position = "none", axis.title = element_text(size = 15), 
        axis.text = element_text(size = 15)) + 
  ylab("-log10(adj p-value)") + xlab("LFC")
p.long


# MA volcano plot selecting the actomyosin genes
res.meta.long.df.sy <- merge(list.genes, res.meta.long.df, all.y = TRUE, 
                             by.x = "ensembl_gene_id", by.y = 0)

# add column indicating if the gene is down, up or not regulated
res.meta.long.df.sy$actomy <- "NO"
res.meta.long.df.sy$actomy[(res.meta.long.df.sy$hgnc_symbol %in% 
                                                          actomyosin)] <- "DEG"
res.meta.long.df.sy$actomy[(res.meta.long.df.sy$hgnc_symbol %in% actomyosin) & 
                             (res.meta.long.df.sy$lfc > 1.5) & 
                             (res.meta.long.df.sy$adj.pval < 0.05)] <- "UP"
res.meta.long.df.sy$actomy[(res.meta.long.df.sy$hgnc_symbol %in% actomyosin) & 
                             (res.meta.long.df.sy$lfc < -1.5) & 
                             (res.meta.long.df.sy$adj.pval < 0.05)] <- "DOWN"

# choose colors
mycolors <- c("darkgray", "blue", "red", "gray40")
names(mycolors) <- c("DOWN", "UP", "NO", "DEG")
res.meta.long.df.sy$actomy <- factor(res.meta.long.df.sy$actomy, 
                                     levels = c("NO", "DOWN", "UP", "DEG"))

# order data
res.meta.long.df.sy <- res.meta.long.df.sy[order(res.meta.long.df.sy$actomy), ]
res.meta.long.df.sy$delabel <- NA
res.meta.long.df.sy$delabel[
  res.meta.long.df.sy$actomy != "NO"] <- res.meta.long.df.sy$hgnc_symbol[
                                             res.meta.long.df.sy$actomy != "NO"]

# plot
p2.long <- ggplot(data = res.meta.long.df.sy, aes(x = lfc, y = -log10(adj.pval),
    col = actomy, label = delabel)) + 
  geom_point(shape = 21, size = 2.2, 
             fill = mycolors[res.meta.long.df.sy$actomy]) +
    theme_light() + scale_colour_manual(values = c("darkgray", "black", "black",
    "black")) + labs(x = "LFC", y = "log10(p-valor ajustado)") + 
  theme(legend.position = "none", axis.title = element_text(size = 15), 
        axis.text = element_text(size = 15)) +
    scale_x_continuous(breaks = seq(-8, 8, by = 1)) + 
  geom_vline(xintercept = c(-1.5, 1.5), col = "gray15", lty = 2) + 
  geom_hline(yintercept = -log10(0.05), col = "gray15", lty = 2) + 
  geom_label_repel(data = subset(res.meta.long.df.sy, adj.pval < 0.05 &
                                  (lfc > 1.5 | lfc < -1.5)), max.overlaps = 100, 
                   position = position_nudge_repel(y = 30), size = 2.7)

p2.long


## forest plots of some genes

# MYL9
forest.meta(genes.meta.lines.long[['ENSG00000101335']], 
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2', 
            studlab = sub('treatment', 'MAPKi', 
                          genes.meta.lines.long[[1]]$data$studylab),
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# MYH7
forest.meta(genes.meta.lines.long[['ENSG00000092054']], 
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2', 
            studlab = sub('treatment', 'MAPKi', 
                          genes.meta.lines.long[[1]]$data$studylab),
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# MYOZ
forest.meta(genes.meta.lines.long[['ENSG00000177791']], 
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2', 
            studlab = sub('treatment', 'MAPKi', 
                          genes.meta.lines.long[[1]]$data$studylab),
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# MYOT
forest.meta(genes.meta.lines.long[['ENSG00000120729']], 
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2', 
            studlab = sub('treatment', 'MAPKi', 
                          genes.meta.lines.long[[1]]$data$studylab),
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# LMNB1
forest.meta(genes.meta.lines.long[['ENSG00000113368']], 
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2', 
            studlab = sub('treatment', 'MAPKi', 
                          genes.meta.lines.long[[1]]$data$studylab),
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# CIT
forest.meta(genes.meta.lines.long[['ENSG00000122966']], 
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2', 
            studlab = sub('treatment', 'MAPKi', 
                          genes.meta.lines.long[[1]]$data$studylab),
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# ITGA9
forest.meta(genes.meta.lines.long[['ENSG00000144668']], 
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2', 
            studlab = sub('treatment', 'MAPKi', 
                          genes.meta.lines.long[[1]]$data$studylab),
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# ITGB3
forest.meta(genes.meta.lines.long[['ENSG00000259207']], 
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2', 
            studlab = sub('treatment', 'MAPKi', 
                          genes.meta.lines.long[[1]]$data$studylab),
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# ACTA2
forest.meta(genes.meta.lines.long[['ENSG00000107796']], 
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2', 
            studlab = sub('treatment', 'MAPKi', 
                          genes.meta.lines.long[[1]]$data$studylab),
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# ACTBL2
forest.meta(genes.meta.lines.long[['ENSG00000169067']], 
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2', 
            studlab = sub('treatment', 'MAPKi', 
                          genes.meta.lines.long[[1]]$data$studylab),
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# TMOD2
forest.meta(genes.meta.lines.long[['ENSG00000128872']], 
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2', 
            studlab = sub('treatment', 'MAPKi', 
                          genes.meta.lines.long[[1]]$data$studylab),
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# ANK3
forest.meta(genes.meta.lines.long[['ENSG00000151150']], 
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2', 
            studlab = sub('treatment', 'MAPKi', 
                          genes.meta.lines.long[[1]]$data$studylab),
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# MYO15B
forest.meta(genes.meta.lines.long[['ENSG00000266714']], 
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2', 
            studlab = sub('treatment', 'MAPKi', 
                          genes.meta.lines.long[[1]]$data$studylab),
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# PPP1R14C
forest.meta(genes.meta.lines.long[['ENSG00000198729']], 
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2', 
            studlab = sub('treatment', 'MAPKi', 
                          genes.meta.lines.long[[1]]$data$studylab),
            print.I2 = FALSE,
            print.pval.Q = FALSE)


## actomyosin network heatmap
lfc.genes.long.sy <- merge(list.genes, lfc.genes.long[, -12], all.y = TRUE, 
                           by.x = "ensembl_gene_id", by.y = 0)
lfc.genes.long.sy.actomy <- lfc.genes.long.sy[
                         which(lfc.genes.long.sy$hgnc_symbol %in% actomyosin), ]

lfc.genes.long.sy.actomy <- lfc.genes.long.sy.actomy[, -1]
dim(lfc.genes.long.sy.actomy)
rownames(lfc.genes.long.sy.actomy) <- lfc.genes.long.sy.actomy$hgnc_symbol
lfc.genes.long.sy.actomy <- lfc.genes.long.sy.actomy[, -1]
colnames(lfc.genes.long.sy.actomy) <- sub("treatment", "MAPKi", 
                                          colnames(lfc.genes.long.sy.actomy))

pheatmap(na.omit(lfc.genes.long.sy.actomy), col = col.hm, show_rownames = TRUE, 
         breaks = seq(-4, 4, length = 50), fontsize_col = 15, border_color = NA)

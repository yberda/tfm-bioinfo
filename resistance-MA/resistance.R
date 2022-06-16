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
library("GSVA")

library("ggplot2")
library("ggfortify")
library("RColorBrewer")
library("pheatmap")
library("corrplot")
library("ggVennDiagram")
library("egg")
library("gridExtra")
library("grid")
library("plotly")
library("ggrepel")


setwd("./tfm-bioinfo/resistance-MA")

# load salmon output (gse objects)
load("./salmon-output-resistance.RData")

# in case of running the code, load the data above and skip the data loading 
# lines ("load transcriptomes" section) up to the definition of gse objects
# (excluding the cell lines and the P/R definition)

# ---------------------------- plot functions ----------------------------------

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


# volcano plot using limma results
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
    p <- ggplot(data = genes, aes(x = logFC, 
                                  y = -log10(adj.P.Val), col = diffexpressed)) +
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
      theme(plot.title = element_text(hjust = 0.5,
        size = 20), axis.title = element_text(size = 17), 
        axis.text = element_text(size = 17))
    return(p3)
}


# This script contains the code used to do the following:

# 1. Data preparation: data are loaded and processed to be used for 
#    metaanalyses.
# 2. Resistance metaanalyses: this section includes the effect sizes 
#    metaanalyses, the exploration of the results (including the RankProd and 
#    RankSum comparison and clusterProfiler, a resistance signature and cell  
#    cycle analysis) and a subgroup analysis.
# 3. Double-drug resistant lines analysis and comparison
# 4. Treatment vs. resistance metaanalyses: this section compares the results 
#    obtained in the metaanalyses of resistance and parental cells treatment.
# 5. Metaanalyses separating subgroups (Subgroup A and B): this section 
#    includes the metaanalyses for each subgrup and the exploration of the  
#    results (including clusterProfiler, GSEA, TFEA.ChIP and actomyosin 
#    network analysis)



################################################################################
#                            1. DATA PREPARATION                               #
################################################################################

# -------------------------- load transcriptomes -------------------------------

## singleton ----

# load data
files.singleton <- read.table("../../files/singleton.txt")
files.singleton <- file.path(unlist(files.singleton))
all(file.exists(files.singleton))
tdata.singleton <- as.data.frame(files.singleton)
colnames(tdata.singleton) <- "files"
cell.lines.singleton <- rep("A375", 6)
tdata.singleton$cell <- factor(cell.lines.singleton)
treat.singleton <- relevel(factor(c(rep("P", 3), rep("R", 3)), 
                                  levels = c("P", "R")),
                           ref = "P")
tdata.singleton$treatment <- treat.singleton
tdata.singleton$names <- paste("Sample", 1:6, sep = "")

# reference transcriptome checksum recognized by tximeta
se.singleton <- tximeta(tdata.singleton)  

# summarize quantifications to gene-level
gse.singleton <- summarizeToGene(se.singleton, 
                                 countsFromAbundance = "lengthScaledTPM")  
y.singleton <- makeDGEList(gse.singleton)

# filter null and low counts
dim(y.singleton)
design.singleton <- model.matrix(~treat.singleton)
keep.singleton <- filterByExpr(y.singleton, design.singleton)
dge.singleton <- y.singleton[keep.singleton, , keep.lib.sizes = FALSE]
dim(dge.singleton)[1]/dim(y.singleton)[1] * 100

# plots
plot.filt(y.singleton, dge.singleton, c(0, 0.8))
plotMDS(dge.singleton, col = as.numeric(as.factor(treat.singleton)))

# differential expression analysis
v.singleton <- voom(dge.singleton, design.singleton, plot = TRUE)

fit.singleton <- lmFit(v.singleton, design.singleton)
fit.singleton <- eBayes(fit.singleton, trend = FALSE)

# A375 genes
genes.singleton <- topTable(fit.singleton, coef = 2, adjust.method = "BH", 
                            sort.by = "p", number = dim(fit.singleton)[1])
sum(genes.singleton$adj.P.Val < 0.05)  
sum(genes.singleton$adj.P.Val < 0.001)  

# volcano plot
volcanoplot.colors.limma(genes = genes.singleton, adjpval = 0.05, logfc = 1.5, 
                         title = "singleton, A375R")


## ho ----

# load data
files.ho <- read.table("../../files/ho.txt")
files.ho <- file.path(unlist(files.ho))
all(file.exists(files.ho))
tdata.ho <- as.data.frame(files.ho)
colnames(tdata.ho) <- "files"
cell.lines.ho <- c(rep("451Lu", 4), rep("A375", 4))
tdata.ho$cell <- factor(cell.lines.ho)
treat.ho <- relevel(factor(rep(c(rep("P", 2), rep("R", 2)), 2), 
                           levels = c("P", "R")), ref = "P")
tdata.ho$treatment <- treat.ho
tdata.ho$names <- paste("Sample", 1:8, sep = "")

# reference transcriptome checksum recognized by tximeta
se.ho <- tximeta(tdata.ho)  

# summarize quantifications to gene-level
gse.ho <- summarizeToGene(se.ho, countsFromAbundance = "lengthScaledTPM")  
y.ho <- makeDGEList(gse.ho)

# filter null and low counts
dim(y.ho)
group.ho <- paste(cell.lines.ho, treat.ho, sep = ".")
design.ho <- model.matrix(~0 + group.ho)
colnames(design.ho) <- c("Lu.P", "Lu.R", "A375.P", "A375.R")
keep.ho <- filterByExpr(y.ho, design.ho)
dge.ho <- y.ho[keep.ho, , keep.lib.sizes = FALSE]
dim(dge.ho)[1]/dim(y.ho)[1] * 100

# plots
plot.filt(y.ho, dge.ho, c(0, 0.85))
plotMDS(dge.ho, col = as.numeric(as.factor(treat.ho)))

# differential expression analysis
contrast.ho <- makeContrasts(Lu = Lu.R - Lu.P, A375 = A375.R - A375.P, 
                             levels = design.ho)

v.ho <- voom(dge.ho, design.ho, plot = TRUE)

fit.ho <- lmFit(v.ho, design.ho)
fit.ho <- contrasts.fit(fit.ho, contrast.ho)
fit.ho <- eBayes(fit.ho, trend = FALSE)

# 451Lu genes
genes.ho.lu <- topTable(fit.ho, coef = "Lu", adjust.method = "BH", 
                        sort.by = "p", number = dim(fit.ho)[1])
sum(genes.ho.lu$adj.P.Val < 0.05)
sum(genes.ho.lu$adj.P.Val < 0.001)

# A375 genes
genes.ho.A375 <- topTable(fit.ho, coef = "A375", adjust.method = "BH", 
                          sort.by = "p", number = dim(fit.ho)[1])
sum(genes.ho.A375$adj.P.Val < 0.05)
sum(genes.ho.A375$adj.P.Val < 0.001)

# volcano plots
volcanoplot.colors.limma(genes = genes.ho.lu, adjpval = 0.05, logfc = 1.5, 
                         title = "Ho, 451LuR")
volcanoplot.colors.limma(genes = genes.ho.A375, adjpval = 0.05, logfc = 1.5, 
                         title = "Ho, A375R")


# berico ----

# load data
files.berico <- read.table("../../files/berico.txt")
files.berico <- file.path(unlist(files.berico))
all(file.exists(files.berico))
tdata.berico <- as.data.frame(files.berico)
colnames(tdata.berico) <- "files"
cell.lines.berico <- rep("MM074", 4)
tdata.berico$cell <- factor(cell.lines.berico)
treat.berico <- relevel(factor(c(rep("P", 2), rep("R", 2)), 
                               levels = c("P", "R")), ref = "P")
tdata.berico$treatment <- treat.berico
tdata.berico$names <- paste("Sample", 1:4, sep = "")

# reference transcriptome checksum recognized by tximeta
se.berico <- tximeta(tdata.berico)  

# summarize quantifications to gene-level
gse.berico <- summarizeToGene(se.berico, 
                              countsFromAbundance = "lengthScaledTPM")  
y.berico <- makeDGEList(gse.berico)

# filter null and low counts
dim(y.berico)
design.berico <- model.matrix(~ treat.berico)
keep.berico <- filterByExpr(y.berico, design.berico)
dge.berico <- y.berico[keep.berico, , keep.lib.sizes = FALSE]
dim(dge.berico)[1]/dim(y.berico)[1] * 100 

# plots
plot.filt(y.berico, dge.berico, c(0, 0.8))
plotMDS(dge.berico, col = as.numeric(as.factor(treat.berico)))

# differential expression analysis
v.berico <- voom(dge.berico, design.berico, plot = TRUE)

fit.berico <- lmFit(v.berico, design.berico)
fit.berico <- eBayes(fit.berico, trend = FALSE)

# MM074R genes
genes.berico <- topTable(fit.berico, coef = 2, adjust.method = "BH", 
                         sort.by = "p", number = dim(fit.berico)[1])
sum(genes.berico$adj.P.Val < 0.05)
sum(genes.berico$adj.P.Val < 0.001)

# volcano plot
volcanoplot.colors.limma(genes = genes.berico, adjpval = 0.05, logfc = 1.5, 
                         title = "berico, MM074R")


## song ----

# load data
files.song <- read.table("../../files/song_resistance.txt")
files.hugo <- read.table("../../files/hugo_resistance.txt")  # same lab
files.song <- file.path(unlist(files.song)[-c(5, 6, 15, 18)])
files.hugo <- file.path(unlist(files.hugo)[-c(1, 6)])
files.song <- c(files.song, files.hugo)  # combine both datasets
all(file.exists(files.song))
tdata.song <- as.data.frame(files.song)
colnames(tdata.song) <- "files"
cell.lines.song <- c(rep("M395", 4), rep("M229", 4), rep("M238", 4), 
                     rep("SKMEL28", 2), rep("M229", 2), rep("M238", 2), 
                     rep("SKMEL28", 2))
tdata.song$cell <- factor(cell.lines.song)
treat.song <- relevel(factor(c(rep(c("P", "R"), 2), 
                               rep(c(rep("P", 2), rep("R", 2)), 2), "P", "R", 
                               rep(c("P", "R"), 3)), 
                             levels = c("P", "R", "DR")), ref = "P")
tdata.song$treatment <- treat.song
tdata.song$names <- paste("Sample", 1:20, sep = "")

# reference transcriptome checksum recognized by tximeta
se.song <- tximeta(tdata.song)  

# summarize quantifications to gene-level
gse.song <- summarizeToGene(se.song, countsFromAbundance = "lengthScaledTPM")  
y.song <- makeDGEList(gse.song)

# filter null and low counts
dim(y.song)
group.song <- paste(cell.lines.song, treat.song, sep = ".")
design.song <- model.matrix(~0 + group.song)
colnames(design.song) <- c("M229.P", "M229.R", "M238.P", "M238.R", "M395.P", 
    "M395.R", "SKMEL28.P", "SKMEL28.R")
keep.song <- filterByExpr(y.song, design.song)
dge.song <- y.song[keep.song, , keep.lib.sizes = FALSE]
dim(dge.song)[1]/dim(y.song)[1] * 100 

# plots
plot.filt(y.song, dge.song, c(0, 0.8))
plotMDS(y.song, col = as.numeric(as.factor(treat.song)), 
        labels = paste(dge.song$samples$cell, dge.song$samples$treatment, 
                                                                 sep = "."), 
        top = dim(dge.song$counts)[1])

# differential expression analysis
contrast.song <- makeContrasts(M229 = M229.R - M229.P, 
                               M238 = M238.R - M238.P, 
                               M395 = M395.R - M395.P, 
                               SKMEL = SKMEL28.R - SKMEL28.P, 
                               levels = design.song)

v.song <- voom(dge.song, design.song, plot = TRUE)

fit.song <- lmFit(v.song, design.song)
fit.song <- contrasts.fit(fit.song, contrast.song)
fit.song <- eBayes(fit.song, trend = FALSE)

# M229 genes
genes.song.m229 <- topTable(fit.song, coef = "M229", adjust.method = "BH",
                            sort.by = "p", number = dim(fit.song)[1])
sum(genes.song.m229$adj.P.Val < 0.05) 
sum(genes.song.m229$adj.P.Val < 0.001)  

# M238 genes
genes.song.m238 <- topTable(fit.song, coef = "M238", adjust.method = "BH", 
                            sort.by = "p", number = dim(fit.song)[1])
sum(genes.song.m238$adj.P.Val < 0.05)  
sum(genes.song.m238$adj.P.Val < 0.001)

# M395 genes
genes.song.m395 <- topTable(fit.song, coef = "M395", adjust.method = "BH", 
                            sort.by = "p", number = dim(fit.song)[1])
sum(genes.song.m395$adj.P.Val < 0.05)  
sum(genes.song.m395$adj.P.Val < 0.001)

# SKMEL28 genes
genes.song.skmel <- topTable(fit.song, coef = "SKMEL", adjust.method = "BH", 
                             sort.by = "p", number = dim(fit.song)[1])
sum(genes.song.skmel$adj.P.Val < 0.05)  
sum(genes.song.skmel$adj.P.Val < 0.001)  


# volcano plots
volcanoplot.colors.limma(genes = genes.song.m229, adjpval = 0.05, logfc = 1.5, 
                         title = "Song, M229R")
volcanoplot.colors.limma(genes = genes.song.m238, adjpval = 0.05, logfc = 1.5, 
                         title = "Song, M238R")
volcanoplot.colors.limma(genes = genes.song.m395, adjpval = 0.05, logfc = 1.5, 
                         title = "Song, M395R")
volcanoplot.colors.limma(genes = genes.song.skmel, adjpval = 0.05, logfc = 1.5, 
                         title = "Song, SKMEL28R")


# ----------------------------- logFC and SE -----------------------------------

# singleton A375
lfcSE.singleton.A375 <- (sqrt(fit.singleton$s2.post) * 
                           fit.singleton$stdev.unscaled)[, "treat.singletonR"]
genes.singleton.A375.se <- cbind(genes.singleton, 
                        lfcSE = lfcSE.singleton.A375[rownames(genes.singleton)])
genes.singleton.A375.lfc <- genes.singleton.A375.se[, c("logFC", "lfcSE")]

# ho 421Lu
lfcSE.ho.lu <- (sqrt(fit.ho$s2.post) * fit.ho$stdev.unscaled)[, "Lu"]
genes.ho.lu.se <- cbind(genes.ho.lu, lfcSE = lfcSE.ho.lu[rownames(genes.ho.lu)])
genes.ho.lu.lfc <- genes.ho.lu.se[, c("logFC", "lfcSE")]

# ho A375
lfcSE.ho.A375 <- (sqrt(fit.ho$s2.post) * fit.ho$stdev.unscaled)[, "A375"]
genes.ho.A375.se <- cbind(genes.ho.A375, 
                          lfcSE = lfcSE.ho.A375[rownames(genes.ho.A375)])
genes.ho.A375.lfc <- genes.ho.A375.se[, c("logFC", "lfcSE")]

# berico MM074
lfcSE.berico.mm <- (sqrt(fit.berico$s2.post) * 
                      fit.berico$stdev.unscaled)[, "treat.bericoR"]
genes.berico.mm.se <- cbind(genes.berico, 
                            lfcSE = lfcSE.berico.mm[rownames(genes.berico)])
genes.berico.mm.lfc <- genes.berico.mm.se[, c("logFC", "lfcSE")]

# song M229
lfcSE.song.m229 <- (sqrt(fit.song$s2.post) * fit.song$stdev.unscaled)[, "M229"]
genes.song.m229.se <- cbind(genes.song.m229, 
                            lfcSE = lfcSE.song.m229[rownames(genes.song.m229)])
genes.song.m229.lfc <- genes.song.m229.se[, c("logFC", "lfcSE")]

# song M238
lfcSE.song.m238 <- (sqrt(fit.song$s2.post) * fit.song$stdev.unscaled)[, "M238"]
genes.song.m238.se <- cbind(genes.song.m238, 
                            lfcSE = lfcSE.song.m238[rownames(genes.song.m238)])
genes.song.m238.lfc <- genes.song.m238.se[, c("logFC", "lfcSE")]

# song M395
lfcSE.song.m395 <- (sqrt(fit.song$s2.post) * fit.song$stdev.unscaled)[, "M395"]
genes.song.m395.se <- cbind(genes.song.m395, 
                            lfcSE = lfcSE.song.m395[rownames(genes.song.m395)])
genes.song.m395.lfc <- genes.song.m395.se[, c("logFC", "lfcSE")]

# song SKMEL28
lfcSE.song.skmel <- (sqrt(fit.song$s2.post) * fit.song$stdev.unscaled)[,"SKMEL"]
genes.song.skmel.se <- cbind(genes.song.skmel, 
                          lfcSE = lfcSE.song.skmel[rownames(genes.song.skmel)])
genes.song.skmel.lfc <- genes.song.skmel.se[, c("logFC", "lfcSE")]


################################################################################
#                      2. RESISTANCE METAANALYSES                              #
################################################################################

# save gene names
all.genes <- rownames(y.berico)

# function to format the results of the differential expression analysis  
# and prepare the data frame for the metaanalyses
search.gene <- function(gene) {
    mat <- matrix(ncol = 6, nrow = 0)
    df <- data.frame(mat)

    # singleton
    if (nrow(genes.singleton.A375.lfc[
                          rownames(genes.singleton.A375.lfc) == gene,]) == 1) {
        a <- genes.singleton.A375.lfc[
              rownames(genes.singleton.A375.lfc) == gene, ][1, 1]
        b <- genes.singleton.A375.lfc[
              rownames(genes.singleton.A375.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Singleton, A375R", "Singleton", "A375", "R"))
    } else {
        df <- rbind(df, c(NA, NA, "Singleton, A375R", "Singleton", "A375", "R"))
    }

    # ho 451Lu
    if (nrow(genes.ho.lu.lfc[rownames(genes.ho.lu.lfc) == gene, ]) == 1) {
        a <- genes.ho.lu.lfc[rownames(genes.ho.lu.lfc) == gene, ][1, 1]
        b <- genes.ho.lu.lfc[rownames(genes.ho.lu.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Ho, 451LuR", "Ho", "451Lu", "R"))
    } else {
        df <- rbind(df, c(NA, NA, "Ho, 451LuR", "Ho", "451Lu", "R"))
    }

    # ho A375
    if (nrow(genes.ho.A375.lfc[rownames(genes.ho.A375.lfc) == gene, ]) == 1) {
        a <- genes.ho.A375.lfc[rownames(genes.ho.A375.lfc) == gene, ][1, 1]
        b <- genes.ho.A375.lfc[rownames(genes.ho.A375.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Ho, A375R", "Ho", "A375", "R"))
    } else {
        df <- rbind(df, c(NA, NA, "Ho, A375R", "Ho", "A375", "R"))
    }

    # berico
    if (nrow(genes.berico.mm.lfc[rownames(genes.berico.mm.lfc) == gene,]) == 1){
        a <- genes.berico.mm.lfc[rownames(genes.berico.mm.lfc) == gene, ][1, 1]
        b <- genes.berico.mm.lfc[rownames(genes.berico.mm.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Berico, MM074R", "Berico", "MM074", "R"))
    } else {
        df <- rbind(df, c(NA, NA, "Berico, MM074R", "Berico", "MM074", "R"))
    }

    # song M229
    if (nrow(genes.song.m229.lfc[rownames(genes.song.m229.lfc) == gene,]) == 1){
        a <- genes.song.m229.lfc[rownames(genes.song.m229.lfc) == gene, ][1, 1]
        b <- genes.song.m229.lfc[rownames(genes.song.m229.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Song, M229R", "Song", "M229", "R"))
    } else {
        df <- rbind(df, c(NA, NA, "Song, M229R", "Song", "M229", "R"))
    }

    # song M238
    if (nrow(genes.song.m238.lfc[rownames(genes.song.m238.lfc) == gene,]) == 1){
        a <- genes.song.m238.lfc[rownames(genes.song.m238.lfc) == gene, ][1, 1]
        b <- genes.song.m238.lfc[rownames(genes.song.m238.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Song, M238R", "Song", "M238", "R"))
    } else {
        df <- rbind(df, c(NA, NA, "Song, M238R", "Song", "M238", "R"))
    }

    # song M395
    if (nrow(genes.song.m395.lfc[rownames(genes.song.m395.lfc) == gene,]) == 1){
        a <- genes.song.m395.lfc[rownames(genes.song.m395.lfc) == gene, ][1, 1]
        b <- genes.song.m395.lfc[rownames(genes.song.m395.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Song, M395R", "Song", "M395", "R"))
    } else {
        df <- rbind(df, c(NA, NA, "Song, M395R", "Song", "M395", "R"))
    }

    # song skmel
    if (nrow(genes.song.skmel.lfc[
                              rownames(genes.song.skmel.lfc) == gene, ]) == 1) {
        a <- genes.song.skmel.lfc[
              rownames(genes.song.skmel.lfc) == gene, ][1, 1]
        b <- genes.song.skmel.lfc[
              rownames(genes.song.skmel.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Song, SKMEL28R", "Song", "SKMEL28", "R"))
    } else {
        df <- rbind(df, c(NA, NA, "Song, SKMEL28R", "Song", "SKMEL28", "R"))
    }

    # end of search
    colnames(df) <- c("logFC", "lfcSE", "studylab", "author", "line", 
                      "phenotype")
    return(df)
}

# example of use
search.gene("ENSG00000000003.15")

# iterate over all the genes using the function
data.meta.resist <- lapply(all.genes, FUN = function(x) {
    search.gene(x)
})
names(data.meta.resist) <- all.genes

length(data.meta.resist)

# remove null results
for (g in all.genes) {
    if (all(is.na(data.meta.resist[[g]][, c(1, 2)]))) {
        data.meta.resist[[g]] <- NULL
    }
}
length(data.meta.resist)

# remove genes with info from only 1 study
for (g in names(data.meta.resist)) {
    if (sum(is.na(data.meta.resist[[g]])) == 14) {
        data.meta.resist[[g]] <- NULL
    }
}
length(data.meta.resist)

# remove digits after dots and save
names(data.meta.resist) <- as.vector(mapply(names(data.meta.resist), 
                                   FUN = function(x) {gsub("\\.\\d+", "", x)}))

# saveRDS(data.meta.resist, file = 'data-meta-resist.rds')

# metaanalyses function
metaanalyses <- function(data.metaanalyses, gene) {
    data.m <- data.metaanalyses[[gene]]
    data.m[, c("logFC", "lfcSE")] <- sapply(data.m[c("logFC", "lfcSE")], 
                                            as.numeric)
    m.gen <- metagen(TE = logFC, seTE = lfcSE, studlab = studylab, 
                     data = data.m, sm = "MD", fixed = FALSE, random = TRUE, 
                     method.tau = "REML", id = author, title = gene, 
                     prediction = TRUE, control = list(optimizer = "nmk"))
    summary(m.gen)
    return(m.gen)
}

# loop to iterate over all the genes using the function
counter <- 0
total <- length(data.meta.resist)
genes.meta.resist <- lapply(names(data.meta.resist), FUN = function(x) {
    counter <<- counter + 1
    print(paste0("iteration number ", counter, "; total ", total))
    metaanalyses(data.meta.resist, x)
})

names(genes.meta.resist) <- names(data.meta.resist)

# remove null results
length(genes.meta.resist)

for (g in names(data.meta.resist)) {
    if (is.null(genes.meta.resist[[g]])) {
        genes.meta.resist[[g]] <- NULL
    }
}
length(genes.meta.resist)

# saveRDS(genes.meta.resist, file='genes-meta-resist.rds') 

# summarize results in a data frame
results.meta <- function(meta, gene) {
    res <- data.frame(num.studies = meta[[gene]]$k.study, 
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
results.meta(genes.meta.resist, "ENSG00000000003")

# loop to iterate over all the  genes using the function
res.meta.resist <- lapply(names(genes.meta.resist), FUN = function(x) {
    results.meta(genes.meta.resist, x)
})
names(res.meta.resist) <- names(genes.meta.resist)
res.meta.resist.df <- do.call(rbind.data.frame, res.meta.resist)
dim(res.meta.resist.df)

# multiple testing correction: BH
res.meta.resist.df$adj.pval <- p.adjust(res.meta.resist.df$pval, method = "BH", 
                                        n = length(res.meta.resist.df$pval))

# multiple testing correction: Bonferroni
res.meta.resist.df$adj.pval.bonfer <- p.adjust(res.meta.resist.df$pval, 
                                               method = "bonferroni",
                                            n = length(res.meta.resist.df$pval))

# saveRDS(res.meta.resist.df, file='res-meta-resist-df.rds')


# ------------------------ metaanalyses results --------------------------------

# metaanalyses and experiments correlation 

## experiments

# save the LFC values of the experiments in a new variable
lfc.genes <- lapply(data.meta.resist, '[[', 'logFC')
lfc.genes <- do.call('rbind.data.frame', lfc.genes)
lfc.genes <- sapply(lfc.genes, as.numeric)
colnames(lfc.genes) <- data.meta.resist[[1]]$studylab
rownames(lfc.genes) <- names(data.meta.resist)
lfc.genes <- lfc.genes[names(genes.meta.resist), ]
dim(lfc.genes)

## metaanalyses

# save the LFC values of the metaanalyses in a new variable
lfc.meta <- lapply(genes.meta.resist, '[[', 'TE.random')
lfc.meta <- do.call('rbind.data.frame', lfc.meta)
lfc.meta <- sapply(lfc.meta, as.numeric)
colnames(lfc.meta) <- 'lfc.meta'
rownames(lfc.meta) <- names(genes.meta.resist)

# merge both
lfc.genes <- merge(lfc.genes, lfc.meta, by = 0)
rownames(lfc.genes) <- lfc.genes$Row.names
lfc.genes <- lfc.genes[ , -1]

# plot correlations using all the genes
lfc.genes.cor <- lfc.genes
lfc.genes.cor[lfc.genes[, "lfc.meta"] < 0, "lfc.meta"] <- -1 * 
                         lfc.genes.cor[lfc.genes[, "lfc.meta"] < 0, "lfc.meta"]

par(mfrow = c(2, 4))
for (i in 1:8) {
    lfc.genes.cor[lfc.genes[, "lfc.meta"] < 0, i] <- -1 * 
                                   lfc.genes.cor[lfc.genes[, "lfc.meta"] < 0, i]
    plot(lfc.genes.cor[, i], lfc.genes.cor[, "lfc.meta"], 
         main = colnames(lfc.genes.cor)[i],
         ylab = "meta lfc", xlab = "individual lfc", pch = 21, col = "blue4", 
         bg = "blue1", cex = 1, lwd = 1)
}
par(mfrow = c(1, 1))

# compute correlations using all the genes
first.cor <- cor.test(lfc.genes.cor[, 1], lfc.genes.cor[, "lfc.meta"], 
                      method = "pearson")
cor.res.all <- data.frame(cor = first.cor$estimate, pval = first.cor$p.value)
for (i in 2:8) {
    cor <- cor.test(lfc.genes.cor[, i], lfc.genes.cor[, "lfc.meta"], 
                    method = "pearson")
    cor.res.all <- rbind(cor.res.all, c(cor$estimate, cor$p.value))
}
cor.res.all

# now take only differentially expressed genes (DEGs) in order to reduce noise 
# in correlations and plots
dim(lfc.genes)
lfc.genes.05 <- lfc.genes[rownames(res.meta.resist.df[
                                      res.meta.resist.df$adj.pval <= 0.05, ]), ]
dim(lfc.genes.05)

# plot DEGs: down in blue and up in red
par(mfrow = c(4, 4))
for (i in 1:8) {
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

# correlations
lfc.genes.05.cor <- lfc.genes.05
for (i in 1:8) {
    lfc.genes.05.cor[lfc.genes.05.cor[, 9] < 0, i] <- -1 * 
                               (lfc.genes.05.cor[lfc.genes.05.cor[, 9] < 0, i])
}

# change downregulated DEGs sign first
lfc.genes.05.cor[lfc.genes.05.cor[, 9] < 0, 9] <- -1 * 
                               (lfc.genes.05.cor[lfc.genes.05.cor[, 9] < 0, 9])

# compute correlations
first.cor <- cor.test(lfc.genes.05.cor[, 1], lfc.genes.05.cor[, "lfc.meta"], 
                      method = "pearson")
cor.res <- data.frame(cor = first.cor$estimate, pval = first.cor$p.value)

for (i in 2:8) {
    cor <- cor.test(lfc.genes.05.cor[, i], lfc.genes.05.cor[, "lfc.meta"], 
                    method = "pearson")
    cor.res <- rbind(cor.res, c(cor$estimate, cor$p.value))
}

dim(cor.res)
cor.res

# plot DEGs
par(mfrow = c(2, 4))
for (i in 1:8) {
    plot(lfc.genes.05.cor[, i], lfc.genes.05.cor[, "lfc.meta"],
         main = colnames(lfc.genes.05.cor)[i], ylab = "meta lfc", 
         xlab = "individual lfc", pch = 21, col = "darkorchid4",
         bg = "darkorchid2", cex = 1, lwd = 1)
}
par(mfrow = c(1, 1))


# correlation between the data two-by-two 
# ! these correlations were computed taking into account all the genes, 
# not only DEGs
lfc.genes.pw.cor <- lfc.genes[,-9]
colnames(lfc.genes.pw.cor) <- c('Si.A375R', 'H.451LuR', 'H.A375R', 'B.MM074R', 
                              'So.M229R', 'So.M238R', 'S.M395R', 'S.SKMEL28R')

cor.matrix.lfc <- cor(lfc.genes.pw.cor, use = "pairwise.complete.obs")
col4 <- colorRampPalette(c("yellow2", "white", "purple4"))(50)
cor.pheat <- pheatmap(cor.matrix.lfc, col = col4, display_numbers = TRUE, 
                      breaks = seq(-1, 1, length = 50))
cor.pheat


# comparison with RankProd and RankSum

# load rankprod and ranksum results
rp.results.tables <- readRDS("data/rp-results-tables.rds")
rs.results.tables <- readRDS("data/rs-results-tables.rds")

RP.genes <- c(rownames(rp.results.tables$Table1), 
              rownames(rp.results.tables$Table2))
RS.genes <- c(rownames(rs.results.tables$Table1), 
              rownames(rs.results.tables$Table2))

# now take DEGs that were obtained with values of expression in all samples 
# in TE-MA
res.meta.resist.df.05.8stud <- res.meta.resist.df[
                         res.meta.resist.df$num.data == 8 &
                           res.meta.resist.df$adj.pval < 0.05, ]

# venn diagram
ggVennDiagram(list(ES_MA = rownames(res.meta.resist.df.05.8stud), 
                   RankProd_MA = RP.genes,
                   RankSum_MA = RS.genes), 
              color = "black", lwd = 0.8, lty = 1, label = "count") +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  scale_color_manual(values = rep("darkgray", 3)) + theme_void() + 
  theme(legend.position = "none")


# now compute correlations selecting only DEGs detected by the 3 MA
degs.8st <- intersect(rownames(res.meta.resist.df[
  res.meta.resist.df$adj.pval <= 0.05, ]), intersect(RP.genes, RS.genes))
lfc.genes.05.cor.8st <- lfc.genes.05.cor[degs.8st, ]
dim(lfc.genes.05.cor.8st)

first.cor <- cor.test(lfc.genes.05.cor.8st[, 1], 
                      lfc.genes.05.cor.8st[, "lfc.meta"], method = "pearson")
cor.res.degs.8st <- data.frame(cor = first.cor$estimate, 
                               pval = first.cor$p.value)

for (i in 2:8) {
    cor <- cor.test(lfc.genes.05.cor.8st[, i], 
                    lfc.genes.05.cor.8st[, "lfc.meta"], method = "pearson")
    cor.res.degs.8st <- rbind(cor.res.degs.8st, c(cor$estimate, cor$p.value))
}

dim(cor.res.degs.8st)
cor.res.degs.8st

# plot DEGs 
par(mfrow = c(2, 4))
for (i in 1:8) {
    plot(lfc.genes.05.cor.8st[, i], lfc.genes.05.cor.8st[, "lfc.meta"], 
         main = colnames(lfc.genes.05.cor.8st)[i], ylab = "meta lfc", 
         xlab = "individual lfc", pch = 21, col = "darkorchid4",
         bg = "darkorchid2", cex = 1, lwd = 1)
}
par(mfrow = c(1, 1))


# metaanalyses and experiments values comparison

## metaanalyses
res.meta.resist.df.05 <- res.meta.resist.df[degs.8st, ]
dim(res.meta.resist.df.05)[1]/dim(res.meta.resist.df)[1] * 100 
dim(res.meta.resist.df.05)

# upregulated genes
sum(res.meta.resist.df.05$lfc > 0)  # n upregulated genes
median(res.meta.resist.df.05[res.meta.resist.df.05$lfc > 0, ]$lfc)  

# downregulated genes
sum(res.meta.resist.df.05$lfc < 0)  # n downregulated genes
median(res.meta.resist.df.05[res.meta.resist.df.05$lfc < 0, ]$lfc)  


## experiments
degs.lfc.indiv <- merge(genes.singleton[
                       genes.singleton$adj.P.Val < 0.05, "logFC", drop = FALSE], 
              genes.ho.lu[genes.ho.lu$adj.P.Val < 0.05, "logFC", drop = FALSE],
              all = TRUE, by = 0) %>%
    merge(genes.ho.A375[genes.ho.A375$adj.P.Val < 0.05, "logFC", drop = FALSE], 
        all = TRUE, by.x = "Row.names", by.y = 0) %>%
    merge(genes.berico[genes.berico$adj.P.Val < 0.05, "logFC", drop = FALSE], 
        all = TRUE, by.x = "Row.names", by.y = 0) %>%
    merge(genes.song.m229[
                       genes.song.m229$adj.P.Val < 0.05, "logFC", drop = FALSE],
        all = TRUE, by.x = "Row.names", by.y = 0) %>%
    merge(genes.song.m238[
                       genes.song.m238$adj.P.Val < 0.05, "logFC", drop = FALSE],
        all = TRUE, by.x = "Row.names", by.y = 0) %>%
    merge(genes.song.m395[
                       genes.song.m395$adj.P.Val < 0.05, "logFC", drop = FALSE],
        all = TRUE, by.x = "Row.names", by.y = 0) %>%
    merge(genes.song.skmel[
                      genes.song.skmel$adj.P.Val < 0.05, "logFC", drop = FALSE],
        all = TRUE, by.x = "Row.names", by.y = 0)

rownames(degs.lfc.indiv) <- degs.lfc.indiv$Row.names
degs.lfc.indiv <- degs.lfc.indiv[, -1]
rownames(degs.lfc.indiv) <- as.vector(mapply(rownames(degs.lfc.indiv), 
                                   FUN = function(x) {gsub("\\.\\d+", "", x)}))
dim(degs.lfc.indiv)

# take the median from each experiment and then compute mean of medians
md.up <- c()
md.down <- c()
for (i in 1:8) {
    md.up <- c(md.up, median(degs.lfc.indiv[degs.lfc.indiv[, i] > 0, i], 
                             na.rm = TRUE))
    md.down <- c(md.down, median(degs.lfc.indiv[degs.lfc.indiv[, i] < 0, i], 
                                 na.rm = TRUE))
}
mean(md.up)
mean(md.down)

# plot distribution
p.degs.df <- rbind(data.frame(LFC = res.meta.resist.df.05$lfc, 
                              COND = rep('MA', 
                                         length(res.meta.resist.df.05$lfc))),
                   data.frame(LFC = degs.lfc.indiv[ , 1], 
                              COND = rep('EXP1', length(degs.lfc.indiv[ , 1]))),
                   data.frame(LFC = degs.lfc.indiv[ , 2], 
                              COND = rep('EXP2', length(degs.lfc.indiv[ , 2]))),
                   data.frame(LFC = degs.lfc.indiv[ , 3], 
                              COND = rep('EXP3', length(degs.lfc.indiv[ , 3]))),
                   data.frame(LFC = degs.lfc.indiv[ , 4], 
                              COND = rep('EXP4', length(degs.lfc.indiv[ , 4]))),
                   data.frame(LFC = degs.lfc.indiv[ , 5], 
                              COND = rep('EXP5', length(degs.lfc.indiv[ , 5]))),
                   data.frame(LFC = degs.lfc.indiv[ , 6], 
                              COND = rep('EXP6', length(degs.lfc.indiv[ , 6]))),
                   data.frame(LFC = degs.lfc.indiv[ , 7], 
                              COND = rep('EXP7', length(degs.lfc.indiv[ , 7]))),
                   data.frame(LFC = degs.lfc.indiv[ , 8], 
                              COND = rep('EXP8', length(degs.lfc.indiv[ , 8]))))

p.degs.df$COND <- factor(p.degs.df$COND, levels = c("EXP1", "EXP2", "EXP3", 
    "EXP4", "EXP5", "EXP6", "EXP7", "EXP8", "MA"), ordered = TRUE)
p.degs.df$EXPRESSION <- NA
p.degs.df$EXPRESSION[p.degs.df$LFC > 0] <- "UP"
p.degs.df$EXPRESSION[p.degs.df$LFC < 0] <- "DOWN"

p.degs <- ggplot(subset(p.degs.df, !is.na(LFC)), 
                 aes(EXPRESSION, LFC, fill = COND)) +
    geom_boxplot(alpha = 1, outlier.shape = NA) + theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, size = 20), 
        axis.title = element_text(size = 8), axis.text = element_text(size = 8),
    legend.title = element_blank()) + labs(y = "LFC", x = "") + ylim(c(-4, 4)) +
    scale_x_discrete(labels = c("DOWN", "UP")) + 
  scale_fill_manual(values = c(rep("gray68", 8), "limegreen")) + 
  geom_hline(yintercept = c(mean(md.up), mean(md.down)), col = "red",
    lty = 5)

p.degs

## heatmap of DEGs
col.hm <- colorRampPalette(c("green", "black", "red"))(50)

p.heatmap.degs <- pheatmap(lfc.genes[degs.8st, -9], col = col.hm, 
                           cluster_rows = TRUE, cluster_cols = TRUE, 
                           breaks = seq(-8, 8, length = 50), 
                           show_rownames = FALSE, fontsize_col = 10)

p.heatmap.degs



# -------------------------- subgroup analysis ---------------------------------

# select genes expressed in all the samples
genes.meta.resist.8st <- genes.meta.resist[res.meta.resist.df$num.data == 8]

# define subgroups
res.mapk <- c("red", "act", "red", "act", "red", "red", "act", "red")

# subgroup analysis function
i <- 0
error.subg.genes <- c()
sub.resistance.fx <- function(x, meta.list) {
    print(paste0("iteration number ", i))
    i <<- i + 1
    res <- try(update.meta(meta.list[[x]], subgroup = res.mapk, 
                           tau.common = TRUE))
    if (class(res) == "try-error") {
        res <- try(update.meta(meta.list[[x]], subgroup = as.factor(res.mapk), 
                               tau.common = TRUE))
        if (class(res) == "try-error") {
            res <- NULL
            error.genes <<- c(error.subg.genes, x)
        }
    }
    return(res)
}

# example of use
sub.resistance.fx("ENSG00000179476", genes.meta.resist.8st)

# loop to iterate over all the genes using the function
subg.meta.res <- lapply(names(genes.meta.resist.8st), FUN = function(x) {
    sub.resistance.fx(x, genes.meta.resist.8st)
})
names(subg.meta.res) <- names(genes.meta.resist.8st)

# remove null results
null.subg <- lapply(names(subg.meta.res), function(x) {
    is.null(subg.meta.res[[x]])
})

null.subg.names <- names(subg.meta.res)[
                           which(do.call(rbind.data.frame, null.subg) == TRUE)]
length(subg.meta.res)
subg.meta.res[null.subg.names] <- NULL
length(subg.meta.res)

# function to retrieve subgroup analysis significances
diff.subg.meta <- function(meta, gene) {
    res <- data.frame(pval.btw = meta[[gene]]$pval.Q.b.random, 
                      pval.wth = meta[[gene]]$pval.Q.w.random)
    return(res)
}

# example of use
diff.subg.meta(subg.meta.res, "ENSG00000000003")

# loop to iterate over all the genes using the function
diff.subg.data <- lapply(names(subg.meta.res), function(x) {
    diff.subg.meta(subg.meta.res, x)
})
names(diff.subg.data) <- names(subg.meta.res)
diff.subg.data <- do.call(rbind.data.frame, diff.subg.data)

# select genes with the same expression within groups but not between them
diff.genes.btw <- diff.subg.data[diff.subg.data$pval.btw <= 0.05 &
                                   diff.subg.data$pval.wth > 0.05, ]

# add symbol column
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
btw.sy <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), 
                filters = "ensembl_gene_id",
                values = rownames(diff.genes.btw), mart = mart)

diff.genes.btw <- merge(btw.sy, diff.genes.btw, by.x = "ensembl_gene_id", 
                        by.y = 0, all = TRUE)
rownames(diff.genes.btw) <- diff.genes.btw$ensembl_gene_id
diff.genes.btw <- diff.genes.btw[, -1]
dim(diff.genes.btw)


# heatmap
lfc.genes.subgdiff <- lfc.genes[rownames(diff.genes.btw), -9]
p.subgdiff <- pheatmap(lfc.genes.subgdiff, col = col.hm, cluster_rows = TRUE, 
                       cluster_cols = TRUE,
    breaks = seq(-2, 2, length = 50), show_rownames = FALSE, fontsize_row = 7)

p.subgdiff


# -------------------------------- TFEA-ChIP -----------------------------------

# retrieve data
tfea.data <- res.meta.resist.df[, c("lfc", "pval", "adj.pval")]
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
CM_list_UP <- contingency_matrix(genes.upreg, genes.control)

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

head(TF_ranking[["TF_ranking"]])

plot_CM(pval_mat_DOWN)  # plot p-values against ORs

tfs.down <- rownames(tf.ranking.down)
names(tfs.down) <- rownames(tf.ranking.down)
col.down <- sample(viridis::turbo(length(tfs.down)), length(tfs.down))

# plot p-values against ORs highlighting indicated TFs
plot_CM(pval_mat_DOWN, specialTF = tfs.down, TF_colors = col.down)  


# --------------------- GSVA - resistance signature ----------------------------

# retrieve counts
counts.df <- merge(y.singleton$counts, y.ho$counts, by = 0, all = TRUE) %>%
    merge(y.berico$counts, by.x = "Row.names", by.y = 0) %>%
    merge(y.song$counts, by.x = "Row.names", by.y = 0)
rownames(counts.df) <- as.vector(mapply(counts.df$Row.names, FUN = function(x) {
    gsub("\\.\\d+", "", x)
}))
counts.df <- counts.df[, -1]

# add colnames
colnames(counts.df) <- c(paste("singleton", y.singleton$samples$cell, 
                               y.singleton$samples$treatment, sep = "."), 
                         paste("ho", y.ho$samples$cell, 
                               y.ho$samples$treatment, sep = "."),
                         paste("berico", y.berico$samples$cell, 
                               y.berico$samples$treatment, sep = "."),
                         paste("song", y.song$samples$cell, 
                               y.song$samples$treatment, sep = "."))

# read resistance signature
res.sig <- read.table("data/resistance-sig.grp")

# get symbols
list.genes <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), 
                    values = rownames(res.meta.resist.df),
                    mart = mart)
res.sig.ensemblid <- list.genes[which(list.genes$hgnc_symbol %in% res.sig$V1), ]
sign <- list(c(res.sig.ensemblid$ensembl_gene_id))

# gsva
gsva.resistance <- gsva(as.matrix(counts.df), sign, verbose = FALSE, 
                        method = "ssgsea")
colfunc <- colorRampPalette(c("yellow", "goldenrod1", "orangered"))
pheatmap(gsva.resistance, col = colfunc(100), display_numbers = TRUE, 
         cluster_rows = FALSE, cluster_cols = TRUE, main = "Resistance")

# repeat the same but computing average between 'technical' replicates first
counts.df.avg <- data.frame(row.names = rownames(counts.df))
counts.df.avg$singleton.A395 <- rowMeans2(
  as.matrix(counts.df[, grepl("singleton.A375.P", colnames(counts.df))]))
counts.df.avg$singleton.A395R <- rowMeans2(
  as.matrix(counts.df[, grepl("singleton.A375.R", colnames(counts.df))]))
counts.df.avg$ho.451Lu <- rowMeans2(
             as.matrix(counts.df[, grepl("ho.451Lu.P", colnames(counts.df))]))
counts.df.avg$ho.451LuR <- rowMeans2(
             as.matrix(counts.df[, grepl("ho.451Lu.R", colnames(counts.df))]))
counts.df.avg$ho.A375 <- rowMeans2(
             as.matrix(counts.df[, grepl("ho.A375.P", colnames(counts.df))]))
counts.df.avg$ho.A375R <- rowMeans2(
             as.matrix(counts.df[, grepl("ho.A375.R", colnames(counts.df))]))
counts.df.avg$berico.MM074 <- rowMeans2(
          as.matrix(counts.df[, grepl("berico.MM074.P", colnames(counts.df))]))
counts.df.avg$berico.MM074R <- rowMeans2(
          as.matrix(counts.df[, grepl("berico.MM074.R", colnames(counts.df))]))
counts.df.avg$song.M395 <- rowMeans2(
          as.matrix(counts.df[, grepl("song.M395.P", colnames(counts.df))]))
counts.df.avg$song.M395R <- rowMeans2(
  as.matrix(counts.df[, grepl("song.M395.R", colnames(counts.df))]))
counts.df.avg$song.M229 <- rowMeans2(
  as.matrix(counts.df[, grepl("song.M229.P", colnames(counts.df))]))
counts.df.avg$song.M229R <- rowMeans2(
  as.matrix(counts.df[, grepl("song.M229.R", colnames(counts.df))]))
counts.df.avg$song.M238 <- rowMeans2(
  as.matrix(counts.df[, grepl("song.M238.P", colnames(counts.df))]))
counts.df.avg$song.M238R <- rowMeans2(
  as.matrix(counts.df[, grepl("song.M238.R", colnames(counts.df))]))
counts.df.avg$song.SKMEL28 <- rowMeans2(
  as.matrix(counts.df[, grepl("song.SKMEL28.P", colnames(counts.df))]))
counts.df.avg$song.SKMEL28R <- rowMeans2(
  as.matrix(counts.df[, grepl("song.SKMEL28.R", colnames(counts.df))]))

# gsva
gsva.resistance.avg <- gsva(as.matrix(counts.df.avg), sign, verbose = FALSE, 
                            method = "ssgsea")
# heatmap gsva enrichment
pheat.gsva.res <- pheatmap(gsva.resistance.avg[, order(gsva.resistance.avg),
                                               drop = FALSE], col = colfunc(50), 
                           display_numbers = TRUE, cluster_rows = FALSE, 
                           cluster_cols = FALSE)
pheat.gsva.res

## experiments LFC heatmap (selecting this signature)

# add a symbol column in the experiments LFC data
lfc.genes.sy <- merge(list.genes, lfc.genes, by.x = 'ensembl_gene_id', by = 0,
                      all.y = TRUE)

# select signature genes
lfc.genes.pheat.resist <- lfc.genes.sy[lfc.genes.sy$ensembl_gene_id %in% 
                                         res.sig.ensemblid$ensembl_gene_id, ]
rownames(lfc.genes.pheat.resist) <- lfc.genes.pheat.resist$ensembl_gene_id
lfc.genes.pheat.resist <- lfc.genes.pheat.resist[, -c(1, 2)]
pheat.resist <- pheatmap(lfc.genes.pheat.resist, col = col.hm, 
                         breaks = seq(-2, 2, length = 50), cluster_rows = TRUE, 
                         cluster_cols = TRUE, fontsize_col = 13, 
                         na_col = "black", display_numbers = FALSE)
pheat.resist

# now select the second cluster of the previous heatmap 
hcluster <- pheat.resist$tree_row
lbl <- cutree(hcluster, 2)  # 2 is the number of gene-groups 
which(lbl == 1)  # find genes corresponding to he 1st group
resist.genes.clean <- names(lbl)[lbl == 1]

# plot
pheat.resist.clean <- pheatmap(lfc.genes[resist.genes.clean, -9], col = col.hm, 
                               show_rownames = FALSE, 
                               breaks = seq(-2, 2, length = 50), 
                               fontsize_row = 7, fontsize = 15)
pheat.resist.clean

# repeat gsva using clean gene set
sign2 <- list(resist.genes.clean)
gsva.resistance.avg.clean <- gsva(as.matrix(counts.df.avg), sign2, 
                                            verbose = FALSE, method = "ssgsea")

# gsva heatmap
pheat.gsva.res.2 <- pheatmap(
                  gsva.resistance.avg.clean[, order(gsva.resistance.avg.clean),
                                                                  drop = FALSE], 
                  col = colfunc(50), display_numbers = TRUE, 
                  cluster_rows = FALSE,  cluster_cols = FALSE, 
                  main = "Resistance")
pheat.gsva.res.2

# paired data plot
paired.plot.df <- data.frame(ES = c(gsva.resistance.avg.clean), 
                             lines = c(rep("A375.singleton", 2), 
                                       rep("451Lu", 2), rep("A375.ho", 2), 
                                       rep("MM074", 2), rep("M395", 2), 
                                       rep("M229", 2), rep("M238", 2), 
                                       rep("SKMEL28", 2)), 
                             group = c(rep(c("P", "R"), 8)))

paired.plot <- ggplot(paired.plot.df, aes(x = group, y = ES, label = lines)) + 
  geom_line(aes(group = lines), color = "darkgrey") + 
  geom_point(size = 2, aes(color = group)) + scale_x_discrete("") + 
  scale_color_manual(values = c("darkblue", "darkblue")) + theme_bw() + 
  theme(legend.position = "none", axis.text = element_text(size = 14), 
        axis.title = element_text(size = 14)) +
  ylab("NES") + ylim(c(1, 2.25)) + 
  geom_text_repel(position = position_nudge_repel())

paired.plot

# compare A375 singleton vs. ho
pheatmap(na.omit(lfc.genes[, grepl("A375", colnames(lfc.genes))]), col = col.hm,
    show_rownames = FALSE, breaks = seq(-2, 2, length = 50), 
    fontsize_row = 7, fontsize = 15)

# PCA using all the data
pca.all <- prcomp(t(na.omit(lfc.genes[, -9])))
ggplotly(autoplot(pca.all, label = FALSE, 
           col = c("coral2", "seagreen4")[factor(c(1, 2, 1, 2, 1, 1, 2, 1))]) + 
           aes(text = colnames(lfc.genes[, -9])), tooltip = "text")

autoplot(pca.all, x = 1, y = 3, label = FALSE, size = 2,
         col = c("red", "darkblue")[factor(c(1, 2, 1, 2, 1, 1, 2, 1))]) + 
  theme_bw() + geom_text_repel(mapping = aes(label = colnames(lfc.genes[, -9])))
  

# ------------------------- cell cycle analysis --------------------------------

kegg.cell.cycle.ensemblid.clean <- read.table('data/cell-cycle-clean.txt')[,1]
pheat.cell.cycle.clean <- pheatmap(na.omit(lfc.genes[, -9][
                                            kegg.cell.cycle.ensemblid.clean, ]), 
                                   col = col.hm, show_rownames = FALSE, 
                                   breaks = seq(-5, 5, length=50), 
                                   fontsize_row = 7, border_color = NA,
                                   fontsize = 15)
pheat.cell.cycle.clean


################################################################################
#         3. DOUBLE-DRUG RESISTANT LINES ANALYSIS AND COMPARISON               #
################################################################################

# load data
files.ddr <- read.table("../output/files/song-part2.txt")[c(10, 11, 13:16), ]
files.ddr2 <- read.table("../output/files/song.txt")[c(20, 21, 40), ]
files.hugo.ddr <- read.table("../output/files/hugo.txt")[c(62, 64, 65, 67, 69),]
files.ddr <- c(file.path(unlist(c(files.ddr, files.ddr2, files.hugo.ddr))))
all(file.exists(files.ddr))
tdata.ddr <- as.data.frame(files.ddr)
colnames(tdata.ddr) <- "files"
cell.lines.ddr <- c("SKMEL28", rep("M249", 3), "SKMEL28", rep("M229", 3), 
    "M249", rep("M229", 2), "M238", rep("SKMEL28", 2))
tdata.ddr$cell <- factor(cell.lines.ddr)
treat.ddr <- relevel(factor(c(rep("P", 2), rep("DR", 4), rep("P", 4), "DR", 
    rep("P", 2), "DR"), levels = c("P", "DR")), ref = "P")
tdata.ddr$treatment <- treat.ddr
tdata.ddr$names <- paste("Sample", 1:14, sep = "")

# reference transcriptome checksum recognized by tximeta
se.ddr <- tximeta(tdata.ddr) 

# summarize quantifications to gene-level
gse.ddr <- summarizeToGene(se.ddr, countsFromAbundance = "lengthScaledTPM") 
y.ddr <- makeDGEList(gse.ddr)
y.ddr$samples

# filter null and low counts
dim(y.ddr)
design.ddr <- model.matrix(~cell.lines.ddr + treat.ddr)
keep.ddr <- filterByExpr(y.ddr, design.ddr)
dge.ddr <- y.ddr[keep.ddr, , keep.lib.sizes = FALSE]
dim(dge.ddr)[1]/dim(y.ddr)[1] * 100  

# plots
plot.filt(y.ddr, dge.ddr, c(0, 1))
plotMDS(dge.ddr, col = as.numeric(as.factor(treat.ddr)))

# differential expression analysis
v.ddr <- voom(dge.ddr, design.ddr, plot = TRUE)

fit.ddr <- lmFit(v.ddr, design.ddr)
fit.ddr <- eBayes(fit.ddr, trend = FALSE)

# DDR genes
genes.ddr <- topTable(fit.ddr, coef = "treat.ddrDR", adjust.method = "BH", 
                      sort.by = "p", number = dim(fit.ddr)[1])
sum(genes.ddr$adj.P.Val < 0.05)
sum(genes.ddr$adj.P.Val < 0.001)

# volcano plot
volcanoplot.colors.limma(genes = genes.ddr, adjpval = 0.05, logfc = 1.5, 
                         title = "DDR")

## enrichment analysis 
# GO BP
# no results
ans.go.ddr <- enrichGO(gene = rownames(genes.ddr[genes.ddr$adj.P.Val < 0.05, ]),
                       keyType = "ENSEMBL", 
                       ont = "BP", 
                       OrgDb = "org.Hs.eg.db", 
                       universe = rownames(genes.ddr),
                       readable = TRUE, 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05)
tab.go.ddr <- as.data.frame(ans.go.ddr)
head(tab.go.ddr)


# GO CC
# no results
ans.go.ddr.cc <- enrichGO(gene = rownames(genes.ddr[
                                                 genes.ddr$adj.P.Val < 0.05, ]),
                          keyType = "ENSEMBL",
                          ont = "CC", 
                          OrgDb = "org.Hs.eg.db", 
                          universe = rownames(genes.ddr),
                          readable = TRUE, 
                          pAdjustMethod = "BH", 
                          pvalueCutoff = 0.05)
tab.go.ddr.cc <- as.data.frame(ans.go.ddr.cc)
head(tab.go.ddr.cc)


# intersection btw treatment and DDR DEGs
rownames(genes.ddr) <- as.vector(mapply(rownames(genes.ddr), 
                                   FUN = function(x) {gsub("\\.\\d+", "", x)}))
# venn diagram
ggVennDiagram(list(SDR = degs.8st, 
                   DDR = rownames(genes.ddr)[genes.ddr$adj.P.Val < 0.05]), 
              label_alpha = 0, label = "count", label_size = 7) + 
  scale_fill_gradient(low = alpha("forestgreen", 0.3), 
                      high = alpha("forestgreen", 0.5)) + 
  scale_color_manual(values = c(DDR_genes = "gray",
                                MA_SDR_genes = "gray")) + 
  theme_void() + theme(legend.position = "none")

# intersection
inter.ddr <- intersect(rownames(genes.ddr)[genes.ddr$adj.P.Val < 0.05], 
                       degs.8st)

inter.ddr.lfc <- merge(res.meta.resist.df[inter.ddr, "lfc", drop = FALSE], 
                       genes.ddr[inter.ddr, "logFC", drop = FALSE], by = 0)
rownames(inter.ddr.lfc) <- paste(inter.ddr.lfc$Row.names, 
                                 list.genes[list.genes$ensembl_gene_id %in%
                                   inter.ddr.lfc$Row.names, "hgnc_symbol"], 
                                 sep = " - ")
inter.ddr.lfc <- inter.ddr.lfc[, -1]
colnames(inter.ddr.lfc) <- c("SDR", "DDR")

# select genes with a LFC > 1.5 or < -1.5 in any of the conditions
inter.ddr.lfc.sel <- inter.ddr.lfc[abs(inter.ddr.lfc$SDR) > 1.5 | 
                                     abs(inter.ddr.lfc$DDR) > 1.5, ]
inter.ddr.pheat <- pheatmap(inter.ddr.lfc.sel, cluster_cols = FALSE, 
    col = col.hm, show_rownames = TRUE, border_color = NA, fontsize_col = 10, 
    breaks = seq(-6, 6, length = 50))

inter.ddr.pheat

# enrichment analysis (intersection)
# GO CC
ans.go.ddr.inter.cc <- enrichGO(gene = inter.ddr,
                                keyType = 'ENSEMBL',
                                ont = "CC",
                                OrgDb ="org.Hs.eg.db",
                                universe = rownames(genes.ddr),
                                readable = TRUE,
                                pAdjustMethod = 'BH',
                                pvalueCutoff = 0.05)
tab.go.ddr.inter.cc <- as.data.frame(ans.go.ddr.inter.cc)
head(tab.go.ddr.inter.cc)

# plot
p.ddr.inter.cc <- barplot(ans.go.ddr.inter.cc, showCategory = 25)
p.ddr.inter.cc
p.ddr.inter.2 <- dotplot(ans.go.ddr.inter.cc, showCategory = 25)
p.ddr.inter.2


################################################################################
#              4. TREATMENT vs. RESISTANCE METAANALYSIS                        #
################################################################################

# venn diagram to compare treatment and resistance DEGs
inter.genes.resist <- degs.8st
inter.genes.long <- read.table("data/inter-genes-treatment.txt")$V1
inter.genes.treat.resist <- intersect(inter.genes.long, inter.genes.resist)

# venn diagram
ggVennDiagram(list(treatment = inter.genes.long, 
                   resistance = inter.genes.resist),
    label_alpha = 0, label = "count", set_size = 6, label_size = 7) + 
  scale_fill_gradient(low = alpha("steelblue3", 0.7), high = "steelblue3") + 
  scale_color_manual(values = c(TreatmentMA = "gray", ResistanceMA = "gray")) + 
  theme_void() + theme(legend.position = "none")

## heatmap

# load "long-term" (24-48 h) treatment data first
res.meta.long.df <- readRDS("data/res-meta-treatment-df.rds")

# merge data to be compared
res.meta.treat.resist <- merge(res.meta.long.df[inter.genes.treat.resist, "lfc",
    drop = FALSE], res.meta.resist.df[inter.genes.treat.resist, "lfc", 
    drop = FALSE], by = 0)
res.meta.treat.resist <- merge(list.genes, res.meta.treat.resist, 
                     by.x = "ensembl_gene_id", by.y = "Row.names", all.y = TRUE)

rownames(res.meta.treat.resist) <- paste(res.meta.treat.resist$ensembl_gene_id, 
                                 res.meta.treat.resist$hgnc_symbol, sep = " - ")
colnames(res.meta.treat.resist)[3:4] <- c("treatment", "resistance")

# select only genes with abs(LFC) > 1.5 in "treatment" or "resistance" and plot
# (heatmap)
pheatmap(res.meta.treat.resist[, 3:4][abs(res.meta.treat.resist[, 3]) > 1.5 | 
                                      abs(res.meta.treat.resist[, 4]) > 1.5, ], 
         color = col.hm, border_color = NA, breaks = seq(-4, 4, length = 50),
         cluster_cols = FALSE, cluster_rows = TRUE, fontsize_col = 15)


## association test: treatment and resistance DEGs

# contingency table
#                ..............................
#                | DEGs_res_yes | DEGs_res_no |
# .............................................
# DEGs_treat_yes |    A = 88    |  C = 3396   |
# .............................................
# DEGs_treat_no  |    B = 171   |  D = 23114  |
# .............................................

# numbers
A <- length(intersect(inter.genes.long, inter.genes.resist))
B <- length(intersect(inter.genes.resist, setdiff(rownames(res.meta.long.df), 
                                                  inter.genes.long)))
C <- length(intersect(inter.genes.long, setdiff(rownames(res.meta.resist.df), 
                                                inter.genes.resist)))
D <- length(intersect(setdiff(rownames(res.meta.long.df), inter.genes.long), 
                    setdiff(rownames(res.meta.resist.df), inter.genes.resist)))

# data frame
treat.res.assoc <- data.frame(DEGs_res_yes = c(A, B), DEGs_res_no = c(C, D), 
                              row.names = c("DEGs_treat_yes", "DEGs_treat_no"), 
                              stringsAsFactors = FALSE)

# fisher test
f.test <- fisher.test(treat.res.assoc)
f.test


## GO terms enrichment (intersection treatment & resistance)
# no results (BP/CC)
ans.go.treat.resist <- enrichGO(gene = inter.genes.treat.resist, 
                                keyType = "ENSEMBL", 
                                ont = "CC", 
                                OrgDb = "org.Hs.eg.db", 
                                universe = rownames(res.meta.resist.df),
                                readable = TRUE, 
                                pAdjustMethod = "BH", 
                                pvalueCutoff = 0.05)
tab.go.treat.resist <- as.data.frame(ans.go.treat.resist)
head(tab.go.treat.resist)


## intersection vs. the cytoskeleton list

# read list containing actomyosin genes
actomyosin <- read.table("data/actomyosin.txt")
actomyosin <- actomyosin$V1

# intersection
intersect(res.meta.treat.resist$hgnc_symbol, actomyosin)  # only CIT


# ## comparison treatment vs. Subgroup A metaanalyses
# 
# inter.genes.treat.subA <- intersect(inter.genes.long, degs.5st)
# 
# ggVennDiagram(list(Tratamiento = inter.genes.long, Resistencia = degs.5st), 
#               label_alpha = 0, label = "count", set_size = 7, label_size = 7) + 
#   scale_fill_gradient(low = alpha("forestgreen", 0.7), 
#                       high = alpha("coral2", 0.7)) + 
#   scale_color_manual(values = c(TreatmentMA = "gray", esistanceMA = "gray")) + 
#   theme_void() + theme(legend.position = "none")
# 
# ## heatmap
# res.meta.treat.subA <- merge(res.meta.long.df[inter.genes.treat.subA, 
#                                                            "lfc", drop = FALSE],
#                               res.meta.resist.subA.df[inter.genes.treat.subA,
#                                                    "lfc", drop = FALSE], by = 0)
# res.meta.treat.subA <- merge(list.genes, res.meta.treat.subA, 
#                               by.x = "ensembl_gene_id", by.y = "Row.names", 
#                               all.y = TRUE)
# rownames(res.meta.treat.subA) <- paste(res.meta.treat.subA$ensembl_gene_id, 
#                                         res.meta.treat.subA$hgnc_symbol, 
#                                         sep = " - ")
# colnames(res.meta.treat.subA)[3:4] <- c("Treatment", "Resistance")
# 
# pheatmap(res.meta.treat.subA[, 3:4][abs(res.meta.treat.subA[, 3]) > 1.5 | 
#     abs(res.meta.treat.subA[, 4]) > 1.5, ], color = col.hm, border_color = NA, 
#     breaks = seq(-4, 4, length = 50), cluster_cols = FALSE, cluster_rows = TRUE)


################################################################################
#                    5. METAANALYSES SEPARATING SUBGROUPS                      #
################################################################################

# ----------- Subgroup A metaanalyses: Song SKMEL28R, Song M229R, --------------
# ----------------- Song M238R, Singleton A375R and Ho A375R -------------------

# prepare the data to metaanalyze
search.gene.subA <- function(gene) {
    mat <- matrix(ncol = 6, nrow = 0)
    df <- data.frame(mat)

    # singleton
    if (nrow(genes.singleton.A375.lfc[
                           rownames(genes.singleton.A375.lfc) == gene,]) == 1) {
        a <- genes.singleton.A375.lfc[
                             rownames(genes.singleton.A375.lfc) == gene, ][1, 1]
        b <- genes.singleton.A375.lfc[
                             rownames(genes.singleton.A375.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Singleton, A375R", "Singleton", "A375", "R"))
    } else {
        df <- rbind(df, c(NA, NA, "Singleton, A375R", "Singleton", "A375", "R"))
    }

    # ho A375
    if (nrow(genes.ho.A375.lfc[rownames(genes.ho.A375.lfc) == gene, ]) == 1) {
        a <- genes.ho.A375.lfc[rownames(genes.ho.A375.lfc) == gene, ][1, 1]
        b <- genes.ho.A375.lfc[rownames(genes.ho.A375.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Ho, A375R", "Ho", "A375", "R"))
    } else {
        df <- rbind(df, c(NA, NA, "Ho, A375R", "Ho", "A375", "R"))
    }

    # song M229
    if (nrow(genes.song.m229.lfc[rownames(genes.song.m229.lfc) == gene,]) == 1){
        a <- genes.song.m229.lfc[rownames(genes.song.m229.lfc) == gene, ][1, 1]
        b <- genes.song.m229.lfc[rownames(genes.song.m229.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Song, M229R", "Song", "M229", "R"))
    } else {
        df <- rbind(df, c(NA, NA, "Song, M229R", "Song", "M229", "R"))
    }

    # song M238
    if (nrow(genes.song.m238.lfc[rownames(genes.song.m238.lfc) == gene,]) == 1){
        a <- genes.song.m238.lfc[rownames(genes.song.m238.lfc) == gene, ][1, 1]
        b <- genes.song.m238.lfc[rownames(genes.song.m238.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Song, M238R", "Song", "M238", "R"))
    } else {
        df <- rbind(df, c(NA, NA, "Song, M238R", "Song", "M238", "R"))
    }

    # song skmel
    if (nrow(genes.song.skmel.lfc[
                              rownames(genes.song.skmel.lfc) == gene, ]) == 1) {
        a <- genes.song.skmel.lfc[rownames(genes.song.skmel.lfc) == gene,][1, 1]
        b <- genes.song.skmel.lfc[rownames(genes.song.skmel.lfc) == gene,][1, 2]
        df <- rbind(df, c(a, b, "Song, SKMEL28R", "Song", "SKMEL28", "R"))
    } else {
        df <- rbind(df, c(NA, NA, "Song, SKMEL28R", "Song", "SKMEL28", "R"))
    }

    # end of search
    colnames(df) <- c("logFC", "lfcSE", "studylab", "author", "line", 
                      "phenotype")
    return(df)
}

# example of use
search.gene.subA("ENSG00000000003.15")

# loop to iterate over all the genes uusing the function
data.meta.resist.subA <- lapply(all.genes, FUN = search.gene.subA)
names(data.meta.resist.subA) <- all.genes

# remove null results
length(data.meta.resist.subA)
for (g in all.genes) {
    if (all(is.na(data.meta.resist.subA[[g]][, c(1, 2)]))) {
        data.meta.resist.subA[[g]] <- NULL
    }
}
length(data.meta.resist.subA)

# remove genes with info from only 1 study
for (g in names(data.meta.resist.subA)) {
    if (sum(is.na(data.meta.resist.subA[[g]])) == 8) {
        data.meta.resist.subA[[g]] <- NULL
    }
}
length(data.meta.resist.subA)

# remove digits after dots and save
names(data.meta.resist.subA) <- as.vector(mapply(names(data.meta.resist.subA),
    FUN = function(x) {
        gsub("\\.\\d+", "", x)
    }))

# saveRDS(data.meta.resist.subA, file = "data-meta-resist-subA.rds")

## metaanalyses
counter <- 0
total <- length(data.meta.resist.subA)
genes.meta.resist.subA <- lapply(names(data.meta.resist.subA), 
                                                          FUN = function(x) {
    counter <<- counter + 1
    print(paste0("iteration number ", counter, "; total ", total))
    metaanalyses(data.meta.resist.subA, x)
})

names(genes.meta.resist.subA) <- names(data.meta.resist.subA)

# remove null results and save
length(genes.meta.resist.subA)

for (g in names(data.meta.resist.subA)) {
    if (is.null(genes.meta.resist.subA[[g]])) {
        genes.meta.resist.subA[[g]] <- NULL
    }
}
length(genes.meta.resist.subA)

# saveRDS(genes.meta.resist.subA, file='genes-meta-resist-subA.rds')


# summarize results in a data frame
res.meta.resist.subA <- lapply(names(genes.meta.resist.subA), 
                                                      FUN = function(x) {
    results.meta(genes.meta.resist.subA, x)
})
names(res.meta.resist.subA) <- names(genes.meta.resist.subA)
res.meta.resist.subA.df <- do.call(rbind.data.frame, res.meta.resist.subA)
dim(res.meta.resist.subA.df)

# multiple testing correction: BH
res.meta.resist.subA.df$adj.pval <- p.adjust(res.meta.resist.subA.df$pval, 
                                              method = "BH",
                                      n = length(res.meta.resist.subA.df$pval))

# multiple testing correction: Bonferroni
res.meta.resist.subA.df$adj.pval.bonfer <- p.adjust(
                   res.meta.resist.subA.df$pval, method = "bonferroni", 
                   n = length(res.meta.resist.subA.df$pval))

# save saveRDS(res.meta.resist.subA.df, file='res-meta-resist-subA-df.rds')


# -------------------------- metaanalyses results ------------------------------

# metaanalyses and experiments correlation 

## experiments

# save the LFC values of the experiments in a new variable
lfc.genes.subA <- lapply(data.meta.resist.subA, '[[', 'logFC')
lfc.genes.subA <- do.call('rbind.data.frame', lfc.genes.subA)
lfc.genes.subA <- sapply(lfc.genes.subA, as.numeric)
colnames(lfc.genes.subA) <- data.meta.resist.subA[[1]]$studylab
rownames(lfc.genes.subA) <- names(data.meta.resist.subA)
lfc.genes.subA <- lfc.genes.subA[names(genes.meta.resist.subA), ]
dim(lfc.genes.subA)


## metaanalyses

# save the LFC values of the metaanalyses in a new variable
lfc.meta.subA <- lapply(genes.meta.resist.subA, '[[', 'TE.random')
lfc.meta.subA <- do.call('rbind.data.frame', lfc.meta.subA)
lfc.meta.subA <- sapply(lfc.meta.subA, as.numeric)
colnames(lfc.meta.subA) <- 'lfc.meta.subA'
rownames(lfc.meta.subA) <- names(genes.meta.resist.subA)

# merge both
lfc.genes.subA <- merge(lfc.genes.subA, lfc.meta.subA, by = 0)
rownames(lfc.genes.subA) <- lfc.genes.subA$Row.names
lfc.genes.subA <- lfc.genes.subA[ , -1]


# now take only DEGs in order to reduce noise in correlations and plots
dim(lfc.genes.subA)
lfc.genes.subA.05 <- lfc.genes.subA[rownames(res.meta.resist.subA.df[
                                res.meta.resist.subA.df$adj.pval <= 0.05, ]), ]
dim(lfc.genes.subA.05)

# correlation
lfc.genes.subA.05.cor <- lfc.genes.subA.05
for (i in 1:5) {
    lfc.genes.subA.05.cor[lfc.genes.subA.05.cor[, 6] < 0, i] <- -1 * 
                    (lfc.genes.subA.05.cor[lfc.genes.subA.05.cor[, 6] < 0, i])
}

# change meta downregulated DEGs sign first
lfc.genes.subA.05.cor[lfc.genes.subA.05.cor[, 6] < 0, 6] <- -1 * 
  (lfc.genes.subA.05.cor[lfc.genes.subA.05.cor[, 6] < 0, 6])

# compute correlations
first.cor <- cor.test(lfc.genes.subA.05.cor[, 1], 
                      lfc.genes.subA.05.cor[, "lfc.meta.subA"], 
                      method = "pearson")
cor.res.subA <- data.frame(cor = first.cor$estimate, pval = first.cor$p.value)

for (i in 2:5) {
    cor <- cor.test(lfc.genes.subA.05.cor[, i], 
                    lfc.genes.subA.05.cor[, "lfc.meta.subA"], 
                    method = "pearson")
    cor.res.subA <- rbind(cor.res.subA, c(cor$estimate, cor$p.value))
}

dim(cor.res.subA)
cbind(studylab = genes.meta.resist.subA[[1]]$data$studylab, cor.res.subA)

# plot DEGs 
par(mfrow = c(2, 3))
for (i in 1:5) {
    plot(lfc.genes.subA.05.cor[, i], 
         lfc.genes.subA.05.cor[, "lfc.meta.subA"],
         main = colnames(lfc.genes.subA.05.cor)[i], 
         ylab = "meta lfc", xlab = "individual lfc",
        pch = 21, col = "darkorchid4", bg = "darkorchid2", cex = 1, lwd = 1)
}
par(mfrow = c(1, 1))


# comparison with RankProd and RankSum

# load rankprod and ranksum results
rp.results.subA.tables <- readRDS("data/rp-results-subA-tables.rds")
rs.results.subA.tables <- readRDS("data/rs-results-subA-tables.rds")

RP.genes.subA <- c(rownames(rp.results.subA.tables$Table1), 
                    rownames(rp.results.subA.tables$Table2))
RS.genes.subA <- c(rownames(rs.results.subA.tables$Table1), 
                    rownames(rs.results.subA.tables$Table2))

# now take DEGs that were obtained with values of expression in all samples 
# in TE-MA
res.meta.subA.df.05.5st <- res.meta.resist.subA.df[
  res.meta.resist.subA.df$num.data == 5 & 
    res.meta.resist.subA.df$adj.pval < 0.05, ]

# venn diagram
ggVennDiagram(list(ES_MA = rownames(res.meta.subA.df.05.5st), 
                   RankProd_MA = RP.genes.subA,
                   RankSum_MA = RS.genes.subA), 
              color = "black", lwd = 0.8, lty = 1, label = "count", 
              label_alpha = 0, 
              category.names = c("TE-MA", "RankProd", "RankSum"), 
              label_size = 7) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  scale_color_manual(values = rep("darkgray", 3)) + 
  theme_void() + theme(legend.position = "none")


# compute correlations using DEGs selected by the 3 MA
degs.5st <- intersect(rownames(res.meta.subA.df.05.5st), 
                      intersect(RP.genes.subA, RS.genes.subA))
lfc.genes.subA.05.cor.5st <- lfc.genes.subA.05.cor[degs.5st, ]
dim(lfc.genes.subA.05.cor.5st)

first.cor <- cor.test(lfc.genes.subA.05.cor.5st[, 1], 
                      lfc.genes.subA.05.cor.5st[, "lfc.meta.subA"], 
                      method = "pearson")
cor.res.degs.5st <- data.frame(cor = first.cor$estimate, 
                               pval = first.cor$p.value)
for (i in 2:5) {
    cor <- cor.test(lfc.genes.subA.05.cor.5st[, i], 
                    lfc.genes.subA.05.cor.5st[, "lfc.meta.subA"], 
                    method = "pearson")
    cor.res.degs.5st <- rbind(cor.res.degs.5st, c(cor$estimate, cor$p.value))
}

dim(cor.res.degs.5st)
cor.res.degs.5st

# plot DEGs: experiment vs. TE-MA LFC
par(mfrow = c(2, 3))
for (i in 1:5) {
    plot(lfc.genes.subA.05.cor.5st[, i], 
         lfc.genes.subA.05.cor.5st[, "lfc.meta.subA"],
        main = colnames(lfc.genes.subA.05.cor.5st)[i], ylab = "meta lfc", 
        xlab = "individual lfc", pch = 21, col = "darkorchid4", 
        bg = "darkorchid2", cex = 1, lwd = 1)
}
par(mfrow = c(1, 1))


# select DEGs in MA and in experiments and then compare medians

## metaanalyses
res.meta.subA.df.05 <- res.meta.resist.subA.df[degs.5st, ]
dim(res.meta.subA.df.05)[1]/dim(res.meta.resist.subA.df)[1] * 100
dim(res.meta.subA.df.05)

# upregulated genes
sum(res.meta.subA.df.05$lfc > 0)  # n upregulated genes
median(res.meta.subA.df.05[res.meta.subA.df.05$lfc > 0, ]$lfc) 

# downregulated genes
sum(res.meta.subA.df.05$lfc < 0)  # n downregulated genes
median(res.meta.subA.df.05[res.meta.subA.df.05$lfc < 0, ]$lfc)

## experiments
degs.lfc.indiv.subA <- merge(genes.singleton[
           genes.singleton$adj.P.Val < 0.05, "logFC", drop = FALSE], 
           genes.ho.A375[genes.ho.A375$adj.P.Val < 0.05, "logFC", drop = FALSE],
           all = TRUE, by = 0) %>%
    merge(genes.song.m229[genes.song.m229$adj.P.Val < 0.05, "logFC", 
                   drop = FALSE], all = TRUE, by.x = "Row.names", by.y = 0) %>%
    merge(genes.song.m238[genes.song.m238$adj.P.Val < 0.05, "logFC", 
                   drop = FALSE], all = TRUE, by.x = "Row.names", by.y = 0) %>%
    merge(genes.song.skmel[genes.song.skmel$adj.P.Val < 0.05, "logFC", 
                        drop = FALSE], all = TRUE, by.x = "Row.names", by.y = 0)

rownames(degs.lfc.indiv.subA) <- degs.lfc.indiv.subA$Row.names
degs.lfc.indiv.subA <- degs.lfc.indiv.subA[, -1]
rownames(degs.lfc.indiv.subA) <- as.vector(mapply(
    rownames(degs.lfc.indiv.subA), FUN = function(x) {gsub("\\.\\d+", "", x)}))
dim(degs.lfc.indiv.subA)

# compute mean of medians
md.up.subA <- c()
md.down.subA <- c()
for (i in 1:5) {
    md.up.subA <- c(md.up.subA, 
                  median(degs.lfc.indiv.subA[degs.lfc.indiv.subA[, i] > 0, i], 
                         na.rm = TRUE))
    md.down.subA <- c(md.down.subA, 
                  median(degs.lfc.indiv.subA[degs.lfc.indiv.subA[, i] < 0, i], 
                         na.rm = TRUE))
}
mean(md.up.subA)
mean(md.down.subA)

# plot distributions
p.degs.subA.df <- rbind(data.frame(LFC = res.meta.subA.df.05$lfc, 
                                    COND = rep('MA', 
                                            length(res.meta.subA.df.05$lfc))),
                         data.frame(LFC = degs.lfc.indiv.subA[ , 1], 
                                    COND = rep('EXP1', 
                                          length(degs.lfc.indiv.subA[ , 1]))),
                         data.frame(LFC = degs.lfc.indiv.subA[ , 2], 
                                    COND = rep('EXP2', 
                                          length(degs.lfc.indiv.subA[ , 2]))),
                         data.frame(LFC = degs.lfc.indiv.subA[ , 3], 
                                    COND = rep('EXP3', 
                                          length(degs.lfc.indiv.subA[ , 3]))),
                         data.frame(LFC = degs.lfc.indiv.subA[ , 4], 
                                    COND = rep('EXP4', 
                                          length(degs.lfc.indiv.subA[ , 4]))),
                         data.frame(LFC = degs.lfc.indiv.subA[ , 5], 
                                    COND = rep('EXP5', 
                                           length(degs.lfc.indiv.subA[ , 5]))))


p.degs.subA.df$COND <- factor(p.degs.subA.df$COND, levels = c("EXP1", "EXP2", 
                                 "EXP3", "EXP4", "EXP5", "MA"), ordered = TRUE)
p.degs.subA.df$EXPRESSION <- NA
p.degs.subA.df$EXPRESSION[p.degs.subA.df$LFC > 0] <- "UP"
p.degs.subA.df$EXPRESSION[p.degs.subA.df$LFC < 0] <- "DOWN"

p.degs.subA <- ggplot(subset(p.degs.subA.df, !is.na(LFC)), 
                       aes(EXPRESSION, LFC, fill = COND)) + 
  geom_boxplot(alpha = 1, outlier.shape = NA) + theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 20), 
          axis.title = element_text(size = 12), 
          axis.text = element_text(size = 12), legend.position = "none") + 
  labs(y = "LFC", x = "") + scale_x_discrete(labels = c("DOWN", "UP")) + 
  scale_fill_manual(values = c(rep("gray68", 5), "limegreen")) + 
  geom_hline(yintercept = c(mean(md.up.subA), mean(md.down.subA)), 
             col = "red", lty = 5) + ylim(c(-6, 6))

p.degs.subA


# ----------------- enrichment analysis - clusterProfiler ----------------------

# GO BP
ans.go.subA <- enrichGO(gene = degs.5st, keyType = "ENSEMBL", ont = "BP", 
                         OrgDb = "org.Hs.eg.db", 
                         universe = rownames(res.meta.resist.subA.df), 
                         readable = TRUE, pAdjustMethod = "BH",
                         pvalueCutoff = 0.05)
tab.go.subA <- as.data.frame(ans.go.subA)
head(tab.go.subA)

# plot
p.subA <- barplot(ans.go.subA, showCategory = 15)
p.subA
p.subA.1.2 <- dotplot(ans.go.subA, showCategory = 15) + ggtitle("GO BP")
p.subA.1.2

# GO CC
ans.go.subA.2 <- enrichGO(gene = degs.5st, keyType = "ENSEMBL", ont = "CC", 
                           OrgDb = "org.Hs.eg.db",
                           universe = rownames(res.meta.resist.subA.df), 
                           readable = TRUE, pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)
tab.go.subA.2 <- as.data.frame(ans.go.subA.2)
head(tab.go.subA.2)

# plot
p.subA.2 <- barplot(ans.go.subA.2, showCategory = 15)
p.subA.2
p.subA.2.1 <- dotplot(ans.go.subA, showCategory = 15) + ggtitle("GO BP")
p.subA.2.1

# ---------------------------- GSEA results ------------------------------------

# plot results tft
tft.gsea.pos.subA <- read_tsv("data/subA-tft-pos.tsv")
tft.gsea.neg.subA <- read_tsv("data/subA-tft-neg.tsv")
tft.gsea.pos.subA <- tft.gsea.pos.subA[1:5, c(1, 4, 6, 8)]
tft.gsea.neg.subA <- tft.gsea.neg.subA[1:5, c(1, 4, 6, 8)]

# positive NES
tft.pos.subA <- ggplot(tft.gsea.pos.subA, 
                        aes(x = NES, y = `FDR q-val`,  
                            color = `FDR q-val` < 0.05, label = NAME)) + 
  geom_point(size = tft.gsea.pos.subA$SIZE/100) + theme_light() +
  scale_color_manual(values = c("black", "red")) + 
  geom_hline(yintercept = 0.05, col = "gray10", lty = 2) + ylab("") + 
  xlab("NES") + ylim(c(-0.02, 0.27)) + 
  theme(legend.position = "none", axis.title = element_text(size = 17), 
        axis.text = element_text(size = 17)) +
    geom_label_repel(position = position_dodge(width = 6), size = 4, 
                     max.overlaps = 100)

tft.pos.subA

# negative NES
tft.neg.subA <- ggplot(tft.gsea.neg.subA, 
                        aes(x = NES, y = `FDR q-val`, 
                            color = `FDR q-val` < 0.05, label = NAME)) + 
  geom_point(size = tft.gsea.neg.subA$SIZE/100) + theme_light() +
    scale_color_manual(values = c("black", "darkblue")) + 
  geom_hline(yintercept = 0.05, col = "gray10", lty = 2) + ylab("FDR") + 
  xlab("NES") + ylim(c(-0.02, 0.27)) + 
  theme(legend.position = "none", axis.title = element_text(size = 17),
        axis.text = element_text(size = 17)) +
    geom_label_repel(position = position_dodge(width = 6), 
                     size = 4, max.overlaps = 100)

tft.neg.subA

# plot
tft.plot.subA <- ggarrange(tft.neg.subA, tft.pos.subA, ncol = 2, nrow = 1)

tft.plot.subA


# plot results hallmarks
hallm.gsea.pos.subA <- read_tsv("data/subA-hallmarks-pos.tsv")
hallm.gsea.neg.subA <- read_tsv("data/subA-hallmarks-neg.tsv")
hallm.gsea.pos.subA <- hallm.gsea.pos.subA[1:15, c(1, 4, 6, 8)]
hallm.gsea.neg.subA <- hallm.gsea.neg.subA[1:7, c(1, 4, 6, 8)]

# positive NES
hallm.pos.subA <- ggplot(hallm.gsea.pos.subA, 
                          aes(x = NES, y = `FDR q-val`, 
                              color = `FDR q-val` < 0.05, label = 1:15)) + 
  geom_point(size = hallm.gsea.pos.subA$SIZE/100) + theme_light() +
  scale_color_manual(values = c("black", "red")) + 
  geom_hline(yintercept = 0.05, col = "gray10", lty = 2) + ylab("") + 
  xlab("NES") + ylim(c(-0.02, 0.175)) + 
  theme(legend.position = "none", axis.title = element_text(size = 17), 
        axis.text = element_text(size = 17)) +
    geom_label_repel(position = position_dodge(width = 6), size = 4, 
                     max.overlaps = 100)

hallm.pos.subA

# negative NES
hallm.neg.subA <- ggplot(hallm.gsea.neg.subA, 
                          aes(x = NES, y = `FDR q-val`, 
                              color = `FDR q-val` < 0.05, 
                              label = sub("HALLMARK_", "", NAME))) + 
  geom_point(size = hallm.gsea.neg.subA$SIZE/100) +
    theme_light() + scale_color_manual(values = c("black", "darkblue")) + 
  geom_hline(yintercept = 0.05, col = "gray10", lty = 2) + ylab("FDR") + 
  xlab("NES") + ylim(c(-0.02, 0.175)) + 
  theme(legend.position = "none", axis.title = element_text(size = 17), 
        axis.text = element_text(size = 17)) +
    geom_label_repel(position = position_dodge(width = 6), size = 4, 
                     max.overlaps = 100)

hallm.neg.subA

# plot
hallm.plot.subA <- ggarrange(hallm.neg.subA, hallm.pos.subA, ncol = 2, 
                              nrow = 1)

hallm.plot.subA

# save terms that are not shown
write.table(paste(1:15, sub("HALLMARK_", "", hallm.gsea.pos.subA$NAME), 
                  sep = ". "), row.names = FALSE, col.names = FALSE, 
            quote = FALSE, file = "plots-v2/gsea-hallmarks-pos-subA-NAME.txt")


# plot results kegg
kegg.gsea.pos.subA <- read_tsv("data/subA-kegg-pos.tsv")
kegg.gsea.neg.subA <- read_tsv("data/subA-kegg-neg.tsv")
kegg.gsea.pos.subA <- kegg.gsea.pos.subA[1:5, c(1, 4, 6, 8)]
kegg.gsea.neg.subA <- kegg.gsea.neg.subA[1:10, c(1, 4, 6, 8)]

# positive NES
kegg.pos.subA <- ggplot(kegg.gsea.pos.subA, 
                         aes(x = NES, y = `FDR q-val`, 
                             color = `FDR q-val` < 0.05, 
                             label = sub("KEGG_", "", NAME))) + 
  geom_point(size = kegg.gsea.pos.subA$SIZE/100) +
    theme_light() + scale_color_manual(values = c("black", "red")) + 
  geom_hline(yintercept = 0.05, col = "gray10", lty = 2) + ylab("") + 
  xlab("NES") + ylim(c(-0.02, 0.225)) + xlim(c(-2, 6)) + 
  theme(legend.position = "none", axis.title = element_text(size = 17), 
        axis.text = element_text(size = 17)) +
  geom_label_repel(position = position_dodge(width = 6), 
                   size = 4, max.overlaps = 100)

kegg.pos.subA

# negative NES
kegg.neg.subA <- ggplot(kegg.gsea.neg.subA, 
                         aes(x = NES, y = `FDR q-val`, 
                             color = `FDR q-val` < 0.05, 
                             label = sub("KEGG_", "", NAME))) + 
  geom_point(size = kegg.gsea.neg.subA$SIZE/100) +
    theme_light() + scale_color_manual(values = c("black", "darkblue")) + 
  geom_hline(yintercept = 0.05, col = "gray10", lty = 2) + ylab("FDR") + 
  xlab("NES") + ylim(c(-0.02, 0.225)) + xlim(c(-6, 2)) + 
  theme(legend.position = "none", axis.title = element_text(size = 17),
        axis.text = element_text(size = 17)) + 
  geom_label_repel(position = position_dodge(width = 6), 
                   size = 4, max.overlaps = 100)
kegg.neg.subA

# plot
kegg.plot.subA <- ggarrange(kegg.neg.subA, kegg.pos.subA, ncol = 2, nrow = 1)

kegg.plot.subA


# -------------------------------- TFEA.ChIP -----------------------------------

# retrieve data
tfea.data.subA <- res.meta.resist.subA.df[, c("lfc", "pval", "adj.pval")]
list.genes.entrez <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id"), 
                           values = rownames(tfea.data.subA), mart = mart)
tfea.data.subA <- merge(list.genes.entrez, tfea.data.subA, 
                         by.x = "ensembl_gene_id", by.y = 0, all.y = TRUE)
tfea.data.subA <- tfea.data.subA[order(tfea.data.subA$lfc,
                                         decreasing = TRUE), ]
tfea.data.subA <- tfea.data.subA[, -1]
colnames(tfea.data.subA) <- c("Genes", "log2FoldChange", "pvalue", "pval.adj")

# UP-regulated genes extract vector with names of upregulated genes
genes.upreg.subA <- Select_genes(tfea.data.subA, min_LFC = 1)

# extract vector with names of non-responsive genes
genes.control.subA <- Select_genes(tfea.data.subA, min_pval = 0.5, 
                                    max_pval = 1, min_LFC = -0.25, 
                                    max_LFC = 0.25)

## association analysis

# generates list of contingency tables, one per dataset
CM_list_UP_subA <- contingency_matrix(genes.upreg.subA)  

# generates list of p-values and OR from association test
pval_mat_UP_subA <- getCMstats(CM_list_UP_subA)  
head(pval_mat_UP_subA)

TF_ranking_up_subA <- rankTFs(pval_mat_UP_subA, rankMethod = "gsea",
                                                               makePlot = TRUE)
TF_ranking_up_subA[["TFranking_plot"]]

head(TF_ranking_up_subA[["TF_ranking"]])
tf.ranking.up.subA <- na.omit(TF_ranking_up_subA[["TF_ranking"]][
                             TF_ranking_up_subA[["TF_ranking"]]$pVal < 0.05, ])

plot_CM(pval_mat_UP_subA)  # plot p-values against ORs

tfs.up.subA <- rownames(tf.ranking.up.subA)
names(tfs.up.subA) <- rownames(tf.ranking.up.subA)
col.up.subA <- sample(viridis::turbo(length(tfs.up.subA)), 
                       length(tfs.up.subA))

# plot p-values against ORs highlighting indicated TFs
plot_CM(pval_mat_UP_subA, specialTF = tfs.up.subA, TF_colors = col.up.subA)  
TF_ranking_up_subA[["TFranking_plot"]]


# DOWN-regulated genes extract vector with names of upregulated genes
genes.downreg.subA <- Select_genes(tfea.data.subA, max_LFC = -1)

# extract vector with names of non-responsive genes
genes.control.subA <- Select_genes(tfea.data.subA, min_pval = 0.5, max_pval = 1,
    min_LFC = -0.25, max_LFC = 0.25)

## association analysis

# generates list of contingency tables, one per dataset
CM_list_DOWN_subA <- contingency_matrix(genes.downreg.subA) 

# generates list of p-values and OR from association test
pval_mat_DOWN_subA <- getCMstats(CM_list_DOWN_subA)  
head(pval_mat_DOWN_subA)

TF_ranking_down_subA <- rankTFs(pval_mat_DOWN_subA, rankMethod = "gsea", 
                                                               makePlot = TRUE)
TF_ranking_down_subA[["TFranking_plot"]]
tf.ranking.down.subA <- na.omit(TF_ranking_down_subA[["TF_ranking"]][
                           TF_ranking_down_subA[["TF_ranking"]]$pVal < 0.05, ])

head(TF_ranking_down_subA[["TF_ranking"]])

plot_CM(pval_mat_DOWN_subA)  # plot p-values against ORs

tfs.down.subA <- rownames(tf.ranking.down.subA)
names(tfs.down.subA) <- rownames(tf.ranking.down.subA)
col.down.subA <- sample(viridis::turbo(length(tfs.down.subA)), 
                                                        length(tfs.down.subA))

# plot p-values against ORs highlighting indicated TFs
plot_CM(pval_mat_DOWN_subA, specialTF = tfs.down.subA, 
                                                    TF_colors = col.down.subA)  
TF_ranking_down_subA[["TFranking_plot"]]


# ----------------------- actomyosin network analysis --------------------------

## MA volcano plot

# add new column indicating if the gene is down, up or not regulated
res.meta.resist.subA.df$diffexp <- "NO"
res.meta.resist.subA.df$diffexp[res.meta.resist.subA.df$lfc > 1.5 & 
                               res.meta.resist.subA.df$adj.pval < 0.05] <- "UP"
res.meta.resist.subA.df$diffexp[res.meta.resist.subA.df$lfc < -1.5 & 
                             res.meta.resist.subA.df$adj.pval < 0.05] <- "DOWN"

# choose colors
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")

# plot
p.subA <- ggplot(data = res.meta.resist.subA.df, 
                  aes(x = lfc, y = -log10(adj.pval), fill = diffexp)) + 
  geom_point(shape = 21, size = 2.2) + theme_light() + 
  geom_vline(xintercept = c(-1.5, 1.5), col = "darkgray", lty = 2) + 
  geom_hline(yintercept = -log10(0.05), col = "darkgray", lty = 2) + 
  scale_fill_manual(values = c("blue", "black", "red")) + 
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 20), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14)) + 
  ylab("-log10(adj. p-value)") + xlab("LFC") +
  scale_x_continuous(breaks = seq(-10, 10, by = 2), limits = c(-10, 10)) + 
  ylim(c(0, 60))

p.subA


## select actomyosin cyotskeleton genes and plot volcano
res.meta.resist.subA.df.sy <- merge(list.genes, res.meta.resist.subA.df, 
                              all.y = TRUE, by.x = "ensembl_gene_id", by.y = 0)

# add a column indicating if the gene is related to the actomyosin network
res.meta.resist.subA.df.sy$actomy <- "NO"
res.meta.resist.subA.df.sy$actomy[
            (res.meta.resist.subA.df.sy$hgnc_symbol %in% actomyosin)] <- "DEG"
res.meta.resist.subA.df.sy$actomy[
                   (res.meta.resist.subA.df.sy$hgnc_symbol %in% actomyosin) & 
                     (res.meta.resist.subA.df.sy$lfc > 1.5) & 
                     (res.meta.resist.subA.df.sy$adj.pval < 0.05)] <- "UP"
res.meta.resist.subA.df.sy$actomy[
  (res.meta.resist.subA.df.sy$hgnc_symbol %in% actomyosin) & 
    (res.meta.resist.subA.df.sy$lfc < -1.5) & 
    (res.meta.resist.subA.df.sy$adj.pval < 0.05)] <- "DOWN"

# choose colors
mycolors <- c("darkgray", "blue", "red", "gray40")
names(mycolors) <- c("DOWN", "UP", "NO", "DEG")
res.meta.resist.subA.df.sy$actomy <- factor(res.meta.resist.subA.df.sy$actomy,
    levels = c("NO", "DOWN", "UP", "DEG"))

# order data
res.meta.resist.subA.df.sy <- res.meta.resist.subA.df.sy[
                                    order(res.meta.resist.subA.df.sy$actomy), ]
res.meta.resist.subA.df.sy$delabel <- NA
res.meta.resist.subA.df.sy$delabel[
                    res.meta.resist.subA.df.sy$actomy != "NO"] <- 
                              res.meta.resist.subA.df.sy$hgnc_symbol[
                                     res.meta.resist.subA.df.sy$actomy != "NO"]

# plot
p2.volcano.subA <- ggplot(data = res.meta.resist.subA.df.sy, 
                           aes(x = lfc, y = -log10(adj.pval), 
                               col = actomy, label = delabel)) + 
  geom_point(shape = 21, size = 2.2, 
             fill = mycolors[res.meta.resist.subA.df.sy$actomy]) + 
  theme_light() + 
  scale_colour_manual(values = c("darkgray", "black", "black", "black")) + 
  labs(x = "LFC", y = "log10(p-valor ajustado)") + 
  theme(legend.position = "none",  
        plot.title = element_text(hjust = 0.5, size = 20), 
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12)) + 
  scale_x_continuous(breaks = seq(-10, 10, by = 2)) + 
  geom_vline(xintercept = c(-1.5, 1.5), col = "gray15", lty = 2) + 
  geom_hline(yintercept = -log10(0.05), col = "gray15", lty = 2) + 
  geom_label_repel(data = subset(res.meta.resist.subA.df.sy, adj.pval < 0.05 & 
                                  (lfc > 1.5 | lfc < -1.5)), max.overlaps = 100, 
                   position = position_nudge_repel(y = 10)) +
    ylim(c(0, 60)) + ggtitle("Subgroup A")

p2.volcano.subA

## heatmap
actomyosin.ensemblid <- actomyosin
names(actomyosin.ensemblid) <- actomyosin
actomyosin.ensemblid <- merge(list.genes, actomyosin.ensemblid, 
                              by.x = "hgnc_symbol", by.y = 1)
lfc.genes.subA.sy.actomy <- merge(actomyosin.ensemblid, lfc.genes.subA, 
                               by.x = "ensembl_gene_id", by.y = 0, all.x = TRUE)

# remove one repeated line
lfc.genes.subA.sy.actomy <- lfc.genes.subA.sy.actomy[-70, ]

# change row names
rownames(lfc.genes.subA.sy.actomy) <- paste(
                                      lfc.genes.subA.sy.actomy$ensembl_gene_id,
                                      lfc.genes.subA.sy.actomy$hgnc_symbol, 
                                      sep = " - ")
lfc.genes.subA.sy.actomy <- lfc.genes.subA.sy.actomy[, -c(1, 2)]

# plot
p.heatmap.acto.subA <- pheatmap(na.omit(lfc.genes.subA.sy.actomy[, -6]), 
                                 col = col.hm, cluster_rows = TRUE, 
                                 cluster_cols = TRUE, show_rownames = TRUE, 
                                 fontsize_col = 10,
    breaks = seq(-4, 4, length = 50))

p.heatmap.acto.subA


## forest plots of some genes
# MYL9
forest.meta(genes.meta.resist.subA[['ENSG00000101335']],
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2',
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# MYH9
forest.meta(genes.meta.resist.subA[['ENSG00000100345']],
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2',
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# MYL7
forest.meta(genes.meta.resist.subA[['ENSG00000106631']],
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2',
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# ITGB3
forest.meta(genes.meta.resist.subA[['ENSG00000259207']],
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2',
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# ITGA9
forest.meta(genes.meta.resist.subA[['ENSG00000144668']],
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2',
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# ACTA2
forest.meta(genes.meta.resist.subA[['ENSG00000107796']],
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2',
            print.I2 = FALSE,
            print.pval.Q = FALSE)

# MYH14
forest.meta(genes.meta.resist.subA[['ENSG00000105357']],
            sortvar = logFC,
            print.tau2 = FALSE,
            leftlabs = c("Study", "LFC", "SE"),
            rightlabs = c("LFC", "95%-CI", "Weight"),
            col.diamond.random = 'blue2',
            print.I2 = FALSE,
            print.pval.Q = FALSE)



# ------------------Subgroup B metaanalyses: Song M395R, -----------------------
# --------------------- Ho 451LuR and Berico MM0074R ---------------------------

# prepare the data to metaanalyze
search.gene.subB <- function(gene) {
    mat <- matrix(ncol = 6, nrow = 0)
    df <- data.frame(mat)

    # ho 451Lu
    if (nrow(genes.ho.lu.lfc[rownames(genes.ho.lu.lfc) == gene, ]) == 1) {
        a <- genes.ho.lu.lfc[rownames(genes.ho.lu.lfc) == gene, ][1, 1]
        b <- genes.ho.lu.lfc[rownames(genes.ho.lu.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Ho, 451LuR", "Ho", "451Lu", "R"))
    } else {
        df <- rbind(df, c(NA, NA, "Ho, 451LuR", "Ho", "451Lu", "R"))
    }

    # berico
    if (nrow(genes.berico.mm.lfc[
                               rownames(genes.berico.mm.lfc) == gene, ]) == 1) {
        a <- genes.berico.mm.lfc[rownames(genes.berico.mm.lfc) == gene, ][1, 1]
        b <- genes.berico.mm.lfc[rownames(genes.berico.mm.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Berico, MM074R", "Berico", "MM074", "R"))
    } else {
        df <- rbind(df, c(NA, NA, "Berico, MM074R", "Berico", "MM074", "R"))
    }

    # song M395
    if (nrow(genes.song.m395.lfc[
                               rownames(genes.song.m395.lfc) == gene, ]) == 1) {
        a <- genes.song.m395.lfc[rownames(genes.song.m395.lfc) == gene, ][1, 1]
        b <- genes.song.m395.lfc[rownames(genes.song.m395.lfc) == gene, ][1, 2]
        df <- rbind(df, c(a, b, "Song, M395R", "Song", "M395", "R"))
    } else {
        df <- rbind(df, c(NA, NA, "Song, M395R", "Song", "M395", "R"))
    }

    # end of search
    colnames(df) <- c("logFC", "lfcSE", "studylab", "author", "line", 
                      "phenotype")
    return(df)
}

# loop to iterate over all the genes using the function
data.meta.resist.subB <- lapply(all.genes, FUN = function(x) {
    search.gene.subB(x)
})
names(data.meta.resist.subB) <- all.genes

length(data.meta.resist.subB)

# remove null results
for (g in all.genes) {
    if (all(is.na(data.meta.resist.subB[[g]][, c(1, 2)]))) {
        data.meta.resist.subB[[g]] <- NULL
    }
}
length(data.meta.resist.subB)

# remove genes with info from only 1 study
for (g in names(data.meta.resist.subB)) {
    if (sum(is.na(data.meta.resist.subB[[g]])) == 4) {
        data.meta.resist.subB[[g]] <- NULL
    }
}
length(data.meta.resist.subB)

# remove digits after dots and save
names(data.meta.resist.subB) <- as.vector(mapply(names(data.meta.resist.subB),
    FUN = function(x) {gsub("\\.\\d+", "", x)}))

# saveRDS(data.meta.resist.subB, file = 'data-meta-resist-subB.rds')


# metaanalyses
counter <- 0
total <- length(data.meta.resist.subB)
genes.meta.resist.subB <- lapply(names(data.meta.resist.subB), 
                                                             FUN = function(x) {
    counter <<- counter + 1
    print(paste0("iteration number ", counter, "; total ", total))
    metaanalyses(data.meta.resist.subB, x)
})

names(genes.meta.resist.subB) <- names(data.meta.resist.subB)

# remove null results
length(genes.meta.resist.subB)

for (g in names(data.meta.resist.subB)) {
    if (is.null(genes.meta.resist.subB[[g]])) {
        genes.meta.resist.subB[[g]] <- NULL
    }
}
length(genes.meta.resist.subB)

# summarize results in a data frame
res.meta.resist.subB <- lapply(names(genes.meta.resist.subB), 
                  FUN = function(x) {results.meta(genes.meta.resist.subB, x)})
names(res.meta.resist.subB) <- names(genes.meta.resist.subB)
res.meta.resist.subB.df <- do.call(rbind.data.frame, res.meta.resist.subB)
dim(res.meta.resist.subB.df)

# multiple testing correction: BH
res.meta.resist.subB.df$adj.pval <- p.adjust(res.meta.resist.subB.df$pval, 
                                              method = "BH",
    n = length(res.meta.resist.subB.df$pval))

# multiple testing correction: Bonferroni
res.meta.resist.subB.df$adj.pval.bonfer <- p.adjust(
                                                  res.meta.resist.subB.df$pval,
                                                  method = "bonferroni", 
                                      n = length(res.meta.resist.subB.df$pval))

# save saveRDS(res.meta.resist.subB.df, file='res-meta-resist-subB-df.rds')


# ------------------------- metaanalyses results -------------------------------

# metaanalyses and experiments correlation

## experiments

# save the LFC values of the experiments in a new variable
lfc.genes.subB <- lapply(data.meta.resist.subB, '[[', 'logFC')
lfc.genes.subB <- do.call('rbind.data.frame', lfc.genes.subB)
lfc.genes.subB <- sapply(lfc.genes.subB, as.numeric)
colnames(lfc.genes.subB) <- data.meta.resist.subB[[1]]$studylab
rownames(lfc.genes.subB) <- names(data.meta.resist.subB)
lfc.genes.subB <- lfc.genes.subB[names(genes.meta.resist.subB), ]
dim(lfc.genes.subB)


## metaanalyses

# save the LFC values of the metaanalyses in a new variable
lfc.meta.subB <- lapply(genes.meta.resist.subB, '[[', 'TE.random')
lfc.meta.subB <- do.call('rbind.data.frame', lfc.meta.subB)
lfc.meta.subB <- sapply(lfc.meta.subB, as.numeric)
colnames(lfc.meta.subB) <- 'lfc.meta.subB'
rownames(lfc.meta.subB) <- names(genes.meta.resist.subB)

# merge both
lfc.genes.subB <- merge(lfc.genes.subB, lfc.meta.subB, by = 0)
rownames(lfc.genes.subB) <- lfc.genes.subB$Row.names
lfc.genes.subB <- lfc.genes.subB[ , -1]


# now take only DEGS in order to reduce noise in correlations and plots
dim(lfc.genes.subB)
lfc.genes.subB.05 <- lfc.genes.subB[rownames(res.meta.resist.subB.df[
                                res.meta.resist.subB.df$adj.pval <= 0.05, ]), ]
dim(lfc.genes.subB.05)

# pearson cor with DEGs change downregulated genes sign
lfc.genes.subB.05.cor <- lfc.genes.subB.05
for (i in 1:3) {
    lfc.genes.subB.05.cor[lfc.genes.subB.05.cor[, 4] < 0, i] <- -1 * 
      (lfc.genes.subB.05.cor[lfc.genes.subB.05.cor[, 4] < 0, i])
}

# change meta downregulated DEGs sign
lfc.genes.subB.05.cor[lfc.genes.subB.05.cor[, 4] < 0, 4] <- -1 * 
  (lfc.genes.subB.05.cor[lfc.genes.subB.05.cor[, 4] < 0, 4])

# compute correlations
first.cor <- cor.test(lfc.genes.subB.05.cor[, 1], 
                      lfc.genes.subB.05.cor[, "lfc.meta.subB"], 
                      method = "pearson")
cor.res.subB <- data.frame(cor = first.cor$estimate, pval = first.cor$p.value)

for (i in 2:3) {
    cor <- cor.test(lfc.genes.subB.05.cor[, i], 
                    lfc.genes.subB.05.cor[, "lfc.meta.subB"],
                    method = "pearson")
    cor.res.subB <- rbind(cor.res.subB, c(cor$estimate, cor$p.value))
}

dim(cor.res.subB)
cbind(studylab = genes.meta.resist.subB[[1]]$data$studylab, cor.res.subB)

# plot DEGs 
par(mfrow = c(1, 3))
for (i in 1:3) {
    plot(lfc.genes.subB.05.cor[, i], 
         lfc.genes.subB.05.cor[, "lfc.meta.subB"],
        main = colnames(lfc.genes.subB.05.cor)[i], ylab = "meta lfc", 
        xlab = "individual lfc", pch = 21, col = "darkorchid4", 
        bg = "darkorchid2", cex = 1, lwd = 1)
}
par(mfrow = c(1, 1))


## comparison with RankProd and RankSum

# load rankprod and ranksum results
rp.results.subB.tables <- readRDS("data/rp-results-subB-tables.rds")
rs.results.subB.tables <- readRDS("data/rs-results-subB-tables.rds")

RP.genes.subB <- c(rownames(rp.results.subB.tables$Table1), 
                    rownames(rp.results.subB.tables$Table2))
RS.genes.subB <- c(rownames(rs.results.subB.tables$Table1), 
                    rownames(rs.results.subB.tables$Table2))

# now take DEGs that were obtained with values of expression in all samples 
# in TE-MA
res.meta.subB.df.05.3st <- res.meta.resist.subB.df[
  res.meta.resist.subB.df$num.data == 3 & 
    res.meta.resist.subB.df$adj.pval < 0.05, ]

# venn diagram
ggVennDiagram(list(ES_MA = rownames(res.meta.subB.df.05.3st), 
                   RankProd_MA = RP.genes.subB,
                   RankSum_MA = RS.genes.subB), 
              color = "black", lwd = 0.8, lty = 1, label = "count",
    label_alpha = 0, category.names = c("TE-MA", "RankProd", "RankSum"),
    label_size = 7) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  scale_color_manual(values = rep("darkgray", 3)) + theme_void() + 
  theme(legend.position = "none")

# now compute correlations selecting only DEGs detected by the 3 MA
degs.3st <- intersect(rownames(res.meta.subB.df.05.3st), 
                      intersect(RP.genes.subB,
                                RS.genes.subB))

lfc.genes.subB.05.cor.3st <- lfc.genes.subB.05.cor[degs.3st, ]
dim(lfc.genes.subB.05.cor.3st)
first.cor <- cor.test(lfc.genes.subB.05.cor.3st[, 1], 
                      lfc.genes.subB.05.cor.3st[, "lfc.meta.subB"], 
                      method = "pearson")
cor.res.degs.3st <- data.frame(cor = first.cor$estimate, 
                               pval = first.cor$p.value)

for (i in 2:3) {
    cor <- cor.test(lfc.genes.subB.05.cor.3st[, i], 
                    lfc.genes.subB.05.cor.3st[, "lfc.meta.subB"], 
                    method = "pearson")
    cor.res.degs.3st <- rbind(cor.res.degs.3st, c(cor$estimate, cor$p.value))
}

dim(cor.res.degs.3st)
cor.res.degs.3st

# plot DEGs: experiment vs. MA LFC
par(mfrow = c(1, 3))
for (i in 1:3) {
    plot(lfc.genes.subB.05.cor.3st[, i], 
         lfc.genes.subB.05.cor.3st[, "lfc.meta.subB"],
         main = colnames(lfc.genes.subB.05.cor.3st)[i], ylab = "meta lfc", 
         xlab = "individual lfc", pch = 21, col = "darkorchid4", 
         bg = "darkorchid2", cex = 1, lwd = 1)
}
par(mfrow = c(1, 1))


# histogram to compare these new correlations vs. what we obtained before
cor.res.compare <- data.frame(cor = cor.res.degs.8st[, 1], group = rep("R"))
cor.res.compare <- rbind(cor.res.compare, 
                         data.frame(cor = cor.res.degs.5st[, 1],
    group = rep("subA")), data.frame(cor = cor.res.degs.3st[, 1], 
                                      group = rep("subB")))

p.hist <- ggplot(data = cor.res.compare, aes(x = cor, fill = group)) + 
  geom_histogram(aes(y = 0.1 * ..density..), color = "black", alpha = 0.5, 
                 position = "identity", binwidth = 0.1) +
    scale_fill_brewer(palette = "Set1", breaks = c("R", "subA", "subB"), 
                      labels = c("All", "Subgroup A", "Subgroup B")) + 
  xlim(c(0, 1)) + theme_light() + 
  theme(plot.title = element_text(hjust = 0.6, size = 20), 
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12), legend.title = element_blank()) + 
  xlab("Pearson correlation coefficient") +
    ylab("Relative Frequency")

p.hist


# select DEGs from metaanalyses and from the experiments and compare medians

## metaanalyses
res.meta.subB.df.05 <- res.meta.resist.subB.df[degs.3st, ]
dim(res.meta.subB.df.05)[1]/dim(res.meta.resist.subB.df)[1] * 100
dim(res.meta.subB.df.05)

# upregulated genes
sum(res.meta.subB.df.05$lfc > 0)  # n upregulated genes
median(res.meta.subB.df.05[res.meta.subB.df.05$lfc > 0, ]$lfc)

# downregulated genes
sum(res.meta.subB.df.05$lfc < 0)  # n downregulated genes
median(res.meta.subB.df.05[res.meta.subB.df.05$lfc < 0, ]$lfc)

## experiments
degs.lfc.indiv.subB <- merge(genes.ho.lu[genes.ho.lu$adj.P.Val < 0.05, 
                                                         "logFC", drop = FALSE], 
                              genes.berico[genes.berico$adj.P.Val < 0.05, 
                                 "logFC", drop = FALSE], all = TRUE, by = 0) %>%
    merge(genes.song.m395[genes.song.m395$adj.P.Val < 0.05, 
                                  "logFC", drop = FALSE], all = TRUE,
                                             by.x = "Row.names", by.y = 0)
rownames(degs.lfc.indiv.subB) <- degs.lfc.indiv.subB$Row.names
degs.lfc.indiv.subB <- degs.lfc.indiv.subB[, -1]
rownames(degs.lfc.indiv.subB) <- as.vector(mapply(rownames(
                                                          degs.lfc.indiv.subB),
    FUN = function(x) {gsub("\\.\\d+", "", x)}))
dim(degs.lfc.indiv.subB)

# compute mean of medians
md.up.subB <- c()
md.down.subB <- c()
for (i in 1:3) {
    md.up.subB <- c(md.up.subB, median(degs.lfc.indiv.subB[
                               degs.lfc.indiv.subB[, i] > 0, i], na.rm = TRUE))
    md.down.subB <- c(md.down.subB, median(degs.lfc.indiv.subB[
                               degs.lfc.indiv.subB[, i] < 0, i], na.rm = TRUE))
}
mean(md.up.subB)
mean(md.down.subB)

# plot distributions
p.degs.subB.df <- rbind(data.frame(LFC = res.meta.subB.df.05$lfc, 
                                    COND = rep('MA', 
                                            length(res.meta.subB.df.05$lfc))),
                         data.frame(LFC = degs.lfc.indiv.subB[ , 1], 
                                    COND = rep('EXP1',
                                          length(degs.lfc.indiv.subB[ , 1]))),
                         data.frame(LFC = degs.lfc.indiv.subB[ , 2], 
                                    COND = rep('EXP2', 
                                          length(degs.lfc.indiv.subB[ , 2]))),
                         data.frame(LFC = degs.lfc.indiv.subB[ , 3], 
                                    COND = rep('EXP3', 
                                          length(degs.lfc.indiv.subB[ , 3]))))

p.degs.subB.df$COND <- factor(p.degs.subB.df$COND, 
                               levels = c("EXP1", "EXP2", "EXP3", "MA"), 
                               ordered = TRUE)
p.degs.subB.df$EXPRESSION <- NA
p.degs.subB.df$EXPRESSION[p.degs.subB.df$LFC > 0] <- "UP"
p.degs.subB.df$EXPRESSION[p.degs.subB.df$LFC < 0] <- "DOWN"

p.degs.subB <- ggplot(subset(p.degs.subB.df, !is.na(LFC)), 
                       aes(EXPRESSION, LFC, fill = COND)) + 
  geom_boxplot(alpha = 1, outlier.shape = NA) + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20), 
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12), legend.position = "none") + 
  labs(y = "LFC", x = "") + ylim(c(-4, 4)) +
  scale_x_discrete(labels = c("DOWN", "UP")) + 
  scale_fill_manual(values = c(rep("gray68", 3), "limegreen")) + 
  geom_hline(yintercept = c(mean(md.up.subB), mean(md.down.subB)), 
             col = "red", lty = 5)

p.degs.subB


# ------------------- enrichment analysis - clusterProfiler --------------------

# plot results tft
tft.gsea.pos.subB <- read_tsv("data/subB-tft-pos.tsv")
tft.gsea.neg.subB <- read_tsv("data/subB-tft-neg.tsv")
tft.gsea.pos.subB <- tft.gsea.pos.subB[1:3, c(1, 4, 6, 8)]
tft.gsea.neg.subB <- tft.gsea.neg.subB[1:5, c(1, 4, 6, 8)]

# positive NES
tft.pos.subB <- ggplot(tft.gsea.pos.subB, 
                        aes(x = NES, y = `FDR q-val`, 
                            color = `FDR q-val` < 0.05, label = NAME)) + 
  geom_point(size = tft.gsea.pos.subB$SIZE/100) + theme_light() +
  scale_color_manual(values = c("black", "red")) + 
  geom_hline(yintercept = 0.05, col = "gray10", lty = 2) + ylab("") + 
  xlab("NES") + ylim(c(-0.02, 0.31)) + 
  theme(legend.position = "none", axis.title = element_text(size = 17),
    axis.text = element_text(size = 17)) + 
  geom_label_repel(position = position_dodge(width = 6), 
                   size = 4, max.overlaps = 100)

tft.pos.subB


# negative NES
tft.neg.subB <- ggplot(tft.gsea.neg.subB, 
                        aes(x = NES, 
                            y = `FDR q-val`, color = `FDR q-val` < 0.05, 
                            label = NAME)) + 
  geom_point(size = tft.gsea.neg.subB$SIZE/100) + theme_light() +
  scale_color_manual(values = c("black", "darkblue")) + 
  geom_hline(yintercept = 0.05, col = "gray10", lty = 2) + ylab("FDR") + 
  xlab("Puntuación normalizada") + ylim(c(-0.02, 0.31)) + 
  theme(legend.position = "none", axis.title = element_text(size = 17), 
        axis.text = element_text(size = 17)) + 
  geom_label_repel(position = position_dodge(width = 6), 
                   size = 4, max.overlaps = 100)
tft.neg.subB

# plot
tft.plot.subB <- ggarrange(tft.neg.subB, tft.pos.subB, ncol = 2, nrow = 1)

tft.plot.subB


# plot results hallmarks
hallm.gsea.pos.subB <- read_tsv("data/subB-hallmarks-pos.tsv")
hallm.gsea.neg.subB <- read_tsv("data/subB-hallmarks-neg.tsv")
hallm.gsea.pos.subB <- hallm.gsea.pos.subB[1:5, c(1, 4, 6, 8)]
hallm.gsea.neg.subB <- hallm.gsea.neg.subB[1:20, c(1, 4, 6, 8)]


# positive NES
hallm.pos.subB <- ggplot(hallm.gsea.pos.subB, 
                          aes(x = NES, y = `FDR q-val`, 
                              color = `FDR q-val` < 0.05, 
                              label = sub("HALLMARK_", "", NAME))) + 
  geom_point(size = hallm.gsea.pos.subB$SIZE/100) + theme_light() + 
  scale_color_manual(values = c("black", "red")) + 
  geom_hline(yintercept = 0.05, col = "gray10", lty = 2) + ylab("") + 
  xlab("NES") + theme(legend.position = "none", 
                      axis.title = element_text(size = 17), 
                      axis.text = element_text(size = 17)) +
    geom_label_repel(position = position_dodge(width = 6), 
                     size = 4, max.overlaps = 100)

hallm.pos.subB


# negative NES
hallm.neg.subB <- ggplot(hallm.gsea.neg.subB, 
                          aes(x = NES, y = `FDR q-val`,
                              color = `FDR q-val` < 0.05, 
                              label = sub("HALLMARK_", "", NAME))) + 
  geom_point(size = hallm.gsea.neg.subB$SIZE/100) + theme_light() + 
  scale_color_manual(values = c("black", "darkblue")) + 
  geom_hline(yintercept = 0.05, col = "gray10", lty = 2) + ylab("FDR") + 
  xlab("NES") + ylim(c(-0.02, 0.08)) + 
  theme(legend.position = "none", axis.title = element_text(size = 17),
        axis.text = element_text(size = 17)) + 
  geom_label_repel(position = position_dodge(width = 6),
                   size = 4, max.overlaps = 100)
hallm.neg.subB

# plot
hallm.plot.subB <- ggarrange(hallm.neg.subB, hallm.pos.subB, 
                              ncol = 2, nrow = 1)
hallm.plot.subB


# plot results kegg
kegg.gsea.pos.subB <- read_tsv("data/subB-kegg-pos.tsv")
kegg.gsea.neg.subB <- read_tsv("data/subB-kegg-neg.tsv")
kegg.gsea.pos.subB <- kegg.gsea.pos.subB[1:5, c(1, 4, 6, 8)]
kegg.gsea.neg.subB <- kegg.gsea.neg.subB[1:5, c(1, 4, 6, 8)]


# positive NES
kegg.pos.subB <- ggplot(kegg.gsea.pos.subB, 
                         aes(x = NES, y = `FDR q-val`, 
                             color = `FDR q-val` < 0.05, label = NAME)) + 
  geom_point(size = kegg.gsea.pos.subB$SIZE/100) + theme_light() +
    scale_color_manual(values = c("black", "red")) + 
  geom_hline(yintercept = 0.05, col = "gray10", lty = 2) + ylab("") + 
  xlab("NES") + theme(legend.position = "none", 
                      axis.title = element_text(size = 17), 
                      axis.text = element_text(size = 17)) +
    geom_label_repel(position = position_dodge(width = 6), 
                     size = 4, max.overlaps = 100)

kegg.pos.subB


# negative NES
kegg.neg.subB <- ggplot(kegg.gsea.neg.subB, 
                         aes(x = NES, y = `FDR q-val`, 
                             color = `FDR q-val` < 0.05, label = NAME)) + 
  geom_point(size = kegg.gsea.neg.subB$SIZE/100) + theme_light() +
    scale_color_manual(values = c("black", "darkblue")) + 
  geom_hline(yintercept = 0.05, col = "gray10", lty = 2) + ylab("FDR") + 
  xlab("NES") + ylim(c(-0.02, 0.08)) + 
  theme(legend.position = "none", axis.title = element_text(size = 17),
    axis.text = element_text(size = 17)) + 
  geom_label_repel(position = position_dodge(width = 6),
    size = 4, max.overlaps = 100)

kegg.neg.subB

# plot
kegg.plot.subB <- ggarrange(kegg.neg.subB, kegg.pos.subB, ncol = 2, nrow = 1)
kegg.plot.subB


# -------------------------------- TFEA.ChIP -----------------------------------

# retrieve data
tfea.data.subB <- res.meta.resist.subB.df[, c("lfc", "pval", "adj.pval")]
list.genes.entrez <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id"), 
                           values = rownames(tfea.data.subB), mart = mart)
tfea.data.subB <- merge(list.genes.entrez, tfea.data.subB, 
                         by.x = "ensembl_gene_id", by.y = 0, all.y = TRUE)

# order data
tfea.data.subB <- tfea.data.subB[
                                order(tfea.data.subB$lfc, decreasing = TRUE), ]
tfea.data.subB <- tfea.data.subB[, -1]
colnames(tfea.data.subB) <- c("Genes", "log2FoldChange", "pvalue", "pval.adj")


# UP-regulated genes extract vector with names of upregulated genes
genes.upreg.subB <- Select_genes(tfea.data.subB, min_LFC = 1)

# extract vector with names of non-responsive genes
genes.control.subB <- Select_genes(tfea.data.subB, min_pval = 0.5, 
    max_pval = 1, min_LFC = -0.25, max_LFC = 0.25)

## association analysis

# generates list of contingency tables, one per dataset
CM_list_UP_subB <- contingency_matrix(genes.upreg.subB)  

# generates list of p-values and OR from association test
pval_mat_UP_subB <- getCMstats(CM_list_UP_subB)  
head(pval_mat_UP_subB)

TF_ranking_up_subB <- rankTFs(pval_mat_UP_subB, 
                                rankMethod = "gsea", makePlot = TRUE)
TF_ranking_up_subB[["TFranking_plot"]]

head(TF_ranking_up_subB[["TF_ranking"]])
tf.ranking.up.subB <- na.omit(TF_ranking_up_subB[["TF_ranking"]][
                             TF_ranking_up_subB[["TF_ranking"]]$pVal < 0.05, ])

plot_CM(pval_mat_UP_subB)  # plot p-values against ORs

tfs.up.subB <- rownames(tf.ranking.up.subB)
names(tfs.up.subB) <- rownames(tf.ranking.up.subB)
col.up.subB <- sample(viridis::turbo(length(tfs.up.subB)), 
                       length(tfs.up.subB))

# plot p-values against ORs highlighting indicated TFs
plot_CM(pval_mat_UP_subB, specialTF = tfs.up.subB, TF_colors = col.up.subB)  
TF_ranking_up_subB[["TFranking_plot"]]


# DOWN-regulated genes extract vector with names of upregulated genes
genes.downreg.subB <- Select_genes(tfea.data.subB, max_LFC = -1)

# extract vector with names of non-responsive genes
genes.control.subB <- Select_genes(tfea.data.subB, min_pval = 0.5, 
                                    max_pval = 1, min_LFC = -0.25, 
                                    max_LFC = 0.25)

## association analysis

# generates list of contingency tables, one per dataset
CM_list_DOWN_subB <- contingency_matrix(genes.downreg.subB)  

# generates list of p-values and OR from association test
pval_mat_DOWN_subB <- getCMstats(CM_list_DOWN_subB)  
head(pval_mat_DOWN_subB)

TF_ranking_down_subB <- rankTFs(pval_mat_DOWN_subB, 
                                 rankMethod = "gsea", makePlot = TRUE)
TF_ranking_down_subB[["TFranking_plot"]]
tf.ranking.down.subB <- na.omit(TF_ranking_down_subB[["TF_ranking"]][
                          TF_ranking_down_subB[["TF_ranking"]]$pVal < 0.05, ])

head(TF_ranking_down_subB[["TF_ranking"]])

plot_CM(pval_mat_DOWN_subB)  # plot p-values against ORs

tfs.down.subB <- rownames(tf.ranking.down.subB)
names(tfs.down.subB) <- rownames(tf.ranking.down.subB)
col.down.subB <- sample(viridis::turbo(length(tfs.down.subB)), 
                         length(tfs.down.subB))

# plot p-values against ORs highlighting indicated TFs
plot_CM(pval_mat_DOWN_subB, specialTF = tfs.down.subB, 
        TF_colors = col.down.subB)  
TF_ranking_down_subB[["TFranking_plot"]]


# ---------------------- actomyosin network analysis ---------------------------

## MA volcano plot

# add new column indicating if the gene is down, up or not regulated
res.meta.resist.subB.df$diffexp <- "NO"
res.meta.resist.subB.df$diffexp[res.meta.resist.subB.df$lfc > 1.5 & 
                              res.meta.resist.subB.df$adj.pval < 0.05] <- "UP"
res.meta.resist.subB.df$diffexp[res.meta.resist.subB.df$lfc < -1.5 & 
                            res.meta.resist.subB.df$adj.pval < 0.05] <- "DOWN"

# choose colors
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")

# plot
p.volcano.subB <- ggplot(data = res.meta.resist.subB.df, 
                          aes(x = lfc, y = -log10(adj.pval), fill = diffexp)) + 
  geom_point(shape = 21, size = 2.2) + theme_light() + 
  geom_vline(xintercept = c(-1.5, 1.5), col = "darkgray", lty = 2) + 
  geom_hline(yintercept = -log10(0.05), col = "darkgray", lty = 2) + 
  scale_fill_manual(values = c("blue", "black", "red")) + 
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 20), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14)) + 
  ylab("-log10(p-valor ajustado)") + xlab("LFC") +
  scale_x_continuous(breaks = seq(-10, 10, by = 2), 
                     limits = c(-10, 10)) + ylim(c(0, 60))
p.volcano.subB


## volcano plot selecting the actomyosin genes

# select actomyosin cytoskeleton genes
res.meta.resist.subB.df.sy <- merge(list.genes, res.meta.resist.subB.df, 
                               all.y = TRUE, by.x = "ensembl_gene_id", by.y = 0)

# add new column indicating if the gene is down, up or not regulated
res.meta.resist.subB.df.sy$actomy <- "NO"
res.meta.resist.subB.df.sy$actomy[
        (res.meta.resist.subB.df.sy$hgnc_symbol %in% actomyosin) & 
                                 (res.meta.resist.subB.df.sy$lfc > 0)] <- "UP"
res.meta.resist.subB.df.sy$actomy[
        (res.meta.resist.subB.df.sy$hgnc_symbol %in% actomyosin) & 
                        (res.meta.resist.subB.df.sy$lfc < 0)] <- "DOWN"

# choose colors
mycolors <- c("darkgray", "blue", "red")
names(mycolors) <- c("DOWN", "UP", "NO")
res.meta.resist.subB.df.sy$actomy <- factor(res.meta.resist.subB.df.sy$actomy,
    levels = c("NO", "DOWN", "UP"))

# order data
res.meta.resist.subB.df.sy <- res.meta.resist.subB.df.sy[
                                    order(res.meta.resist.subB.df.sy$actomy), ]
res.meta.resist.subB.df.sy$delabel <- NA
res.meta.resist.subB.df.sy$delabel[
  res.meta.resist.subB.df.sy$actomy != 
                 "NO"] <- res.meta.resist.subB.df.sy$hgnc_symbol[
                                    res.meta.resist.subB.df.sy$actomy != "NO"]
# plot
p.volcano.subB.2 <- ggplot(data = res.meta.resist.subB.df.sy, 
                            aes(x = lfc, y = -log10(adj.pval),
    col = actomy, label = delabel)) + 
  geom_point(shape = 21, size = 2.2, 
             fill = mycolors[res.meta.resist.subB.df.sy$actomy]) +
  theme_light() + 
  scale_colour_manual(values = c("darkgray", "black", "black")) +
  labs(x = "LFC", y = "log10(p-valor ajustado)") + 
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 20), 
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12)) + 
  scale_x_continuous(breaks = seq(-10, 10, by = 1)) + 
  geom_label_repel(data = subset(
       res.meta.resist.subB.df.sy, adj.pval < 0.05 & (lfc > 1.5 | lfc < -1.5)), 
       max.overlaps = 100) + xlim(c(-9, 9)) + ylim(c(0, 60)) + ggtitle("subB")

p.volcano.subB.2

## actomyosin network heatmap
actomyosin.ensemblid <- actomyosin
names(actomyosin.ensemblid) <- actomyosin
actomyosin.ensemblid <- merge(list.genes, actomyosin.ensemblid, 
                              by.x = "hgnc_symbol", by.y = 1)
lfc.genes.subB.sy.actomy <- merge(actomyosin.ensemblid, lfc.genes.subB,
                                   by.x = "ensembl_gene_id",
                                   by.y = 0, all.x = TRUE)

# remove one repeated line
lfc.genes.subB.sy.actomy <- lfc.genes.subB.sy.actomy[-70, ]

# change rownames
rownames(lfc.genes.subB.sy.actomy) <- paste(
                                     lfc.genes.subB.sy.actomy$ensembl_gene_id,
                            lfc.genes.subB.sy.actomy$hgnc_symbol, sep = " - ")
lfc.genes.subB.sy.actomy <- lfc.genes.subB.sy.actomy[, -c(1, 2)]

# plot
p.heatmap.acto.subB <- pheatmap(na.omit(lfc.genes.subB.sy.actomy[, -6]), 
    col = col.hm, cluster_rows = TRUE, cluster_cols = TRUE, 
    show_rownames = TRUE, fontsize_col = 10, breaks = seq(-4, 4, length = 50))

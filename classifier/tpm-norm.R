# load libraries
library("readr")
library("dplyr")

library("edgeR")


# set working directory
setwd("./tfm-bioinfo/classifier")

# load salmon output: RangedSummarizedExperiment and DGEList objects
load('salmon-output-classifier.RData')


# --------------------------- data normalization -------------------------------

# merge TPMs
data.clas <- merge(assay(gse.singleton, "abundance"), 
                   assay(gse.ho, "abundance"), by = 0) %>%
    merge(assay(gse.berico, "abundance"), by.x = "Row.names", by.y = 0) %>%
    merge(assay(gse.song, "abundance"), by.x = "Row.names", by.y = 0) %>%
    merge(assay(gse.song.clas, "abundance"), by.x = "Row.names", by.y = 0) %>%
    merge(assay(gse.tsoi.clas, "abundance"), by.x = "Row.names", by.y = 0) %>%
    merge(assay(gse.coppe.clas, "abundance"), by.x = "Row.names", by.y = 0) %>%
    merge(assay(gse.waizenegger.clas, "abundance"), 
          by.x = "Row.names", by.y = 0)

rownames(data.clas) <- as.vector(mapply(data.clas$Row.names, FUN = function(x) {
    gsub("\\.\\d+", "", x)
}))
data.clas <- data.clas[, -1]

# change colnames
colnames(data.clas) <- c(
  paste("singleton", y.singleton$samples$names, y.singleton$samples$cell, 
        y.singleton$samples$treatment, sep = "."), 
  paste("ho", y.ho$samples$names, y.ho$samples$cell, y.ho$samples$treatment, 
        sep = "."), 
  paste("berico", y.berico$samples$names, y.berico$samples$cell,
    y.berico$samples$treatment, sep = "."), 
  paste("song", y.song$samples$names, y.song$samples$cell,
    y.song$samples$treatment, sep = "."), 
  paste("song", y.song.clas$samples$names, y.song.clas$samples$cell, 
        y.song.clas$samples$treatment, sep = "."), 
  paste("tsoi", y.tsoi.clas$samples$names, y.tsoi.clas$samples$cell, 
        y.tsoi.clas$samples$treatment, sep = "."), 
  paste("coppe", y.coppe.clas$samples$names, y.coppe.clas$samples$cell, 
        y.coppe.clas$samples$treatment, sep = "."), 
  paste("waizenegger", y.waizenegger.clas$samples$names, 
        y.waizenegger.clas$samples$cell, 
        y.waizenegger.clas$samples$treatment, sep = "."))

# normalization
data.clas.norm <- data.clas
for (col in 1:dim(data.clas)[2]) {
    data.clas.norm[, col] <- 100 * as.numeric(factor(rank(data.clas[, col]))) /
                                 max(as.numeric(factor(rank(data.clas[, col]))))
}
summary(data.clas.norm)


## retrieve P and R cell lines to construct a classifier

# select DEGs in resistance
inter.genes <- read.table('data/inter-genes-resistance.txt')$V1

data.clas.norm.constr <- data.clas.norm[inter.genes, 
                    !grepl("Pt|pdx|DR|DTP|DTPP|BRAF", colnames(data.clas.norm))]

data.clas.norm.constr <- as.data.frame(t(data.clas.norm.constr))
pheno <- c(y.singleton$samples$treatment, y.ho$samples$treatment, 
           y.berico$samples$treatment, y.song$samples$treatment, 
           y.song.clas$samples$treatment, y.tsoi.clas$samples$treatment,
           y.coppe.clas$samples$treatment, y.waizenegger.clas$samples$treatment)
data.clas.norm.constr$phenotype <- pheno[!grepl("Pt|pdx|DR|DTP|DTPP|BRAF",
                                                colnames(data.clas.norm))]

# save table
write.table(data.clas.norm.constr, quote = FALSE, sep = "\t",
            file = "data/norm-data-tpms.txt")


## retrieve waizenegger data to validate the classifier
waizenegger.norm <- data.clas.norm[inter.genes, grepl("waizenegger", 
                                                      colnames(data.clas.norm))]
waizenegger.norm <- as.data.frame(t(waizenegger.norm))
waizenegger.norm$phenotype <- pheno[grepl("waizenegger", 
                                          colnames(data.clas.norm))]
dim(waizenegger.norm)

# save table
write.table(waizenegger.norm, quote = FALSE, sep = "\t", 
            file = "data/waizenegger.txt")


## retrieve coppe data to validate the classifier
coppe.norm <- data.clas.norm[inter.genes, grepl("coppe", 
                                                colnames(data.clas.norm))]
coppe.norm <- as.data.frame(t(coppe.norm))
coppe.norm$phenotype <- pheno[grepl("coppe", colnames(data.clas.norm))]
dim(coppe.norm)

# save table
write.table(coppe.norm, quote = FALSE, sep = "\t", 
            file = "data/coppe.txt")


## retrieve song patients data to validate the classifier
song.pt.norm <- data.clas.norm[inter.genes, 
                               (grepl("song", colnames(data.clas.norm)) &
                                        grepl("Pt", colnames(data.clas.norm)))]
song.pt.norm <- as.data.frame(t(song.pt.norm))
song.pt.norm$phenotype <- pheno[(grepl("song", colnames(data.clas.norm)) & 
                                   grepl("Pt", colnames(data.clas.norm)))]
song.pt.norm$phenotype <- gsub('DR', 'R', song.pt.norm$phenotype )
dim(song.pt.norm)

# save table 
write.table(song.pt.norm, quote = FALSE, sep = "\t", 
            file = "data/song-pt.txt")


# more data to validate the classifier

## tsoi (part 2: resistant and parental cell lines) ----

# load data
treat.tsoi.par <- relevel(factor(c("P", "P", "P", "P", "P", "P", "P", "R", "R", 
                                   "P", "P", "P", "P", "P", "P", "P", "P", "P", 
                                   "P", "P", "P", "P", "R", "P", "R", "P", "P", 
                                   "P", "P", "P", "P", "P", "P", "P", "P", "P", 
                                   "P", "R", "P", "R", "P", "P", "P", "P", "P", 
                                   "P", "P", "P", "P", "P", "P", "P", "P", "R", 
                                   "R", "P", "P", "P"),
                                 levels = c("P", "R")), ref = "P")

# tpms
tsoi.par <- assay(gse.tsoi.par, "abundance")
colnames(tsoi.par) <- c(paste("tsoi", y.tsoi.par$samples$names, 
                              y.tsoi.par$samples$cell,
                              y.tsoi.par$samples$treatment, sep = "."))
rownames(tsoi.par) <- as.vector(mapply(rownames(tsoi.par), 
                                       FUN = function(x) {gsub("\\.\\d+", "", x)}))

# normalization
tsoi.par.norm <- tsoi.par
for (col in 1:dim(tsoi.par)[2]) {
  tsoi.par.norm[, col] <- 100 * 
    as.numeric(factor(rank(tsoi.par[, col]))) /
    max(as.numeric(factor(rank(tsoi.par[, col]))))
}
summary(tsoi.par.norm)
dim(tsoi.par.norm)

tsoi.par.norm <- as.data.frame(t(tsoi.par.norm[inter.genes, ]))
tsoi.par.norm$phenotype <- treat.tsoi.par
dim(tsoi.par.norm)

# save table 
write.table(tsoi.par.norm, quote = FALSE, sep = "\t", 
            file = "data/tsoi-part2.txt")


## biorxiv dataset ----

# load data
treat.biorxiv.clas <- relevel(factor(c("R", "R", "R", "P", "R", "P", "R", "P", 
    "R", "P", "R", "P", "R", "P", "R"), levels = c("P", "R")), ref = "P")

# tpms
biorxiv.norm <- assay(gse.biorxiv.clas, "abundance")
colnames(biorxiv.norm) <- c(paste("biorxiv", y.biorxiv.clas$samples$names, 
                                  y.biorxiv.clas$samples$cell,
                                  y.biorxiv.clas$samples$treatment, sep = "."))

# normalization
biorxiv.norm.clas <- biorxiv.norm
for (col in 1:dim(biorxiv.norm)[2]) {
    biorxiv.norm.clas[, col] <- 100 * 
              as.numeric(factor(rank(biorxiv.norm[, col]))) /
                         max(as.numeric(factor(rank(biorxiv.norm[, col]))))
}
summary(biorxiv.norm.clas)
dim(biorxiv.norm.clas)

rownames(biorxiv.norm.clas) <- as.vector(mapply(rownames(biorxiv.norm.clas), 
                                   FUN = function(x) {gsub("\\.\\d+", "", x)}))
biorxiv.norm.clas <- as.data.frame(t(biorxiv.norm.clas[inter.genes, ]))
biorxiv.norm.clas$phenotype <- treat.biorxiv.clas
dim(biorxiv.norm.clas)

# save table
write.table(biorxiv.norm.clas, quote = FALSE, sep = "\t", 
            file = "data/biorxiv.txt")


## hugo (patients) ----

treat.hugo.pt.clas <- relevel(factor(c("P", "R", "P", "R", "R", "P", "R", "R", 
    "R", "P", "R", "P", "R", "R", "R", "P", "R", "R", "P", "R", "R", "R", "P", 
    "R", "R", "P", "R", "R", "R", "R", "R", "R", "R", "R", "P", "R", "R", "P", 
    "R", "P", "R", "R", "P", "DR", "P", "DR", "DR", "P", "DR", "R", "DR", "DR", 
    "R", "P", "DR", "DR", "DR", "P", "DR", "P", "DR", "DR", "DR"), 
    levels = c("P", "R", "DR")), ref = "P")

# tpms
hugo.pt.clas <- assay(gse.hugo.pt.clas, "abundance")
colnames(hugo.pt.clas) <- c(paste("hugo", y.hugo.pt.clas$samples$names, 
                                  y.hugo.pt.clas$samples$cell, 
                                  y.hugo.pt.clas$samples$treatment, sep = "."))
rownames(hugo.pt.clas) <- as.vector(mapply(rownames(hugo.pt.clas), 
                                    FUN = function(x) {gsub("\\.\\d+", "", x)}))

# patients data transformation
hugo.pt.clas.norm <- hugo.pt.clas
for (col in 1:dim(hugo.pt.clas)[2]) {
    hugo.pt.clas.norm[, col] <- 100 * 
               as.numeric(factor(rank(hugo.pt.clas[, col]))) /
                         max(as.numeric(factor(rank(hugo.pt.clas[, col]))))
}
summary(hugo.pt.clas.norm)
dim(hugo.pt.clas.norm)

hugo.pt.clas.norm <- as.data.frame(t(hugo.pt.clas.norm[inter.genes, ]))
hugo.pt.clas.norm$phenotype <- treat.hugo.pt.clas
dim(hugo.pt.clas.norm)

# save table
write.table(hugo.pt.clas.norm[-c(62, 63),], quote = FALSE, sep = "\t", 
            file = "data/hugo-pt.txt")

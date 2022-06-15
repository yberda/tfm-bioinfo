# load libraries
library(RankProd)
library(dplyr)
library(ggforce)
library(ggplot2)
library(ggVennDiagram)

# load data
setwd('/media/yosra/YOSRA/tfm/resistance/rankprod')
load('resistance-TE-MA.RData')

# merge logCPMs
logcpm <- merge(v.singleton$E, v.ho$E, by = 0, all = TRUE) %>%
  merge(v.berico$E, by.x = 'Row.names', by.y = 0, all = TRUE) %>%
  merge(v.song$E, by.x = 'Row.names', by.y = 0, all = TRUE)
rownames(logcpm) <- as.vector(mapply(logcpm$Row.names, 
                                     FUN = function(x){gsub("\\.\\d+", "", x)}))
logcpm <- logcpm [ , -1]
head(logcpm)

# change colnames
colnames(logcpm) <- c(paste('singleton', y.singleton$samples$names, 
                            y.singleton$samples$cell, 
                            y.singleton$samples$treatment, sep = '.'),
                      paste('ho', y.ho$samples$names, 
                            y.ho$samples$cell, 
                            y.ho$samples$treatment, sep = '.'),
                      paste('berico', y.berico$samples$names, 
                            y.berico$samples$cell, 
                            y.berico$samples$treatment, sep = '.'),
                      paste('song', y.song$samples$names, 
                            y.song$samples$cell, 
                            y.song$samples$treatment, sep = '.'))

# define the origin of each sample
origin <- c(paste('singleton', y.singleton$samples$cell, sep = '.'),
            paste('ho', y.ho$samples$cell, sep = '.'),
            paste('berico', y.berico$samples$cell, sep = '.'),
            paste('song', y.song$samples$cell, sep = '.'))


################################################################################
#                     RANKPROD & RANKSUM METAANALYSIS                          #
################################################################################

# RankProd analysis
rp.results <- RP.advance(data = na.omit(logcpm), 
                         cl = as.numeric(!grepl('.P', colnames(logcpm))), 
                         origin = origin,
                         logged = TRUE, 
                         na.rm = FALSE, 
                         gene.names = rownames(na.omit(logcpm)))

rp.results.tables <- topGene(rp.results, cutoff = 0.05, method = "pfp", 
                             logged = TRUE, logbase = 2, 
                             gene.names = rownames(na.omit(logcpm)))

# saveRDS(rp.results.tables, 'rp-results-tables.rds')


# RankSum analysis
rs.results <- RP.advance(data = na.omit(logcpm), 
                         cl = as.numeric(!grepl('.P', colnames(logcpm))), 
                         origin = origin,
                         logged = TRUE, 
                         na.rm = FALSE, 
                         calculateProduct = FALSE,
                         gene.names = rownames(na.omit(logcpm)))

rs.results.tables <- topGene(rs.results, cutoff = 0.05, method = "pfp", 
                             logged = TRUE, logbase = 2, 
                             gene.names = rownames(na.omit(logcpm)))

# saveRDS(rs.results.tables, 'rs-results-tables.rds')


# merge results in one table
rp.genes <- rbind(rp.results.tables$Table1, rp.results.tables$Table2)
rs.genes <- rbind(rs.results.tables$Table1, rs.results.tables$Table2)

# compare 3 MA
ggVennDiagram(list(meta = rownames(res.meta.resist.df.05.8stud), 
                  RankProd = rownames(rp.genes), 
                  RankSum = rownames(rs.genes)), label = 'count') +
  theme(legend.position = 'none')

# plot

## rankprod
meta.comp.rankprod <- merge(res.meta.resist.df.05.8stud[ , 'lfc', drop = FALSE],
                            rp.genes[ , 3, drop = FALSE], by = 0)
rownames(meta.comp.rankprod) <- meta.comp.rankprod$Row.names
meta.comp.rankprod <- meta.comp.rankprod[ , -1]

# correct LFC since it was computed by dividing P / T and without log
meta.comp.rankprod[ , 2] <- - log2(meta.comp.rankprod[ , 2])


## ranksum
meta.comp.ranksum <- merge(res.meta.resist.df.05.8stud[ , 'lfc', drop = FALSE],
                            rs.genes[ , 3, drop = FALSE], by = 0)
rownames(meta.comp.ranksum) <- meta.comp.ranksum$Row.names
meta.comp.ranksum <- meta.comp.ranksum[ , -1]

# correct LFC since it was computed by dividing P / T and without log
meta.comp.ranksum[ , 2] <- - log2(meta.comp.ranksum[ , 2])


# plots
par(mfrow = c(2, 2))
plot(meta.comp.rankprod[ , 2], meta.comp.rankprod[ , 1], 
     main = 'TE-MA vs.\nRankProd DEGs logFC',
     ylab = 'TE-MA DEGs (LFC)',
     xlab = 'RankProd metaanalysis DEGs (LFC)',
     pch = 21,
     col = "blue4", 
     bg = "blue1",
     ylim = c(0, 5),
     xlim = c(0, 5),
     cex = 1,      
     lwd = 1,
     panel.first = abline(coef = c(0,1), col = 'darkgrey', lwd = 3))
plot(meta.comp.ranksum[ , 2], meta.comp.ranksum[ , 1], 
     main = 'TE-MA vs.\nRankSum DEGs logFC',
     ylab = 'TE-MA DEGs (lfc)',
     xlab = 'RankSum metaanalysis DEGs (LFC)',
     pch = 21,
     col = "blue4", 
     bg = "blue1",
     ylim = c(0, 5),
     xlim = c(0, 5),
     cex = 1,      
     lwd = 1,
     panel.first = abline(coef = c(0,1), col = 'darkgrey', lwd = 3))
plot(meta.comp.rankprod[ , 2], meta.comp.rankprod[ , 1], 
     main = 'TE-MA vs.\nRankProd DEGs logFC',
     ylab = 'TE-MA DEGs (lfc)',
     xlab = 'RankProd metaanalysis DEGs (LFC)',
     pch = 21,
     col = "red4", 
     bg = "red1",
     ylim = c(-5, 0),
     xlim = c(-5, 0),
     cex = 1,      
     lwd = 1,
     panel.first = abline(coef = c(0,1), col = 'darkgrey', lwd = 3))
plot(meta.comp.ranksum[ , 2], meta.comp.ranksum[ , 1], 
     main = 'TE-MA vs.\nRankSum DEGs logFC',
     ylab = 'TE-MA DEGs (lfc)',
     xlab = 'RankSum metaanalysis DEGs (LFC)',
     pch = 21,
     col = "red4", 
     bg = "red1",
     ylim = c(-5, 0),
     xlim = c(-5, 0),
     cex = 1,      
     lwd = 1,
     panel.first = abline(coef = c(0,1), col = 'darkgrey', lwd = 3))
par(mfrow = c(1, 1))

# correlation: effect sizes MA (ES) vs. Rankprod or RankSum
print(paste('Correlation ES - RankProd MA', cor(meta.comp.rankprod[
                                                 meta.comp.rankprod$lfc > 0, 1], 
                                                meta.comp.rankprod[
                                                 meta.comp.rankprod$lfc > 0, 2], 
                                                method = 'pearson') )) 
print(paste('Correlation ES - RankProd MA', cor(meta.comp.rankprod[
                                                 meta.comp.rankprod$lfc < 0, 1], 
                                                meta.comp.rankprod[
                                                 meta.comp.rankprod$lfc < 0, 2], 
                                                method = 'pearson') )) 
print(paste('Correlation ES - RankSum MA', cor(meta.comp.ranksum[
                                                  meta.comp.ranksum$lfc > 0, 1], 
                                               meta.comp.ranksum[
                                                  meta.comp.ranksum$lfc > 0, 2], 
                                               method = 'pearson') )) 
print(paste('Correlation ES - RankSum MA', cor(meta.comp.ranksum[
                                                  meta.comp.ranksum$lfc < 0, 1], 
                                               meta.comp.ranksum[
                                                  meta.comp.ranksum$lfc < 0, 2], 
                                               method = 'pearson') )) 


################################################################################
#             RANKPROD & RANKSUM METAANALYSIS SEPARATING SUBGROUPS             #
################################################################################

# ----------- Subgroup A metaanalysis: Song SKMEL28R, Song M229R, --------------
# ----------------- Song M238R, Singleton A375R and Ho A375R -------------------

# RankProd analysis 
rp.results.subA <- RP.advance(data = na.omit(logcpm[ , c(1:6, 11:14, 23:38)]), 
                              cl = as.numeric(!grepl('.P', colnames(logcpm[ , 
                                                       c(1:6, 11:14, 23:38)]))), 
                              origin = origin[c(1:6, 11:14, 23:38)],
                              logged = TRUE, 
                              na.rm = FALSE, 
                              gene.names = rownames(na.omit(logcpm[ , 
                                                        c(1:6, 11:14, 23:38)])))

rp.results.subA.tables <- topGene(rp.results.subA, cutoff = 0.05, 
                             method = "pfp",  logged = TRUE, logbase = 2, 
                             gene.names = rownames(na.omit(logcpm[ , 
                                                        c(1:6, 11:14, 23:38)])))

# saveRDS(rp.results.subA.tables, 'rp-results-subA-tables.rds')


# RankSum analysis 
rs.results.subA <- RP.advance(data = na.omit(logcpm[ , c(1:6, 11:14, 23:38)]), 
                              cl = as.numeric(!grepl('.P', colnames(logcpm[ , 
                                                       c(1:6, 11:14, 23:38)]))), 
                              origin = origin[c(1:6, 11:14, 23:38)],
                              logged = TRUE, 
                              na.rm = FALSE, 
                              calculateProduct = FALSE,
                              gene.names = rownames(na.omit(logcpm[ , 
                                                        c(1:6, 11:14, 23:38)])))

rs.results.subA.tables <- topGene(rs.results.subA, cutoff = 0.05, 
                             method = "pfp", logged = TRUE, logbase = 2, 
                             gene.names = rownames(na.omit(logcpm[ , 
                                                        c(1:6, 11:14, 23:38)])))

# saveRDS(rs.results.subA.tables, 'rs-results-subA-tables.rds')


# merge results in one table
rp.results.subA.tables <- readRDS('rp-results-subA-tables.rds')
rs.results.subA.tables <- readRDS('rs-results-subA-tables.rds')

rp.subA.genes <- rbind(rp.results.subA.tables$Table1, 
                        rp.results.subA.tables$Table2)
rs.subA.genes <- rbind(rs.results.subA.tables$Table1, 
                        rs.results.subA.tables$Table2)

# compare 3 MA
ggVennDiagram(list(meta = rownames(res.meta.resist.subA.df.5st), 
                   RankProd = rownames(rp.subA.genes), 
                   RankSum = rownames(rs.subA.genes)), label = 'count') +
  theme(legend.position = 'none')


# plot

## rankprod
meta.comp.rankprod.subA <- merge(res.meta.resist.subA.df.5st[
             res.meta.resist.subA.df.5st$adj.pval < 0.05, 'lfc', drop = FALSE],
             rp.subA.genes[ , 3, drop = FALSE], by = 0)
rownames(meta.comp.rankprod.subA) <- meta.comp.rankprod.subA$Row.names
meta.comp.rankprod.subA <- meta.comp.rankprod.subA[ , -1]

# correct LFC since it was computed by dividing P / T and without log
meta.comp.rankprod.subA[ , 2] <- - log2(meta.comp.rankprod.subA[ , 2])


## ranksum
meta.comp.ranksum.subA <- merge(res.meta.resist.subA.df.5st[
             res.meta.resist.subA.df.5st$adj.pval < 0.05, 'lfc', drop = FALSE],
             rs.subA.genes[ , 3, drop = FALSE], by = 0)
rownames(meta.comp.ranksum.subA) <- meta.comp.ranksum.subA$Row.names
meta.comp.ranksum.subA <- meta.comp.ranksum.subA[ , -1]

# correct LFC since it was computed by dividing P / T and without log
meta.comp.ranksum.subA[ , 2] <- - log2(meta.comp.ranksum.subA[ , 2])


# plots
# first merge data
meta.compare.subA <- merge(res.meta.resist.subA.df[
               res.meta.resist.subA.df$adj.pval < 0.05 & 
                 res.meta.resist.subA.df$num.data == 5, 'lfc', drop = FALSE],
                    rp.subA.genes[ , 3, drop = FALSE], by = 0, all = TRUE) %>%
  merge(rs.subA.genes[ , 3, drop = FALSE], by.x = 'Row.names', 
        by.y = 0, all = TRUE)
rownames(meta.compare.subA) <- meta.compare.subA$Row.names
meta.compare.subA <- meta.compare.subA[ , -1]
colnames(meta.compare.subA) <- c('TE-MA', 'RankProd', 'RankSum')

# correct RankProd and RankSum LFC since it was computed by dividing 
# P / T and without log
meta.compare.subA[ , c(2, 3)] <- - log2(meta.compare.subA[ , c(2, 3)])

# separate UP and DOWN genes
meta.compare.subA.UP <- meta.compare.subA[meta.compare.subA$RankProd > 0, ]
meta.compare.subA.DOWN <- meta.compare.subA[meta.compare.subA$RankProd < 0, ]

# and plot
## UP
plot.MA.cor.subA.UP <- ggplot(meta.compare.subA.UP, 
                               aes(x = .panel_x, y = .panel_y)) + 
  geom_point(alpha = 0.2, shape = 16, size = 2, color = 'red') + 
  theme_light() +
  theme(strip.background = element_rect(fill = "gray33", 
                                        size = 1.5, linetype = "solid"),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  facet_matrix(vars(everything())) +
  coord_cartesian(ylim = c(0, 7), xlim = c(0, 7))

plot.MA.cor.subA.UP

## DOWN
plot.MA.cor.subA.DOWN <- ggplot(meta.compare.subA.DOWN, 
                                 aes(x = .panel_x, y = .panel_y)) + 
  geom_point(alpha = 0.2, shape = 16, size = 2, color = 'darkblue') + 
  theme_light() +
  theme(strip.background = element_rect(fill = "gray33", 
                                        size = 1.5, linetype = "solid"),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  facet_matrix(vars(everything())) +
  coord_cartesian(ylim = c(-7, 0), xlim = c(-7, 0))

plot.MA.cor.subA.DOWN


# correlation: effect sizes MA (ES) vs. RankProd or RankSum
print(paste('Correlation ES - RankProd MA', cor(meta.comp.rankprod.subA[
                                           meta.comp.rankprod.subA$lfc > 0, 1], 
                                                meta.comp.rankprod.subA[
                                           meta.comp.rankprod.subA$lfc > 0, 2], 
                                           method = 'pearson') ))
print(paste('Correlation ES - RankProd MA', cor(meta.comp.rankprod.subA[
                                           meta.comp.rankprod.subA$lfc < 0, 1], 
                                                meta.comp.rankprod.subA[
                                           meta.comp.rankprod.subA$lfc < 0, 2], 
                                           method = 'pearson') ))

print(paste('Correlation ES - RankSum MA', cor(meta.comp.ranksum.subA[
                                            meta.comp.ranksum.subA$lfc > 0, 1], 
                                               meta.comp.ranksum.subA[
                                            meta.comp.ranksum.subA$lfc > 0, 2], 
                                            method = 'pearson') )) 
print(paste('Correlation ES - RankSum MA', cor(meta.comp.ranksum.subA[
                                            meta.comp.ranksum.subA$lfc < 0, 1], 
                                               meta.comp.ranksum.subA[
                                            meta.comp.ranksum.subA$lfc < 0, 2], 
                                            method = 'pearson') ))


# ------------------Subgroup B metaanalysis: Song M395R, -----------------------
# --------------------- Ho 451LuR and Berico MM0074R ---------------------------

# RankProd analysis
rp.results.subB <- RP.advance(data = na.omit(logcpm[ , -c(1:6, 11:14, 23:38)]), 
                               cl = as.numeric(!grepl('.P', colnames(logcpm[ , 
                                                      -c(1:6, 11:14, 23:38)]))), 
                               origin = origin[-c(1:6, 11:14, 23:38)],
                               logged = TRUE, 
                               na.rm = FALSE, 
                               gene.names = rownames(na.omit(logcpm[ , 
                                                       -c(1:6, 11:14, 23:38)])))

rp.results.subB.tables <- topGene(rp.results.subB, cutoff = 0.05, 
                                   method = "pfp", logged = TRUE, logbase = 2, 
                                   gene.names = rownames(na.omit(logcpm[ ,
                                                       -c(1:6, 11:14, 23:38)])))

# saveRDS(rp.results.subB.tables, 'rp-results-subB-tables.rds')


# RankSum analysis
rs.results.subB <- RP.advance(data = na.omit(logcpm[ , -c(1:6, 11:14, 23:38)]), 
                               cl = as.numeric(!grepl('.P', colnames(logcpm[ , 
                                                      -c(1:6, 11:14, 23:38)]))), 
                               origin = origin[-c(1:6, 11:14, 23:38)],
                               logged = TRUE, 
                               na.rm = FALSE, 
                               calculateProduct = FALSE,
                               gene.names = rownames(na.omit(logcpm[ , 
                                                       -c(1:6, 11:14, 23:38)])))

rs.results.subB.tables <- topGene(rs.results.subB, cutoff = 0.05, 
                                   method = "pfp", logged = TRUE, logbase = 2, 
                                   gene.names = rownames(na.omit(logcpm[ , 
                                                       -c(1:6, 11:14, 23:38)])))

# saveRDS(rs.results.subB.tables, 'rs-results-subB-tables.rds')


# merge data in one table
rp.results.subB.tables <- readRDS('rp-results-subB-tables.rds')
rs.results.subB.tables <- readRDS('rs-results-subB-tables.rds')

rp.subB.genes <- rbind(rp.results.subB.tables$Table1, 
                        rp.results.subB.tables$Table2)
rs.subB.genes <- rbind(rs.results.subB.tables$Table1, 
                        rs.results.subB.tables$Table2)


# compare 3 MA
ggVennDiagram(list(meta = rownames(res.meta.resist.subB.df.3st), 
                   RankProd = rownames(rp.subB.genes), 
                   RankSum = rownames(rs.subB.genes)), label = 'count') +
  theme(legend.position = 'none')


# plot

## rankprod
meta.comp.rankprod.subB <- merge(res.meta.resist.subB.df.3st[
        res.meta.resist.subB.df.3st$adj.pval < 0.05, 'lfc', drop = FALSE],
                                  rp.subB.genes[ , 3, drop = FALSE], by = 0)
rownames(meta.comp.rankprod.subB) <- meta.comp.rankprod.subB$Row.names
meta.comp.rankprod.subB <- meta.comp.rankprod.subB[ , -1]

# correct LFC since it was computed by dividing P / T and without log
meta.comp.rankprod.subB[ , 2] <- - log2(meta.comp.rankprod.subB[ , 2])


# ranksum
meta.comp.ranksum.subB <- merge(res.meta.resist.subB.df.3st[
       res.meta.resist.subB.df.3st$adj.pval < 0.05, 'lfc', drop = FALSE],
                                 rs.subB.genes[ , 3, drop = FALSE], by = 0)
rownames(meta.comp.ranksum.subB) <- meta.comp.ranksum.subB$Row.names
meta.comp.ranksum.subB <- meta.comp.ranksum.subB[ , -1]

# correct LFC since it was computed by dividing P / T and without log
meta.comp.ranksum.subB[ , 2] <- - log2(meta.comp.ranksum.subB[ , 2])


# plots
# first merge data
meta.compare.subB <- merge(res.meta.resist.subB.df[
                           res.meta.resist.subB.df$adj.pval < 0.05 & 
                res.meta.resist.subB.df$num.data == 3, 'lfc', drop = FALSE],
                    rp.subB.genes[ , 3, drop = FALSE], by = 0, all = TRUE) %>%
  merge(rs.subB.genes[ , 3, drop = FALSE], by.x = 'Row.names', 
        by.y = 0, all = TRUE)
rownames(meta.compare.subB) <- meta.compare.subB$Row.names
meta.compare.subB <- meta.compare.subB[ , -1]
colnames(meta.compare.subB) <- c('TE-MA', 'RankProd', 'RankSum')

# correct RankProd and RankSum LFC since it was computed by ividing 
# P / T and without log
meta.compare.subB[ , c(2, 3)] <- - log2(meta.compare.subB[ , c(2, 3)])

# separate UP and DOWN
meta.compare.subB.UP <- meta.compare.subB[meta.compare.subB$`TE-MA` > 0, ,
                                            drop = FALSE]
meta.compare.subB.DOWN <- meta.compare.subB[meta.compare.subB$`TE-MA` < 0, ]

# and plot
## UP
plot.MA.cor.subB.UP <- ggplot(meta.compare.subB.UP, 
                               aes(x = .panel_x, y = .panel_y)) + 
  geom_point(alpha = 0.2, shape = 16, size = 2, color = 'red') + 
  theme_light() +
  theme(strip.background = element_rect(fill = "gray33", 
                                        size = 1.5, linetype = "solid"),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  facet_matrix(vars(everything())) +
  coord_cartesian(ylim = c(0, 6.5), xlim = c(0, 6.5))

plot.MA.cor.subB.UP

## DOWN
plot.MA.cor.subB.DOWN <- ggplot(meta.compare.subB.DOWN, 
                                 aes(x = .panel_x, y = .panel_y)) + 
  geom_point(alpha = 0.2, shape = 16, size = 2, color = 'darkblue') + 
  theme_light() +
  theme(strip.background = element_rect(fill="gray33", 
                                        size=1.5, linetype="solid"),
        strip.text.x = element_text(size=20),
        strip.text.y = element_text(size=20),
        axis.text = element_text(size = 15)) +
  facet_matrix(vars(everything())) +
  coord_cartesian(ylim = c(-6.5, 0), xlim = c(-6.5, 0))

plot.MA.cor.subB.DOWN

# correlation: effect sizes MA (ES) vs. RankProd or RankSum
print(paste('Correlation ES - RankProd MA', cor(meta.comp.rankprod.subB[
                                           meta.comp.rankprod.subB$lfc > 0, 1], 
                                                meta.comp.rankprod.subB[
                                           meta.comp.rankprod.subB$lfc > 0, 2], 
                                           method = 'pearson') ))
print(paste('Correlation ES - RankProd MA', cor(meta.comp.rankprod.subB[
                                           meta.comp.rankprod.subB$lfc < 0, 1], 
                                                meta.comp.rankprod.subB[
                                           meta.comp.rankprod.subB$lfc < 0, 2], 
                                           method = 'pearson') ))

print(paste('Correlation ES - RankSum MA', cor(meta.comp.ranksum.subB[
                                            meta.comp.ranksum.subB$lfc > 0, 1], 
                                               meta.comp.ranksum.subB[
                                            meta.comp.ranksum.subB$lfc > 0, 2], 
                                            method = 'pearson') )) 
print(paste('Correlation ES - RankSum MA', cor(meta.comp.ranksum.subB[
                                            meta.comp.ranksum.subB$lfc < 0, 1], 
                                               meta.comp.ranksum.subB[
                                            meta.comp.ranksum.subB$lfc < 0, 2], 
                                            method = 'pearson') ))

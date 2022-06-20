# load libraries
library('RankProd')
library('dplyr')

library('ggVennDiagram')
library('ggplot2')
library('ggforce')

# set working directory
setwd('tfm-bioinfo/treatment-MA')

# load transcriptomic data
load('rankprod-data-treatment.RData')

# this script contains the Rankprod and Ranksum metaanalyses and some 
# comparisons with TE-MA

################################################################################
#                              DATA PREPARATION                                #
################################################################################

# merge logCPMs
logcpm <- merge(v.obenauf.lines$E, v.fallahi.lines$E, by = 0, all = TRUE) %>%
  merge(v.corre$E, by.x = 'Row.names', by.y = 0, all = TRUE) %>%
  merge(v.gerosa$E, by.x = 'Row.names', by.y = 0, all = TRUE) %>%
  merge(v.smalley$E, by.x = 'Row.names', by.y = 0, all = TRUE) %>%
  merge(v.reganfendt$E, by.x = 'Row.names', by.y = 0, all = TRUE) %>%
  merge(v.song.lines$E, by.x = 'Row.names', by.y = 0, all = TRUE)
rownames(logcpm) <- as.vector(mapply(logcpm$Row.names, 
                                     FUN = function(x){gsub("\\.\\d+", "", x)}))
logcpm <- logcpm [ , -1]
head(logcpm)

# change colnames
colnames(logcpm) <- c(paste('obenauf', y.obenauf.lines$samples$names, 
                            y.obenauf.lines$samples$cell, 
                            y.obenauf.lines$samples$treatment, sep = '.'),
                      paste('fallahi', y.fallahi.lines$samples$names, 
                            y.fallahi.lines$samples$cell, 
                            y.fallahi.lines$samples$treatment, sep = '.'),
                      paste('corre', y.corre$samples$names, 
                            y.corre$samples$cell, 
                            y.corre$samples$treatment, sep = '.'),
                      paste('gerosa', y.gerosa$samples$names, 
                            y.gerosa$samples$cell, 
                            y.gerosa$samples$treatment, sep = '.'),
                      paste('smalley', y.smalley$samples$names, 
                            y.smalley$samples$cell, 
                            y.smalley$samples$treatment, sep = '.'),
                      paste('reganfendt', y.reganfendt$samples$names, 
                            y.reganfendt$samples$cell, 
                            y.reganfendt$samples$treatment, sep = '.'),
                      paste('song', y.song$samples$names, 
                            y.song$samples$cell, 
                            y.song$samples$treatment, sep = '.'))

# duplicate parental columns that are compared with more than one condition
logcpm <- cbind(logcpm, logcpm[ , 1:3], logcpm[ , 18:19], logcpm[ , 24:25])

# define the origin of each sample
origin <- c(paste('obenauf', y.obenauf.lines$samples$cell, sep = '.'),
            paste('fallahi', y.fallahi.lines$samples$cell, sep = '.'),
            paste('corre', y.corre$samples$cell, sep = '.'),
            paste('gerosa', y.gerosa$samples$cell, sep = '.'),
            paste('smalley', y.smalley$samples$cell, sep = '.'),
            paste('reganfendt', y.reganfendt$samples$cell, sep = '.'),
            paste('song', y.song$samples$cell, sep = '.'))

# add a "2" to repeated columns so the RankProd function does not consider
# different conditions the same
origin <- c(origin, paste(origin[1:3], '2', sep = '.'), 
            paste(origin[18:19], '2', sep = '.'),
            paste(origin[24:25], '2', sep = '.'))

origin[7:9] <- paste(origin[7:9], '2', sep = '.')
origin[22:23] <- paste(origin[22:23], '2', sep = '.')
origin[28:29] <- paste(origin[28:29], '2', sep = '.')


################################################################################
#              RANKPROD METAANALYSIS USING 6-48H TREATMENT DATA                #
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


## analysis of the results

# merge results in one table
rp.genes <- rbind(rp.results.tables$Table1, rp.results.tables$Table2)
rs.genes <- rbind(rs.results.tables$Table1, rs.results.tables$Table2)

# load TE-MA results
res.meta.df.05.14st <- readRDS('data/res-meta-df-05-14st.rds')

# compare 3 MA
ggVennDiagram(list(TE_MA = rownames(res.meta.df.05.14st), 
                   RankProd = rownames(rp.genes), 
                   RankSum = rownames(rs.genes)), label = 'count') +
  theme(legend.position = 'none')

# now plot

## rankprod
meta.comp.rankprod <- merge(res.meta.df.05.14st[ , 'lfc', drop = FALSE],
                            rp.genes[ , 3, drop = FALSE], by = 0)
rownames(meta.comp.rankprod) <- meta.comp.rankprod$Row.names
meta.comp.rankprod <- meta.comp.rankprod[ , -1]

# correct RankProd LFC since it was computed by dividing P / T and without log
meta.comp.rankprod[ , 2] <- - log2(meta.comp.rankprod[ , 2])


## ranksum
meta.comp.ranksum <- merge(res.meta.df.05.14st[ , 'lfc', drop = FALSE],
                           rs.genes[ , 3, drop = FALSE], by = 0)
rownames(meta.comp.ranksum) <- meta.comp.ranksum$Row.names
meta.comp.ranksum <- meta.comp.ranksum[ , -1]

# correct RankSum LFC since it was computed by dividing P / T and without log
meta.comp.ranksum[ , 2] <- - log2(meta.comp.ranksum[ , 2])

# plot
par(mfrow = c(2, 2))
plot(meta.comp.rankprod[meta.comp.rankprod > 0, 2], 
     meta.comp.rankprod[meta.comp.rankprod > 0, 1], 
     main='TE-MA metaanalysis vs.\nRankProd DEGs logFC',
     ylab = 'TE-MA metaanalysis DEGs (LFC)',
     xlab = 'RankProd metaanalysis DEGs (LFC)',
     pch = 21,
     col = "blue4", 
     bg = "blue1",
     ylim = c(0, 6),
     xlim = c(0, 6),
     cex = 1,      
     lwd = 1,
     panel.first = abline(coef = c(0,1), col = 'darkgrey', lwd = 3))
plot(meta.comp.ranksum[meta.comp.ranksum > 0, 2], 
     meta.comp.ranksum[meta.comp.ranksum > 0, 1], 
     main='TE-MA metaanalysis vs.\nRankSum DEGs logFC',
     ylab = 'TE-MA metaanalysis DEGs (LFC)',
     xlab = 'RankSum metaanalysis DEGs (LFC)',
     pch = 21,
     col = "blue4", 
     bg = "blue1",
     ylim = c(0, 6),
     xlim = c(0, 6),
     cex = 1,      
     lwd = 1,
     panel.first = abline(coef = c(0,1), col = 'darkgrey', lwd = 3))
plot(meta.comp.rankprod[meta.comp.rankprod$lfc < 0, 2], 
     meta.comp.rankprod[meta.comp.rankprod$lfc < 0, 1], 
     main='TE-MA metaanalysis vs.\nRankProd DEGs logFC',
     ylab = 'TE-MA metaanalysis DEGs (LFC)',
     xlab = 'RankProd metaanalysis DEGs (LFC)',
     pch = 21,
     col = "red4", 
     bg = "red1",
     ylim = c(-6, 0),
     xlim = c(-6, 0),
     cex = 1,      
     lwd = 1,
     panel.first = abline(coef = c(0,1), col = 'darkgrey', lwd = 3))
plot(meta.comp.ranksum[meta.comp.ranksum$lfc < 0, 2], 
     meta.comp.ranksum[meta.comp.ranksum$lfc < 0, 1], 
     main='TE-MA metaanalysis vs.\nRankSum DEGs logFC',
     ylab = 'TE-MA metaanalysis DEGs (LFC)',
     xlab = 'RankSum metaanalysis DEGs (LFC)',
     pch = 21,
     col = "red4", 
     bg = "red1",
     ylim = c(-6, 0),
     xlim = c(-6, 0),
     cex = 1,      
     lwd = 1,
     panel.first = abline(coef = c(0,1), col = 'darkgrey', lwd = 3))
par(mfrow = c(1, 1))

# correlation : effect sizes MA (TE-MA) vs. RankProd or RankSum
print(paste('Correlation TE-MA - RankProd MA', 
            cor(meta.comp.rankprod[meta.comp.rankprod$lfc > 0, 1], 
                meta.comp.rankprod[meta.comp.rankprod$lfc > 0, 2], 
                method = 'pearson') )) 

print(paste('Correlation TE-MA - RankProd MA', 
            cor(meta.comp.rankprod[meta.comp.rankprod$lfc < 0, 1], 
                meta.comp.rankprod[meta.comp.rankprod$lfc < 0, 2], 
                method = 'pearson') )) 

print(paste('Correlation TE-MA - RankSum MA', 
            cor(meta.comp.ranksum[meta.comp.ranksum$lfc > 0, 1],
                meta.comp.ranksum[meta.comp.ranksum$lfc > 0, 2], 
                method = 'pearson') )) 

print(paste('Correlation TE-MA - RankSum MA', 
            cor(meta.comp.ranksum[meta.comp.ranksum$lfc < 0, 1], 
                meta.comp.ranksum[meta.comp.ranksum$lfc < 0, 2], 
                method = 'pearson') )) 


################################################################################
#             RANKPROD METAANALYSIS USING 24-48H TREATMENT DATA                #
################################################################################

# RankProd analysis
rp.results.long <- RP.advance(data = na.omit(logcpm[ , -c(1:6, 42:53)]), 
                         cl = as.numeric(!grepl('.P', 
                                         colnames(logcpm[ , -c(1:6, 42:53)]))), 
                         origin = origin[-c(1:6, 42:53)],
                         logged = TRUE, 
                         na.rm = FALSE, 
                      gene.names = rownames(na.omit(logcpm[ , -c(1:6, 42:53)])))

rp.results.long.tables <- topGene(rp.results.long, cutoff= 0.05, method = "pfp", 
                             logged = TRUE, logbase = 2, 
                      gene.names = rownames(na.omit(logcpm[ , -c(1:6, 42:53)])))

# saveRDS(rp.results.long.tables, 'rp-results-long-tables.rds')


# RankSum analysis
rs.results.long <- RP.advance(data = na.omit(logcpm[ , -c(1:6, 42:53)]), 
                              cl = as.numeric(!grepl('.P', 
                                          colnames(logcpm[ , -c(1:6, 42:53)]))), 
                              origin = origin[-c(1:6, 42:53)],
                              logged = TRUE, 
                              na.rm = FALSE, 
                              calculateProduct = FALSE,
                      gene.names = rownames(na.omit(logcpm[ , -c(1:6, 42:53)])))

rs.results.long.tables <- topGene(rs.results.long, cutoff= 0.05, method = "pfp", 
                             logged = TRUE, logbase = 2, 
                      gene.names = rownames(na.omit(logcpm[ , -c(1:6, 42:53)])))

# saveRDS(rs.results.long.tables, 'rs-results-long-tables.rds')


# merge data in one table
rp.long.genes <- rbind(rp.results.long.tables$Table1, 
                       rp.results.long.tables$Table2)
rs.long.genes <- rbind(rs.results.long.tables$Table1, 
                       rs.results.long.tables$Table2)

# load TE-MA results
res.meta.long.df.05.11st <- readRDS('data/res-meta-long-df-05-11st.rds')

# compare
ggVennDiagram(list(meta = rownames(res.meta.long.df.05.11st), 
                   RankProd = rownames(rp.long.genes), 
                   RankSum = rownames(rs.long.genes)), label = 'count') +
  theme(legend.position = 'none')


# rankprod vs. effect sizes MA (TE-MA)
meta.comp.rankprod.long <- merge(res.meta.long.df.05.11st[
                 res.meta.long.df.05.11st$adj.pval < 0.05, 'lfc', drop = FALSE],
                                  rp.long.genes[ , 3, drop = FALSE], by = 0)
rownames(meta.comp.rankprod.long) <- meta.comp.rankprod.long$Row.names
meta.comp.rankprod.long <- meta.comp.rankprod.long[ , -1]

# correct LFC since it was computed by dividing P / T and without log
meta.comp.rankprod.long[ , 2] <- - log2(meta.comp.rankprod.long[ , 2])


# ranksum vs. effect sizes MA
meta.comp.ranksum.long <- merge(res.meta.long.df.05.11st[
                 res.meta.long.df.05.11st$adj.pval < 0.05, 'lfc', drop = FALSE],
                                 rs.long.genes[ , 3, drop = FALSE], by = 0)
rownames(meta.comp.ranksum.long) <- meta.comp.ranksum.long$Row.names
meta.comp.ranksum.long <- meta.comp.ranksum.long[ , -1]

# correct LFC since it was computed by dividing P / T and without log
meta.comp.ranksum.long[ , 2] <- - log2(meta.comp.ranksum.long[ , 2])


# rankprod vs. rankprod
rankprod.comp.ranksum.long <- merge(rp.long.genes[ , 3, drop = FALSE],
                                rs.long.genes[ , 3, drop = FALSE], by = 0)
rownames(rankprod.comp.ranksum.long) <- rankprod.comp.ranksum.long$Row.names
rankprod.comp.ranksum.long <- rankprod.comp.ranksum.long[ , -1]

# correct LFC since it was computed dividing P / T and without log
rankprod.comp.ranksum.long <- - log2(rankprod.comp.ranksum.long)


## plots

# first merge data
meta.compare <- merge(res.meta.long.df.05.11st[
                 res.meta.long.df.05.11st$adj.pval < 0.05, 'lfc', drop = FALSE],
                     rp.long.genes[ , 3, drop = FALSE], by = 0, all = TRUE) %>%
                merge(rs.long.genes[ , 3, drop = FALSE], by.x = 'Row.names', 
                                                       by.y = 0, all = TRUE)
rownames(meta.compare) <- meta.compare$Row.names
meta.compare <- meta.compare[ , -1]
colnames(meta.compare) <- c('TE-MA', 'RankProd', 'RankSum')

# correct RankProd and RankSum LFC since it was computed by dividing P / T 
# and without log
meta.compare[ , c(2, 3)] <- - log2(meta.compare[ , c(2, 3)])

# separate UP and DOWN
meta.compare.UP <- meta.compare[meta.compare$RankProd > 0, ]
meta.compare.DOWN <- meta.compare[meta.compare$RankProd < 0, ]

# and plot
## UP
plot.MA.cor.UP <- ggplot(meta.compare.UP, aes(x = .panel_x, y = .panel_y)) + 
  geom_point(alpha = 0.2, shape = 16, size = 2, color = 'red') + 
  theme_light() +
  theme(strip.background = element_rect(fill = "gray33", size = 1.5, 
                                        linetype = "solid"),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  facet_matrix(vars(everything())) +
  coord_cartesian(ylim = c(0, 6.5), xlim = c(0, 6.5))

plot.MA.cor.UP

## DOWN
plot.MA.cor.DOWN <- ggplot(meta.compare.DOWN, aes(x = .panel_x, y = .panel_y)) + 
  geom_point(alpha = 0.2, shape = 16, size = 2, color = 'darkblue') + 
  theme_light() +
  theme(strip.background = element_rect(fill = "gray33", size = 1.5, 
                                        linetype = "solid"),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  facet_matrix(vars(everything())) +
  coord_cartesian(ylim = c(-6.5, 0), xlim = c(-6.5, 0))

plot.MA.cor.DOWN


# correlation : effect sizes MA (ES) vs. RankProd or RankSum
print(paste('Correlation TE-MA - RankProd MA', cor(meta.comp.rankprod.long[
                                            meta.comp.rankprod.long$lfc > 0, 1], 
                                                meta.comp.rankprod.long[
                                            meta.comp.rankprod.long$lfc > 0, 2], 
                                            method = 'pearson')))

print(paste('Correlation TE-MA - RankProd MA', cor(meta.comp.rankprod.long[
                                            meta.comp.rankprod.long$lfc < 0, 1], 
                                                meta.comp.rankprod.long[
                                            meta.comp.rankprod.long$lfc < 0, 2], 
                                            method = 'pearson') ))

print(paste('Correlation TE-MA - RankSum MA', cor(meta.comp.ranksum.long[
                                             meta.comp.ranksum.long$lfc > 0, 1], 
                                               meta.comp.ranksum.long[
                                            meta.comp.ranksum.long$lfc > 0, 2], 
                                            method = 'pearson') )) 

print(paste('Correlation TE-MA - RankSum MA', cor(meta.comp.ranksum.long[
                                             meta.comp.ranksum.long$lfc < 0, 1], 
                                               meta.comp.ranksum.long[
                                             meta.comp.ranksum.long$lfc < 0, 2], 
                                             method = 'pearson') ))

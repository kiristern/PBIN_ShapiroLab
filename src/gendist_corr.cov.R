#https://github.com/benjjneb/dada2/issues/947
#https://github.com/adriaaulaICM/bbmo_niche_sea/blob/master/src/analysis/nucdist_otuclustering.R
library(tibble)
library(DECIPHER)
library(Biostrings)
library(reshape2)

#sparcc values
options(stringsAsFactors=FALSE)
.Options$sig_pval <- 0.05
# # Read in data as matrix
# d.corr_full <- read.table('median_correlation.tsv', sep='\t', header=TRUE, comment.char='', row.names=1, check.names=FALSE)
# d.pval_full <- read.table('pvalues.tsv', sep='\t', header=TRUE, comment.char='', row.names=1, check.names=FALSE)
# # Mask upper triangle with NAs then melt, excluding upper triangle values
# d.corr_full[upper.tri(d.corr_full, diag=TRUE)] <- NA
# d.pval_full[upper.tri(d.pval_full, diag=TRUE)] <- NA
# d.corr <- melt(as.matrix(d.corr_full), na.rm=TRUE, varnames=c('otu_1', 'otu_2'), value.name='correlation')
# d.pval <- melt(as.matrix(d.pval_full), na.rm=TRUE, varnames=c('otu_1', 'otu_2'), value.name='pvalue')
# # Merge correlations and pvalues, apply multiple testing correction
# sparc <- merge(d.corr, d.pval)
# sparc$pvalue_adjusted <- p.adjust(sparc$pvalue, method='BH')
# # Print rows with a significant pvalue, after multiple testing correction
# sparc.signif <- print(sparc[sparc$pvalue_adjusted <= .Options$sig_pval, ], row.names=FALSE)


sparcc.corr <- read.table("corr_sparcc_new.tsv", sep = "\t", header=T)
sparcc.signif <- sparcc.corr[sparcc.corr$pvalue_adjusted <= .Options$sig_pval, ]



# BiocManager::install("DECIPHER")
 library("DECIPHER")
seqs <- readDNAStringSet("ASVs.fa")
seqs.or <- OrientNucleotides(seqs)
aligned <- AlignSeqs(seqs.or)
d.align <- DECIPHER::DistanceMatrix(aligned)
d.align[1:5]

# gendist <- tidyr::gather(as.data.frame(d.align))
# as.data.frame(d.align)[1:5, 1:5]
# head(gendist)
# #write.csv(gendist, "gendist.csv")
# gendist$to <- paste0("ASV_", rownames(gendist))
# names(gendist)[names(gendist) == "key"] <- "from"
# dim(gendist)
vir.cov.tab <- read.csv("spiec_vir-vir.covar.csv", header = T, row.names = 1)
head(vir.cov.tab)
dim(vir.cov.tab)
names(vir.cov.tab)[names(vir.cov.tab) == "weight"] <- "covariance"
# dim(gendist.table)

#create df with genetic distance and covariance between unique viral ASV pairs
indx <- which(upper.tri(d.align, diag = F), arr.ind = TRUE)
nn <- dimnames(d.align)
gendist.tab <- data.frame(from = nn[[1]][indx[, 1]],
           to = nn[[2]][indx[, 2]],
           gen.dist = d.align[indx])
head(gendist.tab, n=20)
dim(gendist.tab)
df3 <- inner_join(gendist.tab, vir.cov.tab, by = c("to", "from"))
df4 <- inner_join(gendist.tab, vir.cov.tab, by = c("to" = "from", "from" = "to"))
gendist.covar <- rbind(df3,df4) %>% unique()
#write.csv(gendist.covar, "gendist.covar.csv")
head(gendist.covar <- read.csv("gendist.covar.csv", header = T, row.names = 1))

df5 <- inner_join(gendist.covar,sparcc.signif, by = c("to" = "otu_1", "from" = "otu_2"))
df6 <- inner_join(gendist.covar, sparcc.signif, by = c("from" = "otu_1", "to" = "otu_2"))
gendist.sparcc <- rbind(df5, df6) %>% unique()
head(gendist.sparcc)
#write.csv(gendist.sparcc, "gendist.sparcc.csv")
head(gendist.sparcc <- read.csv("gendist.sparcc.csv", header = T, row.names = 1))


dim(gendist.covar)
dim(gendist.sparcc)
head(gendist.covar)
head(gendist.sparcc)
#write.csv(gendist.covar, "gendist.covar.csv")


### plot boxplots ###

#### filter limits ####
gc <- gendist.covar %>%
  filter(gen.dist > 0 & gen.dist <0.05)

gc.pos <- gendist.covar %>%
  filter(gen.dist > 0 & gen.dist <0.05) %>%
  filter(covariance > 0)
dim(gc.filt)

gc.neg <- gendist.covar %>%
  filter(gen.dist > 0 & gen.dist <0.05) %>%
  filter(covariance < 0)

fit <- lm(covariance ~ gen.dist , data = gc.pos)
coefs <- coef(fit)
summary(fit)
#get r2
(r2 = signif(summary(fit)$adj.r.squared))
#get p-value
(pval <- signif(summary(fit)$coefficients[2,4]))
ggplot(gc.pos, aes(x=gen.dist, y=covariance, group=gen.dist))+
  geom_violin()+
  geom_jitter()+
  stat_summary(fun=mean, geom="point", shape=23, color="cyan4", fill="cyan4")+
  geom_smooth(method = "lm", se=F, color="red", aes(group=1))+
  labs(title = paste("Adj R2 = ", r2, 
                     #"Intercept =", signif(coefs[1]),
                     "Slope =", signif(coefs[2]),
                     "p-value =", pval))


gs <- gs.pos <- gendist.sparcc %>%
  filter(gen.dist <0.05)

gs.pos <- gendist.sparcc %>%
  filter(correlation > 0 & gen.dist <0.05)

gs.neg <- gendist.sparcc %>%
  filter(correlation < 0 & gen.dist <0.05)

fit <- lm(correlation ~ gen.dist , data = gs.pos)
(coefs <- coef(fit))
#get r2
(r2 = signif(summary(fit)$adj.r.squared))
#get p-value
(pval <- signif(summary(fit)$coefficients[2,4]))
ggplot(gs.pos, aes(x=gen.dist, y=correlation, group=gen.dist))+
  geom_point()+
  geom_violin()+
  geom_jitter()+
  stat_summary(fun=mean, geom="point", shape=23, size = 3,color="cyan4", fill="cyan4")+
  geom_smooth(method = "lm", se=F, color="red", aes(group=1))+
  labs(title = paste("Adj R2 = ", r2,
                     "Slope =", signif(coefs[2]),
                      "p-value =", pval))
  # coord_cartesian(xlim = c(0.0, 0.05)) #zoom in over these coords only without rm other data points 


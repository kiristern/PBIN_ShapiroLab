# No transformations: Alpha and Beta diversity, Mantel, Procruste
setwd("~/Documents/GitHub/PBIN_ShapiroLab/data")

#import scripts
source("../src/functions.R")
source("../src/preprocess.R")


## sequences per sample
sums_Phy <- data.frame(colSums(otu_table(virps3000filt)))
colnames(sums_Phy) <- "Sample_TotalSeqs"
sums_Phy$sample <- row.names(sums_Phy)
sums_Phy <- arrange(sums_Phy, Sample_TotalSeqs)
ggplot(sums_Phy, aes(x=reorder(sample, Sample_TotalSeqs), y = Sample_TotalSeqs)) +
  ylab("Number of Sequences per Sample") +
  geom_bar(stat = "identity", colour="black",fill="cornflowerblue")  + xlab("Sample Name") +
  ggtitle("Total Number of Sequences per Sample") +  #scale_y_continuous(breaks =seq(0, 1000000, 10000))+
  theme(axis.text.x = element_text(colour = "black", size=6, angle=45, hjust = 1, vjust = 1))


#check data
print(virps3000filt)

#https://microbiome.github.io/tutorials/
summarize_phyloseq(virps3000filt)
#sparsity is how populated is the data with zeros.

# separate into pelagic and littoral phyloseq objects
vir_ps_lit <- subset_samples(viral_physeq, Site == "Littoral")
vir_ps_pel <- subset_samples(viral_physeq, Site == "Pelagic")

vir_ps_lit_filt <- subset_samples(virps3000filt, Site == "Littoral")
vir_ps_pel_filt <- subset_samples(virps3000filt, Site == "Pelagic")



###### RELATIVE ABUNDANCE ######

### Top 20 ###
taxa_abun_tab_lit <- getTop20(vir_ps_lit_filt)
taxa_abun_tab_pel <- getTop20(vir_ps_pel_filt)

#get month-day column
md.l <- getMonthDay(taxa_abun_tab_lit)
md.p <- getMonthDay(taxa_abun_tab_pel)

#ensure colours are the same across both plots
dd <- union(taxa_abun_tab_lit$species, taxa_abun_tab_pel$species)

#convert taxa to factor with levels equal to dd variable (to incl all levels from both datasets)
taxa_abun_tab_lit$species <- factor(taxa_abun_tab_lit$species, dd)
taxa_abun_tab_pel$species <- factor(taxa_abun_tab_pel$species, dd)

#generate distinct colours for each asv
library("RColorBrewer")
set.seed(24)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
colpals <- qual_col_pals[c("Set1", "Dark2", "Set3"),]
col_vector = unlist(mapply(brewer.pal, colpals$maxcolors, rownames(colpals)))
dd.col=sample(col_vector, length(dd))
names(dd.col) <- dd
dd.col[names(dd.col) == "Other"] <- "lightgrey"

#plot relative abundance
rel_ab_plot_lit <- plotRelAb(taxa_abun_tab_lit, md.l, "Relative Abundance (littoral)")
rel_ab_plot_pel <- plotRelAb(taxa_abun_tab_pel, md.p, "Relative Abundance (pelagic)")

#combine both plots to one figure
library(ggpubr)
ggpubr::ggarrange(rel_ab_plot_pel, rel_ab_plot_lit, ncol=2, nrow=1, common.legend = T, legend="bottom")


#which asvs are the same between samples
uniq.lit <- unique(taxa_abun_tab_lit$species)
uniq.pel <- unique(taxa_abun_tab_pel$species)

uniq.lit[uniq.lit %in% uniq.pel]
uniq.lit[! uniq.lit %in% uniq.pel]
uniq.pel[! uniq.pel %in% uniq.lit]



#### ALPHA DIV ####
library(breakaway)

ba.dates <- meta %>% dplyr::select(Date)

vir_ps_lit 
vir_ps_pel 

#richness by year
ba.pel <- plot.ba(vir_ps_pel, "Breakaway richness of pelagic samples")
ba.lit <- plot.ba(vir_ps_lit, "Breakaway richness of littoral samples")

ggpubr::ggarrange(ba.pel, ba.lit, ncol=2, nrow=1, common.legend = T, legend="bottom")


# richness by years
ba_year_lit <- box.years(vir_ps_lit, "Observed richness by year (littoral)")
ba_year_lit
ba_year_pel <- box.years(vir_ps_pel, "Observed richness by year (pelagic)")
ba_year_pel
ba_year <- box.years(viral_physeq, "Observed richness by year")
ba_year


#richness by bloom/no-bloom and site
box.bloom <- plotbox(viral_physeq, "bloom2")
box.bloom[[2]]
box.bloom.df <- box.bloom[[1]]

box.site <- plotbox(viral_physeq, "Site")
box.site[[2]]
box.site.df <- box.site[[1]]


#t test
#http://www.sthda.com/english/wiki/unpaired-two-samples-t-test-in-r
#tests to check independent t-test assumptions
#Assumption 1: are the two samples independent? Yes. not taken from the same time 
#Assumption 2: does the data from each of the 2 groups follow a normal distirbution?
#Use Shapiro-Wilk normaility test
#Null hypothesis: the data are normally distributed
#Alternative hypothesis: the data are not normally distributed

# Shapiro-Wilk normality test for Bloom
with(box.site.df, shapiro.test(ba_observed_richness[env.var == "Littoral"])) # p = 0.03722 for lit; p-val = 0.2278 for bloom
# Shapiro-Wilk normality test for No Bloom
with(box.site.df, shapiro.test(ba_observed_richness[env.var == "Pelagic"])) # p = 0.6408 for pel; p=0.9858 for bloom
#the two p-values >0.05, thus fail to reject null; implying that the distribution of the data are not significantly different from the normal distribution. 
#Ie, we can assume the normality. 

#visualize for normality
library(ggpubr)
ggdensity(box.site.df$ba_observed_richness)
ggdensity(box.bloom.df$ba_observed_richness)

ggqqplot(box.site.df$ba_observed_richness)
ggqqplot(box.bloom.df$ba_observed_richness)

#Do the two vars have the same variances?
site.ftest <- var.test(ba_observed_richness ~ env.var, data=box.site.df) # p-val = 0.3056 > 0.05; no significant diff. b/w variances, therefore can use classic t-test
bloom.ftest <- var.test(ba_observed_richness ~ env.var, data=box.bloom.df)# p-val = 0.2611 

#t-test
site.ttest <- t.test(ba_observed_richness ~ env.var, data=box.site.df, var.equal=T)
site.ttest$p.value # 0.09212665 > 0.05; no signif differences b/w groups
bloom.ttest <- t.test(ba_observed_richness ~ env.var, data=box.bloom.df, var.equal=T)
bloom.ttest$p.value #0.4781509 > 0.05; no signif difference b/w sites




### sampling depth (reads per sample)
#https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
summary(sample_sums(viral_physeq)) #large difference in number of reads, min=22; max=130183
sort(sample_sums(viral_physeq))

ba <- breakaway(viral_physeq)

rich_depth <- data.frame("total_reads" =  phyloseq::sample_sums(viral_physeq),
                         "richness" = (ba %>% summary)$estimate)

#finding linear regression:
fit <- lm(richness ~ total_reads , data = rich_depth)
coefs <- coef(fit)
summary(fit)
#get r2
(r2 = signif(summary(fit)$adj.r.squared))
#get p-value
(pval <- signif(summary(fit)$coefficients[2,4]))

ggplot(data = rich_depth,
       aes(x = total_reads, y = richness)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x = "\nTotal Reads", y = "Richness\n",
       title = "Breakaway richness by sampling depth")+
    geom_abline(intercept = coefs[1], slope = coefs[2])+
    annotate(geom="text", x = 2e+04, y=700, size = 3,
             label= paste("Adj R2 = ", r2,
                           "p-val = ", pval))



ggplot(data = data.frame("total_reads" =  phyloseq::sample_sums(viral_physeq),
                         "Years" = viral_physeq %>% sample_data %>% get_variable("Years")),
       aes(x = total_reads, y = Years)) +
  geom_bar(stat = "identity") +
  labs(x = "\nTotal Reads", y = "Year\n")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), #rotate axis labels
        plot.title = element_text(hjust = 0.5))+ #center title
  ggtitle("Number of reads per Year")+
  coord_flip()






#### Shannon diversity ####
shan.lit <- plot.shannon(vir_ps_lit, "Shannon diversity of littoral samples")
shan.lit.plot <- shan.lit[[2]]
shan.lit.df <- shan.lit[[1]]
shan.pel <- plot.shannon(vir_ps_pel, "Shannon diversity of pelagic samples")
shan.pel.plot <- shan.pel[[2]]
shan.pel.df <- shan.pel[[1]]

ggpubr::ggarrange(shan.lit.plot, shan.pel.plot, ncol=2, nrow=1, common.legend = T, legend="bottom")

vir_shannon <- plot.shannon(viral_physeq, "Shannon diversity by sample")
vir_shannon.plot <- vir_shannon[[2]]
vir_shan.df <- vir_shannon[[1]]

#shannon by years with regression
shan.years(vir_shan.df)
shan.years(shan.pel.df)
shan.years(shan.lit.df)


# shannon by bloom/no bloom and site
box.shan.bloom <- box.shannon(viral_physeq, "bloom2")
box.shan.bloom[[2]]
shan.bloom <- box.shan.bloom[[1]]

box.shan.site <- box.shannon(viral_physeq, "Site")
box.shan.site[[2]]
shan.site <- box.shan.site[[1]]


#t test
# Shapiro-Wilk normality test 
with(shan.bloom, shapiro.test(Shannon[env.var == "yes"])) # p = 0.237 for lit; p-val = 0.087 for bloom
# Shapiro-Wilk normality test for No Bloom
with(shan.bloom, shapiro.test(Shannon[env.var == "no"])) # p = 0.0015 for pel; p=0.028 for bloom
#the two p-values >0.05, thus fail to reject null; implying that the distribution of the data are not significantly different from the normal distribution. 
#Ie, we can assume the normality. 
#if the data are not normally distributed, it’s recommended to use the non parametric two-samples Wilcoxon rank test.

#visualize for normality
library(ggpubr)
ggdensity(shan.site$Shannon)
ggdensity(shan.bloom$Shannon)

ggqqplot(shan.site$Shannon)
ggqqplot(shan.bloom$Shannon)

#Do the two vars have the same variances?
site.ftest <- var.test(Shannon ~ env.var, data=shan.site) # p-val = 0.7187 > 0.05; no significant diff. b/w variances, therefore can use classic t-test
bloom.ftest <- var.test(Shannon ~ env.var, data=shan.bloom)# p-val = 0.4085 

# #Wilcoxon test
# # Question: Is there any significant changes in the richness of ASV during and not during bloom?
# wilco.site <- wilcox.test(Shannon ~ env.var, data = shan.site)
# wilco.bloom <- wilcox.test(Shannon ~ env.var, data = shan.bloom)
# wilco.site$p.value# 0.3164 > 0.05 therefore no significant difference between sites
# wilco.bloom$p.value # 0.9461972 no signif difference between bloom

#t-test
site.ttest <- t.test(Shannon ~ env.var, data=shan.site, var.equal=T)
site.ttest$p.value #  0.2216704 > 0.05; no signif differences b/w sites
bloom.ttest <- t.test(Shannon ~ env.var, data=shan.bloom, var.equal=T)
bloom.ttest$p.value #0.7008102 > 0.05; no signif difference b/w groups




#### shannon vs. depth ####
shan <- estimate_richness(viral_physeq, measures="shannon")
head(shan)

viral_df = data.frame("shannon" = shan$Shannon,
                      "sample" = (viral_physeq %>% sample_data)$description,
                      "Years" = viral_physeq %>% sample_data %>% get_variable("Years"),
                      "depth" = sample_sums(viral_physeq))
head(viral_df)
str(viral_df)
ggplot(viral_df, aes(x = forcats::fct_inorder(sample), y = shannon, color = Years))+ #fct_inorder ensures plotting in order of sample date
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5), #rotate axis labels
        plot.title = element_text(hjust = 0.5))+ #center title
  ggtitle("Shannon diversity by sample")+
  scale_x_discrete(labels = viral_physeq %>% sample_data %>% get_variable("Months"), name="Month")#change x-axis sample name to Month

### sampling depth (reads per sample)
summary(sample_sums(viral_physeq)) #large difference in number of reads, min=22; max=130183
sort(sample_sums(viral_physeq))

fit <- lm(shannon ~ depth , data = viral_df)
coefs <- coef(fit)
summary(fit)
#get r2
(r2 = signif(summary(fit)$adj.r.squared))
#get p-value
(pval <- signif(summary(fit)$coefficients[2,4]))

ggplot(data = viral_df, aes(x = depth, y = shannon)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x = "\nTotal Reads", y = "Shannon\n")+
  geom_abline(intercept = fit$coefficients[1], slope =  fit$coefficients[2])+
  annotate(geom="text", x = 2e+04, y=4.75, label= paste("Adj R2 = ", r2,
                                                        "p-val = ", pval))+
  ggtitle("Shannon diversity by depth")


ggplot(data = data.frame("total_reads" =  phyloseq::sample_sums(viral_physeq),
                         "Years" = viral_physeq %>% sample_data %>% get_variable("Years")),
       aes(x = total_reads, y = Years)) +
  geom_bar(stat = "identity") +
  labs(x = "\nTotal Reads", y = "Year\n")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), #rotate axis labels
        plot.title = element_text(hjust = 0.5))+ #center title
  ggtitle("Number of reads per Year")+
  coord_flip()




#### PERMANOVA ####
#http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova
library(vegan)
citation("vegan")
virps3000filt %>% sample_data() %>% head
jsd <- sqrt(phyloseq::distance(virps3000filt, method = "jsd")) #jsd is more robust to sampling depth
sampledf <- data.frame(sample_data(virps3000filt)) #make a df from the sample_data
adonis(jsd ~ Period + Years + Months + Site, data = sampledf)
adonis(jsd ~ Site, data = sampledf)
adonis(jsd ~ Months, data = sampledf)

#homogeneity of dispersion test
betadisp <- betadisper(jsd, sampledf$Period)
permutest(betadisp)



### NMDS ###
nmds=metaMDS(comm = sqrt(phyloseq::distance(virps3000filt, "jsd")), k=3, trymax = 100)
# nmds <- ordinate(physeq = virps3000,
#                  method = "NMDS",
#                  distance = "jsd",
#                  trymax=500)
plot_ordination(physeq = virps3000filt,
                ordination = nmds,
                color = "Years",
                #shape = "Site",
                title = "NMDS of Lake Champlain viral Communities")+
  geom_point(aes(color=Years))+
  scale_color_brewer(palette = "Paired")
 # geom_point(colour="grey90", size=1.5)
nmds$stress #0.1432 good -- convergence. high stress value means that the algorithm had a hard time representing the distances between samples in 2 dimensions (anything <0.2 is considered acceptable)













#https://bioconductor.riken.jp/packages/3.8/bioc/vignettes/microbiome/inst/doc/vignette.html
#Visually-Weighted Regression curve with smoothed error bars
# Estimate Shannon diversity and add it to the phyloseq object
sample_data(filt_virseq)$ShannonDiv <- 
  metadata$ShannonDiv <- filt_virseq %>% otu_table %>% microbiome::alpha() %>% select("diversity_shannon")

#compare year and microbiome shannon diversity
microbiome::plot_regression(ShannonDiv ~ Years, metadata) #doesn't work!


#visualize the core microbiota (set of taxa that are detected in a remarkable fraction of the population above a give abundance threshold)
library(ggplot2, quiet = TRUE)
p <- plot_core(transform(filt_virseq, "compositional"), 
               plot.type = "heatmap", 
               colours = gray(seq(0,1,length=5)),
               prevalences = seq(.05, 1, .05), 
               detections = 10^seq(log10(1e-3), log10(.2), length = 10), 
               horizontal = TRUE) +
  xlab("Detection Threshold (Relative Abundance (%))") 
print(p)    


# #compositional heatmap
# tmp <- plot_composition(filt_virseq, plot.type="heatmap", transform = "hellinger")
# 
# #compositional barplot
# plot_composition(transform(filt_virseq, "compositional"), 
#                  plot.type = "barplot", sample.sort = "neatmap", label=F)











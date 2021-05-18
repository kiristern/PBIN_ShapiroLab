# No transformations: Alpha and Beta diversity, Mantel, Procruste

packageVersion("Phyloseq")
citation("phyloseq")

### Functions ###
##function to set vir_hel to same format as bact_hel
colsamp2date <- function(tab2format){
  colnames(tab2format) <- sub("*._*._*._*._*._*._*._","", colnames(tab2format))
  colnames(tab2format) <- gsub("_", ".", colnames(tab2format))
  return(tab2format)
}


##################################

library(dplyr)
library(stringr)
library(vegan)

setwd("~/Documents/GitHub/PBIN_ShapiroLab/data")

#### UPLOAD DATA ####
#upload viral ASV count table and metadata
ASV_count <- read.table("ASVs_counts_copy.tsv", row.names = 1, header=T)
str(ASV_count)
dim(ASV_count)
range(ASV_count)
apply(ASV_count, 2, median) #get median for each col
median(unlist(ASV_count), na.rm = T) #et median for whole df
summary(ASV_count)

colnames(ASV_count)[colnames(ASV_count) == "FLD0295_15_05_2011_1"] <- "FLD0295_15_05_2011_2" #dates were duplicated therefore need to correct
head(ASV_count, n=2)
meta <- read.csv("meta_cmd.csv", row.names = 1, header = T)
head(meta, n=2)
#meta$Years <- as.factor(meta$Years)
meta$Years <- as.factor(meta$Years)
meta$Date <- as.Date(meta$Date)
str(meta)

#order meta by date
meta <- meta[order(meta$Date),]

#get basic meta data info for methods section of report
length(ASV_count)
nrow(meta)
(amt <- meta %>% group_by(Years) %>% summarise(amount=length(Years))) #view how many samples per year
min(amt[,2])
max(amt[,2])
median(as.numeric(unlist(amt[,2]))) #get median samples per year
meta %>% group_by(Site) %>% summarise(amount=length(Site))


#ensure same samples between ASV_count and meta
asv_count<- ASV_count[,(colnames(ASV_count) %in% rownames(meta))]
length(asv_count)
nrow(meta)



#### CREATE PHYLOSEQ OBJECT ####
library(phyloseq)
library(microbiome)
#VIRAL

#add ASV count table, metadata, virTree to phyloseq table
nozero <- asv_count[rownames(asv_count) %in% names(rowSums(asv_count > 0)),]
length(rowSums(asv_count > 0)) == nrow(asv_count)
count_phy <- otu_table(asv_count, taxa_are_rows=T)
sample_info <- sample_data(meta)
virTree <- read_tree("viral_tree")

fake_taxa <- read.table("fake_viral_tax.txt", header = T, row.names = 1, fill=T)
mock_taxa <- tax_table(fake_taxa)
mock_taxa[,7] <- str_remove(mock_taxa[,7], "s__")
row.names(mock_taxa) <- mock_taxa[,7]
colnames(mock_taxa)[7] <- "species"
head(mock_taxa)
#add to phyloseq object
viral_physeq <- phyloseq(count_phy, sample_info, virTree, mock_taxa)
viral_physeq %>% otu_table( ) %>% dim

# rm reads less than 3000
virps3000 = prune_samples(sample_sums(viral_physeq)>=3000, viral_physeq)

virps3000filt <- filter_taxa(virps3000, function(x) sum(x > 1) > (0.10*length(x)), TRUE)

virfiltotu <- virps3000filt %>% otu_table()
write.csv(virfiltotu, "viralCounts_filtered.csv")


#which rownames are different
virfiltotu

viral_relab_ps <- transform(viral_physeq, "compositional", target="OTU")
viral_relab_otu <- viral_relab_ps %>% otu_table()

viralfilt_relab_otu <- viral_relab_otu[rownames(virfiltotu),]
dim(viralfilt_relab_otu)
dim(virfiltotu)
write.csv(viralfilt_relab_otu, "viralfilt_relab.csv")



which(taxa_sums(virps3000filt) == 0)

#asv_count
asv_count_3000 <- virps3000filt %>% otu_table()
#write.table(asv_count_3000, "asv_count_3000.tsv")

vir_abun <- virps3000 %>% otu_table()

sums_Phy <- data.frame(colSums(otu_table(virps3000)))
colnames(sums_Phy) <- "Sample_TotalSeqs"
sums_Phy$sample <- row.names(sums_Phy)
sums_Phy <- arrange(sums_Phy, Sample_TotalSeqs)
ggplot(sums_Phy, aes(x=reorder(sample, Sample_TotalSeqs), y = Sample_TotalSeqs)) +
  ylab("Number of Sequences per Sample") +
  geom_bar(stat = "identity", colour="black",fill="cornflowerblue")  + xlab("Sample Name") +
  ggtitle("Total Number of Sequences per Sample") +  #scale_y_continuous(breaks =seq(0, 1000000, 10000))+
  theme(axis.text.x = element_text(colour = "black", size=6, angle=45, hjust = 1, vjust = 1))


#check data
print(virps3000)

#https://microbiome.github.io/tutorials/
summarize_phyloseq(virps3000)
#sparsity is how populated is the data with zeros.

# separate into pelagic and littoral phyloseq objects
vir_ps_lit <- subset_samples(viral_physeq, Site == "Littoral")
vir_ps_pel <- subset_samples(viral_physeq, Site == "Pelagic")

vir_ps_lit_filt <- subset_samples(virps3000, Site == "Littoral")
vir_ps_pel_filt <- subset_samples(virps3000, Site == "Pelagic")



###### RELATIVE ABUNDANCE ######

### Top 20 ###
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("gmteunisse/Fantaxtic")
library("fantaxtic")
getTop20 <- function(sampType){
  ps_relab <- transform_sample_counts(sampType, function(OTU) OTU/sum(OTU))
  
  relab_top20 <- get_top_taxa(ps_relab, 20, relative = TRUE, discard_other = F,
                              other_label = "Other")
  taxa_abun_tab <- psmelt(relab_top20)
  return(taxa_abun_tab)
}

taxa_abun_tab_lit <- getTop20(vir_ps_lit_filt)
taxa_abun_tab_pel <- getTop20(vir_ps_pel_filt)

#create date col
library(lubridate)
getMonthDay <- function(table){
  m <- month(table$Date)
  day <- day(table$Date)
  md <- paste(day, m, sep="-")
  return(md)
}
md.l <- getMonthDay(taxa_abun_tab_lit)
md.p <- getMonthDay(taxa_abun_tab_pel)

#ensure colours are the same across both plots (need to run pelagic script below)
dd <- union(taxa_abun_tab_lit$species, taxa_abun_tab_pel$species)

#generate distinct colours for each asv
library("RColorBrewer")
set.seed(24)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
colpals <- qual_col_pals[c("Set1", "Dark2", "Set3"),]
col_vector = unlist(mapply(brewer.pal, colpals$maxcolors, rownames(colpals)))
dd.col=sample(col_vector, length(dd))
names(dd.col) <- dd
dd.col[names(dd.col) == "Other"] <- "lightgrey"

plotRelAb <- function(rel_ab_tab, md, plotTitle){
  rel_ab_plot <- rel_ab_tab %>% 
    ggplot(aes(x =Sample, y = Abundance, fill = species, order = -species)) +
    geom_bar(stat = "identity",position = position_stack(reverse = T)) + #position: Other should be top stack
    scale_fill_manual("ASV", values = dd.col)+
    labs(x = "",
         y = "Relative Abundance",
         title = plotTitle) +
    facet_grid(~ Years, scales = "free") +
    theme(
      axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 12),
      legend.text = element_text(size = 10),
      strip.text = element_text(size = 12)
    )+
    scale_x_discrete(labels = md, name="Sample date")+
    guides(fill=guide_legend(reverse = T)) #match the legend order to the order of the stack bars
return(rel_ab_plot)
}

rel_ab_plot_lit <- plotRelAb(taxa_abun_tab_lit, md.l, "Relative Abundance (littoral)")
rel_ab_plot_pel <- plotRelAb(taxa_abun_tab_pel, md.p, "Relative Abundance (pelagic)")

library(ggpubr)
ggpubr::ggarrange(rel_ab_plot_pel, rel_ab_plot_lit, ncol=2, nrow=1, common.legend = F, legend="bottom")





#### ALPHA DIV ####
library(breakaway)

ba.dates <- meta %>% dplyr::select(Date)

vir_ps_lit 
vir_ps_pel 

#richness by year
library(lubridate)

plot.ba <- function(ps.obj, plot.title){
  ba <- breakaway(ps.obj)
  ymd <- ps.obj %>% sample_data %>% get_variable("Date")
  m <- month(ymd)
  d <- day(ymd)
  md <- paste( d, m, sep="-")
  
  ba_vir_df = data.frame("richness" = (ba %>% summary)$estimate,
                          #"sample" = (ba %>% summary)$sample_names,
                          "error" = (ba %>% summary)$error,
                          "Years" = ps.obj %>% sample_data %>% get_variable("Years"),
                          "Upper" = (ba %>% summary)$upper,
                          "Lower" = (ba %>% summary)$lower,
                         "sample"= ps.obj %>% sample_data %>% get_variable("description"))
  head(ba_vir_df)
  baPlot <- ggplot(ba_vir_df, aes(x = forcats::fct_inorder(sample), y = richness, color = Years))+ #fct_inorder ensures plotting in order of sample date
    geom_point(size=3)+
    geom_errorbar(aes(ymin=richness-abs(richness-Lower), ymax=richness+abs(richness-Upper), width=0.05))+ 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5), #rotate axis labels
          plot.title = element_text(hjust = 0.5))+ #center title
    ggtitle(plot.title)+
    scale_x_discrete(labels = md, name="Sample date")+ #change x-axis sample name to Month-Day
    scale_y_continuous(name="Richness")+
    ylim(0, 650)
  return(baPlot)
}

ba.pel <- plot.ba(vir_ps_pel, "Breakaway richness of pelagic samples")
ba.lit <- plot.ba(vir_ps_lit, "Breakaway richness of littoral samples")

ggpubr::ggarrange(ba.pel, ba.lit, ncol=2, nrow=1, common.legend = T, legend="bottom")


#boxplot years
ba_year = data.frame("ba_observed_richness" = (ba %>% summary)$estimate,
                     "Years" = virps3000 %>% sample_data %>% get_variable("Years"))
(ba_plot <-  ggplot(ba_year, aes(x = Years, y = ba_observed_richness))+
    geom_point()) + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                                 geom="crossbar", width=0.5) + theme_minimal()+
  ggtitle("Observed richness by year")+
  theme(plot.title = element_text(hjust=0.5))+
  scale_y_continuous(name = "Observed richness")


#boxplot years
means <- aggregate(ba_observed_richness ~ Years, ba_year, mean)
means$step <- 1:nrow(means)

#finding linear regression:
#http://r-statistics.co/Linear-Regression.html 
fit <- lm(ba_observed_richness ~ step , data = means)
coefs <- coef(fit)
summary(fit)
#get r2
(r2 = signif(summary(fit)$adj.r.squared))
#get p-value
(pval <- signif(summary(fit)$coefficients[2,4]))

(ba_plot <-  ggplot(ba_year, aes(x = Years, y = ba_observed_richness))+
    geom_point()+
    geom_abline(intercept = coefs[1], slope = coefs[2])+
    annotate(geom="text", x = 3.5, y=620, label= paste("Adj R2 = ", r2,
                                                     "p-val = ", pval))+
    labs(title = "Observed richness by year")+
    stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                       geom="crossbar", width=0.5)+ 
    theme_minimal())




#richness by bloom/no-bloom

ba_virps3000 <- breakaway(virps3000)

#boxplot bloom
ba_bloom = data.frame("ba_observed_richness" = (ba_virps3000 %>% summary)$estimate,
                      "Bloom" = virps3000 %>% sample_data %>% get_variable("bloom2"))
ba_bloom_yn <- na.omit(ba_bloom[1:2])

(ba_plot <-  ggplot(ba_bloom_yn, aes(x = Bloom, y = ba_observed_richness))+
    geom_point()) + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                                 geom="crossbar", width=0.5) + theme_minimal()

#t test
#http://www.sthda.com/english/wiki/unpaired-two-samples-t-test-in-r
group_by(ba_bloom_yn, Bloom) %>%
  summarise(
    # count = n(),
    mean = mean(ba_observed_richness),
    stdev = sd(ba_observed_richness)
  )

#tests to check independent t-test assumptions
#Assumption 1: are the two samples independent? Yes. not taken from the same time 
#Assumption 2: does the data from each of the 2 groups follow a normal distirbution?
#Use Shapiro-Wilk normaility test
#Null hypothesis: the data are normally distributed
#Alternative hypothesis: the data are not normally distributed

# Shapiro-Wilk normality test for Bloom
with(ba_bloom, shapiro.test(ba_observed_richness[Bloom == "yes"])) # p = 0.008 for site; p-val = 0.22 for bloom
# Shapiro-Wilk normality test for No Bloom
with(ba_bloom, shapiro.test(ba_observed_richness[Bloom == "no"])) # p = 0.447 for site; p=0.63 for bloom
#the two p-values are not greater than the significance level 0.05 implying that the distribution of the 
#data are significantly different from the normal distribution. Ie, we cannot assume the normality.
#if the data are not normally distributed, it’s recommended to use the non parametric two-samples Wilcoxon rank test.

#Wilcoxon test
# Question: Is there any significant changes in the richness of ASV during and not during bloom?
wilco <- wilcox.test(ba_observed_richness ~ Bloom, data = ba_bloom)
wilco

#print p-value only
wilco$p.value
# 0.07 > 0.05 therefore no significant difference between sites
# 0.56 no signif difference between bloom




### sampling depth (reads per sample)
#https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
summary(sample_sums(viral_physeq)) #large difference in number of reads, min=22; max=130183
sort(sample_sums(viral_physeq))

ba <- breakaway(virps3000)

rich_depth <- data.frame("total_reads" =  phyloseq::sample_sums(virps3000),
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



ggplot(data = data.frame("total_reads" =  phyloseq::sample_sums(virps3000),
                         "Years" = virps3000 %>% sample_data %>% get_variable("Years")),
       aes(x = total_reads, y = Years)) +
  geom_bar(stat = "identity") +
  labs(x = "\nTotal Reads", y = "Year\n")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), #rotate axis labels
        plot.title = element_text(hjust = 0.5))+ #center title
  ggtitle("Number of reads per Year")+
  coord_flip()










#### Shannon diversity ####
vir_ps_pel
vir_ps_lit

plot.shannon <- function(ps.obj, plotTitle){
  vir_shannon <- estimate_richness(ps.obj, measures="Shannon")
  
  vir_shannon$sample <- rownames(vir_shannon)
  vir_shannon$Years <- sub("^([^_]*.[^_]*.[^_]*.[^_]*).*$",'\\1', rownames(vir_shannon)) #rm everything after 4th _
  vir_shannon$Years <- sub(".*[/_]", "", vir_shannon$Years) #remove everything before 3rd _
  # gsub("^.*\\_","", vir_shannon$Years) #another way to do same thing
  vir_shannon$date <- sub("^([^_]*.[^_]*.[^_]*.[^_]*).*$",'\\1', rownames(vir_shannon))
  vir_shannon$date <- gsub("^[^_]*_", "",vir_shannon$date) #remove sample name -- just keep date
  vir_shannon$date <- as.Date(vir_shannon$date, format="%d_%m_%Y")
  vir_shannon$bloom <- ps.obj %>% sample_data %>% get_variable("bloom2")
    
  head(vir_shannon)
  
  m <- month(vir_shannon$date)
  d <- day(vir_shannon$date)
  md <- paste(d, m, sep="-")
  
  plot.shan <- ggplot(vir_shannon, aes(x = forcats::fct_inorder(sample), y = Shannon, color = Years))+ #fct_inorder ensures plotting in order of sample date
    geom_point(size=3)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5), #rotate axis labels
          plot.title = element_text(hjust = 0.5))+ #center title
    ggtitle(plotTitle)+
    scale_x_discrete(labels = md, name="Sample date")+ #change x-axis sample name to Month-Day
    scale_y_continuous(name="Shannon diversity")+
    ylim(1.5, 5)
    #facet_grid(~ Years, scales = "free")
  return(plot.shan)
}
shan.lit <- plot.shannon(vir_ps_lit, "Shannon diversity of littoral samples")
shan.pel <- plot.shannon(vir_ps_pel, "Shannon diversity of pelagic samples")

ggpubr::ggarrange(shan.lit, shan.pel, ncol=2, nrow=1, common.legend = T, legend="bottom")



#boxplot years
means <- aggregate(Shannon ~ Years, vir_shannon, mean)
means$step <- 1:nrow(means)

#finding linear regression:
#http://r-statistics.co/Linear-Regression.html 
fit <- lm(Shannon ~ step , data = means)
coefs <- coef(fit)
summary(fit)
names(summary(fit)) #see calls you can make

#get r2
(r2 = signif(summary(fit)$adj.r.squared))
#get p-value
(pval <- signif(summary(fit)$coefficients[2,4]))

#plotting with reg
(shannon_plot <-  ggplot(vir_shannon, aes(x = Years, y = Shannon))+
    geom_point()) + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                                 geom="crossbar", width=0.5) + theme_minimal()+
  theme(plot.title = element_text(hjust=0.5))+
  scale_y_continuous(name = "Shannon diversity")+
  geom_abline(intercept = fit$coefficients[1], slope =  fit$coefficients[2])+
  annotate(geom="text", x = 3, y=4.5, label= paste("Adj R2 = ", r2,
                                                  "p-val = ", pval))+
  labs(title = "Shannon diversity by year")

#### shannon vs. depth ####
shan <- estimate_richness(virps3000, measures="shannon")
head(shan)

viral_df = data.frame("shannon" = shan$Shannon,
                      "sample" = (virps3000 %>% sample_data)$description,
                      "Years" = virps3000 %>% sample_data %>% get_variable("Years"),
                      "depth" = sample_sums(virps3000))
head(viral_df)
str(viral_df)
ggplot(viral_df, aes(x = forcats::fct_inorder(sample), y = shannon, color = Years))+ #fct_inorder ensures plotting in order of sample date
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5), #rotate axis labels
        plot.title = element_text(hjust = 0.5))+ #center title
  ggtitle("Shannon diversity by sample")+
  scale_x_discrete(labels = virps3000 %>% sample_data %>% get_variable("Months"), name="Month")#change x-axis sample name to Month

### sampling depth (reads per sample)
summary(sample_sums(virps3000)) #large difference in number of reads, min=22; max=130183
sort(sample_sums(virps3000))

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


ggplot(data = data.frame("total_reads" =  phyloseq::sample_sums(virps3000),
                         "Years" = virps3000 %>% sample_data %>% get_variable("Years")),
       aes(x = total_reads, y = Years)) +
  geom_bar(stat = "identity") +
  labs(x = "\nTotal Reads", y = "Year\n")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), #rotate axis labels
        plot.title = element_text(hjust = 0.5))+ #center title
  ggtitle("Number of reads per Year")+
  coord_flip()

# Shannon bloom
#boxplot bloom
shannon_bloom = data.frame("Shannon" = vir_shannon$Shannon,
                      "Bloom" = virps3000 %>% sample_data %>% get_variable("bloom2"))
ba_bloom_yn <- na.omit(shannon_bloom[1:2])

(ba_plot <-  ggplot(ba_bloom_yn, aes(x = Bloom, y = Shannon))+
    geom_point()) + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                                 geom="crossbar", width=0.5) + theme_minimal()







#### PERMANOVA ####
#http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova
library(vegan)
package.version("vegan")
citation("vegan")
virps3000 %>% sample_data() %>% head
jsd <- sqrt(phyloseq::distance(virps3000, method = "jsd")) #jsd is more robust to sampling depth
sampledf <- data.frame(sample_data(virps3000)) #make a df from the sample_data
adonis(jsd ~ Period, data = sampledf)

#homogeneity of dispersion test
betadisp <- betadisper(jsd, sampledf$Years)
permutest(betadisp)



### NMDS ###
nmds=metaMDS(comm = sqrt(phyloseq::distance(virps3000, "jsd")), k=3, trymax = 100)
# nmds <- ordinate(physeq = virps3000,
#                  method = "NMDS",
#                  distance = "jsd",
#                  trymax=500)
plot_ordination(physeq = virps3000,
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











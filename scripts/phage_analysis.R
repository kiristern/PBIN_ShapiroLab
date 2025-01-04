# make sure saving in correct conda env path
.libPaths()

# No transformations: Alpha and Beta diversity, Mantel, Procrust
# setwd("/media/nico/MyBook/Virus_amp_seq/2024/scripts/")


#import scripts
source("arxiv/PBIN_ShapiroLab/scripts/functions.R")
source("1_preprocess.R")
library(ggpubr)

## sequences per sample
sums_Phy <- data.frame(colSums(otu_table(viral_physeq)))
colnames(sums_Phy) <- "Sample_TotalSeqs"
sums_Phy$sample <- row.names(sums_Phy)
sums_Phy <- arrange(sums_Phy, Sample_TotalSeqs)
ggplot(sums_Phy, aes(x=reorder(sample, Sample_TotalSeqs), y = Sample_TotalSeqs)) +
  ylab("Number of Sequences per Sample") +
  geom_bar(stat = "identity", colour="black",fill="cornflowerblue")  + xlab("Sample Name") +
  ggtitle("Total Number of Sequences per Sample") +  #scale_y_continuous(breaks =seq(0, 1000000, 10000))+
  theme(axis.text.x = element_text(colour = "black", size=6, angle=45, hjust = 1, vjust = 1))

#Samples that were removed
sums_Phy_filt <- sums_Phy %>% filter(Sample_TotalSeqs < 3000)

#check data
# print(virps3000filt2)

#https://microbiome.github.io/tutorials/
summarize_phyloseq(virps3000filt2)
#sparsity is how populated is the data with zeros.

# sort by date
sample_data(virps3000filt2) <-sample_data(virps3000filt2)[order(sample_data(virps3000filt2)$Date),]

# add new column to sample_data (combine day and month)
sample_data(virps3000filt2)$day_month <- as.factor(paste(sample_data(virps3000filt2)$Day, sample_data(virps3000filt2)$Month, sep = "."))
# drop everything after last "." from description

## not needed -- day_month working with alluvial now.
# split_desc <- str_split(sample_data(virps3000filt2)$description, "\\.", simplify = T)[,1:3]
# # join into single string (save as factor)
# # sample_data(virps3000filt2)$desc_date <- as.factor(apply(split_desc[, 1:3], 1, paste, collapse = "."))
# # another way to do the same thing (save as str)
# sample_data(virps3000filt2)$desc_date <- apply(split_desc[, 1:3], 1, function(x) paste(x, collapse = "."))
# # sample_data(virps3000filt2)$desc_date2 <- NULL # rm accidental added column


# separate into pelagic and littoral phyloseq objects
vir_ps_lit <- subset_samples(virps3000, Site == "Littoral")
vir_ps_pel <- subset_samples(virps3000, Site == "Pelagic")

vir_ps_lit_filt <- subset_samples(virps3000filt2, Site == "Littoral")
vir_ps_pel_filt <- subset_samples(virps3000filt2, Site == "Pelagic")


#Using MicroEco pakcage (Nico)
library(microeco)
library(file2meco)
# littoral
litt_filt_micro <- phyloseq2meco(vir_ps_lit_filt)
t1 <- trans_abund$new(dataset = litt_filt_micro, taxrank = "Species", ntaxa = 20)
# pelagic
pel_filt_micro <- phyloseq2meco(vir_ps_pel_filt)
t2 <- trans_abund$new(dataset = pel_filt_micro, taxrank = "Species", ntaxa = 20)
# filtered
filt_micro <- phyloseq2meco(virps3000filt2)
t3 <- trans_abund$new(dataset = filt_micro, taxrank = "Species", ntaxa = 20)

#### Kiri
# # custom color palette -- doesn't work if want extra customization
# library(ggplot2)
# library(RColorBrewer)
# ## ensure same colors for same species in both plots
# # concat lists and get unique species
# uniq_top <- unique(c(t1$data_taxanames, t2$data_taxanames, t3$data_taxanames)) 
# # define custom color palette to accomodate n_unique species
# custom_colors <- colorRampPalette(
#   brewer.pal(n=9, name="Set1"), 
#   space = "Lab", # helps when colors don't form a natural sequence,
#   )(length(uniq_top)) 
# # make colour dict to ensure same colours for same species in both plots
# color_dict <- setNames(custom_colors, uniq_top)


# documentation: https://www.quantargo.com/help/r/latest/packages/microeco/0.3.3/trans_abund 
# plot top 20 littoral
rel_ab_plot_lit<-t1$plot_bar(
  others_color = "grey90", # http://sape.inf.usi.ch/quick-reference/ggplot2/colour
  # facet = c("Years", "month.numeric"),
  facet="Years",
  # facet2="month.numeric",
  xtext_keep = TRUE,
  legend_text_italic = FALSE,
  # color_values = color_dict,
  # x_axis_name = "Date", # change x-axis name -- specify manually with scale_x_discrete
  barwidth=1,
  use_alluvium = TRUE,
  clustering=TRUE,
) + theme(
  axis.text.y.left = element_text(size = 13),
  axis.text.x = element_text(size = 8, angle =90, vjust = 0.5, hjust = 1, margin = margin(l=10, r=10, t=0, b=0)),
  axis.title.y = element_text(size = 13),
  title = element_text(size = 12),
  legend.position = 'bottom',
) + ggtitle(
  "Top 20 - Littoral"
  ) + guides( # plot legend at bottom in 2 rows
            colour=guide_legend(
                                nrow=2,
                                byrow=TRUE,
                              ),
            fill=guide_legend(
                                nrow=2,
                                byrow=TRUE,
                              )
)
# ) + scale_x_discrete(guide = guide_axis(n.dodge=2) # stagger overlapped text
                            # guide_axis(check.overlap = TRUE) # rm overlapped text
# ) 
plt_lit_ab <- rel_ab_plot_lit + scale_x_discrete(labels = c(sample_data(vir_ps_lit_filt)$day_month))


# plot top 20 pelagic
rel_ab_plot_pel <- t2$plot_bar(
  others_color = "grey90",
  # facet = c("Years", "month.numeric"),
  facet="Years",
  # facet2="month.numeric",
  xtext_keep = TRUE,
  legend_text_italic = FALSE,
  # color_values = color_dict,
  # x_axis_name = "day_month", # change x-axis name
  barwidth=1,
  use_alluvium = TRUE,
  clustering=TRUE,
) + theme(
  axis.text.y.left = element_text(size = 13),
  axis.text.x = element_text(size = 8, angle =90, vjust = 0.5, hjust = 1, margin = margin(l=10, r=10, t=0, b=0)),
  axis.title.y = element_text(size = 13),
  title = element_text(size = 12),
  legend.position = 'bottom',
) + ggtitle(
  "Top 20 - Pelagic"
  ) 
# ) + scale_x_discrete(guide = guide_axis(n.dodge=2) # stagger overlapped text
                            # guide_axis(check.overlap = TRUE) # rm overlapped text
# ) 
# plot legend at bottom in 2 rows
plt_rel_ab_pel <- rel_ab_plot_pel + guides(
                                           colour=guide_legend(
                                                               nrow=2,
                                                               byrow=TRUE,
                                                              ),
                                           fill=guide_legend(
                                                               nrow=2,
                                                               byrow=TRUE,
                                                              )
)
plt_pel_ab <- plt_rel_ab_pel + scale_x_discrete(labels = c(sample_data(vir_ps_pel_filt)$day_month))

#combine both plots to one figure
ggpubr::ggarrange(plt_lit_ab, plt_pel_ab, ncol=1, nrow=2, common.legend = T, legend="bottom")


# require package ggh4x
rel_ab_filt <- t3$plot_bar(others_color = "grey90", 
            facet = c("Site","Years"), 
            xtext_keep = T, 
            legend_text_italic = FALSE,
            color_values = color_dict,
            x_axis_name = "day_month", # change x-axis name
            barwidth=1,
            use_alluvium = TRUE,
            clustering=TRUE,
            ) + 
  theme(axis.text.y.left = element_text(size = 13),
        axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 13),
        title = element_text(size = 12),
        legend.position = 'bottom',
  ) + ggtitle("Top 20 - Littoral & Pelagic") + 
  guides( # plot legend at bottom in 2 rows
            colour=guide_legend(
                                nrow=2,
                                byrow=TRUE,
                              ),
            fill=guide_legend(
                                nrow=2,
                                byrow=TRUE,
                              )
  )
rel_ab_filt


#### nico
#visualize the core microbiota (set of taxa that are detected in a remarkable 
#fraction of the population above a give abundance threshold)
prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-4), log10(.2), length = 10)
virps3000filt2_rel <- microbiome::transform(virps3000filt2, "compositional")
# Also define gray color palette
gray <- gray(seq(0,1,length=5))
p1 <- plot_core(virps3000filt2_rel, 
                plot.type = "heatmap", 
                colours = gray,
                prevalences = prevalences, 
                detections = detections, min.prevalence = .5) +
  xlab("Detection Threshold (Relative Abundance (%))")
p1 <- p1 + theme_bw() + ylab("ASVs")
p1

 


#which asvs are the same between samples
uniq.lit <- unique(taxa_abun_tab_lit$Species)
uniq.pel <- unique(taxa_abun_tab_pel$Species)

uniq.lit[uniq.lit %in% uniq.pel]
uniq.lit[! uniq.lit %in% uniq.pel]
uniq.pel[! uniq.pel %in% uniq.lit]

### Kiri
# t1$data_abund %>% head()
unique_lit <- unique(litt_filt_micro$tax_table$Species)
unique_pel <- unique(pel_filt_micro$tax_table$Species)
unique_lit[unique_lit %in% unique_pel]
unique_lit[! unique_lit %in% unique_pel]
unique_pel[! unique_pel %in% unique_lit]


#### nico
#Rare taxa
library(vegan)
library(TSA)
source("CRT_Functions_v1.1.R")
#virps3000filt <- filter_taxa(virps3000, function(x) sum(x > 1) > (0.10*length(x)), TRUE)
#For diversity analysis and compositional plot
#(virps3000filt2 <- phyloseq::prune_taxa(phyloseq::taxa_sums(virps3000) > 10, virps3000))

virps3000_df <- data.frame(otu_table(virps3000))
write.table(virps3000_df, "virps3000_df.txt", quote = FALSE, sep = "\t")


source("CRT_Functions_v1.1.R")
CRT_analysis<-SimpleRareToPrev.f(otu_fp="virps3000_df.txt", abund_thresh=0.0001, 
                   abund_thresh_ALL=FALSE, b_thresh=0.90, rdp_lastcol=FALSE)

# 2948    6
# 2145    6
# "No. conditionally rare OTUs"
# 2145
# "No. total OTUs"
# 4189
#"Proportion conditional rare / total OTUs"
# 0.5120554
# "No singleton OTUs"
# 4189
# "Proportion conditionally rare / non-singletonOTUs"
# 0.5120554





#### Kiri
#### Shannon with regression line
# https://chiliubio.github.io/microeco_tutorial/diversity-based-class.html#trans_alpha-class
library(magrittr)
filtered_virps <- phyloseq2meco(virps3000filt2) # convert phyloseq to meco
filtered_virps$sample_table$Site %<>% factor(., levels = unique(.)) # groupby Site
t1_alpha <- trans_alpha$new(dataset = filtered_virps, group = "Site", by_group="Years")

# test the differences among groups 
## Kruskal-Wallis Rank Sum Test (overall test when groups > 2) 
## Wilcoxon Rank Sum Tests (for paired groups)
## KW_dunn: Dunn’s Kruskal-Wallis Multiple Comparisons (for paired groups when groups > 2) 
## anova with multiple comparisons 
## Scheirer Ray Hare (nonparametric test) suitable for a two-way factorial experiment
## Linear mixed-effects model--> conditional R2 is the total variance explained by fixed and random effects, and marginal R2 is the variance explained by fixed effects
t1_alpha$cal_diff(method = "wilcox", # KW_dunn | KW | wilcox | t.test | anova | scheirerRayHare | lme
  measure = "Shannon",
  # formula = "Site+Years", # multi-factor analysis of variance, eg, two-way anova
  KW_dunn_letter = TRUE
  )
View(t1_alpha$res_diff)


plot_alpha <- t1_alpha$plot_alpha(
  plot_type = "ggboxplot", # ggboxplot | ggdotplot | ggviolin | ggstripchart | ggerrorplot | errorbar | barerrorbar
  plot_SE = TRUE, 
  measure = "Shannon", 
  add_line = TRUE, 
  line_type = 2,
  add="jitter",
  )
# plot_alpha # boxplot without reg line
# plot_alpha$data %>% head()

### get ggplot params to reuse below
# extract colours used
shannon_palette <- unique(ggplot_build(plot_alpha)$data[[1]]$colour)
# get legend text coords
legend_info <- ggplot_build(plot_alpha)$layout$panel_params
legend_x_coords <- legend_info[[1]]$x.range
legend_y_coords <- legend_info[[1]]$y.range


# test simple abline -- verify overlaying regression line works
# plot_alpha + geom_abline(intercept = 0, slope = 1, color = "red") # works

# plotting Shannon regression
# https://chiliubio.github.io/microeco_tutorial/meconetcomp-package.html#compare-phylogenetic-distances-of-paired-nodes-in-edges
tmp <- t1_alpha$data_stat %>% base::subset(Measure == "Shannon") # subset for Shannon metric
# new mircoeco object
tmp_t1 <- trans_env$new(dataset=NULL, add_data = tmp)
tmp_t1$dataset$sample_table <- tmp_t1$data_env
shannon_reg <- tmp_t1$plot_scatterfit(
  x = "Years", 
  y = "Mean", 
  type = "cor", 
  group = "Site"
  ) + 
    xlab("Years") + 
    ylab("Shannon") + 
    theme(axis.title = element_text(size = 15))
shannon_reg # scatter plot with regression lines
# get the data
reg_lines <- shannon_reg$data

# subset per group
lit_reg <- reg_lines[reg_lines$Group == "Littoral",]
pel_reg <- reg_lines[reg_lines$Group == "Pelagic",]
lm_lit <- lm(y~x, data = lit_reg)
lm_pel <- lm(y~x, data = pel_reg)
# get the coefficients
lm_lit_summ<-summary(lm_lit)
lm_pel_summ<-summary(lm_pel)
# get R-squared and p-value
lm_lit_rsq<-lm_lit_summ$r.squared # R-squared -- fig displays: R-val==sqrt(lm_pel_rsq)=0.562 
lm_lit_pval<-lm_lit_summ$coefficients[2,4] # p-value
lm_pel_rsq<-lm_pel_summ$r.squared # R-squared -- fig displays: R-val==sqrt(lm_lit_rsq)=0.715
lm_pel_pval<-lm_pel_summ$coefficients[2,4] # p-value

## exploration -- see y-intercept to set geom_abline (below)
# lm_lit_df <- fortify(lm_lit)
# lm_lit_df$Year <- factor(lm_lit_df$x, levels = c(seq(2006, 2016, 1))) # convert to factor
# lm_lit_df$reg_fit <- predict(lm_lit, lm_lit_df) # same as .fitted in fortify
# lm_lit$fitted.values


# add regression lines to boxplot
# legend_x_coords # x-coords for legend (extracted above from microeco plot)
plt_shannon <- plot_alpha + 
  # geom_line(data=fortify(lm_lit), aes(x=factor(x), y=.fitted, group=1), color = shannon_palette[1], size=0.8) + 
  # geom_line(data=fortify(lm_pel), aes(x=factor(x), y=.fitted, group=1), color = shannon_palette[2], size=0.8)
  geom_abline(intercept = 3.230468, slope = lm_lit$coefficients[2], 
    color = shannon_palette[1], linewidth=0.5, linetype = "dashed",) +
  geom_abline(intercept = 3.234964, slope = lm_pel$coefficients[2], 
    color = shannon_palette[2], linewidth=0.5, linetype = "dashed",) +
    theme(legend.position=c(0.5, 0.99), # fix legend position
          legend.direction = "horizontal") # horizontal legend

plt_shannon +
  annotate(
    "text", 
    # x=mean(as.numeric(range(levels(plot_alpha$data$Years)))), # center of x-axis
    # y = max(plot_alpha$data$Value)+0.1*diff(range(plot_alpha$data$Value)), # slightly above max y 
    x=(mean(legend_x_coords)-0.5), # center of x-axis
    y=legend_y_coords[2]+0.05*diff(range(legend_y_coords)), # slightly above max y
    label = paste0("R² = ", round(lm_lit_rsq, 3), "\nP = ", signif(lm_lit_pval, 3)),
    hjust = 0.1, vjust=1.0, # place right below legend
    size = 3, color = shannon_palette[1]
    ) +
  annotate(
    "text", 
    x=mean(legend_x_coords)+0.25, # center of x-axis
    y=legend_y_coords[2]+0.05*diff(range(legend_y_coords)), # slightly above max y 
    label = paste0("R² = ", round(lm_pel_rsq, 3), "\nP = ", signif(lm_pel_pval, 3)),
    hjust = 0, vjust=1.0, # place right below legend
    size = 3, color = shannon_palette[2]
  )








#### ALPHA DIV ####
#nico

div<-microbiome::alpha(virps3000filt2, index = c("observed","diversity_shannon","evenness_pielou"))
in.meta <- meta(virps3000filt2)
# Add the rownames as a new colum for easy integration later.
in.meta$sam_name <- rownames(in.meta)
# Add the rownames to diversity table
div$sam_name <- rownames(div)
# merge these two data frames into one
div.df <- merge(div,in.meta, by = "sam_name")
# check the tables
colnames(div.df)
div.df$Years = factor(div.df$Years)
library("ggpubr")
p <- ggboxplot(div.df, 
               x = "Years", 
               y = "diversity_shannon",
               fill = "Years", 
               palette = "jco", add = "jitter", shape="Years" )

p <- p + rotate_x_text()

give.n <- function(x){
  return(c(y = median(x) + 0.1, label = length(x)))
}

p2<-p+ scale_shape_manual(values = 0:10)+     stat_summary(fun.data = give.n, geom = "text", fun = median, position =  position_stack(vjust = .5))
p2+abline(lm(div.df$diversity_shannon~div.df$Years))
print(p2)

div.df2 <- div.df %>% 
  mutate(decimaldat = decimal_date(as.Date(Date, format="%Y-%m-%d")))

ggplot(div.df2, aes(x=decimaldat, y=diversity_shannon)) + 
  geom_point()+
  ylim(0,6)



ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "cyan") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}
fit1 <- lm(diversity_shannon~decimaldat, data = div.df2)
ggplotRegression(fit1) + theme_bw()



#Test for bloom
div.dfnoNA<-div.df[!is.na(div.df$Bloom), ]
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
ggplot(div.dfnoNA, aes(x=Bloom,y=diversity_shannon)) + geom_jitter(position=position_jitter(0.2)) +  stat_summary(fun.data=data_summary, color="blue") + theme_bw()+stat_compare_means(method = "wilcox")



sdiv_lit_fit <- lm(div.df2$decimaldat ~ div.df2$diversity_shannon)
sdiv_lit_coef <- coef(sdiv_lit_fit)
summary(sdiv_lit_fit)
#get r2
(r2 = signif(summary(sdiv_lit_fit)$adj.r.squared))
#get p-value
(pval <- signif(summary(sdiv_lit_fit)$coefficients[2,4]))

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

div.dfnoNA<-div.df[!is.na(div.df$Bloom), ]
p<-boxplot_alpha(virps3000,x_var="Bloom", index = "shannon")
ggplot(div.dfnoNA, aes(x=Bloom,y=diversity_shannon)) + geom_jitter(position=position_jitter(0.2)) +  stat_summary(fun.data=data_summary, color="blue") + theme_bw()+stat_compare_means(method = "wilcox")
wilcox.test(diversity_shannon ~ Bloom, data = div.dfnoNA, exact = FALSE, conf.int = TRUE)


##############################Kiri and Jin
#### Shannon diversity ####
#Shannon diversity by sample
#littoral 
vir_shannon <- estimate_richness(vir_ps_lit_filt, measures="Shannon")

vir_shannon$sample <- rownames(vir_shannon)
vir_shannon$Years <- sub("^([^_]*.[^_]*.[^_]*.[^_]*).*$",'\\1', rownames(vir_shannon)) #rm everything after 4th _
vir_shannon$Years <- sub(".*[/_]", "", vir_shannon$Years) #remove everything before 3rd _
# gsub("^.*\\_","", vir_shannon$Years) #another way to do same thing
vir_shannon$date <- sub("^([^_]*.[^_]*.[^_]*.[^_]*).*$",'\\1', rownames(vir_shannon))
vir_shannon$date <- gsub("^[^_]*_", "",vir_shannon$date) #remove sample name -- just keep date
vir_shannon$date <- as.Date(vir_shannon$date, format="%d_%m_%Y")
vir_shannon$bloom <- vir_ps_lit %>% sample_data %>% get_variable("Bloom")

m <- month(vir_shannon$date)
d <- day(vir_shannon$date)
md <- paste(d, m, sep="-")

ggplot(vir_shannon, aes(x = forcats::fct_inorder(sample), y = Shannon, color = Years))+ #fct_inorder ensures plotting in order of sample date
  geom_point(size=3)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5), #rotate axis labels
        plot.title = element_text(hjust = 0.5))+ #center title
  ggtitle("Shannon diversity of littoral samples")+
  scale_x_discrete(labels = md, name="Sample date")+ #change x-axis sample name to Month-Day
  scale_y_continuous(name="Shannon diversity")+
  ylim(1.5, 5)
#facet_grid(~ Years, scales = "free")
shan <- list(vir_shannon, plot.shan)

length(vir_shannon$sample)
num_samp2 <- c(1:90)
num_samp2
baseplot_sdiv_R2<-data.frame(num_samp2,vir_shannon$Shannon)
baseplot_sdiv_R2

ggplot(baseplot_sdiv_R2, aes(x=num_samp2, y=vir_shannon.Shannon)) + 
  geom_point()+
  ylim(0,6)

sdiv_lit_fit <- lm(num_samp2 ~ vir_shannon$Shannon)
sdiv_lit_coef <- coef(sdiv_lit_fit)
summary(sdiv_lit_fit)
#get r2
(r2 = signif(summary(sdiv_lit_fit)$adj.r.squared))
#get p-value
(pval <- signif(summary(sdiv_lit_fit)$coefficients[2,4]))

#adj.r.square
#[1] 0.0578563
#p-value
#[1] 0.0127477

#Pelagic 
vir_shannon2 <- estimate_richness(vir_ps_pel_filt, measures="Shannon")

vir_shannon2$sample <- rownames(vir_shannon2)
vir_shannon2$Years <- sub("^([^_]*.[^_]*.[^_]*.[^_]*).*$",'\\1', rownames(vir_shannon2)) #rm everything after 4th _
vir_shannon2$Years <- sub(".*[/_]", "", vir_shannon2$Years) #remove everything before 3rd _
# gsub("^.*\\_","", vir_shannon$Years) #another way to do same thing
vir_shannon2$date <- sub("^([^_]*.[^_]*.[^_]*.[^_]*).*$",'\\1', rownames(vir_shannon2))
vir_shannon2$date <- gsub("^[^_]*_", "",vir_shannon2$date) #remove sample name -- just keep date
vir_shannon2$date <- as.Date(vir_shannon2$date, format="%d_%m_%Y")
vir_shannon2$bloom <- vir_ps_pel %>% sample_data %>% get_variable("Bloom")

m <- month(vir_shannon2$date)
d <- day(vir_shannon2$date)
md <- paste(d, m, sep="-")

ggplot(vir_shannon2, aes(x = forcats::fct_inorder(sample), y = Shannon, color = Years))+ #fct_inorder ensures plotting in order of sample date
  geom_point(size=3)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5), #rotate axis labels
        plot.title = element_text(hjust = 0.5))+ #center title
  ggtitle("Shannon diversity of Pelagic samples")+
  scale_x_discrete(labels = md, name="Sample date")+ #change x-axis sample name to Month-Day
  scale_y_continuous(name="Shannon diversity")+
  ylim(1.5, 5)

length(vir_shannon2$sample)
num_samp <- c(1:69)
num_samp
baseplot_sdiv_R2pel<-data.frame(num_samp,vir_shannon2$Shannon)
baseplot_sdiv_R2pel

ggplot(baseplot_sdiv_R2pel, aes(x=num_samp, y=vir_shannon2.Shannon)) + 
  geom_point()+
  ylim(0,6)

sdiv_pel_fit <- lm(num_samp ~ vir_shannon2$Shannon)
sdiv_pel_coef <- coef(sdiv_pel_fit)
summary(sdiv_pel_fit)
#get r2
(r2 = signif(summary(sdiv_pel_fit)$adj.r.squared))
#get p-value
(pval <- signif(summary(sdiv_pel_fit)$coefficients[2,4]))
#adjR2 = -0.00116177
#p-value = 0.340639

#combined plot sdiv by sample date 
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
box.shan.bloom <- box.shannon(virps3000, "Bloom")
box.shan.bloom[[2]]
shan.bloom <- box.shan.bloom[[1]]

box.shan.site <- box.shannon(virps3000, "Site")
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
site.ttest$p.value #  0.2161447 > 0.05; no signif differences b/w sites
bloom.ttest <- t.test(Shannon ~ env.var, data=shan.bloom, var.equal=T)
bloom.ttest$p.value #0.6891384 > 0.05; no signif difference b/w groups




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
  geom_point()+theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 3), #rotate axis labels
        plot.title = element_text(hjust = 0.5))+ #center title
  ggtitle("Shannon diversity by sample")+
  scale_x_discrete(labels = virps3000 %>% sample_data %>% get_variable("Month"), name="Month")#change x-axis sample name to Month

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
  theme_light()+
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
















#### PERMANOVA ####
#http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova
library(vegan)
citation("vegan")
virps3000 %>% sample_data() %>% head
sample_data(virps3000)$Years <- as.factor(sample_data(virps3000)$Years)
sample_data(virps3000)$Month <- as.factor(sample_data(virps3000)$Month)
sample_data(virps3000)$week <- as.factor(sample_data(virps3000)$week)
sample_data(virps3000)$Day <- as.factor(sample_data(virps3000)$Day)

jsd <- sqrt(phyloseq::distance(virps3000, method = "jsd")) #jsd is more robust to sampling depth
sampledf <- data.frame(sample_data(virps3000)) #make a df from the sample_data
adonis2(jsd ~ Period + Site + Years + Month + week + Day, data = sampledf)
adonis2(jsd ~ Site, data = sampledf)
adonis2(jsd ~ Month, data = sampledf)
adonis2(jsd ~ Years, data = sampledf)
adonis2(jsd ~ Day, data = sampledf)
adonis2(jsd ~ week, data = sampledf)
adonis2(jsd ~ Period, data = sampledf)
adonis2(jsd ~ Bloom, data = sampledf,na.action = na.omit)

#Dispersion
#homogeneity of dispersion test
betadisp <- betadisper(jsd, sampledf$Site)
permutest(betadisp)

betadisp <- betadisper(jsd, sampledf$Bloom)
permutest(betadisp,na.action = na.omit)

betadisp <- betadisper(jsd, sampledf$Day)
permutest(betadisp)

betadisp <- betadisper(jsd, sampledf$week)
permutest(betadisp)

betadisp <- betadisper(jsd, sampledf$Month)
permutest(betadisp)

betadisp <- betadisper(jsd, sampledf$Period)
permutest(betadisp)

betadisp <- betadisper(jsd, sampledf$Years)
permutest(betadisp)


##PCoA
sample_data(virps3000)$Years <- as.factor(sample_data(virps3000)$Years)
pcoa=ordinate(virps3000, "PCoA", distance=jsd)
plot_ordination(virps3000, pcoa, color  = "Years") + theme_bw() + scale_colour_manual(values = c("red","blue", "green","brown","purple","yellow","black","grey","pink","orange")) + geom_point(size = 2) + scale_shape_manual(values=c(8, 16, 6)) + theme(axis.text.x  = element_text(vjust=0.5, size=12), axis.text.y  = element_text(vjust=0.5, size=12), axis.title.x = element_text(size = 15, face="bold", color="black"),axis.title.y = element_text(size=15,face="bold",color="black"))


### NMDS ###
nmds=metaMDS(comm = sqrt(phyloseq::distance(virps3000, "jsd")), k=3, trymax = 100)
# nmds <- ordinate(physeq = virps3000,
#                  method = "NMDS",
#                  distance = "jsd",
#                  trymax=500)
plot_ordination(physeq = virps3000,
                ordination = nmds,
                color = as.factor("Years"),
                #shape = "Site",
                title = "NMDS of Lake Champlain viral Communities") + 
  geom_point(aes(color=as.factor(Years)))+
  scale_color_brewer(palette = "Paired")


# geom_point(colour="grey90", size=1.5)
nmds$stress #0.1432 good -- convergence. high stress value means that the algorithm had a hard time representing the distances between samples in 2 dimensions (anything <0.2 is considered acceptable)


#Mantel test
library(readr)
library(dplyr)
library(phyloseq)
library(vegan)
library(ade4)
library(ggplot2)

virus <- as.data.frame(t(otu_table(virps_filt)))
bacteria <- as.data.frame(t(otu_table(bact_physeq)))

# Calculate distance using Bray-Curtis
d_virus <- sqrt(vegdist(virus,"bray"))
d_bacteria <- sqrt(vegdist(bacteria, "bray"))

# Mantel test (using three different methods)
mantel(d_virus, d_bacteria, method = "spearman", permutations = 999) 

plot(d_virus, d_bacteria, type = "p")

procrust.data <- protest(virus, bacteria)
summary(procrust.data)
plot(procrust.data)















library(breakaway)

##error here yet this object doesn't seems to be used in the code
#change to from meta to metadata
ba.dates <- metadata %>% dplyr::select(Date)

vir_ps_pel 
vir_ps_lit 

#richness by year
#pelagic
  ba <- breakaway(vir_ps_pel)
  ymd <- vir_ps_pel %>% sample_data %>% get_variable("Date")
  m <- month(ymd)
  d <- day(ymd)
  md <- paste( d, m, sep="-")
  
  ba_vir_df = data.frame("richness" = (ba %>% summary)$estimate,
                         #"sample" = (ba %>% summary)$sample_names,
                         "error" = (ba %>% summary)$error,
                         "Years" = vir_ps_pel %>% sample_data %>% get_variable("Years"),
                         "Upper" = (ba %>% summary)$upper,
                         "Lower" = (ba %>% summary)$lower,
                         "sample"= vir_ps_pel %>% sample_data %>% get_variable("description"))
  head(ba_vir_df)

  ggplot(ba_vir_df, aes(x = forcats::fct_inorder(sample), y = richness, color = Years))+ #fct_inorder ensures plotting in order of sample date
    geom_point(size=3)+
    theme_classic()+
    geom_errorbar(aes(ymin=richness-abs(richness-Lower), ymax=richness+abs(richness-Upper), width=0.05))+ 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5), #rotate axis labels
          plot.title = element_text(hjust = 0.5))+ #center title
    ggtitle("Breakaway richness of pelagic samples")+
    scale_x_discrete(labels = md, name="Sample date")+ #change x-axis sample name to Month-Day
    scale_y_continuous(name="Richness")+
    ylim(0, 650)+
    geom_smooth(aes(x = as.numeric(forcats::fct_inorder(sample)), y=richness), method = "lm", formula = y~x)+
    annotate("text",x=1,y=1,label=paste0(""))
  
  length(ba_vir_df$sample)
  num_samp <- c(1:72)
  num_samp
  verifyba<-data.frame(num_samp,ba_vir_df$richness)
  verifyba
  
  ggplot(verifyba, aes(x=num_samp, y=ba_vir_df.richness)) + 
    geom_point()+
    ylim(0,600)
  plot(num_samp ~ ba_vir_df$richness)
  fit <- lm(num_samp ~ ba_vir_df$richness)
  coefs <- coef(fit)
  summary(fit)
  #get r2
  (r2 = signif(summary(fit)$adj.r.squared))
  #get p-value
  (pval <- signif(summary(fit)$coefficients[2,4]))
  
#Littoral
  ba <- breakaway(vir_ps_lit)
  ymd <- vir_ps_lit %>% sample_data %>% get_variable("Date")
  m <- month(ymd)
  d <- day(ymd)
  md <- paste( d, m, sep="-")
  
  ba_vir_df2 = data.frame("richness" = (ba %>% summary)$estimate,
                         #"sample" = (ba %>% summary)$sample_names,
                         "error" = (ba %>% summary)$error,
                         "Years" = vir_ps_lit %>% sample_data %>% get_variable("Years"),
                         "Upper" = (ba %>% summary)$upper,
                         "Lower" = (ba %>% summary)$lower,
                         "sample"= vir_ps_lit %>% sample_data %>% get_variable("description"))
  head(ba_vir_df2)
  ba.model<- 
    ggplot(ba_vir_df2, aes(x = forcats::fct_inorder(sample), y = richness, color = Years))+ #fct_inorder ensures plotting in order of sample date
    geom_point(size=3)+
    theme_classic()+
    geom_errorbar(aes(ymin=richness-abs(richness-Lower), ymax=richness+abs(richness-Upper), width=0.05))+ 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5), #rotate axis labels
          plot.title = element_text(hjust = 0.5))+ #center title
    ggtitle("Breakaway richness of littoral samples")+
    scale_x_discrete(labels = md, name="Sample date")+ #change x-axis sample name to Month-Day
    scale_y_continuous(name="Richness")+
    ylim(0, 650)+
    geom_smooth(aes(x = as.numeric(forcats::fct_inorder(sample)), y=richness), method = "lm", formula = y~x)+
    annotate("text",x=1,y=1,label=paste0(""))
  
  length(ba_vir_df2$sample)
  num_samp2 <- c(1:90)
  num_samp2
  verifyba2<-data.frame(num_samp2,ba_vir_df2$richness)
  verifyba2
  
  ggplot(verifyba2, aes(x=num_samp2, y=ba_vir_df2.richness)) + 
    geom_point()+
    ylim(0,700)
  
  
  fit2 <- lm(num_samp2 ~ ba_vir_df2$richness)
  coefs2 <- coef(fit2)
  summary(fit2)
  #get r2
  (r2 = signif(summary(fit2)$adj.r.squared))
  #get p-value
  (pval <- signif(summary(fit2)$coefficients[2,4]))
  
ba.pel <- plot.ba(vir_ps_pel, "Breakaway richness of pelagic samples")
ba.lit <- plot.ba(vir_ps_lit, "Breakaway richness of littoral samples")

ba.pel
ba.lit 
ggpubr::ggarrange(ba.pel, ba.lit, ncol=2, nrow=1, common.legend = T, legend="bottom")

##scatter.ba(vir_ps_pel)

# richness by years
ba_year_lit <- box.years(vir_ps_lit, "Observed richness by year (littoral)")
ba_year_lit
ba_year_pel <- box.years(vir_ps_pel, "Observed richness by year (pelagic)")
ba_year_pel
ba_year <- box.years(virps3000, "Observed richness by year")
ba_year


#richness by bloom/no-bloom and site
box.bloom <- plotbox(virps3000, "Bloom")
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
site.ttest$p.value # 0.0968814 > 0.05; no signif differences b/w groups
bloom.ttest <- t.test(ba_observed_richness ~ env.var, data=box.bloom.df, var.equal=T)
bloom.ttest$p.value # 0.506417 > 0.05; no signif difference b/w sites




### sampling depth (reads per sample)
#https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
summary(sample_sums(library(ggpubr)
ggdensity(box.site.df$ba_observed_richness)
ggdensity(box.bloom.df$ba_observed_richness))) #large difference in number of reads, min=22; max=130183
sort(sample_sums(virps3000))

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
                           "p-val = ", pval))+
  theme_light()



ggplot(data = data.frame("total_reads" =  phyloseq::sample_sums(virps3000),
                         "Years" = virps3000 %>% sample_data %>% get_variable("Years")),
       aes(x = total_reads, y = Years)) +
  geom_bar(stat = "identity") +
  labs(x = "\nTotal Reads", y = "Year\n")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), #rotate axis labels
        plot.title = element_text(hjust = 0.5))+ #center title
  ggtitle("Number of reads per Year")+
  ylim(2005,2020)+
  coord_flip()




#https://bioconductor.riken.jp/packages/3.8/bioc/vignettes/microbiome/inst/doc/vignette.html
#Visually-Weighted Regression curve with smoothed error bars
# Estimate Shannon diversity and add it to the phyloseq object
sample_data(virps3000filt)$ShannonDiv <- metadata$ShannonDiv <- virps3000filt %>% otu_table %>% microbiome::alpha() %>% select("diversity_shannon")

#compare year and microbiome shannon diversity
microbiome::plot_regression(ShannonDiv ~ Years, metadata) #doesn't work!



# #compositional heatmap
# tmp <- plot_composition(filt_virseq, plot.type="heatmap", transform = "hellinger")
# 
# #compositional barplot
# plot_composition(transform(filt_virseq, "compositional"), 
#                  plot.type = "barplot", sample.sort = "neatmap", label=F)

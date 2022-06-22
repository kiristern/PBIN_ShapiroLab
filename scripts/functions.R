### Functions ###
library(lubridate)

## set vir_hel to same format as bact_hel
colsamp2date <- function(tab2format){
  colnames(tab2format) <- sub("*._*._*._*._*._*._*._","", colnames(tab2format))
  colnames(tab2format) <- gsub("_", ".", colnames(tab2format))
  return(tab2format)
}

## Relative abundance##
# Get top 20
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

#create date (month-day) col
library(lubridate)
getMonthDay <- function(table){
  m <- month(table$Date)
  day <- day(table$Date)
  md <- paste(day, m, sep="-")
  return(md)
}

#plot relative abundance ### X-axis by season not date.
plotRelAb <- function(rel_ab_tab, md, plotTitle){
  rel_ab_plot <- rel_ab_tab %>% 
    ggplot(aes(x =Sample, y = Abundance, fill = species, order = -species)) +
    geom_bar(stat = "identity")+
    scale_fill_manual("ASV", values = dd.col, drop = F)+ #drop=F to prevent dropping of unused factor levels
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

#breakaway richness plot
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

#boxplot richness by years
box.years <- function(ps.obj, plotTitle){
  ba <- breakaway(ps.obj)
  ba_year <- data.frame("ba_observed_richness" = (ba %>% summary)$estimate,
                        "Years" = ps.obj %>% sample_data %>% get_variable("Years"))
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
      labs(title = plotTitle)+
      stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                   geom="crossbar", width=0.5)+ 
      theme_minimal())
  return(ba_plot)
}

# richness by bloom/no bloom and site
plotbox <- function(ps.obj, env.var){
  ba <- breakaway(ps.obj)
  ba.df = data.frame("ba_observed_richness" = (ba %>% summary)$estimate,
                     env.var = ps.obj %>% sample_data %>% get_variable(env.var))
  ba_bloom_yn <- na.omit(ba.df[1:2])
  
  ba_plot <-  ggplot(ba_bloom_yn, aes(x = env.var, y = ba_observed_richness))+
    geom_point() + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                                geom="crossbar", width=0.5) + theme_minimal()+
    scale_x_discrete(name=env.var)
  ba.box.list <- list(ba_bloom_yn, ba_plot)
  return(ba.box.list)
}


#Shannon diversity by sample
plot.shannon <- function(ps.obj, plotTitle){
  vir_shannon <- estimate_richness(ps.obj, measures="Shannon")
  
  vir_shannon$sample <- rownames(vir_shannon)
  vir_shannon$Years <- sub("^([^_]*.[^_]*.[^_]*.[^_]*).*$",'\\1', rownames(vir_shannon)) #rm everything after 4th _
  vir_shannon$Years <- sub(".*[/_]", "", vir_shannon$Years) #remove everything before 3rd _
  # gsub("^.*\\_","", vir_shannon$Years) #another way to do same thing
  vir_shannon$date <- sub("^([^_]*.[^_]*.[^_]*.[^_]*).*$",'\\1', rownames(vir_shannon))
  vir_shannon$date <- gsub("^[^_]*_", "",vir_shannon$date) #remove sample name -- just keep date
  vir_shannon$date <- as.Date(vir_shannon$date, format="%d_%m_%Y")
  vir_shannon$bloom <- ps.obj %>% sample_data %>% get_variable("Bloom")
  
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
  shan <- list(vir_shannon, plot.shan)
  return(shan)
}

#Shannon by years
shan.years <- function(df){
  means <- aggregate(Shannon ~ Years, df, mean)
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
  shannon_plot <-  ggplot(df, aes(x = Years, y = Shannon))+
      geom_point() + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                                   geom="crossbar", width=0.5) + theme_minimal()+
    theme(plot.title = element_text(hjust=0.5))+
    scale_y_continuous(name = "Shannon diversity")+
    geom_abline(intercept = fit$coefficients[1], slope =  fit$coefficients[2])+
    annotate(geom="text", x = 3, y=4.5, label= paste("Adj R2 = ", r2,
                                                     "p-val = ", pval))+
    labs(title = "Shannon diversity by year")
  return(shannon_plot)
}

#boxplot shannon bloom/no bloom and site
box.shannon <- function(ps.obj, env.var){
  shannon_bloom = data.frame("Shannon" = vir_shannon[[1]]$Shannon,
                             env.var = ps.obj %>% sample_data %>% get_variable(env.var))
  shan_yn <- na.omit(shannon_bloom[1:2])
  
  shan_plot <-  ggplot(shan_yn, aes(x = env.var, y = Shannon))+
    geom_point() + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                                geom="crossbar", width=0.5) + theme_minimal()+
    scale_x_discrete(name=env.var)
  shan.list <- list(shan_yn, shan_plot)
  return(shan.list)
}

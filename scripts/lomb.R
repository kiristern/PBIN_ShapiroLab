# source("lomb_fnc.R")

library(tidyverse)
library(phyloseq)
library(lubridate)
library(lomb)
library(ggplot2)

# load seasonal data from nico
bact_szn <- read.csv("lomb_seasonality_bact.csv")
vir_szn <- read.csv("lomb_seasonality_virus.csv")

vir_szn %>% head()

############################
# Import from .RData
############################
# load from backup -- faster -- use if only need phyloseq data
lomb_env <- new.env() # Create a new environment
# load("arxiv/Large_data/lomb backup.RData", envir = lomb_env) # Load the file into the new environment
load("lomb2025.RData", envir = lomb_env) # Load the file into the new environment
ls(lomb_env)  # List the variables in the new environment
### lomb.sea.02V is the filtered output
### lomb.02v is the unfiltered output of lomb

lomb_filt <- lomb_env$lomb.sea.02V
dplyr::glimpse(lomb_filt) # compact view of data

virps <- lomb_env$virps3000filt
bactps <- lomb_env$bact3000filt

##############
# Modify taxonomy table -- simplify
# tax.raw <- bactps %>% tax_table() %>% as.data.frame() 

# # subset taxa
# tax.genus <- c("g__Microcystis", "g__Dolichospermum", "g__Synechococcus", "g__Fimbriimonas", "g__Pseudomonas","g__Steroidobacter", "g__Nitrospira", "g__Pseudanabaena")

# # Access the taxonomy table
# # Extract the taxonomy table
# tax_table_df <- as.data.frame(tax_table(bactps))

# # Step 2: Identify taxa to keep and modify others
# tax_table_df <- tax_table_df %>%
#   mutate(
#     Kingdom = ifelse(Genus %in% tax.genus, Kingdom, "Other"),
#     Phylum = ifelse(Genus %in% tax.genus, Phylum, "Other"),
#     Class = ifelse(Genus %in% tax.genus, Class, "Other"),
#     Order = ifelse(Genus %in% tax.genus, Order, "Other"),
#     Family = ifelse(Genus %in% tax.genus, Family, "Other"),
#     Genus = ifelse(Genus %in% tax.genus, Genus, "Other")
#   )
# # Step 3: Update the taxonomy table in phyloseq object
# tax_table(bactps) <- as.matrix(tax_table_df)

# # check unique taxa
for (i in colnames(tax_table(virps))) {
  print(i) # Print the taxonomic rank (column name)
  # Extract the column and get unique values
  if (i != "Species"){
    unique_values <- unique(tax_table(virps)[, i])
    print(unique_values) # Print unique values for the rank
  }
}
##########################

relab.bact <- transform_sample_counts(bactps, function(x) x / sum(x))
relab.virps <- transform_sample_counts(virps, function(x) x / sum(x))
taxa_sums(relab.virps) %>% head()

### check to make sure ASV abundances are not all zero
# dat <- psmelt(virps)
# dat %>% head()
# # make decimal date
# dat$decimaldat <- lubridate::decimal_date(as.Date(dat$Date))

# dat$OTU %>% head()
# dat %>% select(OTU, decimaldat, Abundance) %>%
# group_by(OTU) %>%
# summarize( howmany0 = sum(Abundance == 0)) %>%
# arrange(-howmany0)

# ASVs0 <- dat %>%
# select(OTU, decimaldat, Abundance) %>%
# group_by(OTU) %>%
# summarize( howmany0 = sum(Abundance == 0)) %>%
# arrange(-howmany0) %>%
# filter(howmany0 == 20) %>%
# pull(OTU)

# # seasonal niche 
# lomb_object <- dat %>%
# filter(! OTU %in% ASVs0) %>%
# split(.$OTU) %>%
# map(~randlsp( x =.x$Abundance,
# times = .x$decimaldat,
# type = 'period',
# plot = F))
####################


source("bbmo_timeseries_fraction/src/sourcefiles/main_functions.R")
# psmelt.raw <- psmelt_dplyr(relab.virps) 
# psmelt.raw %>% head()

# all <- psmelt.raw %>% 
#   filter( Abundance > 0) %>% 
#   mutate(seasonal = case_when(
#     OTU %in% vir_szn$asv & Site == 'Littoral' ~ 'seasonal',
#     OTU %in% vir_szn$asv & Site == 'Pelagic' ~ 'seasonal',
#     TRUE ~ 'no seasonal'
#   )) %>% 
#   group_by(Site,OTU,seasonal) %>% 
#   dplyr::summarize(Abundance = sum(Abundance),
#                    count = unique(OTU) %>% length()) 
# all %>% head()


# all.perfect <- all %>%
#   group_by(Site, OTU) %>% 
#   dplyr::mutate(relab.ab = Abundance / sum(Abundance),
#                    relab.count = count / sum(count)) %>% 
#   ungroup() %>% 
#   mutate(fraction = if_else(Site == "Littoral",
#                                    "Littoral" , "Pelagic")) #%>% 
#   # italizyce_synes()
# all.perfect %>% head()

# ggplot(all.perfect) +
#   geom_bar(aes( x = OTU,
#                 y = relab.ab, fill = seasonal), stat='identity') +
#   geom_point(data = filter(all.perfect, seasonal == "seasonal"),
#              aes( x = OTU, y = relab.count, color = seasonal) ) +
#   geom_text(data = filter(all.perfect, seasonal == "seasonal"),
#             aes(x = OTU, y = 1.05, label = count ),
#             size = 2.5, 
#             # color = economist_pal()(1)
#             ) + 
#   geom_hline(yintercept = 0.5, linewidth=.2, linetype="dotted", color = "grey") +
#   coord_flip() +
#   facet_wrap(~Site) +
#   scale_fill_manual(values = c('#A9D0D0', '#60ABAE'), name = "") + 
#   # scale_colour_economist() +
#   theme_minimal() +
#   # lil.strip +
#   theme( axis.text.y = ggtext::element_markdown()) + 
#   # scale_x_discrete(limits = rev(tax.order.palette)) +
#   scale_y_continuous( labels = scales::percent) + 
#   guides(color = FALSE) +
#   ylab("Relative abundance") +
#   xlab( "Taxonomy" ) 

# ggsave("./figs250109/summary-seasonality-lomb-bact.png",
#        width = 9, height = 7)


#### Distributino along year
sea.asvs <- unique(vir_szn$asv)
taxV <- as(tax_table(virps), 'matrix') %>% as_tibble(rownames = 'asv')

sea.category <-  data.frame( asv = taxV$asv) %>% 
  mutate( seasonal = case_when( asv %in% sea.asvs ~ 'seasonal',
                                TRUE ~ 'no seasonal' ), 
          seasonal.type = case_when(
            asv %in% vir_szn$asv ~ 'seasonal',
            TRUE ~ 'no seasonal'),
          seasonal.type = factor(seasonal.type, 
                                 levels = c('seasonal',
                                            'no seasonal'))
  ) 

psmelt.relab <- psmelt_dplyr(relab.virps) 

psmelt.relab %>% 
  left_join(sea.category, by = c('OTU' = 'asv')) %>% 
  group_by(Site, Sample, Month, seasonal.type
    ) %>% 
  summarize(relab = sum(Abundance)) %>% 
  ggplot( aes(Month, relab, color = seasonal.type
    )
  ) + 
  geom_boxplot(aes(fill = seasonal.type), 
               position = position_dodge(width = 0.7),
               alpha = 0.7, outlier.colour = 'transparent') + 
  geom_point(position= position_dodge(width = 0.7))  + 
  facet_wrap(~ str_c('Site ', Site), ncol = 1) + 
  # lil.strip + 
  # scale_x_month + 
  # scale_fill_pander() + 
  # scale_color_pander() + 
  scale_y_continuous(labels = scales::percent, name = 'Relative abundance (%)')

ggsave('./figs250109/vir_str_typeseasonality.png',
       width = 9, height = 6)


# Occurence
psmelt.relab %>% 
  left_join(sea.category, by = c('OTU' = 'asv')) %>% 
  filter( Abundance > 0) %>% 
  group_by(Site, Sample, Month,  seasonal.type) %>% 
  summarize(count = n()) %>% 
  ggplot( aes(Month, count, color = seasonal.type)) + 
  geom_boxplot(aes(fill = seasonal.type), 
               position = position_dodge(width = 0.7),
               alpha = 0.7, outlier.colour = 'transparent') + 
  geom_point(position= position_dodge(width = 0.7))  + 
  facet_wrap(~ str_c('Site ', Site), ncol = 1) + 
  # lil.strip + 
  # scale_x_month + 
  # scale_fill_pander() + 
  # scale_color_pander() + 
  ylab('Ocurrence (n. ASVs)')
ggsave('./figs250109/vir_nOccurence_typeseasonal.png',
       width = 9, height = 6)

# Polar plot
polar_plot <- function(df, yaxis = 'peak', shape = 'Site', link = 'OTU'){
  
  df %>% 
  ggplot( aes_string("Month", yaxis)) + 
  geom_vline(data =  data.frame( values = seq(0.5, 12.5, by = 1)),
             aes(xintercept = values), color = 'gray') + 
  geom_jitter(aes_string(fill = "OTU", shape = shape),
              size = 2.4,
              width = 0.3)  +
  coord_polar(theta = 'x', start = 12.3) + 
  guides( fill = guide_legend( override.aes = list(size = 3, shape = 21)))  + 
  guides( shape = guide_legend( override.aes = list(size = 3)))  + 
  # lil.strip + 
  # bac.fillScale +
  # scale_x_month + 
  cowplot::theme_minimal_hgrid() 
}

# power.lomb <- bind_rows(
#   list( `0.2` = lomb_filt, 
#         `3` =  lomb_filt),
#   .id = 'Site')

maxima.median <- psmelt.relab %>% 
  left_join(sea.category, by = c('OTU' = 'asv')) %>% 
  group_by(Site,OTU, Month) %>% 
  summarize(median.relab = median(Abundance)) %>% 
  group_by(Site,OTU) %>% 
  filter(median.relab == max(median.relab)) %>% 
  arrange(OTU)

szn_vir_lomb <- left_join(maxima.median, lomb_filt, 
          by = c('OTU' = 'asv')) %>% 
  left_join(sea.category, by = c('OTU' = 'asv')) %>% 
  select(Site, OTU, Month, peak, seasonal.type) %>% 
  left_join(taxV, by = c('OTU' = 'asv'))  
  # italizyce_synes() %>% 

szn_vir_lomb %>% 
  polar_plot(shape = "seasonal.type") + 
  ylab('Strength recurrence') + 
  facet_wrap(~Site, ncol = 2)  + 
  # italics.strip.legend + 
  scale_shape_manual(values = 21:25) +
  theme(legend.position = "none")

ggsave('./figs250109/polar_plot_general.png', width = 9, height = 9)






# Bps <- lomb_env$bact_physeq # unfiltered
bact_ps <- lomb_env$bact3000filt
vir_ps <- lomb_env$virps3000filt

# add day of year col to metadata
sample_data(bact_ps)$day_of_year <- as.numeric(format(as.Date(sample_data(bact_ps)$Date), format = "%j"))

tax_table(bact_ps) %>% head()
tax_table(vir_ps) %>% head()


b.phy.relab <- transform_sample_counts(bact_ps, function(x) x / sum(x))
taxa_sums(b.phy.relab) %>% head()

#########################
# taxa abundance prevalence
#########################
source('CRT_fncs.R')
# abundance prevalence DF
abunprevB <- calculate_abund_prevalen.df(b.phy.relab)

# incorporate conditionally rare taxa
abunprev.crt <- incorporate_crt_results(abunprevB, bact_ps, string = 'fl') %>% 
  mutate( behavior = ifelse(crt, 'CRT', behavior), 
          presence = factor(presence),
          behavior = factor(behavior,
                            levels = c('CRT', 'Narrow', 'Intermediate', 'Broad')))

taxB <- as(tax_table(bact_ps), 'matrix') %>% as_tibble(rownames = 'asv')
taxB %>% head()

# distribution of behaviours
abprevtax <- abunprev.crt %>%
    left_join(taxB, by='asv') %>%
    mutate(seasonal = ifelse(asv %in% bact_szn$asv, TRUE, FALSE))

abprevsea <- abprevtax %>%
    filter(asv %in% bact_szn$asv) %>%
    select(asv, behavior, presence)

abprevsea %>% 
  select(behavior, presence) %>% 
  table()

abprevsea %>% filter(behavior == 'Broad') 
abunprev.crt %>% filter(behavior == 'Broad') %>% 
  left_join(taxB, by = 'asv')

abprevtax %>% filter(behavior == 'Narrow') %>% pull(Class) %>% table()
abprevtax %>% filter(behavior == 'Broad') %>% pull(Class) %>% table()


#######################
# # plotting
lil.strip <- theme(strip.background = element_blank(),
                   strip.text.x =element_text(margin = margin(.05, 0, .1, 0, "cm")))

geoMean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

b.phy.clr <- bact_ps %>% 
  transform_sample_counts(., function(x) x+1) %>%
  transform_sample_counts(function(x) log(x/geoMean(x)))
taxa_sums(b.phy.clr) %>% head()


# plot
# # cumulative day numbers for the start of each month
num.days.mnt <- c(0,31,28,31,30,31,30,31,31,30,31,30)
cumnum <- cumsum(num.days.mnt)
print(cumnum)
# month order
month.order <- c('january', 'february', 'march', 'april', 'may', 'june', 'july', 'august', 'september', 'october', 'november', 'december')


# plot - CLR
sea.asvs <- psmelt(b.phy.clr) %>% 
  filter(OTU %in% c('ASV_121', 'ASV_138', 'ASV_17')) %>% 
  ggplot(aes(day_of_year, Abundance)) + 
  geom_jitter(aes(color = OTU), alpha = 0.6, show.legend = F) + 
  stat_smooth(aes(x = day_of_year,
                  group = OTU,
                  color = OTU), # continuous x-axis
              method = "gam",
              formula = y ~ s(x, k =12, bs = 'cc'),
              se = F, size = 2, show.legend = F) + 
  facet_wrap(~str_c(Class, ', ', str_to_upper(OTU))) + 
  lil.strip + 
  ylab('Centered log ratio (log10( ASV read count / geoMean))') + 
  scale_x_continuous(breaks = cumnum,
                     name = 'Month',
                     labels = str_to_title(month.order) %>% str_sub(1,3)
                    )
sea.asvs

# relab
sea.asvs <- psmelt(b.phy.relab) %>% 
  filter(OTU %in% c('ASV_121', 'ASV_138', 'ASV_17')) %>% 
  ggplot(aes(day_of_year, Abundance)) + 
  geom_jitter(aes(color = OTU), alpha = 0.6, show.legend = F) + 
  stat_smooth(aes(x = day_of_year,
                  group = OTU,
                  color = OTU), # continuous x-axis
              method = "gam",
              formula = y ~ s(x, k =12, bs = 'cc'),
              se = F, size = 2, show.legend = F) + 
  facet_wrap(~str_c(Class, ', ', str_to_upper(OTU))) + 
  lil.strip + 
  ylab('Relative Abundance') + 
  scale_x_continuous(breaks = cumnum,
                     name = 'Month',
                     labels = str_to_title(month.order) %>% str_sub(1,3)
                    )
sea.asvs



#######################
# Comparison trends ASVs
#######################
tax_table(bactps) %>% head()

vals <- abprevtax %>% 
  pull(Genus) %>% 
  table() %>% sort() 
vals

# compare g__Dolichospermum and g__Microcystis
genera <- vals %>% 
           names() %>% tail(n=5)

asvs <- abprevtax %>% 
  filter(Genus %in% c(genera, 'g__Microcystis')) %>% 
  pull(asv)




# plot
b.phy.asinh <- transform_sample_counts(bact_ps, function(x) asinh(x))
taxa_sums(b.phy.asinh) %>% head()

psmelt(b.phy.asinh) %>%
    filter(OTU %in% asvs) %>%
    ggplot(aes(Month, Abundance)) + 
    geom_jitter(alpha = 0.6) + 
    stat_smooth(aes(x = Month,
                    group = OTU,
                    color = Family),
                method = "gam",
                formula = y ~ s(x, k =12, bs = 'cc'),
                se = F, size = 2,
                show.legend = F, alpha = 0.7) + 
    facet_wrap(~Genus)




Goverview <- taxB %>% 
  filter(Genus %in% c('g__Microcystis', 'g__Dolichospermum')) %>% #| order %in% 'HIMB59' ) %>% 
  pull(asv)

psm.asinh %>% 
  filter(OTU %in% Goverview) %>% 
  mutate(seasonal = ifelse(OTU %in% bact_szn$asv, TRUE, FALSE)) %>% 
  ggplot(aes(Month, Abundance)) + 
    geom_jitter(alpha = 0.6) + 
    stat_smooth(aes(x = Month,
                    group = OTU,
                    color = Genus),
                method = "gam",
                formula = y ~ s(x, k =12, bs = 'cc'),
                se = F, size = 2,
                show.legend = F, alpha = 0.7) + 
    facet_wrap(~Genus + seasonal)







# ##### Get the data from the new environment
# data_lomb <- lomb_env$lomb.02
# # data_lomb <- lomb_env[["lomb.sea.02"]] # alternative way to get the data
# dplyr::glimpse(data_lomb) # compact view of data

# # # select the data (change threshold b/c saved $lomb.sea.02 is empty)
# lomb.sea.02 <- tibble( asv = names(data_lomb),
#                        pval = map_dbl(data_lomb, ~.x$p.value),
#                        peak = map_dbl(data_lomb, ~.x$peak),
#                        interval  = map(data_lomb, ~.x$peak.at),
#                        int.min = map_dbl(interval, ~.[[2]]),
#                        int.max = map_dbl(interval, ~.[[1]])) %>% 
#   mutate( qval = fdrtool::fdrtool(pval, statistic = 'pvalue')$qval) %>% 
# #   filter(qval <= 0.01, peak >= 0.150, int.max <= 2)
#     mutate(seasonality = ifelse(pval < 0.01 & peak >= 0.150 & int.max <= 2,
#                                    'seasonal', 'non_seasonal')) # add column indicating if ASV is seasonal or not
# lomb.sea.02 %>% dplyr::glimpse() # compact view of data

# # View data
# lomb.sea.02 %>% 
#   mutate( index = str_replace(asv, 'asv', '') %>% as.integer()) %>% 
#   arrange(index) %>% View()

# # filter for seasonal and non-seasonal ASVs
# seasonalASV <- lomb.sea.02 %>% filter(seasonality == 'seasonal') # filter for seasonal ASVs
# nonseasonalASV <- lomb.sea.02 %>% filter(seasonality == 'non_seasonal') # filter for non-seasonal ASVs

# ##### also in saved data (call lomb_env$results.lomb02)
# results.lomb02_seasonal <- data_lomb[seasonalASV %>% pull(asv)]
# dplyr::glimpse(results.lomb02_seasonal) # compact view of data

# # plotting
# asv_periods <- map(results.lomb02_seasonal, ~tibble( scanned = .x$scanned,
#                                              power = .x$power)) %>% 
#                                             bind_rows(.id = 'asv') %>% 
#                                             split(.$asv)
# asv_periods %>% head() # check data


# periodoplots <- asv_periods %>% 
#   map(~ggplot(.x, aes(scanned, power)) + 
#         geom_line(aes(group = asv)) + 
#         facet_wrap(~asv) + 
#         lil.strip
#         )
# periodoplots$ASV_30 # example plot for ASV_30



##############
# check unique taxa
for (i in colnames(tax_table(bact_ps))) {
  print(i) # Print the taxonomic rank (column name)
  # Extract the column and get unique values
  if (i != "ASV"){
    unique_values <- unique(tax_table(bact_ps)[, i])
    print(unique_values) # Print unique values for the rank
  }
}
# only have Virus (kingdom) & Myoviridae (family) in this virps
# bact_ps has more taxonomic ranks

############## 
bact_szn$asv 
vir_szn$asv
# seasonal ASVs bact: 1069, 1149
# seasonal ASVs vir: 10, 100
asv.sel <- str_c('ASV_', c(221, 25, 19, 145, 161, 61,66,71,68, 93, 153, 203, 34, 40, 47))
# asv.sel <- bact_szn$asv #[1:10] # select first 5 ASVs for testing



# run for virps and bact_ps separately
ps <- bact_ps
# ps <- vir_ps.copy
taxa_sums(ps) %>% head() # sum of taxa (ASVs x taxa) -- raw counts
ps <- transform_sample_counts(ps, function(x) x / sum(x))
taxa_sums(ps) %>% head() # sum of taxa (ASVs x taxa) -- transformed relab

# check phyloseq object
tax_table(ps) %>% head() # taxonomy table (ASVs x taxa)
otu_table(ps) %>% head() # abundance counts (ASVs x samples) -- transformed relab (line 54)
sample_data(ps) %>% head() # sample metadata (samples x metadata)

# just subset for select ASV
ps_specific <- psmelt(ps) %>%
    as_tibble() %>%
    filter(OTU %in% asv.sel) # filter for specific ASVs



###################
# plot seasonality
###################
# https://github.com/adriaaulaICM/bbmo_niche_sea/blob/6cef1b004e75a88a007975f6c5ebc37a40d32b0e/src/figures/sea_explanation.R#L45
gam.gg <- ggplot(data = ps_specific, aes(day_of_year,Abundance)) + 
  geom_jitter(aes(color = OTU), alpha = 0.4) + 
  stat_smooth(aes(x = day_of_year,
                  group = OTU,
                  color = OTU
                  ),
              method = "gam",
              formula = y ~ s(x, k =12, bs = 'cc'),
              se = TRUE, 
              size = 0.5,
              show.legend = TRUE, 
              alpha = 0.3
              ) + 
  # split into separate facets (plots) for each OTU
#   facet_wrap(#~str_c(OTU),
#             ~str_c(str_to_upper(OTU), Class, sep = ', '),
#             scales = 'free_y'
#             ) + 
  # display y-axis as percentage
#   scale_y_continuous(labels = scales::percent_format(accuracy = 2L)) +
  scale_x_continuous(
                    breaks = cumnum,
                     # Display month abbrv first 3 letters
                     labels = str_to_title(month.order) %>% str_sub(1,3),
                     name = 'Month',
                     ) +
#   guides(color = "none") + 
  ylab('Relative abundance') + 
  lil.strip + 
  theme(
        # legend.position = 'bottom',
        axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels
  )+
  ggtitle('Seasonal ASVs relative abundance (bacteria)')
gam.gg

# periodo.gg <- map(data_lomb, ~tibble( scanned = .x$scanned,
#                                     power = .x$power)) %>%
#   bind_rows(.id = 'asv') %>% 
#   filter(asv %in% asv.sel) %>% 
#   mutate( desc = ifelse(asv == 'ASV_10', 'Strong peak with a periodicity of 1 year',
#                         'No patterns, random walk') %>% as.factor() %>% fct_inorder()) %>% 
#   ggplot(aes(scanned, power)) + 
#   geom_line( aes(group = asv)) + 
#   ylab('Strength recurrence') + 
#   # TODO: have multiple years... is strength reoccurrence computed across all years?
#   xlab('Periods checked (0.083 ~ 1/12 months)') + 
#   facet_wrap(~desc) +
#   lil.strip
# periodo.gg
# # gam.gg + periodo.gg + plot_layout(ncol = 1)

# # save plot
# figpath <- 'figs250104' # today's date
# ggsave(filename = str_c(figpath, 'example_seasonality.png'), plot=gam.gg,
#        width = 11, height = 7)


###############################
# seasonality contribution to total relative abundance
################################
# https://github.com/adriaaulaICM/bbmo_niche_sea/blob/6cef1b004e75a88a007975f6c5ebc37a40d32b0e/src/figures/seasonal_numbers.Rmd#L388
ls(lomb_env)
# lomb.sea.02 %>% dplyr::glimpse()

bact_ps <- transform_sample_counts(bact_ps, function(x) x / sum(x))
vir_ps <- transform_sample_counts(vir_ps, function(x) x / sum(x))

# seasonal ASVs bacteria
season_relab_bact <- taxa_sums(bact_ps)[bact_szn$asv] %>% sum(na.rm = TRUE)
total.relab_bact <- taxa_sums(bact_ps) %>% sum()
proportion_szn_of_tot_relab_bact <- season_relab_bact * 100 / total.relab_bact
print(proportion_szn_of_tot_relab_bact)

# seasonal ASVs virus
season_relab_vir <- taxa_sums(vir_ps)[vir_szn$asv] %>% sum(na.rm = TRUE)
total.relab_vir <- taxa_sums(vir_ps) %>% sum()
proportion_szn_of_tot_relab_vir <- season_relab_vir * 100 / total.relab_vir
print(proportion_szn_of_tot_relab_vir)


####################
# contribution of each behaviour
####################
sum(taxa_sums(b.phy.relab)[abprevsea %>% filter(behavior == 'Broad') %>% pull(asv)]) * 100 / total.relab_bact
sum(taxa_sums(b.phy.relab)[abprevsea %>% filter(behavior == 'Narrow') %>% pull(asv)]) * 100 / total.relab_bact
sum(taxa_sums(b.phy.relab)[abprevsea %>% filter(behavior %in% c('Other', 'CRT')) %>% pull(asv)]) * 100 / total.relab_bact



###############################
# plotting
################################
# bacterial seasonal
szn_bact_otu <- prune_taxa(taxa_names(bact_ps) %in% bact_szn$asv, bact_ps)
sam.sea_bact <- otu_table(szn_bact_otu) %>% colSums()
total.sea_bact <- otu_table(bact_ps) %>% colSums()

relabs_bact <- sam.sea_bact / total.sea_bact 

mean(relabs_bact)
sd(relabs_bact)

tmpB <-  data.frame(season_relab_bact = relabs_bact, 
           sample_data(bact_ps))

tmpB %>% 
  ggplot( aes(season_relab_bact)) + 
  geom_histogram()

tmpB %>% head()
tmpB$Date <- as.Date(tmpB$Date)
# arrange by date
tmpB <- tmpB %>%
  arrange(Date)
tmpB$Year <- format(tmpB$Date, "%Y")  # Extract year as a character or integer

# fix ordering by date
ggplot(tmpB, 
       aes(x = Date, y = season_relab_bact, color=Month)) + 
  geom_jitter(alpha = 0.7) +
  ggtitle('Seasonal bacterial ASVs relative abundance by month') +
  # axis labels by month
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Extract Month from Date and convert to a properly ordered factor
tmpB$Month <- factor(format(tmpB$Date, "%B"), 
                      levels = month.name)

ggplot(tmpB, 
       aes(x = Month, y = season_relab_bact, color=Year)) + 
  geom_jitter(alpha = 0.7) +
  ggtitle('Seasonal bacterial ASVs relative abundance by month') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# viral seasonal
szn_vir_otu <- prune_taxa(taxa_names(vir_ps) %in% vir_szn$asv, vir_ps)
sam.sea_vir <- otu_table(szn_vir_otu) %>% colSums()
total.sea_vir <- otu_table(vir_ps) %>% colSums()

relabs_vir <- sam.sea_vir / total.sea_vir 

mean(relabs_vir)
sd(relabs_vir)

tmpV <-  data.frame(season_relab_vir = relabs_vir, 
           sample_data(vir_ps))

tmpV %>% head()
tmpV$Date <- as.Date(tmpV$Date)
# arrange by date
tmpV <- tmpV %>%
  arrange(Date)

tmpV %>% 
  ggplot( aes(season_relab_vir)) + 
  geom_histogram()

tmpV %>% head()

# Extract Month from Date and convert to a properly ordered factor
tmpV$Month <- factor(format(tmpV$Date, "%B"), 
                      levels = month.name)
# Extract year as a character
tmpV$Year <- format(tmpV$Date, "%Y")  

# fix ordering by date
ggplot(tmpV, 
       aes(x = Date, y = season_relab_bact, color=Month)) + 
  geom_jitter(alpha = 0.7) +
  ggtitle('Seasonal viral ASVs relative abundance by month') +
  # axis labels by month
#   scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# by Month coloured by Year
ggplot(tmpV, 
       aes(x = Month, y = season_relab_vir, color=Year)) + 
  geom_jitter(alpha = 0.7) +
  ggtitle('Seasonal viral ASVs relative abundance by month') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))








####################
library(microeco)
# https://chiliubio.github.io/microeco_tutorial/composition-based-class.html






# https://github.com/adriaaulaICM/bbmo_niche_sea/blob/6cef1b004e75a88a007975f6c5ebc37a40d32b0e/src/analysis/seasonality_taxranks.R#L125

# create df with ASV and taxonomic ranks (Kingdom and Family only)
rawdat <- as(tax_table(ps), "matrix") %>%  
  as_tibble(rownames = 'asv') %>% 
  filter(!is.na(ta5)) %>% select(-species) %>%
#   filter(class %in% c('Alphaproteobacteria',
#                       'Gammaproteobacteria',
#                       'Bacteroidia')) %>% 
  pivot_longer(names_to = 'rank', values_to = 'value', cols = -asv) %>%
  filter(rank != 'ta2', rank != 'ta3', rank != 'ta4', rank!='ta6') %>%
  filter(!is.na(value))

# filter for seasonal ASVs
seasonallevels <- rawdat %>% 
  filter(asv %in% lomb.sea.02$asv)

# check length
length(unique(rawdat$asv)) # 585
length(unique(lomb.sea.02$asv)) # 132
length(unique(seasonallevels$asv)) # 125


# ranks_more10 <-  rawdat  %>%
#   group_by(rank,value) %>%
#   filter(n() >= 10) %>%
#   filter(value %in% seasonallevels$value)  %>% 
#   distinct(rank,value)

# # note, increasing iteration (5) may be long
# distributions <- ranks_more10 %>% 
#   mutate( values = map2(.x = value, .y = rank,
#                         ~check_at_rank_level(ps, 5, .x, .y)))

# without20max <- ranks_more10 %>% 
#   mutate( values = map2(.x = value, .y = rank,
#                         ~check_situation_max(ps, .x, .y)))

# normal <- ranks_more10 %>% 
#   mutate( values = map2(.x = value, .y = rank,
#                         ~check_normal(ps, .x, .y)))
# # look at outputs
# distributions %>%
#     unnest(values)

# without20max %>%
#     unnest(values)

# normal %>%
#     unnest(values)







# # https://github.com/adriaaulaICM/bbmo_niche_sea/blob/6cef1b004e75a88a007975f6c5ebc37a40d32b0e/src/figures/cohesivity_taxrank.R#L48C1-L56C22
# dist <- distributions %>% 
#   ungroup() %>%
#   unnest(cols = values) #%>% 
# #     #Some changes in the nomenclature since we have ids repeated 
# #   mutate( value = case_when(
# #     rank == 'family' & value == 'D2472' ~ 'D2472_f',
# #     rank == 'genus' & value == 'D2472' ~ 'D2472_g',
# #     TRUE ~ value)) %>% 
# #   valtax_to_italics()




# tax <- as(tax_table(ps), 'matrix') %>% as_tibble(rownames = 'asv')

# corr <- tax %>% 
#   select(ta1, ta5)  %>% 
# #   filter(class %in% c('Alphaproteobacteria',
# #                       'Gammaproteobacteria',
# #                       'Bacteroidia')) %>% 
#   pivot_longer(names_to = 'level', values_to = 'value', cols = -ta5) %>% 
# #   mutate( value = case_when(
# #     level == 'family' & value == 'D2472' ~ 'D2472_f',
# #     level == 'genus' & value == 'D2472' ~ 'D2472_g',
# #     TRUE ~ value)) %>% 
#   select(-level)

# dist.hist <- distributions %>% 
#   unnest(values) %>%
#   mutate(rank = factor(rank, levels = c('ta1', 'ta5')),
#          significance = ifelse(pval < 0.01 & peak >= 0.154 & int.max <= 2,
#                                'seasonal', 'non seasonal')) %>% 
#   left_join(corr, by = 'value') %>% 
#   mutate(ta5 = ifelse(is.na(ta5), value, ta5))

# global <- dist.hist %>% 
#   ggplot( aes(peak)) + 
#   geom_histogram(aes(fill = rank)) +
#   facet_wrap(~significance, scales = 'free_y', ncol = 1) + 
#   theme_bw() + 
# #   lil.strip + 
# #   scale_fill_ptol() + 
#   xlab('Peak power of recurrence (strength of signal)') + 
#   ylab('Count of statistic') #+ 
# #   leg.bottom  
# global


# seas.label <- data.frame(label = c('', 'seasonal'),
#                          y = c('ta5', 'ta5'),
#                          x = c(7.3,16.9))
  
# global_res <- ggplot(dist.hist %>%
#                        mutate(ta5 = factor(ta5,
#                                              levels = c('k__Virus', 'f__Myoviridae'))),
#                      aes(x = peak,
#                          y = rank,
#                          fill = rank)) + 
#   geom_vline(xintercept = 10, linetype = 2, color ='tomato3') + 
# #   geom_density_ridges(show.legend = F) + 
#   geom_text(data = seas.label,
#             aes(label = label, x = x, y = y),
#             inherit.aes = FALSE,
#             size = 3,
#             nudge_y = -0.3, color = 'tomato3') + 
#   facet_wrap(~ta5) + 
#   scale_y_discrete(label = str_to_title) + 
# #   scale_fill_ptol() + 
# #   lil.strip + 
#   xlab('Peak power of recurrence (strength of signal)') + 
#   ylab('Rank level') 

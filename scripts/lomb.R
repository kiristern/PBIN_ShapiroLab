# source("lomb_fnc.R")

library(tidyverse)
library(phyloseq)
library(lubridate)
library(lomb)

# load seasonal data from nico
bact_szn <- read.csv("lomb_seasonality_bact.csv")
vir_szn <- read.csv("lomb_seasonality_virus.csv")

bact_szn %>% head()


# load from backup
lomb_env <- new.env() # Create a new environment
load("arxiv/Large_data/lomb backup.RData", envir = lomb_env) # Load the file into the new environment
ls(lomb_env)  # List the variables in the new environment

bact_ps <- lomb_env$bact3000filt
vir_ps <- lomb_env$virps3000filt





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

# # plotting
lil.strip <- theme(strip.background = element_blank(),
                   strip.text.x =element_text(margin = margin(.05, 0, .1, 0, "cm")))
# periodoplots <- asv_periods %>% 
#   map(~ggplot(.x, aes(scanned, power)) + 
#         geom_line(aes(group = asv)) + 
#         facet_wrap(~asv) + 
#         lil.strip
#         )
# periodoplots$ASV_30 # example plot for ASV_30





tax_table(bact_ps) %>% head()
tax_table(vir_ps) %>% head()


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
# asv.sel <- str_c('ASV_', c(1069, 1149))
asv.sel <- bact_szn$asv[1:5] # select first 5 ASVs for testing



# run for virps and bact_ps separately
ps <- bact_ps
# ps <- vir_ps.copy
ps <- transform_sample_counts(ps, function(x) x / sum(x))

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
# day of year col
ps_specific$day_of_year <- as.numeric(format(as.Date(ps_specific$Date), format = "%j"))

# # cumulative day numbers for the start of each month
num.days.mnt <- c(0,31,28,31,30,31,30,31,31,30,31,30)
cumnum <- cumsum(num.days.mnt)
print(cumnum)
# month order
date_order <- c('january', 'february', 'march', 'april', 'may', 'june', 'july', 'august', 'september', 'october', 'november', 'december')

# https://github.com/adriaaulaICM/bbmo_niche_sea/blob/6cef1b004e75a88a007975f6c5ebc37a40d32b0e/src/figures/sea_explanation.R#L45
gam.gg <- ggplot(data = ps_specific, aes(day_of_year,Abundance)) + 
  geom_jitter(aes(color = OTU), alpha = 0.4) + 
  stat_smooth(aes(x = day_of_year,
                  group = OTU,
                  color = OTU),
              method = "gam",
              formula = y ~ s(x, k =12, bs = 'cc'),
              se = F, size = 1,
              show.legend = F, alpha = 0.7) + 
  # split into separate facets (plots) for each OTU
  facet_wrap(~str_c(str_to_upper(OTU), Class, sep = ', '), scales = 'free_y') + 
  # display y-axis as percentage
  scale_y_continuous(labels = scales::percent_format(accuracy = 2L)) +
  scale_x_continuous(
                    breaks = cumnum,
                     # Display month abbrv first 3 letters
                     labels = str_to_title(date_order) %>% str_sub(1,3),
                     name = 'Month',
                     ) +
  guides(color = "none") + 
  ylab('Relative abundance') + 
  lil.strip + 
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels
  )
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
season_relab_vir <- taxa_sums(vir_ps)[bact_szn$asv] %>% sum(na.rm = TRUE)
total.relab_vir <- taxa_sums(vir_ps) %>% sum()
proportion_szn_of_tot_relab_vir <- season_relab_vir * 100 / total.relab_vir
print(proportion_szn_of_tot_relab_vir)


szn_bact_otu <- prune_taxa(taxa_names(bact_ps) %in% bact_szn$asv, bact_ps)
sam.sea_bact <- otu_table(szn_bact_otu) %>% colSums()
total.sea_bact <- otu_table(bact_ps) %>% colSums()

relabs_bact <- sam.sea_bact / total.sea_bact 

mean(relabs_bact)
sd(relabs_bact)

therel <-  data.frame(season_relab_bact = relabs_bact, 
           sample_data(bact_ps))

therel %>% 
  ggplot( aes(season_relab_bact)) + 
  geom_histogram()

therel %>% head()

ggplot(therel, 
       aes(Month, season_relab_bact)) + 
  geom_jitter(alpha = 0.7) +
  ggtitle('Seasonal bacterial ASVs relative abundance by month')






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

#https://github.com/adriaaulaICM/bbmo_niche_sea/blob/master/src/analysis/seasonality_asvs.R

library(tidyverse)
library(phyloseq)
library(lubridate)
library(lomb)

packageVersion("lomb")
citation("lomb")

bl.phy <- viral_physeq

# data_transformation <- function(physeq){
#   require(tidyverse)
#   require(phyloseq)
#   
#   
#   
#   # Relative abundance of OTUs
#   phy.relab = transform_sample_counts(physeq, function(x) x / sum(x))
#   
#   # Log10 abundance of OTUs
#   phy.log = transform_sample_counts(physeq, function(x) log10(x+1))
#   
#   # Rarefied dataset
#   phy.raref = rarefy_even_depth(physeq,rngseed = 42)
#   
#   # Log centered ratio (Aitchison)
#   geoMean = function(x, na.rm=TRUE){
#     exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
#   }
#   
#   phy.clr <- transform_sample_counts(physeq, function(x) x+1) %>%
#     transform_sample_counts(function(x) log(x/geoMean(x)))
#   
#   print('if this is not close to zero or zero, something is wrong! Beware')
#   print(sum(otu_table(phy.clr)[1,]))
#   
#   print("Starting a psmelt of a large dataset, this could take a while!")
#   
#   return(list( relab = phy.relab, 
#                log = phy.log,
#                raref = phy.raref,
#                clr = phy.clr))
# }
# 
# physeq.list <- data_transformation(bl.phy)
# 
# bl.phy.relab <- physeq.list$relab
# bl.phy.raref <- physeq.list$raref
# bl.phy.log <- physeq.list$log
# bl.phy.clr <- physeq.list$clr
# (virps_helli <- transform(viral_physeq, transform = "hellinger", target = "OTU"))
# bl.phy.asinh <- transform_sample_counts(bl.phy, function(x) asinh(x))
# physeq.list <- NULL


psmelt_dplyr = function(physeq) {
  #Implementation of the psmelt from phyloseq with dplyr
  # It is indeed faster, without further complications or differences. 
  # (well, in fact, its way more prone to give errors if the output is not well established, 
  # it doesn't check anything)
  
  sd = data.frame(sample_data(physeq)) %>% 
    rownames_to_column("Sample")
  TT = data.frame(as(tax_table(physeq),'matrix')) %>%
    rownames_to_column("OTU")
  if(taxa_are_rows(physeq)){
    otu.table = data.frame(as(otu_table(physeq),"matrix"),
                           check.names = FALSE) %>%
      rownames_to_column("OTU")
  } else {
    otu.table = data.frame(t(as(otu_table(physeq),"matrix")),
                           check.names = FALSE) %>%
      rownames_to_column("OTU")
  }
  
  all <- otu.table %>%
    left_join(TT, by = "OTU") %>%
    gather_("Sample", "Abundance", setdiff(colnames(otu.table), "OTU")) %>%
    left_join(sd, by = "Sample") %>% 
    select(Sample, Abundance, OTU, everything())
  
  return(all)
}

tsdf.0.2 <- bl.phy %>% 
  psmelt_dplyr() %>% 
  mutate(decimaldat = decimal_date(Date)) 
head(tsdf.0.2)
tail(tsdf.0.2)


#troubleshooting why lomb was not working
tsdf.0.2$OTU %>% head()
tsdf.0.2 %>% select(OTU, decimaldat, Abundance) %>%
  group_by(OTU) %>%
  summarize( howmany0 = sum(Abundance == 0)) %>%
  arrange(-howmany0)
# ~4000+ ASVs that have abundance of 0 in 166 samples

ASVs0 <- tsdf.0.2 %>%
  select(OTU, decimaldat, Abundance) %>%
  group_by(OTU) %>%
  summarize( howmany0 = sum(Abundance == 0)) %>%
  arrange(-howmany0) %>%
  filter(howmany0 == 166) %>%
  pull(OTU)


lomb.02 <- tsdf.0.2 %>% 
  filter(! OTU %in% ASVs0) %>% 
  split(.$OTU) %>% 
  map(~randlsp( x =.x$Abundance,
                times = .x$decimaldat,
                type = 'period',
                plot = F))

par(mar=c(1,1,1,1)) #fix error: Error in plot.new() : figure margins too large

lomb.sea.02 <- tibble( asv = names(lomb.02),
                       pval = map_dbl(lomb.02, ~.x$p.value),
                       peak = map_dbl(lomb.02, ~.x$peak),
                       interval  = map(lomb.02, ~.x$peak.at),
                       int.min = map_dbl(interval, ~.[[2]]),
                       int.max = map_dbl(interval, ~.[[1]])) %>% 
  mutate( qval = fdrtool::fdrtool(pval, statistic = 'pvalue')$qval) %>% 
  filter(qval <= 0.01, peak >= 10, int.max <= 2)
head(lomb.sea.02)
head(lomb.season <- as.data.frame(lomb.sea.02[,c(1, 2,7)]))
write.csv(lomb.season, "lomb_seasonality.csv")

write_rds(lomb.02, 'data/analysis/lomball.rds')
write_rds(lomb.sea.02, 'data/analysis/lombsea.rds')

results.lomb02 <- lomb.02[lomb.sea.02 %>% pull(asv)]

# Strip problems
lil.strip <- theme(strip.background = element_blank(),
                   strip.text.x =element_text(margin = margin(.05, 0, .1, 0, "cm")))

periodoplots <- map(results.lomb02, ~tibble( scanned = .x$scanned,
                                             power = .x$power)) %>% 
  bind_rows(.id = 'asv') %>% 
  split(.$asv) %>% 
  map(~ggplot(.x, aes(scanned, power)) + 
        geom_line(aes(group = asv)) + 
        facet_wrap(~asv) + 
        lil.strip)

periodoplots$ASV_302




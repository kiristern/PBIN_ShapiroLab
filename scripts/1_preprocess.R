packageVersion("Phyloseq")
citation("phyloseq")


##################################

library(tidyverse)
library(stringr)
library(vegan)
library(lubridate)


#### UPLOAD DATA ####
#upload viral ASV count table and vir_metaa
ASV_count <- read.table("data/ASVs_counts_copy.tsv", row.names = 1, header=T)
str(ASV_count)
dim(ASV_count)
range(ASV_count)
apply(ASV_count, 2, median) #get median for each col
median(unlist(ASV_count), na.rm = T) #et median for whole df
summary(ASV_count)

colnames(ASV_count)[colnames(ASV_count) == "FLD0295_15_05_2011_1"] <- "FLD0295_15_05_2011_2" #dates were duplicated therefore need to correct

head(ASV_count, n=2)

source("scripts/metadata.R")
head(vir_meta, n=2)
vir_meta$Years <- as.factor(vir_meta$Years)
vir_meta$Date <- as.Date(vir_meta$Date)
vir_meta$Day <- format(as.Date(vir_meta$Date, format="%Y-%m-%d"), "%d")
vir_meta$month.numeric <- format(as.Date(vir_meta$Date, format="%Y-%m-%d"), "%m")
vir_meta$week <- lubridate::week(ymd(vir_meta$Date))
str(vir_meta)

#order meta by date
vir_meta <- vir_meta[order(vir_meta$Date),]

#get basic vir_meta data info for methods section of report
#check to make sure ASV data and vir_metaa dates are the same
if (length(ASV_count) == nrow(vir_meta)){
  print("ASV and vir_metaa have same sampling date lengths -- GREAT")
  } else {
    print("ASV and vir_metaa have different sampling date lengths -- NOT GREAT")
}


(amt <- vir_meta %>% group_by(Years) %>% summarise(amount=length(Years))) #view how many samples per year
min(amt[,2])
max(amt[,2])
median(as.numeric(unlist(amt[,2]))) #get median samples per year
vir_meta %>% group_by(Site) %>% summarise(amount=length(Site))

# #ensure same samples between ASV_count and meta
# names(ASV_count) %in% rownames(vir_meta)
# #find which one is not the same
# which(!names(ASV_count) %in% rownames(vir_meta))
# names(ASV_count[100])

#fix error above (commented out)
rownames(vir_meta)[rownames(vir_meta) == "FLD0295_15_05_2011_1"] <- "FLD0295_15_05_2011_2" #dates were duplicated therefore need to correct
vir_meta$description[rownames(vir_meta) == 'FLD0295_15_05_2011_2'] <- '15.05.2011.2' #update description to proper date

# asv_count <- ASV_count[,(colnames(ASV_count) %in% rownames(vir_meta))]
# length(asv_count)
# nrow(vir_meta)

asv_count <- ASV_count



#### CREATE PHYLOSEQ OBJECT ####
library(phyloseq)
library(microbiome)
#VIRAL

#add ASV count table, vir_metaa, virTree to phyloseq table
nrow(asv_count)
nozero <- asv_count[rownames(asv_count) %in% names(rowSums(asv_count > 0)),] #only keep ASVs that are present in data
nrow(nozero)
length(rowSums(asv_count > 0)) == nrow(asv_count)

count_phy <- otu_table(asv_count, taxa_are_rows=T)
dim(count_phy)
sample_info <- sample_data(vir_meta)
dim(sample_info)
#virTree <- read_tree("data/viral_tree")

fake_taxa <- read.table("data/fake_viral_tax.txt", header = T, row.names = 1, fill=T)
mock_taxa <- tax_table(fake_taxa)
mock_taxa[,7] <- stringr::str_remove(mock_taxa[,7], "s__")
row.names(mock_taxa) <- mock_taxa[,7]
colnames(mock_taxa)[7] <- "species"
head(mock_taxa)
dim(mock_taxa)
#add to phyloseq object
viral_physeq <- phyloseq(count_phy, sample_info, mock_taxa) #, virTree) 
viral_physeq %>% otu_table( ) %>% dim

print("it appears the viral phylogenetic tree removes some ASVs. Consider removing virTree from ps object")

# rm reads fewer than 3000 -- filtering samples
virps3000 <- prune_samples(sample_sums(viral_physeq)>=3000, viral_physeq)
virps3000

# filter taxa not seen more than once in 10% of samples -- filtering taxa
virps3000filt <- filter_taxa(virps3000, function(x) sum(x > 1) > (0.10*length(x)), TRUE)
virps3000filt

virfiltotu <- virps3000filt %>% otu_table()
#write.csv(virfiltotu, "viralCounts_filtered.csv")

#Data for Dawson
# #which rownames are different
# virfiltotu
# 
# viral_relab_ps <- transform(viral_physeq, "compositional", target="OTU")
# viral_relab_otu <- viral_relab_ps %>% otu_table()
# 
# viralfilt_relab_otu <- viral_relab_otu[rownames(virfiltotu),]
# dim(viralfilt_relab_otu)
# dim(virfiltotu)
# #write.csv(viralfilt_relab_otu, "viralfilt_relab.csv")

#which(taxa_sums(virps3000filt) == 0)

#asv_count
# asv_filt <- virps3000filt %>% otu_table()
#write.table(asv_filt, "asv_filt.tsv")




### Bacterial preprocessing -- for SpiecEasi ###
source("scripts/bacteria_analysis.R")
#subset by taxa
cyano_ps <- subset_taxa(bact_physeq, Phylum == "p__Cyanobacteria")
doli_ps <- subset_taxa(bact_physeq, Genus == "g__Dolichospermum")
micro_ps <- subset_taxa(bact_physeq, Genus == "g__Microcystis")

cyano_ps_relab <- subset_taxa(bact_relab_ps, Phylum == "p__Cyanobacteria")
doli_ps_relab <- subset_taxa(bact_relab_ps, Genus == "g__Dolichospermum")
micro_ps_relab <- subset_taxa(bact_relab_ps, Genus == "g__Microcystis")

cyano_ps_helli <- subset_taxa(bact_helli_ps, Phylum == "p__Cyanobacteria")
doli_ps_helli <- subset_taxa(bact_helli_ps, Genus == "g__Dolichospermum")
micro_ps_helli <- subset_taxa(bact_helli_ps, Genus == "g__Microcystis")

# cyano_ps_clr <- subset_taxa(bact_clr_ps, Phylum == "p__Cyanobacteria")
# doli_ps_clr <- subset_taxa(bact_clr_ps, Genus == "g__Dolichospermum")
# micro_ps_clr <- subset_taxa(bact_clr_ps, Genus == "g__Microcystis")

#replace name so don't have to edit whole script
#cyano_ps <- micro_ps

# #ensure viral ps has same samples as cyano_ps 
# virbact_meta

#put bacterial counts per sample in vir_metaa
virbact_meta$micro.sum <- micro_ps %>% otu_table() %>% colSums()
virbact_meta$doli.sum <- doli_ps %>% otu_table() %>% colSums()
virbact_meta$cyano.sum <- cyano_ps %>% otu_table() %>% colSums()

virbact_meta$micro.relab.sum <- micro_ps_relab %>% otu_table() %>% colSums()
virbact_meta$doli.relab.sum <- doli_ps_relab %>% otu_table() %>% colSums()
virbact_meta$cyano.relab.sum <- cyano_ps_relab %>% otu_table() %>% colSums()

virbact_meta$micro.helli.sum <- micro_ps_helli %>% otu_table() %>% colSums()
virbact_meta$doli.helli.sum <- doli_ps_helli %>% otu_table() %>% colSums()
virbact_meta$cyano.helli.sum <- cyano_ps_helli %>% otu_table() %>% colSums()

# virbact_meta$micro.clr.sum <- micro_ps_clr %>% otu_table() %>% colSums()
# virbact_meta$doli.clr.sum <- doli_ps_clr %>% otu_table() %>% colSums()
# virbact_meta$cyano.clr.sum <- cyano_ps_clr %>% otu_table() %>% colSums()

head(virbact_meta)
dim(meta_all)
#write.csv(virbact_meta, "~/Desktop/meta_w_dolimicroSums.csv")




virps3000_samemeta <- virps3000filt

sample_data(virps3000_samemeta) <- sample_data(virps3000filt)[get_variable(virps3000filt, "description") %in% virbact_meta$description]

sample_names(virps3000_samemeta) <- sample_data(virps3000_samemeta)$description
virps3000_samemeta %>% sample_data() %>% head()


#taxa_names(viral_physeq) <- paste0("vir_", taxa_names(viral_physeq))
taxa_names(doli_ps) <- paste0("doli_", taxa_names(doli_ps))
taxa_names(micro_ps) <- paste0("micro_", taxa_names(micro_ps))
taxa_names(virps3000_samemeta) <- paste0("vir_", taxa_names(virps3000_samemeta))


doli.ps <- prune_samples(rownames(sample_data(doli_ps)) %in% rownames(sample_data(virps3000_samemeta)), doli_ps)
micro.ps <- prune_samples(rownames(sample_data(micro_ps)) %in% rownames(sample_data(virps3000_samemeta)), micro_ps)

#reorder phyloseq by chronological date
map <- sample_data(virps3000_samemeta)[order(sample_data(virps3000_samemeta)$Date),]
toorder <- rownames(map)

otu_table(virps3000_samemeta) <- otu_table(virps3000_samemeta)[,toorder]
otu_table(doli.ps) <- otu_table(doli.ps)[,toorder]
otu_table(micro.ps) <- otu_table(micro.ps)[,toorder]
otu_table(bact_physeq) <- otu_table(bact_physeq)[,toorder]
otu_table(cyano_ps) <- otu_table(cyano_ps)[,toorder]

bactnoCyan <- subset_taxa(bact_physeq, !Phylum == "p__Cyanobacteria")
bactnoCyan_filt <- filter_taxa(bactnoCyan, function(x) sum(x > 1) > (0.10*length(x)), TRUE)

virps_filt <- filter_taxa(virps3000_samemeta, function(x) sum(x > 1) > (0.10*length(x)), TRUE)

cyano.ps_filt <- filter_taxa(cyano_ps, function(x) sum(x > 1) > (0.10*length(x)), TRUE)
# doli_filt <- filter_taxa(doli_ps, function(x) sum(x > 1) > (0.10*length(x)), TRUE)
# micro_filt <- filter_taxa(micro_ps, function(x) sum(x > 1) > (0.10*length(x)), TRUE)


virps_filt
cyano.ps_filt 





packageVersion("Phyloseq")
citation("phyloseq")


##################################

library(tidyverse)
library(stringr)
library(vegan)
library(lubridate)
source("scripts/functions.R")

##### GET raw data #####


#### UPLOAD DATA ####
#upload viral ASV count table and metadataa
ASV_count <- read.table("data/ASVs_counts_copy.tsv", row.names = 1, header=T)
str(ASV_count)
dim(ASV_count)
range(ASV_count)
apply(ASV_count, 2, median) #get median for each col
median(unlist(ASV_count), na.rm = T) #et median for whole df
summary(ASV_count)

colnames(ASV_count)[colnames(ASV_count) == "FLD0295_15_05_2011_1"] <- "FLD0295_15_05_2011_2" #dates were duplicated therefore need to correct

head(ASV_count, n=2)

# source("scripts/metadata.R")
# head(metadata, n=2)
# metadata$Years <- as.factor(metadata$Years)
# metadata$Date <- as.Date(metadata$Date)
# metadata$Day <- format(as.Date(metadata$Date, format="%Y-%m-%d"), "%d")
# metadata$month.numeric <- format(as.Date(metadata$Date, format="%Y-%m-%d"), "%m")
# metadata$week <- lubridate::week(ymd(metadata$Date))
# str(metadata)

metadata <- read.csv('data/PBIN_metadata - METAFINAL.csv')

#order meta by date
metadata <- metadata[order(metadata$Date),]

#get basic metadata data info for methods section of report
(amt <- metadata %>% group_by(Years) %>% summarise(amount=length(Years))) #view how many samples per year
min(amt[,2])
max(amt[,2])
median(as.numeric(unlist(amt[,2]))) #get median samples per year
metadata %>% group_by(Site) %>% summarise(amount=length(Site))


#check to make sure ASV data and metadataa dates are the same
if (length(ASV_count) == nrow(metadata)){
  print("ASV and metadata have same sampling date lengths -- GREAT")
} else {
  print("ASV and metadataa have different sampling date lengths -- NOT GREAT")
}

# see what is not the same between ASV samples and metadata
names(ASV_count) %in% metadata$sampleID
#find which one is not the same
names(ASV_count[which(!names(ASV_count) %in% metadata$sampleID)])

## another way to do above:
# metadata[!(metadata$sampleID %in% colnames(ASV_count)),] # all ASV samples in metadata
# head(ASV_count[,!(colnames(ASV_count) %in% metadata$sampleID)]) # FLD0235_02_06_2007 FLD0247_08_09_2008_3d FLD0295_15_05_2011_2 missing from metadata


#fix error above (commented out)
# rownames(metadata)[rownames(metadata) == "FLD0295_15_05_2011_1"] <- "FLD0295_15_05_2011_2" #dates were duplicated therefore need to correct
metadata$description[metadata$sampleID == 'FLD0295_15_05_2011_2'] <- '15.05.2011.2' #update description to proper date

asv_count <- ASV_count[,(colnames(ASV_count) %in% metadata$sampleID)]
length(asv_count)
nrow(metadata)



#### CREATE PHYLOSEQ OBJECT ####
library(phyloseq)
library(microbiome)
#VIRAL

#add ASV count table, metadataa, virTree to phyloseq table
nrow(asv_count)
nozero <- asv_count[rownames(asv_count) %in% names(rowSums(asv_count > 0)),] #only keep ASVs that are present in data
nrow(nozero)
length(rowSums(asv_count > 0)) == nrow(asv_count)

count_phy <- otu_table(asv_count, taxa_are_rows=T)
dim(count_phy)
rownames(metadata) <- metadata$sampleID
sample_info <- sample_data(metadata)
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

vp <- viral_physeq %>% otu_table( )
write.csv(vp, "data/raw data/viral_phyloseq.csv")

print("it appears the viral phylogenetic tree removes some ASVs. Consider removing virTree from ps object")

# rm reads fewer than 3000 -- filtering samples
virps3000 <- prune_samples(sample_sums(viral_physeq)>=3000, viral_physeq)

vp3000 <- virps3000 %>% otu_table()
write.csv(vp3000, 'data/raw data/viral_physeq_3000.csv')

# filter taxa not seen more than once in 10% of samples -- filtering taxa
virps3000filt <- filter_taxa(virps3000, function(x) sum(x > 1) > (0.10*length(x)), TRUE)

vpfilt3000 <- virps3000filt %>% otu_table()
write.csv(vpfilt3000, 'data/raw data/viral_physeq_filt3000.csv')


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


bact_physeq %>% sample_data() %>% head()
head(virbact_meta)
meta_all <- metadata
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
bnoC <- bactnoCyan %>% otu_table()
write.csv(bnoC, 'data/raw data/bact_count_noCyano.csv')

bactnoCyan_filt <- filter_taxa(bactnoCyan, function(x) sum(x > 1) > (0.10*length(x)), TRUE)
bnoC_filt <- bactnoCyan_filt %>% otu_table()
write.csv(bnoC_filt, 'data/raw data/bact_count_noCyanofilt.csv')

virps_filt <- filter_taxa(virps3000_samemeta, function(x) sum(x > 1) > (0.10*length(x)), TRUE)

cyano.ps_filt <- filter_taxa(cyano_ps, function(x) sum(x > 1) > (0.10*length(x)), TRUE)
# doli_filt <- filter_taxa(doli_ps, function(x) sum(x > 1) > (0.10*length(x)), TRUE)
# micro_filt <- filter_taxa(micro_ps, function(x) sum(x > 1) > (0.10*length(x)), TRUE)


virps_filt %>% sample_data() %>% head()
cyano.ps_filt 





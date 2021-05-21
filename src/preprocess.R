packageVersion("Phyloseq")
citation("phyloseq")


##################################

library(dplyr)
library(stringr)
library(vegan)


#### UPLOAD DATA ####
#upload viral ASV count table and metadata
ASV_count <- read.table("data/ASVs_counts_copy.tsv", row.names = 1, header=T)
str(ASV_count)
dim(ASV_count)
range(ASV_count)
apply(ASV_count, 2, median) #get median for each col
median(unlist(ASV_count), na.rm = T) #et median for whole df
summary(ASV_count)

colnames(ASV_count)[colnames(ASV_count) == "FLD0295_15_05_2011_1"] <- "FLD0295_15_05_2011_2" #dates were duplicated therefore need to correct
head(ASV_count, n=2)
metadat <- read.csv("data/meta_cmd.csv", row.names = 1, header = T)
head(metadat, n=2)
#metadat$Years <- as.factor(metadat$Years)
metadat$Years <- as.factor(metadat$Years)
metadat$Date <- as.Date(metadat$Date)
str(metadat)

#order meta by date
metadat <- metadat[order(metadat$Date),]

#get basic metadat data info for methods section of report
length(ASV_count)
nrow(metadat)
(amt <- metadat %>% group_by(Years) %>% summarise(amount=length(Years))) #view how many samples per year
min(amt[,2])
max(amt[,2])
median(as.numeric(unlist(amt[,2]))) #get median samples per year
metadat %>% group_by(Site) %>% summarise(amount=length(Site))

#ensure same samples between ASV_count and meta
asv_count<- ASV_count[,(colnames(ASV_count) %in% rownames(metadat))]
length(asv_count)
nrow(metadat)


#### CREATE PHYLOSEQ OBJECT ####
library(phyloseq)
library(microbiome)
#VIRAL

#add ASV count table, metadata, virTree to phyloseq table
nozero <- asv_count[rownames(asv_count) %in% names(rowSums(asv_count > 0)),]
length(rowSums(asv_count > 0)) == nrow(asv_count)
count_phy <- otu_table(asv_count, taxa_are_rows=T)
sample_info <- sample_data(metadat)
virTree <- read_tree("data/viral_tree")

fake_taxa <- read.table("data/fake_viral_tax.txt", header = T, row.names = 1, fill=T)
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

# filter taxa not seen more than once in 10% of samples
virps3000filt <- filter_taxa(virps3000, function(x) sum(x > 1) > (0.10*length(x)), TRUE)

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
asv_filt <- virps3000filt %>% otu_table()
#write.table(asv_filt, "asv_filt.tsv")

vir_abun_filt <- virps3000 %>% otu_table()



### SpiecEasi ###
source("src/bacteria_analysis.R")
#subset by taxa
cyano_ps <- subset_taxa(bact_physeq, Phylum == "p__Cyanobacteria")
doli_ps <- subset_taxa(bact_physeq, Genus == "g__Dolichospermum")
micro_ps <- subset_taxa(bact_physeq, Genus == "g__Microcystis")

#replace name so don't have to edit whole script
#cyano_ps <- micro_ps

#ensure viral ps has same samples as cyano_ps 
meta2

virps3000_samemeta <- virps3000filt

sample_data(virps3000_samemeta) <- sample_data(virps3000filt)[get_variable(virps3000filt, "description") %in% meta2$description]

sample_names(virps3000_samemeta) <- sample_data(virps3000_samemeta)$description

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

cyanops_filt <- filter_taxa(cyano_ps, function(x) sum(x > 1) > (0.10*length(x)), TRUE)
# doli_filt <- filter_taxa(doli_ps, function(x) sum(x > 1) > (0.10*length(x)), TRUE)
# micro_filt <- filter_taxa(micro_ps, function(x) sum(x > 1) > (0.10*length(x)), TRUE)


virps_filt
cyanops_filt 





packageVersion("Phyloseq")
citation("phyloseq")


##################################

library(dplyr)
library(stringr)
library(vegan)


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
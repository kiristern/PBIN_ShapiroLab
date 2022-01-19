#upload cyano ASV data
cyano_counts <- read.table("data/Champ_ASVs_counts.txt", header = TRUE, row.names = 1)
head(cyano_counts)
cyano_taxa <- read.csv("data/ASVs_taxonomy_Champ_Greengenes.csv", header = T, row.names = 1, fill=T)
head(cyano_taxa)

nrow(metadat)
length(cyano_counts)

cyano_taxa$Species <- rownames(cyano_taxa)

colnames(cyano_counts)
#remove X at beginning of date
colnames(cyano_counts)[1:135] <- substring(colnames(cyano_counts)[1:135], 2)
#select dates cols only
cyano_counts <- cyano_counts[1:135]

#match sample dates
meta2 <- metadat

#make sure meta matches cyano samples
nrow(meta2)
length(cyano_counts)

# #relative abundance before rm samples to match meta
# bac_count <- otu_table(cyano_counts, taxa_are_rows = T)
# bact_physeq1 <- phyloseq(bac_count)
# bact_relab <- transform_sample_counts(bact_physeq1, function(x) x / sum(x))
# bact_relab_otu <- bact_relab %>% otu_table()
# 
# cyano_counts <- bact_relab_otu
# go to line 63

#select cols that match dates
bact_counts <- cyano_counts[,(colnames(cyano_counts) %in% meta2$description)]
length(bact_counts)

meta2 <- meta2[meta2$description %in% colnames(bact_counts),]
nrow(meta2)
rownames(meta2) <- meta2$description #need rownames to match colnames of bact_count for phyloseq object

# which(duplicated(meta2$description)) #fixed in preprocess.R
# meta2$description[90]

#Phyloseq
library(phyloseq)

bac_count <- otu_table(bact_counts, taxa_are_rows = T)
dim(bac_count)
bact_tax_tab <- tax_table(cyano_taxa)
rownames(bact_tax_tab) <- rownames(cyano_taxa)

#add to phyloseq object
sample_info_cyano <- sample_data(meta2)
dim(sample_info_cyano)

bact_physeq <- phyloseq(bac_count, bact_tax_tab, sample_info_cyano)
print(bact_physeq)

#rename cols
colnames(tax_table(bact_physeq)) <- c("Kingdom", "Phylum", "Class",
                                      "Order", "Family", "Genus", "ASV")

# #relab of doli and micro
# doli_ps_relab <- subset_taxa(bact_physeq, Genus == "g__Dolichospermum")
# micro_ps_relab <- subset_taxa(bact_physeq, Genus == "g__Microcystis")
# 
# samp_Data <- doli_ps_relab %>% sample_data()
# doli_relab <- doli_ps_relab %>% otu_table()
# samp_Data$doli.relab <- colSums(doli_relab)
# micro_relab <- micro_ps_relab %>% otu_table()
# samp_Data$micro.relab <- colSums(micro_relab)
# 
# write.csv(samp_Data, "~/Desktop/meta_relab_dolimicro.csv")


#quick check
bact_physeq %>% tax_table %>% head()
bact_taxatab <- bact_physeq %>% tax_table
bact_taxatab[82,]

#reorder phyloseq by chronological date
all(colnames(bact_counts) %in% rownames(meta2))
map <- meta2[order(meta2$Date),]
toorder <- rownames(map)
otu_table(bact_physeq) <- otu_table(bact_physeq)[,toorder]
bact_physeq %>% otu_table()


# rm reads less than 3000
bact3000 <- prune_samples(sample_sums(bact_physeq)>=3000, bact_physeq)
print("doesn't appear to have removed anything. Try with a greater filter?")

# filter taxa not seen more than once in 10% of samples
bact3000filt <- filter_taxa(bact3000, function(x) sum(x > 1) > (0.10*length(x)), TRUE)

bactfiltotu <- bact3000filt %>% otu_table()
#write.csv(bactfiltotu, "bacteriaCounts_filtered.csv")


bactotutab <- bact_physeq %>% otu_table()

bact_relab_ps <- transform(bact_physeq, "compositional", target="OTU") #'compositional' abundances are returned as relative abundances in [0, 1] (convert to percentages by multiplying with a factor of 100).
bact_helli_ps <- transform(bact_physeq, "hellinger", target="OTU") #Hellinger transform is square root of the relative abundance but instead given at the scale [0,1]
bact_clr_ps <- transform(bact_physeq, "clr", target="OTU") #CLR transform applies a pseudocount of min(relative abundance)/2 to exact zero relative abundance entries in OTU table before taking logs.


bact_relab_otu <- bact_relab_ps %>% otu_table()

bactfilt_relab_otu <- bact_relab_otu[which(rownames(bactfiltotu) %in% rownames(bactotutab)),]
dim(bactfilt_relab_otu)
dim(bactfiltotu)
#write.csv(bactfilt_relab_otu, "bactfilt_clr.csv")



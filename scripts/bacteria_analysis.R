#upload cyano ASV data
bact_counts <- read.table("data/Champ_ASVs_counts.txt", header = TRUE, row.names = 1)
head(bact_counts)
cyano_taxa <- read.csv("data/ASVs_taxonomy_Champ_Greengenes.csv", header = T, row.names = 1, fill=T)
head(cyano_taxa)

cyano_taxa$Species <- rownames(cyano_taxa)

colnames(bact_counts)
#remove X at beginning of date
colnames(bact_counts)[1:135] <- substring(colnames(bact_counts)[1:135], 2)
#select dates cols only
bact_counts <- bact_counts[1:135]

# get all samples 
bact_meta <- as.data.frame(colnames(bact_counts))
colnames(bact_meta) <- 'sampleID'
bact_meta$Date <- sub('^([^.]+.[^.]+.[^.]+).*', '\\1', bact_meta$sampleID) #keep only dates (rm everything after third .)
bact_meta$Date <- as.Date(format(dmy(bact_meta$Date), '%Y-%m-%d'))

#format dates
bact_meta$Month <- format_month(bact_meta)
bact_meta$Years <- format(as.Date(bact_meta$Date, format='%Y-%m-%d'), "%Y") #get year

bact_meta$Period <- getSeason(bact_meta$Date)

#add new columns to meta df 
for (i in 1:length(bact_meta$Date)){
  bact_meta$Mean_temp_t0_t7[i] <- get_mean_temp(bact_meta$Date[i])
}

#get mean prec for each sample
for (i in 1:length(bact_meta$Date)){
  bact_meta$Cumul_precip_t1_t7_mm[i] <- get_cumul_prec(bact_meta$Date[i])
}

head(bact_meta)



#### Match bacteria metadata dimensions to viral metadata dimensions ####
#match sample dates
virbact_meta <- vir_meta

#make sure meta matches cyano samples
nrow(virbact_meta)
length(bact_counts)

# #relative abundance before rm samples to match meta
# bac_count <- otu_table(bact_counts, taxa_are_rows = T)
# bact_physeq1 <- phyloseq(bac_count)
# bact_relab <- transform_sample_counts(bact_physeq1, function(x) x / sum(x))
# bact_relab_otu <- bact_relab %>% otu_table()
# 
# bact_counts <- bact_relab_otu
# go to line 63

#select cols that match dates
# bact_counts <- bact_counts[,( virbact_meta$description %in% colnames(bact_counts))]
# length(bact_counts)

virbact_meta <- virbact_meta[virbact_meta$description %in% colnames(bact_counts),]
nrow(virbact_meta)
rownames(virbact_meta) <- virbact_meta$description #need rownames to match colnames of bact_count for phyloseq object

# which(duplicated(virbact_meta$description)) #fixed in preprocess.R
# virbact_meta$description[90]

#Phyloseq
library(phyloseq)

bac_count <- otu_table(bact_counts, taxa_are_rows = T)
dim(bac_count)
bact_tax_tab <- tax_table(cyano_taxa)
rownames(bact_tax_tab) <- rownames(cyano_taxa)

#add to phyloseq object
sample_info_cyano <- sample_data(virbact_meta)
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
all(colnames(bact_counts) %in% rownames(virbact_meta))
map <- virbact_meta[order(virbact_meta$Date),]
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

bact_relab_ps <- microbiome::transform(bact_physeq, "compositional", target="OTU") #'compositional' abundances are returned as relative abundances in [0, 1] (convert to percentages by multiplying with a factor of 100).
bact_helli_ps <-microbiome::transform(bact_physeq, "hellinger", target="OTU") #Hellinger transform is square root of the relative abundance but instead given at the scale [0,1]
bact_clr_ps <- microbiome::transform(bact_physeq, "clr", target="OTU") #CLR transform applies a pseudocount of min(relative abundance)/2 to exact zero relative abundance entries in OTU table before taking logs.


bact_relab_otu <- bact_relab_ps %>% otu_table()

bactfilt_relab_otu <- bact_relab_otu[which(rownames(bactfiltotu) %in% rownames(bactotutab)),]
dim(bactfilt_relab_otu)
dim(bactfiltotu)
#write.csv(bactfilt_relab_otu, "bactfilt_clr.csv")



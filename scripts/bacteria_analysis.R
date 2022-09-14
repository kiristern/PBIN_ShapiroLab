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
bact_meta$Month <- lubridate::month(ymd(bact_meta$Date), label=T)
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
virbact_meta <- metadata

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

#find which one is not the same
names(bact_counts[which(!colnames(bact_counts) %in% metadata$description)])

#select cols that match dates
bact_counts <- bact_counts[,( virbact_meta$description %in% colnames(bact_counts))]
length(bact_counts)

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

bp <- bact_physeq %>% otu_table()
write.csv(bp, 'data/bact_phyloseq.csv')

# rm reads less than 3000
bact3000 <- prune_samples(sample_sums(bact_physeq)>=3000, bact_physeq)
print("doesn't appear to have removed anything. Try with a greater filter?")

bp3000 <- bact3000 %>% otu_table()
write.csv(bp3000, 'data/bact3000.csv')

# filter taxa not seen more than once in 10% of samples
bact3000filt <- filter_taxa(bact3000, function(x) sum(x > 1) > (0.10*length(x)), TRUE)

bactfiltotu <- bact3000filt %>% otu_table()
write.csv(bactfiltotu, "data/bact3000filt.csv")


bactotutab <- bact_physeq %>% otu_table()

bact_relab_ps <- microbiome::transform(bact_physeq, "compositional", target="OTU") #'compositional' abundances are returned as relative abundances in [0, 1] (convert to percentages by multiplying with a factor of 100).
bact_helli_ps <-microbiome::transform(bact_physeq, "hellinger", target="OTU") #Hellinger transform is square root of the relative abundance but instead given at the scale [0,1]
bact_clr_ps <- microbiome::transform(bact_physeq, "clr", target="OTU") #CLR transform applies a pseudocount of min(relative abundance)/2 to exact zero relative abundance entries in OTU table before taking logs.

bact_relab_otu <- bact_relab_ps %>% otu_table()

bactfilt_relab_otu <- bact_relab_otu[which(rownames(bactfiltotu) %in% rownames(bactotutab)),]
dim(bactfilt_relab_otu)
dim(bactfiltotu)


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

write.csv(virbact_meta, 'data/bactmeta.csv')



###### ADD BACTERIAL META INFO INTO META TABLE #####
meta_virbact_FINAL <- metadata
meta_virbact_FINAL <- merge(metadata, virbact_meta[c('description', 'Day', 'month.numeric', 'week', 'micro.sum', 'doli.sum', 'cyano.sum',
                                                     'micro.relab.sum', 'doli.relab.sum', 'cyano.relab.sum', 'micro.helli.sum',
                                                     'doli.helli.sum', 'cyano.helli.sum')], 
                            # by.x = c('description', 'Date', 'Years', 'Month',' Period', 'Mean_temp_t0_t7', 
                            #          'Cumul_precip_t1_t7_mm','Site', 'Bloom', 'Tot.P_ug', 'Tot.N_mg', 'Dissolved.P',
                            #           'Dissolved.N', 'Microcystin','P_range', 'N_range','Temp.water', 
                            #          'profondeur_secchi_cm', 'microcystin_ug_L'),
                            # by.y=c('description', 'Date', 'Years', 'Month', 'Period',
                            #        'Mean_temp_t0_t7', 'Cumul_precip_t1_t7_mm', 'Site', 'Bloom', 'Tot.P_ug',
                            #        'Tot.N_mg', 'Dissolved.P', 'Dissolved.N', 'Microcystin', 'P_range', 'N_range',
                            #        'Temp.water', 'profondeur_secchi_cm', 'microcystin_ug_L'), 
                            by='description',
                            all.x=T)
write.csv(meta_virbact_FINAL, 'data/meta_w_bact-abund.csv')
# virbact_meta$micro.clr.sum <- micro_ps_clr %>% otu_table() %>% colSums()
# virbact_meta$doli.clr.sum <- doli_ps_clr %>% otu_table() %>% colSums()
# virbact_meta$cyano.clr.sum <- cyano_ps_clr %>% otu_table() %>% colSums()


#update phyloseq object with new virbact_meta data
sample_info_cyano <- sample_data(virbact_meta)
bact_physeq <- phyloseq(bac_count, bact_tax_tab, sample_info_cyano)

#rename cols
colnames(tax_table(bact_physeq)) <- c("Kingdom", "Phylum", "Class",
                                      "Order", "Family", "Genus", "ASV")

bact_physeq %>% sample_data() %>% head()

 #write.csv(bactfilt_relab_otu, "bactfilt_clr.csv")



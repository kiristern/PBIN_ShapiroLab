.libPaths()

library(tidyverse)
library(phyloseq)
library(lubridate)
library(lomb)
library(ggplot2)

# load seasonal data from nico
cyano_phage <- read.csv("outputs_Alex/cyanovir_module.csv")
cyano_phage.bloom <- read.csv("outputs_Alex/cyanovirB_modules.csv")
cyano_phage.nobloom <- read.csv("outputs_Alex/cyanovirNB_modules.csv")
cyano_phage.mesotrophic <- read.csv("outputs_Alex/cyanovirMesotroph_modules.csv")
cyano_phage.eutrophic <- read.csv("outputs_Alex/cyanovirEutroph_modules.csv")

modularity.phage_cyano <- read.csv("outputs_Alex/Modularity_Phage-Cyano_AM.csv")
nestedness.phage_cyano <- read.csv("outputs_Alex/Nestedness_Phage-Cyano_AM.csv")

# split data to long df
to_long_df <- function(df){
  df <- df %>%
    separate_rows(PhageTaxa, sep = ", ") # Separate the PhageTaxa column into multiple rows
  return(df)
}

cyano_phage <- to_long_df(cyano_phage)
max(cyano_phage$ModuleLabel) # largest module number
# add pval column
# cyano_phage$pval <- 0.05

cyano_phage.bloom <- to_long_df(cyano_phage.bloom)
cyano_phage.nobloom <- to_long_df(cyano_phage.nobloom)
cyano_phage.mesotrophic <- to_long_df(cyano_phage.mesotrophic)
cyano_phage.eutrophic <- to_long_df(cyano_phage.eutrophic)

# cyano_phage %>% head()


############################
# Import from .RData
############################
# load from backup -- faster -- use if only need phyloseq data
lomb_env <- new.env() # Create a new environment
# load("arxiv/Large_data/lomb backup.RData", envir = lomb_env) # Load the file into the new environment
load("lomb2025.RData", envir = lomb_env) # Load the file into the new environment
ls(lomb_env)  # List the variables in the new environment

virps <- lomb_env$virps3000filt2

############################
# bit of preprocessing
# make sure as.Date
sample_data(virps)$Date <- as.Date(sample_data(virps)$Date, format = "%Y-%m-%d")

# sort by date
sample_data(virps) <-sample_data(virps)[order(sample_data(virps)$Date),]
# Create an increasing numeric sequence
sample_data(virps)$Date_numeric <- seq_along(sample_data(virps)$Date)

# add new column to sample_data (combine day and month)
sample_data(virps)$day_month <- as.factor(paste(sample_data(virps)$Day, sample_data(virps)$Month, sep = "."))
# drop everything after last "." from description

taxa_sums(virps) %>% head()

virpsfilt_relab <- transform_sample_counts(virps, function(x) x / sum(x))
taxa_sums(virpsfilt_relab) %>% head()


library(microeco)
library(file2meco)
library(magrittr)
meco_virps <- phyloseq2meco(virps)

library(igraph)
library(dplyr)
library(tidyr)

# https://yunranchen.github.io/intro-net-r/igraph.html#built-networks-from-external-sources-basic-visualization-and-more-on-network-descriptions
igraph_to_meco <- function(df_module, meco_ps){
  net <- df_module %>% select(PhageTaxa, CyanoTaxa, ModuleLabel) %>% 
  group_by(PhageTaxa, CyanoTaxa) 
  # create a graph object
  g <- graph_from_data_frame(net, directed = FALSE)
  # add +ve or -ve edges -- needed for cal_sum_links
  E(g)$label <- ifelse(E(g)$ModuleLabel > 0, '+', '-') 
  # assign group membership to node attribute
  V(g)$group <- components(g)$membership

  new_network <- clone(meco_ps)

  test1 <- trans_network$new(dataset=new_network, 
    graph=g,
    cor_method=NULL,
    taxa_level="OTU",
    )
  test1$cal_network(network_method=NULL)
  test1$tax_table <- new_network$tax_table
  # add graph to network
  test1$res_network <- g
  # test1$save_network("cyano_phage")
  # test1$plot_network()

  # add predefined modules
  test1$res_node_table <- data.frame(name = names(V(g)), module = V(g)$group)

  test1$cal_network_attr()
  test1$res_network_attr

  test1$cal_module(method="cluster_fast_greedy") 
  test1$cal_network_attr()
  # test1$res_network_attr
  test1$get_node_table(node_roles = TRUE)
  # head(test1$res_node_table)
  # test1$plot_taxa_roles()
  test1$get_edge_table()
  # edge_res <- test1$res_edge_table
  test1$get_adjacency_matrix()
  return(test1)
}

meconet_bloom <- igraph_to_meco(cyano_phage.bloom, meco_virps)
meconet_nobloom <- igraph_to_meco(cyano_phage.nobloom, meco_virps)
meconet_cp <- igraph_to_meco(cyano_phage, meco_virps)
meconet_meso <- igraph_to_meco(cyano_phage.mesotrophic, meco_virps)
meconet_eutro <- igraph_to_meco(cyano_phage.eutrophic, meco_virps)

# nice plot
meconet_nobloom$plot_network(edge.arrow.size=.5, 
  vertex.color="gold", 
  vertex.size=3, 
  vertex.frame.color="gray", 
  vertex.label.color="black", 
  vertex.label.cex=.5, 
  vertex.label.dist=2, 
  edge.curved=0.5,
  layout=layout_with_lgl
)

# colour by module
meconet_nobloom$plot_network(
  edge.arrow.size=.5, 
  vertex.color=meconet_nobloom$res_node_table$module, 
  vertex.size=3, 
  vertex.frame.color="gray", 
  vertex.label.color="black", 
  vertex.label.cex=.5, 
  vertex.label.dist=2, 
  edge.curved=0.5,
  # layout=layout_with_lgl
)


# verify same with GT
meconet_cp$res_edge_table %>% 
  select(node1, node2, ModuleLabel) %>%
  filter(ModuleLabel==4)
cyano_phage %>% 
  filter(ModuleLabel==4) %>% 
  select(PhageTaxa, CyanoTaxa, ModuleLabel)

# verify same with GT
meconet_bloom$res_edge_table %>% 
  select(node1, node2, ModuleLabel) %>%
  filter(ModuleLabel==4)
cyano_phage.bloom %>% 
  filter(ModuleLabel==4) %>% 
  select(PhageTaxa, CyanoTaxa, ModuleLabel)



# sample_data(virps) %>% head()
# sample_data(virps)$N_range %>% unique() # suppose this col ?
sample_data(virps)$P_range %>% unique() # Alexis specified it is this col


library(meconetcomp)
# https://chiliubio.github.io/microeco_tutorial/meconetcomp-package.html
# create lists for correlation networks
cp_net <- list()
# select samples of bloom
# clone a deep copy
tmp <- clone(meco_virps)
# change sample_table directly
meco_virps$sample_table$Bloom %>% unique()
tmp$sample_table %<>% subset(Bloom=="yes")
# trim all files in the object
tmp$tidy_dataset()
# use filter_thres parameter to filter the feature with low relative abundance
tmp$cal_network(network_method=NULL)
tmp$tax_table <- tmp$tax_table
# put the network into the list
cp_net$Bloom <- meconet_bloom

# select samples of "no bloom" group
tmp <- clone(meco_virps)
tmp$sample_table %<>% subset(Bloom == "no")
tmp$tidy_dataset()
# tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
# tmp$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
tmp$cal_network(network_method=NULL)
tmp$tax_table <- tmp$tax_table
cp_net$noBloom <- meconet_nobloom

# select samples of "Meso" group
tmp <- clone(meco_virps)
tmp$sample_table %<>% subset(P_range == "Mesotrophic")
tmp$tidy_dataset()
# tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
# tmp$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
tmp$cal_network(network_method=NULL)
tmp$tax_table <- tmp$tax_table
cp_net$Mesotrophic <- meconet_meso

# select samples of "Eutroph" group
tmp <- clone(meco_virps)
tmp$sample_table %<>% subset(P_range == "Eutrophic")
tmp$tidy_dataset()
# tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
# tmp$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
tmp$cal_network(network_method=NULL)
tmp$tax_table <- tmp$tax_table
cp_net$Eutrophic <- meconet_eutro

# select samples of "full" group
tmp <- clone(meco_virps)
# tmp$sample_table %<>% subset(N_range == "Eutrophic")
tmp$tidy_dataset()
# tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
# tmp$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
tmp$cal_network(network_method=NULL)
tmp$tax_table <- tmp$tax_table
cp_net$CyanoPhageFull <- meconet_cp

# Network modularity for all networks -- partition modules for all the networks in the list
cp_net %<>% cal_module(undirected_method="cluster_fast_greedy")

# Network topoligical attributes for all networks -- exact all the res_network_attr tables in the networks and merge them into one final table 
cp_res_net_attr <- cal_network_attr(cp_net)
cp_res_net_attr

# Node and edge properties extraction for all networks
cp_net %<>% get_node_table(node_roles = TRUE) %>% get_edge_table
net_bloom_node_tab <- cp_net$Bloom$res_node_table
net_bloom_edge_tab <- cp_net$Bloom$res_edge_table



# compare nodes across networks
# obtain the node distributions by searching the res_node_table in the object
tmp <- node_comp(cp_net, property = "name")
# obtain nodes intersection
tmp1 <- trans_venn$new(tmp, ratio = "numratio")
g1 <- tmp1$plot_venn(fill_color = FALSE)
g1 <- g1 + ggtitle("Node overlap networks")
g1
ggsave("figs250113/all-node_venn.png", g1, width = 7, height = 6)
# calculate jaccard distance to reflect the overall differences of networks
tmp$cal_betadiv(method = "jaccard")
tmp$beta_diversity$jaccard

# compare edges across networks -- studying edges overlap
# get the edge distributions across networks
tmp <- edge_comp(cp_net)
# obtain edges intersection
tmp1 <- trans_venn$new(tmp, ratio = "numratio")
g1 <- tmp1$plot_venn(fill_color = FALSE)
g1 <- g1 + ggtitle("Edge overlap networks")
g1
ggsave("figs250113/all_edge-overlap.png", g1, width = 7, height = 6)
# calculate jaccard distance
tmp$cal_betadiv(method = "jaccard")
tmp$beta_diversity$jaccard



# extract the subset of edges according to the intersections of edges across networks
# first obtain edges distribution and intersection
tmp <- edge_comp(cp_net)
tmp1 <- trans_venn$new(tmp)
tmp1$data_summary
tmp1$data_details %>% head()

# convert intersection result to a microtable object
tmp2 <- tmp1$trans_comm()

# https://github.com/ChiLiubio/microeco/issues/321#issuecomment-1940357378
# t1 <- clone(meconet_cp)
# t1$get_edge_table()
# t1$get_node_table(node_roles = TRUE)
# t1$cal_network_attr()
# t1$res_network_attr 
# meco_module <- t1$trans_comm(use_col="module") 
# View(meco_module$sample_table)
# meco_module$otu_table %>% head()


# extract the intersection of all the three networks ("Bloom", "noBloom")
# please use colnames(tmp2$otu_table) to find the required name
# tmp2$otu_table %>% head()
names(tmp2$otu_table)

# throws error:
Intersec_all <- subset_network(cp_net, venn = tmp2, name="Bloom&CyanoPhageFull")
# Intersec_all is a trans_network object
# for example, save Intersec_all as gexf format
# Intersec_all$save_network("Intersec_all.gexf")
tmp1$data_details %>% head()


# get overlapping edges (eg. Bloom and CyanoPhageFull)
edge_table_filter <- tmp$otu_table %>% filter(Bloom == 1, CyanoPhageFull == 1)
nodes_keep <- edge_table_filter %>% rownames() %>% 
  str_split(" -- ") %>% # split the string by " -- "
  # map_chr(1) %>% 
  unlist() %>% # flatten list
  unique()

subset_cpFull_bloom <- meconet_cp$subset_network(node=nodes_keep)

# https://github.com/ChiLiubio/microeco/issues/404#issuecomment-2327762965
# igraph to trans_network object
new_obj <- clone(meconet_cp)
new_obj$res_network <- subset_cpFull_bloom
# for example, save Intersec_all as gexf format
new_obj$save_network("intersect_cpFull-Bloom.gexf")


# follow rest of docs: https://chiliubio.github.io/microeco_tutorial/model-based-class.html#trans_network-class

# calculate network attributes
new_obj$cal_network_attr()
new_obj$res_network_attr

# get node properties
new_obj$get_node_table(node_roles = TRUE)
new_obj$res_node_table

# get edge properties
new_obj$get_edge_table()
new_obj$res_edge_table

# get adjacency matrix
new_obj$get_adjacency_matrix()
new_obj$res_adjacency_matrix

# add_label = TRUE can be used to directly add text label for points
new_obj$plot_taxa_roles(use_type = 1)
# plot node roles with phylum information
new_obj$plot_taxa_roles(use_type = 2)

# The eigengene of a module, i.e. the first principal component of PCA, represents the main variance of the abundance in the species of the module.
new_obj$cal_eigen()
new_obj$res_eigen

# create trans_env object
t2 <- trans_env$new(dataset = meco_virps, add_data = sample_data(virps)[, 2:33])
# calculate correlations
t2$cal_cor(add_abund_table = new_obj$res_eigen)
# plot the correlation heatmap
t2$plot_cor()

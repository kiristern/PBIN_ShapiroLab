
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

cyano_phage %>% head()


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
  test1$res_network_attr
  test1$get_node_table(node_roles = TRUE)
  head(test1$res_node_table)
  # test1$plot_taxa_roles()
  test1$get_edge_table()
  # edge_res <- test1$res_edge_table
  test1$get_adjacency_matrix()
  return(test1)
}

meconet_bloom <- igraph_to_meco(cyano_phage.bloom, meco_virps)
meconet_nobloom <- igraph_to_meco(cyano_phage.nobloom, meco_virps)
meconet_cp <- igraph_to_meco(cyano_phage, meco_virps)

meconet_nobloom$plot_network()









cp <- cyano_phage %>% select(PhageTaxa, CyanoTaxa, ModuleLabel) %>% 
  group_by(PhageTaxa, CyanoTaxa) #%>%
  # expand(edge=c(1:ModuleLabel)) %>%
  # select(-edge)
g_cp <- graph_from_data_frame(cp, directed = FALSE)
g_cp
plot(g_cp,edge.arrow.size=.5, vertex.color="gold", vertex.size=3, 
     vertex.frame.color="gray", vertex.label.color="black", 
     vertex.label.cex=.5, vertex.label.dist=2, edge.curved=0.5,
    #  layout=layout_with_lgl
)

# add +ve or -ve edges -- needed for cal_sum_links
# https://github.com/ChiLiubio/microeco/issues/281#issuecomment-1766136594
E(g_cp)$label <- ifelse(E(g_cp)$ModuleLabel > 0, '+', '-') 


# check with GT
cyano_phage %>% select(PhageTaxa, CyanoTaxa, ModuleLabel) %>% 
  filter(ModuleLabel==1) 
# chech membership igraph  
components(g_cp)

# assign group membership to node attribute
# https://stackoverflow.com/questions/61219754/adding-group-membership-as-a-node-attribute-after-simulating-a-network-in-igraph
V(g_cp)$group <- components(g_cp)$membership

plot(g_cp, vertex.color=V(g_cp)$group, vertex.size=3, 
     vertex.frame.color="gray", vertex.label.color="black", 
     vertex.label.cex=.5, vertex.label.dist=2, edge.curved=0.5,
    # layout=layout_with_lgl
)


# https://yunranchen.github.io/intro-net-r/igraph.html#paths-communitites-and-related-visualization
# ground truth modules
GT_mods <- make_clusters(g_cp, V(g_cp)$group)

#non-negative values: different labels; negative values: no labels
# initial=rep(-1,vcount(g_cp))
initial<- V(g_cp)$group
fixed=rep(TRUE,vcount(g_cp))
#need to have names
names(initial)=names(fixed)=V(g_cp)$name 
# initial['Mr Hi']=1
# initial['John A']=2
# fixed['Mr Hi']=fixed['John A']=TRUE
lab=cluster_label_prop(g_cp,initial = initial,fixed = fixed)
set.seed(1)
plot(lab, g_cp, vertex.color=V(g_cp)$group, vertex.size=3, 
     vertex.frame.color="gray", vertex.label.color="black", 
     vertex.label.cex=.5, vertex.label.dist=2, edge.curved=0.5,
    # layout=layout_with_lgl
)

lab






# # https://r.igraph.org/reference/as_membership.html
# modules <- V(g_cp)$group %>% igraph::as_membership()

# modularity(g_cp, membership=modules)

# fc <- cluster_fast_greedy(g_cp, weights=NULL, merges=TRUE, modularity=TRUE, membership=modules)

# compare(fc, modules)
# compare(modules, membership(fc))

# # add modules to the network
# modularity_matrix(g_cp, 
#   membership=modules,
#   weights=NULL,
#   resolution=1,
#   directed=FALSE,
# )
# membership(g_cp)
# g_cp

# https://github.com/ChiLiubio/microeco/issues/404#issuecomment-2327762965
# convert igraph to trans_network object
new_network <- clone(meco_virps)

test1 <- trans_network$new(dataset=new_network, 
  graph=g_cp,
  cor_method=NULL,
  taxa_level="OTU",
  )
test1$cal_network(network_method=NULL)
test1$tax_table <- new_network$tax_table
# add graph to network
test1$res_network <- g_cp
# test1$save_network("cyano_phage")
# test1$plot_network()

# add predefined modules
test1$res_node_table <- data.frame(name = names(V(g_cp)), module = V(g_cp)$group)

test1$cal_network_attr()
test1$res_network_attr

# https://github.com/ChiLiubio/microeco/issues/262
# test2 <- clone(test1)
# test2$cal_module(method="cluster_fast_greedy")

# test2$get_node_table(node_roles = FALSE)
# test2$res_node_table

# test2$get_edge_table()
# test2$res_edge_table

# new_modules <- test2[rownames(test1), "module"]

# test1 <- igraph::set_vertex_attr(g_cp, "module", value = V(g_cp)$group)
# test1$res_network <- test1

# new_modules <- test2$res_node_table[rownames(test1$res_node_table), "module"]



# test1$cal_module(method="cluster_fast_greedy")
# nodes <- V(test1$res_network)$name 
# add rel abundances
# https://github.com/ChiLiubio/microeco/issues/152#issuecomment-1268576783
# V(test$res_network)$RelativeAbundance <- apply(test1$data_abund[, nodes], 2, sum) * 100 / sum(new_network$otu_table)


# then calculate network
# test1$cal_network(COR_p_thres = 0.05, COR_cut = 0.7)
# partition the network into modules ivoking the igraph function
test1$cal_module(method="cluster_fast_greedy") 
test1$cal_network_attr()
test1$res_network_attr
test1$get_node_table(node_roles = TRUE)
head(test1$res_node_table)
test1$plot_taxa_roles()

node_res <- test1$res_node_table

test1$get_edge_table()
edge_res <- test1$res_edge_table

test1$get_adjacency_matrix()
adja_res <- test1$res_adjacency_matrix

# verify same with GT
edge_res %>% select(node1, node2, ModuleLabel) %>%
  filter(ModuleLabel==1)
cyano_phage

# plot the node classification in terms of the within-module connectivity and among-module connectivity.
# test1$plot_taxa_roles(use_type = 1, add_label = TRUE)
# # plot node roles with phylum information
# test1$plot_taxa_roles(use_type = 2)
# test1$cal_sum_links()



















# # https://github.com/ChiLiubio/microeco/issues/373#issuecomment-2167158701
# # convert long format to symmetrical matrix
# # The first and second columns must be names
# vec2mat <- function(data, value_var, value_NA = 0, value_diag = 0){
# 	library(magrittr)
#   datatable <- data %>% as.data.frame
# 	datatable[, value_var] %<>% as.numeric 
# 	use_table <- datatable[, c(2:3)] # using PhageTaxa and CyanoTaxa pairs
# 	use_table[, value_var] <- datatable[, value_var]
# 	colnames(use_table) <- c("t1", "t2", "value")
# 	res_table <- rbind.data.frame(
# 		use_table, 
# 		data.frame(t1 = use_table[, 2], t2 = use_table[, 1], value = use_table[, 3])
# 	)
# 	res_table <- reshape2::dcast(res_table, t1~t2, value.var = "value") %>% 
# 		`row.names<-`(.[,1]) %>% 
# 		.[, -1] %>% 
# 		as.matrix
# 	res_table[is.na(res_table)] <- value_NA
# 	diag(res_table) <- value_diag
# 	return(res_table)
# }

# cor_matrix <- vec2mat(cyano_phage, value_var = "ModuleLabel", value_NA = 0., value_diag = 0.)
# max(cor_matrix) # largest module number

# p_matrix <- vec2mat(cyano_phage, value_var = "pval", value_NA = 1, value_diag = 1)

# # dataset NULL for customized correlation
# test1 <- trans_network$new(dataset = NULL)

# # use list to add these two matrix objects to test1
# test1$res_cor_p <- list(cor = cor_matrix, p = p_matrix)
# # then calculate network
# test1$cal_network(COR_p_thres = 0.05, COR_cut = 0.7, COR_p_adjust="none")
# test1$cal_module()
# test1$cal_network_attr()
# test1$get_node_table()
# head(test1$res_node_table)
# test1$plot_taxa_roles()

# # meco object
# cp_net <- trans_network$new(dataset = NULL)
# cp_net$res_cor_p <- list(cor = cor_matrix)







# # https://github.com/ChiLiubio/microeco/issues/274#issue-1891932090
# # Convert adjacency matrices into long-format data 
# convertCorrelationMatrix <- function(cor_matrix) {
#   cor_matrix[lower.tri(cor_matrix)] <- 0
#   df <- melt(cor_matrix, varnames = c("PhageTaxa", "CyanoTaxa"), value.name = "weight")
#   df <- subset(df, weight != 0 & PhageTaxa != CyanoTaxa)
#   return(df)
# }

# # Calculate the degree of the connected node and distinguish which module the degree of the node belongs to
# calculate_degree_to_each_module <- function(node_table, adja_mtx) {
#   connectivity_edge <- convertCorrelationMatrix(adja_mtx)
#   df_merged <- merge(connectivity_edge, node_table, by.x = "PhageTaxa", by.y = "name", all.x = TRUE)
#   df_merged <- df_merged[, c("PhageTaxa", "CyanoTaxa", "ModuleLabel")]
#   colnames(df_merged)[3] <- "PhageTaxa_ModuleLabel"
#   df_merged <- merge(df_merged, node_table, by.x = "CyanoTaxa", by.y = "name", all.x = TRUE)
#   df_merged <- df_merged[, c("PhageTaxa", "CyanoTaxa", "PhageTaxa_ModuleLabel", "ModuleLabel")]
#   colnames(df_merged)[4] <- "CyanoTaxa_ModuleLabel"
#   edge_module <- df_merged[c("PhageTaxa", "CyanoTaxa", "PhageTaxa_ModuleLabel", "CyanoTaxa_ModuleLabel")]
#   node_to_module <- rbind(edge_module %>% select(node = PhageTaxa, to_mod = CyanoTaxa_module), edge_module %>% select(node = CyanoTaxa, to_mod = PhageTaxa_module))
#   count_node_to_module <- node_to_module %>%
#     group_by(node, to_mod) %>%
#     summarize(count = n())
#   degree_df <- matrix(0, nrow = length(unique(count_node_to_module$node)),
#                       ncol = length(unique(count_node_to_module$to_mod)),
#                       dimnames = list(unique(count_node_to_module$node), sort(unique(count_node_to_module$to_mod)))) %>% as.data.frame()
#   for (i in 1:nrow(count_node_to_module)) {
#     row_name <- as.character(count_node_to_module$node[i])
#     col_name <- as.character(count_node_to_module$to_mod[i])
#     value <- count_node_to_module$count[i]
#     degree_df[row_name, col_name] <- value
#   }
#   degree_df$name <- rownames(degree_df)
#   degree_df
# }

# # Calculate the z-score and participation coeffiecient of the node
# calculate_z_p <- function(node_table, adja_mtx) {
#   degree_df <- calculate_degree_to_each_module(node_table, adja_mtx)
#   tmp_df <- merge(node_table, degree_df, by = 'name', all.x = T)
#   tmp_df[is.na(tmp_df)] <- 0
#   tmp_df$within_degree <- 0
#   tmp_df$total_degree <- 0
#   tmp_df$p_coeff <- 0
  
#   if (0 %in% unique(tmp_df$module)) {
#     module_num <- length(unique(tmp_df$module)) - 1
#   } else {
#     module_num <- length(unique(tmp_df$module))
#   }
  
#   for (i in 1:nrow(tmp_df)) {
#     this_module <- tmp_df$module[i] %>% as.character()
#     if (this_module != '0') {
#       value <- tmp_df[[this_module]][i]
#     } else {
#       value <- 0
#     }
#     tmp_df$within_degree[i] <- value
#     # p
#     mod.degree <- tmp_df[i, paste0('M',1:module_num)]
#     tmp_df$total_degree[i] <- mod.degree %>% sum()
#     tmp_df$p_coeff[i] <- 1 - (sum((mod.degree)**2) / tmp_df$total_degree[i]**2)
#   }
  
#   # z
#   tmp_df %<>%
#     group_by(module) %>%
#     mutate(z_score = (within_degree - mean(within_degree)) / sd(within_degree)) %>%
#     ungroup()
#   # In some cases, z is NaN, which is set to 0 (for example, module with only one node or standard deviation 0).
#   # p is NaN when ki = 0
#   tmp_df$z_score[is.na(tmp_df$z_score)] <- 0
#   tmp_df$p_coeff[is.na(tmp_df$p_coeff)] <- 0
#   tmp_df %>% select(-paste0('M',1:module_num))
# }

# # use node_res and adja_res  can recalculate p and z
# new_node_res <- calculate_z_p(node_res, adja_res)









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
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
# COR_p_thres represents the p value threshold
# COR_cut denotes the correlation coefficient threshold
tmp$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
# put the network into the list
cp_net$Bloom <- tmp

# select samples of "no bloom" group
tmp <- clone(meco_virps)
tmp$sample_table %<>% subset(Bloom == "no")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
tmp$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
cp_net$noBloom <- tmp

# select samples of "mesotrophic" group
tmp <- clone(meco_virps)
meco_virps$sample_table %>% head()



# Network modularity for all networks -- partition modules for all the networks in the list
cp_net %<>% cal_module(undirected_method="cluster_fast_greedy")

# Network topoligiocal attributes for all networks -- exact all the res_network_attr tables in the networks and merge them into one final table 
cp_net %<>% cal_network_attr(cp_net)




# load predefined modules
cyano.phage <- trans_network$new(dataset=cyano_phage)

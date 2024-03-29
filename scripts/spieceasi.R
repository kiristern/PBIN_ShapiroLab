#https://github.com/zdk123/SpiecEasi 
source("scripts/1_preprocess.R")
library(SpiecEasi)
library(igraph)

#spiec easi
SE_viral_cyano <- spiec.easi(list(virps_filt, cyano.ps_filt), method='mb', nlambda=100,
                   lambda.min.ratio=1e-3, pulsar.params = list(thresh = 0.05,
                                                               subsample.ratio=0.8,
                                                               seed = 1234,
                                                               ncores=4))

getStability(SE_viral_cyano)

SE_vir_dol_mic <- spiec.easi(list(virps_filt, doli.ps, micro.ps), method='mb', nlambda=75,
                    lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05,
                                                                subsample.ratio=0.8,
                                                                seed = 1234,
                                                                ncores=4))
getStability(SE_vir_dol_mic)

SE_vir_bactnoCyan <- spiec.easi(list(virps_filt, bactnoCyan_filt), method='mb', nlambda=100,
                                lambda.min.ratio=1e-4, pulsar.params = list(thresh = 0.05,
                                                                            subsample.ratio=0.8,
                                                                            seed = 1234,
                                                                            ncores=4))
getStability(SE_vir_bactnoCyan)


#spiec easi
vir.spie2 <- spiec.easi(virps_filt, method='mb', nlambda=100,
                        lambda.min.ratio=1e-3, pulsar.params = list(thresh = 0.05,
                                                                    subsample.ratio=0.8,
                                                                    seed = 1234,
                                                                    ncores=4))

vir.spie2$select$stars$summary #if coming up with empty network: b/c max value of the StARS summary statistic never crosses the default threshold (0.05). fix by lowering lambda.min.ratio to explore denser networks
getStability(vir.spie2)
sum(getRefit(vir.spie2))/2






#http://psbweb05.psb.ugent.be/conet/microbialnetworks/spieceasi.php
betaMatsym <- as.matrix(symBeta(getOptBeta(vir.spie2)))

#get weights
bm2 <- symBeta(getOptBeta(vir.spie2), mode="maxabs")
diag(bm2) <- 0
weights2 <- Matrix::summary(t(bm2))[,3]
FG.ig.vir <- adj2igraph(Matrix::drop0(getRefit(vir.spie2)),
                        edge.attr=list(weight=weights2),
                        vertex.attr=list(name=taxa_names(virps_filt)))

#plot with weights
#plot_network(FG.ig, list(virps_filt, cyanops_filt))

vir.corr.tab <- igraph::as_data_frame(FG.ig.vir, what="edges")
write.csv(vir.corr.tab, 'filtered_vir-vir_correlation.csv')

#plot vircyn connections with weights only
virplot <- graph_from_data_frame(vir.corr.tab, directed = TRUE, vertices = NULL)

#https://ramellose.github.io/networktutorials/workshop_MDA.html
#Network centrality: degree centrality (ie. degree = number of connections a node has)
virspiec.deg <- igraph::degree(virplot)
hist(virspiec.deg)
range(virspiec.deg)

library(ggplot2)
library(ggnet)
ggnet2(vir.corr.tab,
       alpha=0.75,
       #shape = factor(dtype),
       #shape.legend = "Type",
       node.size = virspiec.deg,
       size.legend = "Degree of Centrality",
       size.cut = 4,
       edge.size = abs(vir.corr.tab[,3]), edge.alpha = 0.5, edge.lty = ifelse(vir.corr.tab$weight > 0, 1, 2),
       label = colnames(vir.spie2$est$data), label.size = 1)+
  ggtitle("Viral correlation network")
# guides(color=FALSE)


#Check which OTUs are part of different modules.
#https://users.dimi.uniud.it/~massimo.franceschet/R/communities.html

#GREEDY COMMUNITY DETECTION
clusters.vir<- cluster_fast_greedy(as.undirected(virplot), weights = abs(E(virplot)$weight))

modularity(clusters.vir)

#membership of nodes
membership(clusters.vir)
#number of communities
length(clusters.vir)
#size of communities
sizes(clusters.vir)
#crossing edges
crossing(clusters.vir, virplot)

#see which edge connects two different communities
which(crossing(clusters.vir, virplot) == T)
length(which(crossing(clusters.vir, virplot) == T)) #number of cross community interactions


#plot communities without shaded regions

ggnet2(vir.corr.tab,
       color = membership(clusters.vir),
       alpha=0.75,
       node.size = virspiec.deg,
       size.legend = "Degree of Centrality",
       size.cut = 8,
       edge.size = abs(vir.corr.tab[,3]), edge.alpha = 0.5, edge.lty = ifelse(vir.corr.tab$weight > 0, 1, 2),
       label = colnames(vir.spie2$est$data), label.size = 1)+
  ggtitle("Viral correlation network by clusters")
#guides(size=FALSE)




#### VIRAL CYANO ###
#http://psbweb05.psb.ugent.be/conet/microbialnetworks/spieceasi.php
betaMatsym <- as.matrix(symBeta(getOptBeta(SE_viral_cyano)))

dim(betaMatsym)

#select for cyano - viral connections only
list(name=c(taxa_names(virps_filt), taxa_names(doli.ps), taxa_names(micro.ps)))
length(taxa_names(virps_filt))
length(c(taxa_names(virps_filt), taxa_names(doli.ps), taxa_names(micro.ps)))


#get weights
#https://github.com/zdk123/SpiecEasi/issues/81 
bm <- symBeta(getOptBeta(SE_viral_cyano), mode="maxabs")
diag(bm) <- 0
weights <- Matrix::summary(t(bm))[,3]
FG.ig.vir.cyn <- adj2igraph(Matrix::drop0(getRefit(SE_viral_cyano)),
                            edge.attr=list(weight=weights),
                            vertex.attr = list(name=c(taxa_names(virps_filt), taxa_names(cyano.ps_filt))))

covar.vir.cyn <- igraph::as_data_frame(FG.ig.vir.cyn, what="edges")

#isolate for viral-cyano interactions only
library(dplyr)
vircyn <- covar.vir.cyn %>% 
  filter(across(to, ~ !grepl('vir_', .))) %>%
  filter(across(from, ~grepl('vir_', .))) 
head(vircyn)
write.csv(vircyn, "filtered_viral-cyano_correlation.csv")

#plot vircyn connections with weights only
vircyan.plot <- graph_from_data_frame(vircyn, directed = TRUE, vertices = NULL)

# get dtype for cyano
dtype.cyan <- as.factor(c(rep("Phage", length(unique(vircyn[,1]))), rep("Cyanobacteria", length(unique(vircyn[,2])))))
otu.id.cyan <- colnames(SE_viral_cyano$est$data)

#https://ramellose.github.io/networktutorials/workshop_MDA.html
#Network centrality: degree centrality (ie. degree = number of connections a node has)
degree.cyan <- igraph::degree(vircyan.plot)
hist(degree.cyan)
range(degree.cyan)

library(ggplot2)
library(ggnet)
ggnet2(vircyn,
       color = dtype.cyan, palette = c("Phage" = "#E1AF00", "Cyanobacteria" = "steelblue"), 
       alpha=0.75,
       #shape = factor(dtype),
       #shape.legend = "Type",
       node.size = degree.cyan,
       size.legend = "Degree of Centrality",
       size.cut = 6,
       edge.size = abs(vircyn[,3]), edge.alpha = 0.5, edge.lty = ifelse(vircyn$weight > 0, 1, 2),
       label = otu.id.cyan, label.size = 1)+
  ggtitle("Viral and Cyanobacterial correlation network")


#Check which OTUs are part of different modules.
#https://users.dimi.uniud.it/~massimo.franceschet/R/communities.html

#GREEDY COMMUNITY DETECTION
clust.cyan<- cluster_fast_greedy(as.undirected(vircyan.plot), weights = abs(E(vircyan.plot)$weight))
modularity(clust.cyan)

#modularity matrix
B.cyan = modularity_matrix(vircyan.plot, membership(clust.cyan))

ggnet2(vircyn,
       color = membership(clust.cyan),
       alpha=0.75,
       shape = factor(dtype.cyan),
       shape.legend = "Type",
       node.size = degree.cyan,
       size.legend = "Degree of Centrality",
       size.cut = 8,
       edge.size = abs(vircyn[,3]), edge.alpha = 0.5, edge.lty = ifelse(vircyn$weight > 0, 1, 2),
       label = otu.id.cyan, label.size = 1)+
  ggtitle("Viral with Cyanobacterial correlation network by clusters")
#guides(size=FALSE)






betaMatsym2 <- as.matrix(symBeta(getOptBeta(SE_vir_dol_mic)))

bm2 <- symBeta(getOptBeta(SE_vir_dol_mic), mode="maxabs")
diag(bm2) <- 0
weights2 <- Matrix::summary(t(bm2))[,3]
FG.ig.vdm <- adj2igraph(Matrix::drop0(getRefit(SE_vir_dol_mic)),
                        edge.attr=list(weight=weights2),
                        vertex.attr = list(name=c(taxa_names(virps_filt), taxa_names(doli.ps), taxa_names(micro.ps))))

covar.vdm <- igraph::as_data_frame(FG.ig.vdm, what="edges")
head(covar.vdm)

vdm <- covar.vdm %>% 
  filter(across(to, ~ !grepl('vir_', .))) %>%
  filter(across(from, ~grepl('vir_', .)))

vdm.plot <- graph_from_data_frame(vdm, directed = TRUE, vertices = NULL)

write.csv(vdm, 'data/vir-cyano_correlation.csv')

#get dtype for doli-micro
library(stringr)
which(as.data.frame(str_count(vdm$to, "micro_"))=="1", arr.ind=T) #see which positions micro_ are in in col 2 of df

uniq.vdm <- unique(vdm$to)

repdm <- list()
for (j in str_count(uniq.vdm, "doli_")){
  if (j == 1){
    repdm[length(repdm)+1] <- print("Dolichospermum")
  } else {
    repdm[length(repdm)+1] <- print("Microcystis")
  }
}
repdm <- unlist(repdm)
head(repdm, n=8)

dtype.vdm <- as.factor(c(rep("Phage", length(unique(vdm[,1]))), repdm))
otu.id.vdm <- colnames(SE_vir_dol_mic$est$data)

degree.vdm <- igraph::degree(vdm.plot)

ggnet2(vdm,
       color = dtype.vdm, palette = c("Phage" = "#E1AF00", "Dolichospermum" = "red", "Microcystis" = "steelblue"), 
       alpha=0.75,
       #shape = factor(dtype),
       #shape.legend = "Type",
       node.size = degree.vdm,
       size.legend = "Degree of Centrality",
       size.cut = 6,
       edge.size = abs(vdm[,3]), edge.alpha = 0.5, edge.lty = ifelse(vdm$weight > 0, 1, 2),
       label = otu.id.vdm, label.size = 1)+
  ggtitle("Viral and Microcystis/Dolichospermum correlation network")

clust.vdm <- cluster_fast_greedy(as.undirected(vdm.plot), weights = abs(E(vdm.plot)$weight))

#membership of nodes
memb.vdm <- membership(clust.vdm)
memb.vdm
#modularity
modularity(clust.vdm, memb.vdm)
B4 = modularity_matrix(vdm.plot, membership(clust.vdm))
#number of communities
length(clust.vdm)
#size of communities
sizes(clust.vdm)
mean(sizes(clust.vdm))
min(sizes(clust.vdm))
max(sizes(clust.vdm))
median(sizes(clust.vdm))

#crossing edges
crossE <- crossing(clust.vdm, vdm.plot)
length(which(crossE == T))


#plot communities without shaded regions
ggnet2(vdm,
       color = membership(clust.vdm),
       alpha=0.75,
       shape = factor(dtype.vdm),
       shape.legend = "Type",
       node.size = degree.vdm,
       size.legend = "Degree of Centrality",
       size.cut = 8,
       edge.size = abs(vdm[,3]), edge.alpha = 0.5, edge.lty = ifelse(vdm$weight > 0, 1, 2),
       label = otu.id.vdm, label.size = 1)+
  ggtitle("Viral with Microcystis and Dolichospermum correlation network by clusters")

# Mean degree (number of edges per ASV)
edges.asv <- as.data.frame(table(covar.vdm$to))
mean(edges.asv$Freq)
#get centrality for bacterial nodes only
doli.E.deg <- edges.asv %>%
  filter(across(Var1, ~grepl("doli_", .)))
micro.E.deg <- edges.asv %>%
  filter(across(Var1, ~grepl("micro_", .)))
md.E.deg <- rbind(doli.E.deg, micro.E.deg)
mean(md.E.deg$Freq)

# Mean centrality
degree.vdm.df <- as.data.frame(degree.vdm)
mean(degree.vdm.df$degree.vdm)
#get centrality for bacterial nodes only
degree.vdm.df$asv <- rownames(degree.vdm.df)
doli.centrality <- degree.vdm.df %>%
  filter(across(asv, ~grepl("doli_", .))) 
micro.centrality <- degree.vdm.df %>%
  filter(across(asv, ~grepl("micro_", .))) 
md.centrality <- rbind(doli.centrality, micro.centrality)
mean(md.centrality$degree.vdm)





### VIRAL - BACTNOCYAN ####
#get weights
#https://github.com/zdk123/SpiecEasi/issues/81 
bm <- symBeta(getOptBeta(SE_vir_bactnoCyan), mode="maxabs")
diag(bm) <- 0
weights <- Matrix::summary(t(bm))[,3]
FG.ig.bactnoCyan <- adj2igraph(Matrix::drop0(getRefit(SE_vir_bactnoCyan)),
                            edge.attr=list(weight=weights),
                            vertex.attr = list(name=c(taxa_names(virps_filt), taxa_names(bactnoCyan_filt))))

covar.vir.bactnoCyan <- igraph::as_data_frame(FG.ig.bactnoCyan, what="edges")

#isolate for viral-cyano interactions only
library(dplyr)
BnoC <- covar.vir.bactnoCyan %>% 
  filter(across(to, ~ !grepl('vir_', .))) %>%
  filter(across(from, ~grepl('vir_', .)))
head(BnoC)
write.csv(BnoC, "filtered_virbactnoCyan_correlation.csv")

#plot vircyn connections with weights only
virbactnocyan.plot <- graph_from_data_frame(BnoC, directed = TRUE, vertices = NULL)

# get dtype for cyano
dtype.cyan <- as.factor(c(rep("Phage", length(unique(BnoC[,1]))), rep("Bacteria", length(unique(BnoC[,2])))))
otu.id.cyan <- colnames(SE_vir_bactnoCyan$est$data)

#https://ramellose.github.io/networktutorials/workshop_MDA.html
#Network centrality: degree centrality (ie. degree = number of connections a node has)
degree.cyan <- igraph::degree(virbactnocyan.plot)
hist(degree.cyan)
range(degree.cyan)

library(ggplot2)
library(ggnet)
ggnet2(BnoC,
       color = dtype.cyan, palette = c("Phage" = "#E1AF00", "Bacteria" = "steelblue"), 
       alpha=0.75,
       #shape = factor(dtype),
       #shape.legend = "Type",
       node.size = degree.cyan,
       size.legend = "Degree of Centrality",
       size.cut = 6,
       edge.size = abs(BnoC[,3]), edge.alpha = 0.5, edge.lty = ifelse(BnoC$weight > 0, 1, 2),
       label = otu.id.cyan, label.size = 1)+
  ggtitle("Viral and bacterial correlation network")


#Check which OTUs are part of different modules.
#https://users.dimi.uniud.it/~massimo.franceschet/R/communities.html

#GREEDY COMMUNITY DETECTION
clust.cyan<- cluster_fast_greedy(as.undirected(virbactnocyan.plot), weights = abs(E(virbactnocyan.plot)$weight))
modularity(clust.cyan)

#modularity matrix
B.cyan = modularity_matrix(virbactnocyan.plot, membership(clust.cyan))

ggnet2(BnoC,
       color = membership(clust.cyan),
       alpha=0.75,
       shape = factor(dtype.cyan),
       shape.legend = "Type",
       node.size = degree.cyan,
       size.legend = "Degree of Centrality",
       size.cut = 8,
       edge.size = abs(BnoC[,3]), edge.alpha = 0.5, edge.lty = ifelse(BnoC$weight > 0, 1, 2),
       label = otu.id.cyan, label.size = 1)+
  ggtitle("Viral with bacterial correlation network by clusters")
#guides(size=FALSE)





#if the degree distribution of a network follows a power law, that network is scale-free
plaw.fit <- fit_power_law(degree.cyan.pos) #The fit_power_law functions fits a power law to the degree distribution of the network.
plaw.fit
#The values for the fit are compared to expected values with the Kolmogorov-Smirnov test. 
#The null hypothesis for this test is that the degree distribution is drawn from a reference distribution. 
#In this case, that reference distribution is generated from a power law.

#The null hypothesis can only be rejected if the p-value of the test is below 0.05. 
#Here, the p-value is 0.321. Therefore, we cannot conclude that the degree distribution is drawn from a different distribution than the power-law distribution.
#Scale-free networks are networks with a degree distribution that follows a power law. 
#Our result indicates that the network may be scale-free and contains nodes with a degree far larger than the average degree. 
#While there is not that much known about the effect of scale-freeness on microbial networks, studies (e.g., Cohen et al 2001, https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.86.3682) 
#indicate that scale-freeness decreases the network’s sensitivity to random attacks. 
#However, we still do not know to what extent biological networks follow a power law as we have few true biological networks.
#Lima-Mendez and van Helden (2009) (https://pubs.rsc.org/en/content/articlehtml/2009/mb/b908681a) discuss some of the weaknesses of this theory.



#plot dendogram
plot_dendrogram(clust.cyan.pos)

#Check which OTUs are part of different modules.
clust.vdm.OneIndices <- which(clust.vdm$membership > 1)
clust.vdm.OneOtus <- clust.vdm$names[clust.vdm.OneIndices]
clust.vdm.OneOtus

clust.vdm[1]


#see which edge connects two different communities
com <- as.data.frame(which(crossing(clust.vdm, vdm.plot) == T))
com$link <- row.names(com)
length(which(crossing(clust.vdm, vdm.plot) == T)) #number of cross community interactions
com <- data.frame(do.call('rbind', strsplit(as.character(com$link),'|',fixed=TRUE)))
com

links <- vdm[vdm$from %in% com$X1,]
node_name <- unique(links$from)

not_node_indices <- which(E(vdm.plot)$start != node_name) 
not_joined_edges <- E(vdm.plot)[not_node_indices] 
n <- delete_edges(vdm.pos.plot, not_joined_edges) 


n1 <- make_ego_graph(vdm.pos.plot, order=1, nodes=node_name)
n2 <- do.call(union, n1)

subgraph(vdm.pos.plot, n2)

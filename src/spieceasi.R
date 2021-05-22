#https://github.com/zdk123/SpiecEasi 
source("src/preprocess.R")
library(SpiecEasi)
library(igraph)

#spiec easi
SE_viral_cyano <- spiec.easi(list(virps_filt, cyanops_filt), method='mb', nlambda=100,
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

# vir.corr.tab.4pos <- vir.corr.tab %>% 
#   filter(weight > 0.4)
# vir.corr.tab.1neg <- vir.corr.tab %>% 
#   filter(weight < -0.1)
# strong.vir.corr.tab <- rbind(vir.corr.tab.4pos, vir.corr.tab.1neg)
# head(strong.vir.corr.tab)

#write.csv(vir.corr.tab, "vir-vir.covar.csv")

#plot vircyn connections with weights only
virplot <- graph_from_data_frame(vir.corr.tab, directed = TRUE, vertices = NULL)

#https://ramellose.github.io/networktutorials/workshop_MDA.html
#Network centrality: degree centrality (ie. degree = number of connections a node has)
virspiec.deg <- igraph::degree(virplot)
hist(virspiec.deg)
range(virspiec.deg)

library(ggplot2)
library(ggnet)
ggnet2(virplot,
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

ggnet2(virplot,
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

# vir.cyan <- betaMatsym[1:576, 577:614]
# 
# #check positive, negative, and total edges (divide by 2 because an edge is represented by 2 entries in the matrix)
# (p.edge =length(betaMatsym[betaMatsym>0])/2)
# (n.edge =length(betaMatsym[betaMatsym<0])/2)
# (tot.edge =length(betaMatsym[betaMatsym!=0])/2)
# 
# #viral-cyano connections only
# (p.edge =length(vir.cyan[vir.cyan>0]))
# (n.edge =length(vir.cyan[vir.cyan<0]))
# (tot.edge =length(vir.cyan[vir.cyan!=0]))


#get weights
#https://github.com/zdk123/SpiecEasi/issues/81 
bm <- symBeta(getOptBeta(SE_viral_cyano), mode="maxabs")
diag(bm) <- 0
weights <- Matrix::summary(t(bm))[,3]
FG.ig.vir.cyn <- adj2igraph(Matrix::drop0(getRefit(SE_viral_cyano)),
                            edge.attr=list(weight=weights),
                            vertex.attr = list(name=c(taxa_names(virps_filt), taxa_names(cyanops_filt))))

covar.vir.cyn <- igraph::as_data_frame(FG.ig.vir.cyn, what="edges")

#isolate for viral-cyano interactions only
library(dplyr)
vircyn <- covar.vir.cyn %>% 
  filter(across(to, ~ !grepl('vir_', .))) %>%
  filter(across(from, ~grepl('vir_', .))) 
head(vircyn)
write.csv(vircyn, "vircyan.cov.csv")

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
ggnet2(vircyan.plot,
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

ggnet2(vircyan.plot,
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

#### POSITIVE COVARIANCE BETWEEN CYANO AND VIRAL ####
# vircyn.pos <- covar.vir.cyn %>% 
#   filter(across(to, ~ !grepl('vir_', .))) %>%
#   filter(across(from, ~grepl('vir_', .))) %>%
#   #rename(weight = V3) %>%
#   filter(weight > 0)
# head(vircyn.pos)

# vircyan.pos.plot <- graph_from_data_frame(vircyn.pos, directed = TRUE, vertices = NULL)
# plot_network(vircyan.pos.plot)

# dtype.cyan.pos <- as.factor(c(rep("Phage", length(unique(vircyn.pos[,1]))), rep("Cyanobacteria", length(unique(vircyn.pos[,2])))))
# otu.id <- colnames(SE_viral_cyano$est$data)

#degree.cyan.pos <- igraph::degree(vircyan.pos.plot)

# dd <- degree.distribution(vircyan.pos.plot)
# plot(0:(length(dd)-1), dd, ylim=c(0,1), type='b', 
#      ylab="Frequency", xlab="Degree", main="Degree Distributions")

# ggnet2(vircyan.pos.plot,
#        color = dtype.cyan.pos, palette = c("Phage" = "#E1AF00", "Cyanobacteria" = "steelblue"), 
#        alpha=0.75,
#        #shape = factor(dtype.cyan.pos),
#        #shape.legend = "Type",
#        node.size = degree.cyan.pos,
#        size.legend = "Degree of Centrality",
#        size.cut = 6,
#        edge.size = vircyn.pos[,3], edge.alpha = 0.5,
#        label = otu.id, label.size = 1)+
#   ggtitle("Viral and Cyanobacterial correlation network")

#clust.cyan.pos <-cluster_fast_greedy(as.undirected(vircyan.pos.plot))
#clust.cyan.pos
#modularity(clust.cyan.pos)
#B.cyan.pos = modularity_matrix(vircyan.pos.plot, membership(clust.cyan.pos))

# ggnet2(vircyan.pos.plot,
#        color = membership(clust.cyan.pos),
#        alpha=0.75,
#        shape = factor(dtype.cyan.pos),
#        shape.legend = "Type",
#        node.size = degree.cyan.pos,
#        size.legend = "Degree of Centrality",
#        size.cut = 8,
#        edge.size = vircyn.pos[,3], edge.alpha = 0.5,
#        label = otu.id, label.size = 1)+
#   ggtitle("Viral with Cyanobacterial correlation network by clusters")
# #guides(size=FALSE)

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

ggnet2(vdm.plot,
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

clust.vdm<- cluster_fast_greedy(as.undirected(vdm.plot), weights = abs(E(vdm.plot)$weight))

modularity(clust.vdm)
# B4 = modularity_matrix(vdm.plot, membership(clust.vdm))
#membership of nodes
membership(clust.vdm)
#number of communities
length(clust.vdm)
#size of communities
sizes(clust.vdm)
#crossing edges
crossing(clust.vdm, vdm.plot)

#plot communities without shaded regions
ggnet2(vdm.plot,
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

#### Viral - Doli Micro Positive ####
# #postive weights only:
# weights.pos <- (1-Matrix::summary(t(bm2))[,3])/2
# FG.ig.vdm.pos <- adj2igraph(Matrix::drop0(getRefit(SE_vir_dol_mic)),
#                             edge.attr=list(weight=weights.pos),
#                             vertex.attr = list(name=c(taxa_names(virps_filt), taxa_names(doli.ps), taxa_names(micro.ps))))
# #plot_network(FG.ig.vdm.pos, list(virps_filt, cyanops_filt))

# vdm.pos <- covar.vdm %>% 
#   filter(across(to, ~ !grepl('vir_', .))) %>%
#   filter(across(from, ~grepl('vir_', .))) %>%
#   #rename(weight = V3) %>%
#   filter(weight > 0)
# head(vdm.pos)

# vdm.pos.plot <- graph_from_data_frame(vdm.pos, directed = TRUE, vertices = NULL)
# plot_network(vdm.pos.plot)

#uniq.vdm.pos <- unique(vdm.pos$to)
# repdm <- list()
# for (j in str_count(uniq.vdm, "doli_")){
#   if (j == 1){
#     repdm[length(repdm)+1] <- print("Dolichospermum")
#   } else {
#     repdm[length(repdm)+1] <- print("Microcystis")
#   }
# }
# repdm <- unlist(repdm)

# dtype.vdm.pos <- as.factor(c(rep("Phage", length(unique(vdm.pos[,1]))), repdm))
# otu.id.vdm.pos <- c(as.character(vdm.pos[,1]), as.character(vdm.pos[,2]))

# degree.vdm.pos <- igraph::degree(vdm.pos.plot)
# hist(degree.vdm.pos)
# range(degree.vdm.pos)

# ggnet2(vdm.pos.plot,
#        color = dtype.vdm.pos, palette = c("Phage" = "#E1AF00", "Dolichospermum" = "red", "Microcystis" = "steelblue"), 
#        alpha=0.75,
#        #shape = factor(dtype),
#        #shape.legend = "Type",
#        node.size = degree.vdm.pos,
#        size.legend = "Degree of Centrality",
#        size.cut = 4,
#        edge.size = vdm.pos[,3], edge.alpha = 0.5,
#        label = otu.id.vdm.pos, label.size = 1)+
#   ggtitle("Viral and Microcystis/Dolichospermum correlation network")
# # guides(color=FALSE)

#clust.vdm.pos<-cluster_fast_greedy(as.undirected(vdm.pos.plot))
#modularity(clust.vdm.pos)
# B2 = modularity_matrix(vdm.pos.plot, membership(clust.vdm.pos))
# round(B2[1,],5)

# ggnet2(vdm.pos.plot,
#        color = membership(clust.vdm.pos),
#        alpha=0.75,
#        shape = factor(dtype.vdm.pos),
#        shape.legend = "Type",
#        node.size = degree.vdm.pos,
#        size.legend = "Degree of Centrality",
#        #size.cut = 6,
#        edge.size = vdm.pos[,3], edge.alpha = 0.5,
#        label = otu.id.vdm.pos, label.size = 1)+
#   ggtitle("Viral with Microcystis and Dolichospermum correlation network by clusters")+
#   guides(size=FALSE)

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
write.csv(BnoC, "virbactnoCyan.cov.csv")

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
ggnet2(virbactnocyan.plot,
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

ggnet2(virbactnocyan.plot,
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
clust.cyan.posOneIndices=which(clust.cyan.pos$membership==1)
clust.cyan.posOneOtus=clust.cyan.pos$names[clust.cyan.posOneIndices]
clust.cyan.posOneOtus
#OR
clust.cyan.pos[2]


names(clust.vdm.pos$membership)[clust.vdm.pos$membership > 1]
clust.vdm.pos




#see which edge connects two different communities
com <- as.data.frame(which(crossing(clust.vdm.pos, vdm.pos.plot) == T))
com$link <- row.names(com)
length(which(crossing(clust.vdm.pos, vdm.pos.plot) == T)) #number of cross community interactions
com <- data.frame(do.call('rbind', strsplit(as.character(com$link),'|',fixed=TRUE)))
com

links <- vdm[vdm$from %in% com$X1,]
node_name <- unique(links$from)

not_node_indices <- which(E(vdm.pos.plot)$start != node_name) 
not_joined_edges <- E(vdm.pos.plot)[not_node_indices] 
n <-delete_edges(vdm.pos.plot, not_joined_edges) 


n1 <- make_ego_graph(vdm.pos.plot, order=1, nodes=node_name)
n2 <- do.call(union, n1)

subgraph(vdm.pos.plot, n2)



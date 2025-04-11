########## SETUP ############

rm(list=ls())
library(sRACIPE)#, lib.loc = "/Users/danramirez/localR/4.2.2-arm")
library(ggplot2)
library(FNN)
library(ClusterR)
library(nnet)
library(reshape2)
library(dplyr)
library(tidyr)
library(prodlim)
library(circlize)
library(viridis)
library(ComplexHeatmap)
library(igraph)
library(keyplayer)
source("R/utils.R")
source("R/utils_clamping.R")
source("R/scratch.R")

# set up directories
topoName <- "emt_bhtopo_26node_CLAMP"
topoDir <- file.path(getwd(),topoName)
plotDir <- file.path(topoDir,"plots_apr2025")
dataDir <- file.path(topoDir,"data")

if(!dir.exists(topoDir)) {
  dir.create(topoDir)
}
if(!dir.exists(dataDir)) {
  dir.create(dataDir)
}
if(!dir.exists(plotDir)) {
  dir.create(plotDir)
}


# load topology
topo <- loadTopo(topoName)
nGenes <- length(unique(c(topo$Source, topo$Target)))
genes_reordered <- c("Foxc2","Zeb1","Klf8","Cdh1","miR101", "Zeb2", "Snai1", "miR141",
                     "Tgfbeta","miR200a","miR200b","miR200c","miR205","miR30c","Snai2",
                     "miR34a","Twist2","miR9","Vim","Twist1","Tcf3","Gsc", "Ovol2", "Grhl2",  "Np63a", "Cldn7")


set.seed(1234)
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(3,2,4:8)]


# print topology
#plotNetwork(topo, topoName)

########## UNTREATED CONDITION ############

numModels <- 1000
numICs <- 100
simTime <- 200


# simulate WT with RACIPE
racipe_fname <- file.path(topoDir, "racipe_200IC.Rds")
if(!file.exists(racipe_fname)) {
  racipe <- sRACIPE::sracipeSimulate(circuit = topo, numModels = numModels, 
                                     nIC = numICs, simulationTime = simTime, initialNoise = 0, nNoise = 0)  
  
  
  # need to extend simulation time to get to convergence, ugh
  #sracipeIC(racipe) <- assay(racipe)
  #racipe <- sracipeSimulate(racipe, genIC = F, genParams = F, simulationTime = 1)
  
  
  saveRDS(racipe, racipe_fname)
} else {
  racipe <- readRDS(racipe_fname)
}




# get moments for normalizing other cases later
genes <- rownames(racipe)
unnormData <- t(assay(racipe))
simExp <- assay(racipe, 1)
simExp <- log2(1+simExp)
tmpMeans <- rowMeans(simExp)
tmpSds <- apply(simExp,1,sd)

racipeNorm <- sracipeNormalize(racipe)
racipeData <- as.data.frame(t(assay(racipeNorm)))



########## IDENTIFY UNIQUE STATES ############

## For each model, id unique states (doesn't need major precision)
summary_df_fname <- file.path(dataDir,"state_summary_df.Rds")
ss_unique_fname <- file.path(dataDir,"ss_unique_df.Rds")
if(!file.exists(summary_df_fname)) {
  # Find unique steady states per model
  ss_rounded <- round(as.data.frame(unnormData), 1)
  ss_rounded$Model <- rep(1:numModels, each=numICs)
  ss_unique <- ss_rounded %>%
    #group_by(Model) %>%
    distinct()
  
  
  clust_fname <- file.path(dataDir,"clust.Rds") 
  if(!file.exists(clust_fname)) {
    
    gmm = GMM(ss_unique[,1:nGenes], 2, dist_mode = "eucl_dist", seed_mode = "random_subset", km_iter = 10,
              em_iter = 10, verbose = F)          
    
    # predict centroids, covariance matrix and weights
    clust = predict(gmm, newdata = ss_unique[,1:nGenes])
    
    # Im reversing it so the initial state is cluster 1
    revClust <- clust
    revClust <- ifelse(clust == 1, 2, ifelse(clust == 2, 1, clust))
    
    #mean(ss_unique$Cdh1[which(clust == 2)])
    #ggplot(pca_df, aes(x=PC1, y=PC2, color=as.factor(revClust))) + geom_point()
    
    saveRDS(revClust, file = clust_fname)
  } else {
    clust <- readRDS(clust_fname)
  }
  ss_unique$Cluster <- clust
  
  
  summary_df <- ss_unique %>%
    select(all_of(c("Model", "Cluster"))) %>%
    group_by(Model) %>%
    summarise(
      NumStates = n(),
      Stability = case_when(
        all(Cluster == 1) ~ '1',
        all(Cluster == 2) ~ '2',
        any(Cluster == 1) & any(Cluster == 2) ~ 'bistable'
      )
    )
  
  
  # save summary
  saveRDS(summary_df, summary_df_fname)
  saveRDS(ss_unique, ss_unique_fname)
} else {
  summary_df <- readRDS(summary_df_fname)
  ss_unique <- readRDS(ss_unique_fname)
}


racipe_bistable_fname <- file.path(dataDir,"racipe_bistable_2025.Rds")
racipe_bistable_indices_fname <- file.path(dataDir,"racipe_bistable_indices_2025.Rds")
models_selected_fname <- file.path(dataDir,"racipe_bistable_indexMap_2025.Rds")
if(!file.exists(racipe_bistable_fname)) {
  # filter for models with <10 states, only bistable
  # save new racipe object w/ only bistable
  models_selected <- unlist(summary_df[which(summary_df$NumStates == 2 & 
                                               summary_df$Stability == "bistable"),"Model"])[1:500]
  
  keepIdx <- c()
  for(model in models_selected) {
    # Select only epithelial state
    #addIdx <- sample((numICs*(model-1)+1):(numICs*model), 1)
    modelStates <- ss_unique[which(ss_unique$Model == model & ss_unique$Cluster == 1),]
    addIdx <- sample(which(row.match(
                              as.data.frame(t(round(assay(racipe), 1)))[(numICs*(model-1)+1):(numICs*model),], 
                              modelStates[,genes], 
                              nomatch = NA) == 1), 1)
    addIdx <- numICs*(model-1)+addIdx # bring back index for original racipe object
    keepIdx <- c(keepIdx, addIdx)
  }
  
  racipe_bistable <- racipe[,keepIdx]
  sracipeParams(racipe_bistable) <- sracipeParams(racipe)[as.numeric(unname(models_selected)),] # parameter numbering is messssssed up in sRACIPE, not my fault
  racipe_bistable@metadata$config$simParams[["numModels"]] <- length(keepIdx)
  
  saveRDS(racipe_bistable, racipe_bistable_fname)
  saveRDS(models_selected, models_selected_fname)
  saveRDS(keepIdx, racipe_bistable_indices_fname)
  racipe_bistable_indices <- keepIdx
} else {
  racipe_bistable <- readRDS(racipe_bistable_fname)
  models_selected <- readRDS(models_selected_fname)
  racipe_bistable_indices <- readRDS(racipe_bistable_indices_fname)
}

racipe_bistable_raw <- racipe_bistable
racipe_bistable_raw@metadata$config$simParams["nIC"] <- 1
unnormData <- t(assay(racipe_bistable_raw))
racipe_bistable <- sracipeNormalize(racipe_bistable)
exprMat <- as.data.frame(t(assay(racipe_bistable)))
exprMat_norm <- as.data.frame(t(assay(racipeNorm)))


########## PCA & CLUSTERING ############

# PCA-based "cluster" assignment
pca_fname <- file.path(dataDir,"pca_2025.Rds")
if(!file.exists(pca_fname)) {
  # PCA on full data
  pca <- prcomp(exprMat_norm[,1:nGenes])
  pca_df_full <- as.data.frame(pca$x)
  pca_df <- pca_df_full[as.character(racipe_bistable_indices),]
  
  ggplot(pca_df_full, aes(x=PC1, y=PC2, color=as.data.frame(t(assay(racipeNorm)))$Cdh1)) + geom_point()
  ggplot(pca_df, aes(x=PC1, y=PC2, color=exprMat$Cdh1)) + geom_point()
  
  #### Old version on bistable models only, doesn't work now 2025 w all E starts
  # pca <- prcomp(exprMat[,1:nGenes])
  # #pca$rotation[,"PC1"] <- -1*pca$rotation[,"PC1"]
  # #pca$x[,1] <- -1*pca$x[,1]
  # 
  # pca_df <- as.data.frame(pca$x)
  # 
  # ggplot(pca_df, aes(x=PC1, y=PC2, color=exprMat$Cdh1)) + geom_point()
  # 
  saveRDS(pca, pca_fname)
} else {
  pca <- readRDS(pca_fname)
  pca_df_full <- as.data.frame(pca$x)
  pca_df <- as.data.frame(pca$x)[as.character(racipe_bistable_indices),]
}

# Clustering bistable models
clust_fname <- file.path(dataDir,"clust_bistable_2025.Rds")
clust_all_fname <- file.path(dataDir,"clust_all_2025.Rds") 
if(!file.exists(clust_fname)) {
  
  
  gmm = GMM(exprMat_norm[,1:nGenes], 2, dist_mode = "eucl_dist", seed_mode = "random_subset", km_iter = 10,
            em_iter = 10, verbose = F)
  
  # predict centroids, covariance matrix and weights
  clust_full = predict(gmm, newdata = exprMat_norm[,1:nGenes])
  clust <- clust_full[racipe_bistable_indices]
  
  # Im reversing it so the initial state is cluster 1
  #revClust <- clust
  #revClust <- ifelse(clust == 1, 2, ifelse(clust == 2, 1, clust))
  
  ggplot(pca_df_full, aes(x=PC1, y=PC2, color=as.factor(clust_full))) + geom_point()
  ggplot(pca_df_full, aes(x=PC1, y=PC2, color=exprMat_norm$Cdh1)) + geom_point()
  ggplot(pca_df, aes(x=PC1, y=PC2, color=as.factor(clust))) + geom_point()
  
  
  
  
  ## Old version which expects models of both clusters
  # gmm = GMM(exprMat[,1:nGenes], 2, dist_mode = "eucl_dist", seed_mode = "random_subset", km_iter = 10,
  #           em_iter = 10, verbose = F)
  # 
  # # predict centroids, covariance matrix and weights
  # clust = predict(gmm, newdata = exprMat[,1:nGenes])
  # 
  # # Im reversing it so the initial state is cluster 1
  # #revClust <- clust
  # #revClust <- ifelse(clust == 1, 2, ifelse(clust == 2, 1, clust))
  # 
  # ggplot(pca_df, aes(x=PC1, y=PC2, color=as.factor(clust))) + geom_point()

  saveRDS(clust, file = clust_fname)
  saveRDS(clust_full, file = clust_all_fname)
} else {
  clust <- readRDS(clust_fname)
  clust_full <- readRDS(clust_all_fname)
}
exprMat$Cluster <- clust


# Plot WT as baseline
#pca_df <- as.data.frame(pca$x)
label_df <- data.frame(x=c(-3,3),y=c(3,3),text=NA)
label_df$text <- c(unname(table(clust)[1]), unname(table(clust)[2]))

image <- ggplot() +
  geom_point(data=pca_df, aes(x=PC1, y=PC2, color=as.factor(clust))) +
  ggtitle(paste0("WT")) +
  guides(color=guide_legend("Cluster")) +
  geom_text(data=label_df, aes(x=x,y=y,label=text))
image

image <- ggplot() +
  geom_point(data=pca_df, aes(x=PC1, y=PC2, color=exprMat$Cdh1)) +
  ggtitle(paste0("WT")) +
  guides(color=guide_legend("CDH1")) +
  geom_text(data=label_df, aes(x=x,y=y,label=text))
image




########## IDENTIFY CLAMP VALUES ############
## Target: dataframe with model id, gene, expression, cluster (very long, maybe a wider format)
clamp_df_fname <- file.path(dataDir,"clamp_values_2025.Rds")
if(!file.exists(clamp_df_fname)) {
  keepIdx <- c()
  for(model in models_selected) {
    # add steady states for cluster 1 and 2
    addIdx <- which(ss_unique$Model == model)
    keepIdx <- c(keepIdx, addIdx)
    
  }
  clamp_df <- pivot_longer(ss_unique[keepIdx,], cols = all_of(genes), 
                           names_to = "Gene", values_to = "Expression")
  clamp_df$ModelIndex <- as.numeric(factor(clamp_df$Model))
  saveRDS(clamp_df, clamp_df_fname)
  
} else {
  clamp_df <- readRDS(clamp_df_fname)
}



image <- ggplot(data = clamp_df, aes(x = Gene, y = Expression, fill = as.factor(Cluster))) + 
  geom_boxplot() +
  labs(title = "Gene Expression by Cluster", 
       x = "Gene", 
       y = "Expression",
       fill = "Cluster") +
  theme_minimal() +
  #scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust text angle for better visibility
image

clamp_boxplot_fname <- file.path(plotDir,"clamp_val_boxplot.pdf")
pdf(clamp_boxplot_fname, height = 10, width = 10)
print(image)
dev.off()

########## BISTABLE WT VISUALIZATIONS ############

# Plot WT distribution
pc1_weight <- round(100*summary(pca)$importance[2,1],2)
pc2_weight <- round(100*summary(pca)$importance[2,2],2)
plot_xlab <- paste("PC1 (",pc1_weight,"%)",sep="")
plot_ylab <- paste("PC2 (",pc2_weight,"%)",sep="")

image <- ggplot(pca_df, aes(x=PC1, y=PC2, color=as.factor(clust))) +
  geom_point(size=3) +
  guides(color=guide_legend(title = "Cluster")) +
  scale_color_manual(values=cbPalette) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  theme_sticcc() +
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color="black"))
image

pca_plot_fname <- file.path(plotDir,"pca_wt_states_preSignaling.pdf")
pdf(pca_plot_fname, height = 10, width = 10)
print(image)
dev.off()

 
# Plot WT gene expression dist
ha_df <- data.frame(Cluster=clust)

# Create an annotation object for the columns
column_annotation <- HeatmapAnnotation(df = ha_df, 
                                       col=list(Cluster=c("1"=unname(cbPalette[1]),"2"=unname(cbPalette[2]),
                                                          "3"=unname(cbPalette[3]),"4"=unname(cbPalette[4]))))
# Create the heatmap with annotation
image <- Heatmap(as.matrix(t(exprMat[,1:nGenes])), 
                 name = "Expression", 
                 top_annotation = column_annotation,
                 row_names_gp=gpar(fontsize=12))
image

wt_hmap_fname <- file.path(plotDir,"wt_hmap_preSignaling.pdf")
pdf(wt_hmap_fname, height = 10, width = 10)
print(image)
dev.off()


########## FULL WT VISUALIZATIONS ############

## Project full data to PCA
ss_unique_norm <- ss_unique[,1:nGenes]
ss_unique_norm[,1:length(genes)] <- log2(1+ss_unique_norm[,1:length(genes)]) # Log transform
ss_unique_norm[,1:length(genes)] <- sweep(ss_unique_norm[,1:length(genes)], 2, tmpMeans, FUN = "-") # scale
ss_unique_norm[,1:length(genes)] <- sweep(ss_unique_norm[,1:length(genes)], 2, tmpSds, FUN = "/") # scale
ss_unique_pca <- as.data.frame(predict(pca, ss_unique_norm))
ss_unique_pca$Cluster <- ss_unique$Cluster

pca_plot_df <- pca_df
pca_plot_df$Cluster <- as.factor(clust)

image <- ggplot(ss_unique_pca, aes(x=PC1, y=PC2, color=as.factor(Cluster))) +
  geom_point(size=3, alpha=0.7) +
  geom_point(data=pca_plot_df, aes(x=PC1, y=PC2, fill=Cluster), color="black", pch=21, size=3) +
  guides(color=guide_legend(title = "Cluster")) +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(values=cbPalette) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  theme_sticcc() +
  guides(fill="none") +
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color="black"))
image
pca_plot_fname <- file.path(plotDir,"fig4a_pca_wt_dist_hilites.pdf")
pdf(pca_plot_fname, height = 10, width = 10)
print(image)
dev.off()



# Plot WT gene expression dist
ha_df <- data.frame(Cluster=ss_unique$Cluster)

# Create an annotation object for the columns
column_annotation <- HeatmapAnnotation(df = ha_df, 
                                       col=list(Cluster=c("1"=unname(cbPalette[1]),"2"=unname(cbPalette[2]),
                                                          "3"=unname(cbPalette[3]),"4"=unname(cbPalette[4]))))
# Create the heatmap with annotation
image <- Heatmap(as.matrix(t(ss_unique_norm[,1:nGenes])), 
                 name = "Expression", 
                 top_annotation = column_annotation,
                 row_names_gp=gpar(fontsize=16))

wt_hmap_fname <- file.path(plotDir,"fig4b_wt_hmap_ALL.pdf")
pdf(wt_hmap_fname, height = 10, width = 10)
print(image)
dev.off()



## Multistability histogram
summary_hist_df_in <- read.table(file.path("/Users/danramirez/Desktop/NEU/projects/emt_state_transitions",
                                           "pyracipe","emt_bhtopo_26node_5k200ic_rep2",
                                           "emt_bhtopo_26node_5k200ic_rep2_pyracipe_summary.csv"),
                                 header = T, row.names = 1, sep = ",")

# Ensure summary has the expected structure
count_df <- summary_hist_df_in %>%
  count(NO_STATES, StateIdentity, name = "Count") %>%  # Count occurrences safely
  arrange(NO_STATES)
count_long <- count_df[which(count_df$NO_STATES > 0),]

# Define colors
count_long$StateIdentity <- factor(count_long$StateIdentity, levels=c("E", "Bistable", "M"))
colors <- c("E" = "blue", "Bistable" = "purple", "M" = "red")

max_bin <- 7
count_long <- count_long %>%
  mutate(NO_STATES = ifelse(NO_STATES >= max_bin, paste0(max_bin, "+"), as.character(NO_STATES)))

# Convert NumStates to factor for proper ordering
count_long$NO_STATES <- factor(count_long$NO_STATES, levels = c(as.character(1:max_bin), paste0(max_bin, "+")))

# Create stacked bar plot
image <- ggplot(count_long, aes(x = factor(NO_STATES), y = Count, fill = StateIdentity)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  labs(x = "States",
       y = "Count",
       fill = "Cluster Identity") +
  theme_sticcc() +
  theme(axis.text = element_text(size=24, inherit.blank = FALSE))
image

stab_hist_fname <- file.path(plotDir,"fig4c_stability_hist.pdf")
pdf(stab_hist_fname, height = 10, width = 10)
print(image)
dev.off()


########## SIMULATIONS WITH CLAMPING ############


# necessary inputs:
# racipe object
# clamp df
# noise level, tcorr
# num genes
# sim time & relaxation time

signal_simTime <- 500
signal_relaxTime <- 50
signal_nGenes <- c(1,2)
signal_noise <- c(0, 0.04, 0.2) #0.04
signal_tcorr <- 10
initClust <- 1
tgtClust <- 2


expName <- paste0("bhtopo_t=",signal_simTime,"_relax_OUnoise=",paste0(signal_noise, collapse = "."),
                  "_tau=",signal_tcorr,"_genes=",paste0(signal_nGenes,collapse = "."),"_CLAMPS_2025")
pca_st <- pca
pca_st$x <- pca_df

# undebug(optimizeST_Clamp)
# paramSets <- optimizeST_Clamp(racipe_bistable_raw, # NON NORMALIZED
#            pca_st, # same size as racipe
#            clust, # vector of length nrow(pca$x)
#            initialClust = 1, # should correspond to clust
#            targetClust = 2, # should correspond to clust
#            nSigGenes = signal_nGenes, # single value or vector of choices
#            clamp_df = clamp_df,
#            outDir = file.path(topoDir),
#            expName = expName,
#            plot=F,
#            noise = signal_noise, # single value or vector of choices
#            forceRerun = F, # whether to rerun signal simulations
#            forceRecompute = F, # whether to recompute final scores
#            checkpointSize=25, # how often to report progress to user
#            totalPerturbations = 500, # Only considered if randomParams is TRUE
#            nPerturbations=30,
#            relax=T, # whether to continue simulations for some time without a signal
#            simTime = signal_simTime,
#            simTimeRelax = signal_relaxTime,
#            onlyParams = F,
#            noise_tcorr = signal_tcorr)

undebug(optimizeST_parallel)
paramSets <- optimizeST_parallel(racipe_bistable_raw,
                                 pca_st,
                                 clust,
                                 totalPerturbations = 500,
                                 initialClust=1,
                                 targetClust=2,
                                 nSigGenes=signal_nGenes,
                                 outDir = file.path(topoDir),
                                 expName = expName,
                                 paramType="G",
                                 nPerturbations = 351,
                                 plot=F,
                                 fcVals = 1,
                                 noise = signal_noise,
                                 checkpointSize = 5,
                                 forceRerun = F,
                                 forceRecompute=F,
                                 anneal = F,
                                 randomParams = F,
                                 relax=T,
                                 simTime = signal_simTime,
                                 simTimeRelax = signal_relaxTime,
                                 nCores = 6,
                                 saveTrajectory = F,
                                 printInterval = 10,
                                 noise_tcorr = signal_tcorr,
                                 clamp=T,
                                 clamp_df=clamp_df,
                                 tmpMeans=tmpMeans,
                                 tmpSds=tmpSds)



# debug(plotCondition)
# plotCondition(racipe_bistable_raw,
#                           pca_st,
#                           clust,
#                           initialClust=1,
#                           targetClust=2,
#                           expName,
#                           setName="Snai1_Zeb1_noise=0.04",
#                           plotDir=plotDir,
#                           expDir=NA,
#                           suffix="test",
#                           outDir=topoName,
#               tmpMeans = tmpMeans,
#               tmpSds = tmpSds
# ) 
# racipe_relax_fname <- file.path(topoName, expName, paste0("racipe_","Snai1_Zeb1_noise=0.04","_relaxed.Rds"))
# racipe2 <- readRDS(racipe_relax_fname)



######## AGGREGATE ANALYSIS #####
resultSet_fname <- file.path(topoDir,expName,"result_summary.Rds")
resultSet_full <- readRDS(resultSet_fname)
resultSet <- resultSet_full
rs_det <- resultSet_full[which(resultSet_full$Noise == 0),]
## this breaks sequence, must run after rs_list - only here as a reminder
#resultSet_full <- rs_list_full[[2]]
#resultSet_full$ConversionPct <- resultSet_full$ConvertingModels / resultSet_full$startingInitPopulation
#saveRDS(resultSet_full, file.path(topoDir,expName,"result_summary.Rds"))

if(!"ConversionPct" %in% colnames(resultSet_full)) {
  rs_full_list <- genes_x_transitions(resultSet_full, # dataframe w/ columns: ModelNo, SetName
                                      topoName = topoName, # string
                                      collectionName = expName, # string
                                      initialClust = 1, # int
                                      targetClust = 2, # int
                                      wt_data = exprMat_norm[,genes], # matrix/dataframe of original data
                                      clust = clust, # vector of length numSamples containing integers
                                      clust_all = clust_full, # full cluster labels matching wt_data
                                      tmpMeans = tmpMeans,
                                      tmpSds = tmpSds
  )
  resultSet_full <- rs_full_list[[2]]
  resultSet_full$ConversionPct <- resultSet_full$ConvertingModels / resultSet_full$startingInitPopulation
  saveRDS(resultSet_full, resultSet_fname)
} 


# ggplot(resultSet_full[which(resultSet_full$Noise %in% c(0,0.04)),], aes(x=factor(Noise), y=ConversionPct, fill=factor(Noise))) +
#   #geom_violin() +
#   geom_boxplot() +
#   xlab("Noise") +
#   ylab("Models Undergoing EMT (%)") +
#   scale_fill_discrete(name="Noise") +
#   theme_sticcc() +
#   theme(axis.line = element_line(color="black"), axis.title = element_text(size=20))


## Basic stats
# overall conversion pcts
summary(resultSet_full[which(resultSet_full$Noise == 0),"ConvertingModels"])
summary(resultSet_full[which(resultSet_full$Noise == 0),"ConversionPct"])

# signals above 50% eff
length(which(resultSet_full[which(resultSet_full$Noise == 0),"ConversionPct"] >= 0.5))
summary(resultSet_full[which(resultSet_full$Noise == 0 & 
                               resultSet_full$ConversionPct > 0.5),"ConversionPct"])
summary(resultSet_full[which(resultSet_full$Noise == 0 &
                               (resultSet_full$`Species 1` == "Zeb1" | 
                                  resultSet_full$`Species 2` == "Zeb1")),"ConversionPct"])
resultSet_full[which(resultSet_full$Noise == 0 & resultSet_full$ConversionPct > 0.5),"SetName"]

# Zeb1 alone deterministic eff
resultSet_full[which(resultSet_full$Noise == 0 & 
                       resultSet_full$`Species 1` == "Zeb1" & 
                       resultSet_full$NumGenes == 1),"ConversionPct"]

# signals above 5% eff
length(which(resultSet_full[which(resultSet_full$Noise == 0),"ConversionPct"] >= 0.05))

# 1-gene best signals
View(resultSet_full[which(resultSet_full$Noise == 0 & resultSet_full$NumGenes == 1),
               c("SetName", "ConversionPct")])

# 2-gene best signals
View(resultSet_full[which(resultSet_full$Noise == 0 & resultSet_full$NumGenes == 2),
               c("SetName", "ConversionPct")])

## 0.04 noise level:
# overall conversion pcts
summary(resultSet_full[which(resultSet_full$Noise == 0.04),"ConvertingModels"])
summary(resultSet_full[which(resultSet_full$Noise == 0.04),"ConversionPct"])

# signals above 50% eff
length(which(resultSet_full[which(resultSet_full$Noise == 0.04),"ConversionPct"] >= 0.5))
summary(resultSet_full[which(resultSet_full$Noise == 0.04 &
                               (resultSet_full$`Species 1` == "Zeb1" | 
                                  resultSet_full$`Species 2` == "Zeb1")),"ConversionPct"])
#resultSet_full[which(resultSet_full$Noise == 0.04 & resultSet_full$ConversionPct > 0.5),"SetName"]

# Zeb1 alone noisy eff
resultSet_full[which(resultSet_full$Noise == 0.04 & 
                       resultSet_full$`Species 1` == "Zeb1" & 
                       resultSet_full$NumGenes == 1),"ConversionPct"]


# Zeb1 vs Snai1 in noisy experiments
summary(resultSet_full[which(resultSet_full$Noise == 0.04 & 
                       (resultSet_full$`Species 1` == "Zeb1" | 
                         resultSet_full$`Species 2` == "Zeb1")),"ConversionPct"])
summary(resultSet_full[which(resultSet_full$Noise == 0.04 & 
                               (resultSet_full$`Species 1` == "Snai1" | 
                               resultSet_full$`Species 2` == "Snai1")),"ConversionPct"])

# signals above 5% eff
length(which(resultSet_full[which(resultSet_full$Noise == 0.04),"ConversionPct"] >= 0.05))


# comparison w older, smaller run
test_rs <- readRDS("/Users/danramirez/Desktop/NEU/projects/emt_state_transitions/EMTDriveR/emt_bhtopo_26node_CLAMP/bhtopo_t=500_relax_OUnoise=0.0.04.0.2_tau=10_genes=1.2_CLAMPS_v2/result_summary.Rds")
length(which(test_rs[which(test_rs$Noise == 0),"ConversionPct"] >= 0.5))
length(which(test_rs[which(test_rs$Noise == 0),"ConversionPct"] >= 0.05))
summary(test_rs[which(test_rs$Noise == 0 & 
                        (test_rs$`Species 1` == "Zeb1" | 
                           test_rs$`Species 2` == "Zeb1")),"ConversionPct"])
summary(test_rs[which(test_rs$Noise == 0.04 & 
                        (test_rs$`Species 1` == "Zeb1" | 
                           test_rs$`Species 2` == "Zeb1")),"ConversionPct"])

#selectedNoise <- c(0, 0.04, 0.2) 
#selectedNoise <- c(0, 0.04)
selectedNoise <- c(0.2)
#selectedNoise <- c(0)
resultSet <- resultSet_full[which(resultSet_full$Noise %in% selectedNoise),]

#debug(genes_x_transitions)
rs_list <- genes_x_transitions(resultSet, # dataframe w/ columns: ModelNo, SetName
                               topoName = topoName, # string
                               collectionName = expName, # string
                               initialClust = 1, # int
                               targetClust = 2, # int
                               wt_data = exprMat_norm[,genes], # matrix/dataframe of size numSamples x numFeatures
                               clust = clust, # vector of length numSamples containing integers
                               clust_all = clust_full,
                               tmpMeans = tmpMeans,
                               tmpSds = tmpSds
)


rsMatrix <- rs_list[[1]]
resultSet <- rs_list[[2]]
resultSet$ConversionPct <- resultSet$ConvertingModels / resultSet$startingInitPopulation
rownames(rsMatrix) <- resultSet$SetName[which(!is.na(resultSet$ConversionPct))]
summary(resultSet$ConversionPct)



# Plot 2D heatmap at different noise levels
eff_hmap_0.04 <- resultSet_full[which(resultSet_full$Noise == selectedNoise),
                                c("Species 1", "Species 2", "ConversionPct")]

eff_hmap_0.04_matrix <- matrix(nrow = 26, ncol = 26)
rownames(eff_hmap_0.04_matrix) <- genes_reordered
colnames(eff_hmap_0.04_matrix) <- genes_reordered
for(gene1 in genes_reordered) {
  for (gene2 in genes_reordered) {
    if(gene1 == gene2) {
      eff_val <- eff_hmap_0.04[which(eff_hmap_0.04$`Species 1` == gene1 & 
                                       is.na(eff_hmap_0.04$`Species 2`)), "ConversionPct"]
    } else if(length(which(eff_hmap_0.04$`Species 1` == gene1 & 
                           eff_hmap_0.04$`Species 2` == gene2)) == 1) {
      eff_val <- eff_hmap_0.04[which(eff_hmap_0.04$`Species 1` == gene1 & 
                                       eff_hmap_0.04$`Species 2` == gene2), "ConversionPct"]
    } else {
      eff_val <- eff_hmap_0.04[which(eff_hmap_0.04$`Species 2` == gene1 & 
                                       eff_hmap_0.04$`Species 1` == gene2), "ConversionPct"]
    }
    eff_hmap_0.04_matrix[gene1, gene2] <- eff_val
    eff_hmap_0.04_matrix[gene2, gene1] <- eff_val
  }
}


image(eff_hmap_0.04_matrix)
eff_hmap_0.04_matrix[lower.tri(eff_hmap_0.04_matrix)] <- NA


# Define the color mapping using viridis
col_fun <- colorRamp2(c(min(eff_hmap_0.04_matrix, na.rm = TRUE), 
                        mean(eff_hmap_0.04_matrix, na.rm = TRUE), 
                        max(eff_hmap_0.04_matrix, na.rm = TRUE)), 
                      viridis(3))  # Generates 3 colors from viridis

# Generate the heatmap
image <- Heatmap(eff_hmap_0.04_matrix, 
                 name = "% EMT", 
                 col = col_fun, 
                 cluster_rows = FALSE, 
                 cluster_columns = FALSE, 
                 show_row_names = TRUE, 
                 show_column_names = TRUE,
                 row_names_gp = gpar(fontsize=20),
                 column_names_gp = gpar(fontsize=20),
                 column_names_rot = 55,
                 column_names_side = "top",
                 row_title = "Gene 1",
                 row_title_gp = gpar(fontsize=24),
                 column_title = "Gene 2",
                 column_title_gp = gpar(fontsize=24),
                 column_title_side = "bottom")
image

plot_fname <- file.path(plotDir,paste0("figxx_2gene_sigEff_heatmap_noise=",paste0(selectedNoise,collapse = ","),".pdf"))
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()


######### Signal Characteristics #########
# Import network into igraph
g <- igraph::graph_from_edgelist(as.matrix(topo[,c("Source","Target")]))

w <- as.matrix(as_adjacency_matrix(g))
gene_id_list <- rownames(w)
colnames(w) <- 1:ncol(w)
rownames(w) <- 1:nrow(w)

for(i in rownames(resultSet)) {
  
  # Identify nodes involved in signal
  genesInSignal <-  c(resultSet[i,"Species 1"], resultSet[i,"Species 2"])
  genesInSignal <- genesInSignal[which(!is.na(genesInSignal))]
  
  genesInSignalIDs <- which(gene_id_list %in% genesInSignal)
  
  # Compute group betweenness centrality
  bet_cent <- kpcent(w, genesInSignalIDs, type="betweenness")
  
  
  # Compute group closeness centrality
  close_cent <- kpcent(w, genesInSignalIDs, type="closeness")
  
  
  # Add to resultSet
  resultSet[i,"GroupBetweenCentrality"] <- bet_cent
  resultSet[i,"GroupClosenessCentrality"] <- close_cent
  
}

if(!"GroupBetweenCentrality" %in% colnames(resultSet_full)) {
  for(i in rownames(resultSet_full)) {
    
    # Identify nodes involved in signal
    genesInSignal <-  c(resultSet_full[i,"Species 1"], resultSet_full[i,"Species 2"])
    genesInSignal <- genesInSignal[which(!is.na(genesInSignal))]
    
    genesInSignalIDs <- which(gene_id_list %in% genesInSignal)
    
    # Compute group betweenness centrality
    bet_cent <- kpcent(w, genesInSignalIDs, type="betweenness")
    
    
    # Compute group closeness centrality
    close_cent <- kpcent(w, genesInSignalIDs, type="closeness")
    
    
    # Add to resultSet
    resultSet_full[i,"GroupBetweenCentrality"] <- bet_cent
    resultSet_full[i,"GroupClosenessCentrality"] <- close_cent
    
  }
}




image <- ggplot(resultSet, aes(x=GroupBetweenCentrality, y=ConversionPct)) +
  geom_point(size=3) +
  theme_sticcc() +
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color="black"),
        axis.title = element_text(size=22)) +
  xlab("Signal Group Betweenness Centrality") +
  ylab("Models Undergoing EMT (%)")

plot_fname <- file.path(plotDir,paste0("eff_vs_betweenness_noise=",
                                       paste0(selectedNoise, collapse = ","),".pdf"))
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()


image <- ggplot(resultSet, aes(x=GroupClosenessCentrality, y=ConversionPct)) +
  geom_point(size=3) +
  theme_sticcc() +
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color="black"),
        axis.title = element_text(size=22)) +
  xlab("Signal Group Closeness Centrality") +
  ylab("Models Undergoing EMT (%)")

plot_fname <- file.path(plotDir,paste0("eff_vs_closeness_noise=",
                                       paste0(selectedNoise, collapse = ","),".pdf"))
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()


image <- ggplot(resultSet, aes(x=TotalOutDegree, y=ConversionPct)) +
  geom_point(size=3) +
  theme_sticcc() +
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color="black")) +
  xlab("Signal Total Out-Degree") +
  ylab("Models Undergoing EMT (%)")

plot_fname <- file.path(plotDir,paste0("eff_vs_outDegree_",
                                       paste0(selectedNoise, collapse = ","),".pdf"))
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()



######### Effect of Noise by Signal #########
resultSet_full$SigName_Short <- paste0(resultSet_full[,c("Species 1")], "_",
                                       resultSet_full[,c("Species 2")])


if(all(selectedNoise == c(0, 0.04, 0.2))) {
  signal_summary_df <- data.frame(Signal=unique(paste0(resultSet_full[,c("Species 1")], "_",
                                                       resultSet_full[,c("Species 2")])),
                                  ConversionPct_D0=NA,
                                  ConvertingModels_D0=NA,
                                  ConversionPct_D0.04=NA,
                                  ConvertingModels_D0.04=NA,
                                  ConversionPct_D0.2=NA,
                                  ConvertingModels_D0.2=NA)
  noises <- c(0, 0.04, 0.2)
  for(sig in signal_summary_df$Signal) {
    signal_summary_df[which(signal_summary_df$Signal==sig),"GroupClosenessCentrality"] <- 
      unique(resultSet_full[which(resultSet_full$SigName_Short == sig),"GroupClosenessCentrality"])
    signal_summary_df[which(signal_summary_df$Signal==sig),"GroupBetweenCentrality"] <- 
      unique(resultSet_full[which(resultSet_full$SigName_Short == sig),"GroupBetweenCentrality"])
    #signal_summary_df[which(signal_summary_df$Signal==sig),"OutDegree"] <- 
    #  unique(resultSet_full[which(resultSet_full$SigName_Short == sig),"TotalOutDegree"])
    for(noise in noises) {
      signal_summary_df[which(signal_summary_df$Signal==sig),paste0("ConversionPct_D",noise)] <- 
        resultSet_full[which(resultSet_full$Noise == noise & 
                               resultSet_full$SigName_Short == sig),"ConversionPct"]
      signal_summary_df[which(signal_summary_df$Signal==sig),paste0("ConvertingModels_D",noise)] <- 
        resultSet_full[which(resultSet_full$Noise == noise & 
                               resultSet_full$SigName_Short == sig),"ConvertingModels"]
    }
  }
  signal_summary_df$D0.04_Gain_Pct <- signal_summary_df$ConversionPct_D0.04 - signal_summary_df$ConversionPct_D0
  signal_summary_df$D0.04_Gain_Pct_ofRemaining <- (signal_summary_df$ConversionPct_D0.04 - signal_summary_df$ConversionPct_D0) / (1-signal_summary_df$ConversionPct_D0)
  
  ggplot(signal_summary_df) +
    geom_histogram(aes(x=ConvertingModels_D0.04-ConvertingModels_D0))
  
  ggplot(signal_summary_df) +
    geom_point(aes(x=GroupClosenessCentrality, y=D0.04_Gain_Pct))
  
  ggplot(signal_summary_df) +
    geom_point(aes(x=GroupBetweenCentrality, y=D0.04_Gain_Pct))
  
  ggplot(signal_summary_df) +
    geom_point(aes(x=GroupBetweenCentrality, y=D0.04_Gain_Pct_ofRemaining)) +
    theme_sticcc() +
    xlab("Betweenness Centrality") +
    ylab("Signal Efficacy, 0.04 vs 0")
  
  # ggplot(signal_summary_df) +
  #   geom_point(aes(x=OutDegree, y=D0.04_Gain_Pct_ofRemaining)) +
  #   theme_sticcc() +
  #   xlab("Out-Degree") +
  #   ylab("Signal Efficacy, 0.04 vs 0")
  
  
  # Heatmap showing gains from noise by signal
  
  eff_gains_0.04_matrix <- matrix(nrow = 26, ncol = 26)
  rownames(eff_gains_0.04_matrix) <- genes_reordered
  colnames(eff_gains_0.04_matrix) <- genes_reordered
  for(gene1 in genes_reordered) {
    for (gene2 in genes_reordered) {
      comb1 <- paste0(gene1, "_", gene2)
      comb2 <- paste0(gene2, "_", gene1)
      if(gene1 == gene2) {
        eff_val <- signal_summary_df[which(signal_summary_df$Signal == paste0(gene1,"_NA")), "D0.04_Gain_Pct_ofRemaining"]
      } else if(length(which(signal_summary_df$Signal == comb1)) == 1) {
        eff_val <- signal_summary_df[which(signal_summary_df$Signal == comb1), "D0.04_Gain_Pct_ofRemaining"]
      } else {
        eff_val <- signal_summary_df[which(signal_summary_df$Signal == comb2), "D0.04_Gain_Pct_ofRemaining"]
      }
      eff_gains_0.04_matrix[gene1, gene2] <- eff_val
      eff_gains_0.04_matrix[gene2, gene1] <- eff_val
    }
  }
  
  
  
  eff_gains_0.04_matrix[lower.tri(eff_gains_0.04_matrix)] <- NA
  
  # Define the color mapping
  library(circlize)
  library(viridis)
  # Define the color mapping using viridis
  col_fun <- colorRamp2(c(min(eff_gains_0.04_matrix, na.rm = TRUE), 
                          mean(eff_gains_0.04_matrix, na.rm = TRUE), 
                          max(eff_gains_0.04_matrix, na.rm = TRUE)), 
                        viridis(3))  # Generates 3 colors from viridis
  
  
  # Generate the heatmap
  image <- Heatmap(eff_gains_0.04_matrix, 
                   name = "Eff. Increase", 
                   col = col_fun, 
                   cluster_rows = FALSE, 
                   cluster_columns = FALSE, 
                   show_row_names = TRUE, 
                   show_column_names = TRUE,
                   row_names_gp = gpar(fontsize=20),
                   column_names_gp = gpar(fontsize=20),
                   column_names_rot = 55,
                   column_names_side = "top",
                   row_title = "Gene 1",
                   row_title_gp = gpar(fontsize=24),
                   column_title = "Gene 2",
                   column_title_gp = gpar(fontsize=24),
                   column_title_side = "bottom")
  image
  
  plot_fname <- file.path(plotDir,paste0("figS1c_2gene_sigEffFromNoise_heatmap_noise=",paste0(selectedNoise,collapse = ","),".pdf"))
  pdf(plot_fname, height = 10, width = 10)
  print(image)
  dev.off()
  
}


######## EFFICACY ANALYSIS ############
#debug(cells_x_signals)
cell_signal_df_fname <- file.path(topoDir,expName,
                                  paste0("cell_signal_df_noise=", 
                                         paste0(selectedNoise, collapse = ","),".Rds")) 
if(!file.exists(cell_signal_df_fname)) {
  #undebug(cells_x_signals)
  cell_signal_df <- cells_x_signals(paramSets = resultSet, # dataframe w/ columns: ModelNo, SetName
                                    topoName = topoName, # string
                                    collectionName = expName, # string
                                    initialClust = 1, # int
                                    targetClust = 2, # int
                                    wt_data = exprMat_norm[,1:nGenes], # matrix/dataframe of size numSamples x numFeatures
                                    clust = clust, # vector of length numSamples containing integers
                                    clust_all = clust_full,
                                    tmpMeans = tmpMeans,
                                    tmpSds = tmpSds,
                                    InitRawStates = unnormData,
                                    useDiffs = F)
  saveRDS(cell_signal_df, cell_signal_df_fname)
  
  
} else {
  cell_signal_df <- readRDS(cell_signal_df_fname)
  #cell_signal_df <- cell_signal_df[,which(resultSet_full$Noise %in% selectedNoise)]
}




numToPlot <- 30

numeric_df <- matrix(as.numeric(factor(as.matrix(cell_signal_df), 
                                       levels = c("Rebellious", "Target->Target", "Init->Init", "Init->Target"))),
                     nrow = nrow(cell_signal_df), ncol = ncol(cell_signal_df))
colnames(numeric_df) <- colnames(cell_signal_df)


# select one-gene signals
filterIdx <- order(resultSet$ConversionPct, decreasing = T)[1:numToPlot]
filter <- resultSet[filterIdx,"SetName"]

filter <- filter[which(filter %in% colnames(numeric_df))]
filterIdx <- filterIdx[which(filter %in% colnames(numeric_df))]
numeric_df <- numeric_df[,filter]


# Convert the categories to factors and then back to numerics to ensure unique coding
numeric_df_factor <- as.numeric(as.factor(numeric_df))

# Create a matrix from the factorized data
matrix_data <- matrix(numeric_df_factor, nrow = nrow(numeric_df), ncol = ncol(numeric_df))
colnames(matrix_data) <- gsub("_noise=0.04", "", colnames(cell_signal_df)[filterIdx])

# Define a color palette with as many colors as there are unique categories
# Replace with the actual number of categories and preferred colors
color_palette <- c("red", "lightblue", "pink", "darkblue") # adjust as needed
color_palette <- c("pink", "darkblue") 





# Perform hierarchical clustering on the columns (Conditions)
model_clusters <- as.factor(cutree(hclust(dist(matrix_data)), k=4))
ha_df <- data.frame(Cluster=model_clusters)

# Create an annotation object for the columns
# Color palette
cbPalette <- palette.colors(palette = "Okabe-Ito")[2:9]
names(cbPalette) <- c(1:8)
column_annotation <- HeatmapAnnotation(df = ha_df, 
                                       col=list(Cluster=c("1"=unname(cbPalette[1]),"2"=unname(cbPalette[2]),
                                                          "3"=unname(cbPalette[3]),"4"=unname(cbPalette[4]))))

# Define a color palette for the heatmap (same as before)
color_palette <- c("red", "lightblue", "pink", "darkblue") # adjust as needed
color_palette <- c("pink", "darkblue")

# Row annotation
colnames(resultSet)
ra_df <- resultSet[filterIdx,c("TotalOutDegree", "GroupBetweenCentrality")]
colnames(ra_df) <- c("OutDeg","Centrality")
#ra_df$Noise <- factor(ra_df$Noise)
row_annotation <- rowAnnotation(df = ra_df)

# Create the heatmap with annotation
image <- Heatmap(t(matrix_data), 
        name = "EMT Result", 
        col = colorRamp2(c(min(matrix_data):max(matrix_data)), color_palette), # using categorical colors
        #top_annotation = column_annotation,
        right_annotation = row_annotation,
        row_names_gp=gpar(fontsize=18),
        show_row_names = T,
        show_column_names = F,
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2")


hmap_fname <- file.path(plotDir,paste0("noise&sig_efficacyHmap_all_noise=",
                                       paste0(selectedNoise, collapse=","),"_top30.pdf"))
pdf(hmap_fname, height = 12, width = 12)
print(image)
dev.off()

dev.off()
lgd = Legend(labels = c("MET","M \U2192 M", "E \U2192 E", "EMT"), 
             title = "EMT Result", 
             legend_gp = gpar(fill = color_palette),
             labels_gp = gpar(fontsize = 16),
             title_gp = gpar(fontsize = 16, fontface="bold"),
             size=unit(3,"cm"))
cairo_pdf(file.path(plotDir,"hmap_legend_manual.pdf"), height=6, width=6, family="Arial")
image <- draw(lgd)
dev.off()


# Plot PCA colored by model clusters
image <- ggplot(pca_df, aes(x=PC1,y=PC2, color=model_clusters)) +
  scale_color_manual(values=cbPalette[1:4]) +
  geom_point(size=3) +
  guides(color=guide_legend(title = "Model Cluster")) + 
  theme_sticcc() +
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color="black"))

plot_fname <- file.path(plotDir,
                        paste0("pca_by_transitionLikelihood_noise=", 
                               paste0(selectedNoise, collapse = ","),".pdf"))
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()




######### Small efficacy hmap #########
### COMMENTED OUT FOR A REASON! Variables get reassigned here, breaking the flow. NB if running this section
# 
# numeric_df <- matrix(as.numeric(factor(as.matrix(cell_signal_df), 
#                                        levels = c("Rebellious", "Target->Target", "Init->Init", "Init->Target"))),
#                      nrow = nrow(cell_signal_df), ncol = ncol(cell_signal_df))
# colnames(numeric_df) <- colnames(cell_signal_df)
# 
# 
# #filterIdx <- order(resultSet$ConversionPct, decreasing = T)[1:50]
# filterIdx <- sample(1:351, 50)
# 
# filter <- resultSet[filterIdx,"SetName"]
# filter <- filter[which(filter %in% colnames(numeric_df))]
# filterIdx <- filterIdx[which(filter %in% colnames(numeric_df))]
# numeric_df <- numeric_df[,filter]
# 
# # Convert the categories to factors and then back to numerics to ensure unique coding
# numeric_df_factor <- as.numeric(as.factor(numeric_df))
# 
# # Create a matrix from the factorized data
# matrix_data <- matrix(numeric_df_factor, nrow = nrow(numeric_df), ncol = ncol(numeric_df))
# #colnames(matrix_data) <- colnames(cell_signal_df)[filterIdx] # old way, keeps the noise
# colnames(matrix_data) <- sub("_noise.*$", "", colnames(cell_signal_df)[filterIdx]) # new way, just the signal name
# 
# # Define a color palette with as many colors as there are unique categories
# # Replace with the actual number of categories and preferred colors
# color_palette <- c("red", "lightblue", "pink", "darkblue") # adjust as needed
# 
# library(circlize)
# library(ComplexHeatmap)
# 
# # Perform hierarchical clustering on the columns (Conditions)
# model_clusters <- as.factor(cutree(hclust(dist(matrix_data)), k=4))
# ha_df <- data.frame(Cluster=model_clusters)
# 
# # Create an annotation object for the columns
# cbPalette <- palette.colors(palette = "Okabe-Ito")[2:9]
# names(cbPalette) <- c(1:8)
# column_annotation <- HeatmapAnnotation(df = ha_df, 
#                                        col=list(Cluster=c("1"=unname(cbPalette[1]),"2"=unname(cbPalette[2]),
#                                                           "3"=unname(cbPalette[3]),"4"=unname(cbPalette[4]))))
# 
# # Define a color palette for the heatmap (same as before)
# color_palette <- c("red", "lightblue", "pink", "darkblue") # adjust as needed
# 
# # Create the heatmap with annotation
# image <- Heatmap(t(matrix_data), 
#                  name = "EMT Result", 
#                  col = colorRamp2(c(min(numeric_df):max(numeric_df)), color_palette), # using categorical colors
#                  top_annotation = column_annotation,
#                  row_names_gp=gpar(fontsize=12))
# 
# 
# hmap_fname <- file.path(plotDir,paste0("fig3a_noise&sig_efficacyHmap_50random_noise=",
#                                        paste0(selectedNoise, collapse=","),".pdf"))
# pdf(hmap_fname, height = 10, width = 10)
# print(image)
# dev.off()
# 
# #dev.off()
# #lgd = Legend(labels = c("MET","M > M", "E > E", "EMT"), title = "EMT Result", legend_gp = gpar(fill = color_palette))
# #image <- draw(lgd)





######### Critical Nodes #########

# Perform hierarchical clustering on the columns (Conditions)
model_clusters <- as.factor(cutree(hclust(dist(matrix_data)), k=4))
ha_df <- data.frame(Cluster=model_clusters)


# Row annotation
ra_df <- as.data.frame(resultSet[which(resultSet$ConversionPct > 0),
                   c("ConversionPct"#, "TotalOutDegree", "GroupClosenessCentrality",
                     #,"Species 1", "Species 2"
                     )])
colnames(ra_df) <- c("% EMT")
plotGene <- "Gsc"
#ra_df$GeneIndicator <- resultSet[which(resultSet$ConversionPct > 0),]$`Species 1` == plotGene | 
#  resultSet[which(resultSet$ConversionPct > 0),]$`Species 2` == plotGene
#row_annotation <- rowAnnotation(df = ra_df,
#                                col=list(GeneIndicator=c("TRUE"="blue", "FALSE"="grey")))
row_annotation <- rowAnnotation(df=ra_df,
                                annotation_name_gp=gpar(fontsize=18))

# optional: try ordering by species 1 & 2
new_order <- order(resultSet[, "Species 1"], resultSet[, "Species 2"])
keep_idx <- intersect(new_order, which(!is.na(rowMeans(rsMatrix))))
rsMatrix_ordered <- rsMatrix[keep_idx,]

rs_df_2 <- resultSet[keep_idx,
                     c("ConversionPct", "TotalOutDegree", "GroupClosenessCentrality", 
                       "Species 1")]
ra_ordered <- rowAnnotation(df = rs_df_2[,],
                            col=list(GeneIndicator=c("TRUE"="blue", "FALSE"="grey")))

library(ComplexHeatmap)
# Create the heatmap (with annotation?)
image <- Heatmap(#na.omit(rsMatrix_ordered),
                 na.omit(rsMatrix),
                 name = "Log FC, EMT vs E>E",
                 show_row_names = F,
                 row_names_gp=gpar(fontsize=4),
                 column_names_gp = gpar(fontsize=18),
                 right_annotation = row_annotation)
                 #cluster_rows = F,
                 #right_annotation = ra_ordered) 
image
##col = colorRamp2(c(min(numeric_df):max(numeric_df)), color_palette), # using categorical colors)


hmap_fname <- file.path(plotDir,paste0("fig3b_criticalNodes_noise=",
                                       paste0(selectedNoise, collapse=","),".pdf"))
pdf(hmap_fname, height = 10, width = 10)
print(image)
dev.off()




######### Critical Parameters #########
crit_param_df_fname <- file.path(topoDir,expName,paste0("crit_params_noise=",
                                                        paste0(selectedNoise, collapse=","),".Rds"))
if(!file.exists(crit_param_df_fname)) {
  crit_param_df <- params_x_transitions(resultSet, # dataframe w/ columns: ModelNo, SetName
                                        topoName = topoName, # string
                                        collectionName = expName, # string
                                        initialClust = 1, # int
                                        targetClust = 2, # int
                                        wt_data = exprMat[,genes], # matrix/dataframe of size numSamples x numFeatures
                                        clust = clust, # vector of length numSamples containing integers
                                        tmpMeans = tmpMeans,
                                        tmpSds = tmpSds
  )
  
  saveRDS(crit_param_df, crit_param_df_fname)
} else {
  crit_param_df <- readRDS(crit_param_df_fname)
}



ra_df <- resultSet[which(resultSet$ConversionPct > 0),
                   c("ConversionPct", "TotalOutDegree"#, "GroupClosenessCentrality"
                     )]
colnames(ra_df) <- c("% EMT", "OutDeg")
plotGene <- "Gsc"
#ra_df$GeneIndicator <- resultSet[which(resultSet$ConversionPct > 0),]$`Species 1` == plotGene | 
#  resultSet[which(resultSet$ConversionPct > 0),]$`Species 2` == plotGene
row_annotation <- rowAnnotation(df = ra_df,
                                annotation_name_gp=gpar(fontsize=18))


# Filter for significant parameters (and lose weak blocks?)
rownames(crit_param_df) <- resultSet$SetName
paramMeans <- colMeans(na.omit(crit_param_df))
minMag <- quantile(paramMeans, probs=seq(0,1,0.1))[8]


param_types <- unlist(lapply(colnames(crit_param_df[,which(paramMeans >= minMag)]), function(x) strsplit(x, "_")[[1]][1]))
ha_df <- data.frame(Type=factor(param_types))
column_annotation <- HeatmapAnnotation(df = ha_df)



image <- Heatmap(na.omit(crit_param_df[,which(paramMeans >= minMag)]), 
                 name = "FC, EMT vs E>E",
                 show_row_names = F,
                 show_column_names = F,
                 right_annotation = row_annotation,
                 top_annotation = column_annotation)
image

# image <- Heatmap(na.omit(crit_param_df[which(resultSet$ConversionPct > 0.6),which(paramMeans >= minMag)]), 
#                  name = "FC, EMT vs E>E",
#                  show_row_names = T,
#                  show_column_names = T)
#col = colorRamp2(c(min(numeric_df):max(numeric_df)), color_palette), # using categorical colors)


hmap_fname <- file.path(plotDir,paste0("fig3c_criticalParams",
                                       paste0(selectedNoise, collapse=","),".pdf"))
pdf(hmap_fname, height = 10, width = 10)
print(image)
dev.off()



######## SYNERGY ANALYSIS ############

# Compute total # transitions (colSums equal to Init->Target)
numTransitions_1gene <- apply(cell_signal_df[,which(resultSet$NumGenes == 1)], 2, 
                              function(x) sum(x == "Init->Target"))
allTransitions_1gene <- apply(cell_signal_df[,which(resultSet$NumGenes == 1)], 2, 
                              function(x) which(x == "Init->Target"))


# Initialize an n x n matrix of zeros
n <- length(numTransitions_1gene)
synergy_hmap_1gene <- matrix(0, nrow = n, ncol = n)

# Loop through the matrix to fill in the sums
for (i in 1:n) {
  for (j in 1:n) {
    #synergy_hmap_1gene[i, j] <- numTransitions_1gene[i] + numTransitions_1gene[j]
    synergy_hmap_1gene[i, j] <- length(union(unlist(allTransitions_1gene[i]), 
                                             unlist(allTransitions_1gene[j])))
  }
}

#synergy_hmap_1gene <- as.data.frame(replicate(72, numTransitions_1gene))
rownames(synergy_hmap_1gene) <- unlist(lapply(names(numTransitions_1gene), function(x) strsplit(x, "_")[[1]][1]))
image(synergy_hmap_1gene)

# Compute total # transitions (colSums equal to Init->Target)
numTransitions_2gene <- apply(cell_signal_df[,which(resultSet$NumGenes == 2)], 2, 
                              function(x) sum(x == "Init->Target"))
allTransitions_2gene <- apply(cell_signal_df[,which(resultSet$NumGenes == 2)], 2, 
                              function(x) which(x == "Init->Target"))

#genes <- unlist(lapply(top50$SetName, function(x) strsplit(x, "_")[[1]][1]))
synergy_hmap_2gene <- matrix(NA, nrow=length(genes), ncol=length(genes))
rownames(synergy_hmap_2gene) <- genes
colnames(synergy_hmap_2gene) <- genes

synergy_hmap_2gene <- synergy_hmap_2gene[rownames(synergy_hmap_1gene), rownames(synergy_hmap_1gene)]
synergy_hmap_2gene_norm <- synergy_hmap_2gene

# Parse individual signals
# sig_df_2gene <- data.frame(Signal=top50_2$SetName,
#                            Gene1=NA,
#                            Gene2=NA,
#                            NumTransitions=numTransitions_2gene)


oneGeneIdxList <- sub("_.*", "", names(numTransitions_1gene))
twoGeneIdx <- which(resultSet$NumGenes == 2)
for(i in seq_along(twoGeneIdx)) {
  subs <- strsplit(resultSet$SetName[twoGeneIdx[i]], "_")
  g1 <- subs[[1]][1]
  g2 <- subs[[1]][2]
  
  synergy_hmap_2gene[g1, g2] <- numTransitions_2gene[i]
  
  
  maxOneGeneEff <- max(numTransitions_1gene[which(oneGeneIdxList == g1)], 
                       numTransitions_1gene[which(oneGeneIdxList == g2)])
  synergy_hmap_2gene_norm[g1, g2] <- (numTransitions_2gene[i] - maxOneGeneEff) / (324-maxOneGeneEff)

    
    
  
  
}

# Show 72x72 matrix
image(synergy_hmap_2gene_norm)


diff_hmap <- synergy_hmap_2gene - synergy_hmap_1gene
image(diff_hmap)

fc_hmap <- synergy_hmap_2gene / synergy_hmap_1gene
image(fc_hmap)

# Melt the matrix into a long format
long_mat <- melt(diff_hmap, na.rm = FALSE)
long_mat$Var1 <- factor(long_mat$Var1, levels=genes)
long_mat$Var2 <- factor(long_mat$Var2, levels=genes)

# Create the heatmap
ggplot(long_mat, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red", na.value = "gray") +
  theme_minimal() +
  labs(fill = "Value", x = "Column", y = "Row") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Function to move values above the diagonal to below the diagonal
adjust_matrix <- function(m) {
  n <- nrow(m)
  c <- ncol(m)
  for (i in 1:n) {
    for (j in 1:c) {
      if (!is.na(m[i, j])) {
        m[j, i] <- m[i, j]
      } else if(!is.na(m[j,i])) {
        m[i, j] <- m[j, i]
      }
    }
  }
  m
}

# Apply the function to the matrix
#fc_hmap_adj <- adjust_matrix(diff_hmap)
fc_hmap_adj <- adjust_matrix(synergy_hmap_2gene_norm)

# Melt the matrix into a long format
long_fc_mat <- melt(fc_hmap_adj, na.rm = FALSE)
long_fc_mat[which(long_fc_mat$Var1 == long_fc_mat$Var2), "Diagonal"] <- T
long_fc_mat$Var1 <- factor(long_fc_mat$Var1, levels=genes)
long_fc_mat$Var2 <- factor(long_fc_mat$Var2, levels=genes)

# Compute the indices for the matrix positions and use these for comparison
long_fc_mat$RowIndex <- as.numeric(long_fc_mat$Var1)
long_fc_mat$ColIndex <- as.numeric(long_fc_mat$Var2)

# Add a condition to keep cells below the diagonal filled
long_fc_mat$keep <- long_fc_mat$RowIndex >= long_fc_mat$ColIndex


# Create the heatmap
image <- ggplot(long_fc_mat, aes(x = Var2, y = Var1, fill = value)) +
  #geom_tile(color = "white") +
  geom_tile(data = subset(long_fc_mat, keep), color = "white") +
  #
  scale_fill_gradient(low = "blue", high = "red", na.value = "gray") +
  geom_tile(data=long_fc_mat[which(long_fc_mat$Diagonal == T),], aes(x = Var2, y = Var1), fill = "black", color="white") +
  theme_minimal() +
  labs(fill = "Synergy Score", x = "Gene 1", y = "Gene 2") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20))

plot_fname <- file.path(plotDir,"fig3e_synergy_hmap_NORMALIZED_D=0.04.pdf")
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()


## Also plot absolute transition numbers
# Apply the function to the matrix
synergy_hmap_2gene_adj <- adjust_matrix(synergy_hmap_2gene)

# Melt the matrix into a long format
long_2gene_mat <- melt(synergy_hmap_2gene_adj, na.rm = FALSE)
long_2gene_mat[which(long_2gene_mat$Var1 == long_2gene_mat$Var2), "Diagonal"] <- T

# Compute the indices for the matrix positions and use these for comparison
long_2gene_mat$Var1 <- factor(long_2gene_mat$Var1, levels=genes)
long_2gene_mat$Var2 <- factor(long_2gene_mat$Var2, levels=genes)
long_2gene_mat$RowIndex <- as.numeric(long_2gene_mat$Var1)
long_2gene_mat$ColIndex <- as.numeric(long_2gene_mat$Var2)

# Add a condition to keep cells below the diagonal filled
long_2gene_mat$keep <- long_2gene_mat$RowIndex >= long_2gene_mat$ColIndex


# Create the heatmap
image <- ggplot(long_2gene_mat, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(data = subset(long_2gene_mat, keep), color = "white") +
  #geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red", na.value = "gray") +
  geom_tile(data=long_2gene_mat[which(long_2gene_mat$Diagonal == T),], aes(x = Var2, y = Var1), fill = "black", color="white") +
  theme_minimal() +
  labs(fill = "Sig. Eff.", x = "Gene 1", y = "Gene 2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size=12))

plot_fname <- file.path(plotDir,"fig3e_synergy_hmap_counts_D=0.04.pdf")
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()


# Plot example synergistic conditions
plotSets_comb <- resultSet[which.max(resultSet$EndTargetStates),"SetName"]

plotCondition(racipe = racipe_bistable_raw,
              pca = pca,
              clust = clust,
              initialClust = initClust,
              targetClust = tgtClust,
              expName = expName,
              setName = plotSets_comb[1],
              suffix="_01052025_synergy_ex_OUnoise",
              outDir = file.path(topoDir))





########## RANKING ANALYSIS ############


## Correlation of 1-gene signals
sigEffs_1gene_racipe <- resultSet[which(resultSet$NumGenes == 1 & resultSet$Noise == selectedNoise), c("Species 1", "ConversionPct")]
# manually putting in signal data from spin model
sigEffs_1gene_spin <- data.frame(Species=sigEffs_1gene_racipe$`Species 1`, ConversionPct=NA)
sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "Cdh1"),"ConversionPct"] <- 0.65
sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "miR200b"),"ConversionPct"] <- 0.45
sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "miR200c"),"ConversionPct"] <- 0.25
sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "miR34a"),"ConversionPct"] <- 0.30
sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "Ovol2"),"ConversionPct"] <- 0.45
sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "Grhl2"),"ConversionPct"] <- 0.05
sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "Foxc2"),"ConversionPct"] <- 0.35
sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "Zeb1"),"ConversionPct"] <- 1
sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "Zeb2"),"ConversionPct"] <- 1
sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "Snai1"),"ConversionPct"] <- 1
sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "Vim"),"ConversionPct"] <- 0.30
sigEffs_1gene_spin[which(sigEffs_1gene_spin$Species == "Twist1"),"ConversionPct"] <- 0.85


sigEffs_1gene_racipe$ConversionPct_Spin <- sigEffs_1gene_spin$ConversionPct

ggplot() +
  geom_point(data=sigEffs_1gene_racipe, aes(x=1:26, y=ConversionPct), color="red") +
  geom_point(data=sigEffs_1gene_spin, aes(x=1:26, y=ConversionPct), color="blue")

diff_df <- data.frame(Signal=sigEffs_1gene_racipe$`Species 1`,
                      Diff=sigEffs_1gene_racipe$ConversionPct_Spin - sigEffs_1gene_racipe$ConversionPct)
ggplot(data=diff_df[which(!is.na(diff_df$Diff)),]) +
  geom_bar(stat="identity", aes(x=Signal, fill=Signal, y=Diff)) +
  ylab("Conversion %, Spin - RACIPE") +
  xlab("Signal") +
  theme_sticcc() +
  theme(axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle=90))

ggplot() +
  geom_point(data=sigEffs_1gene_racipe, aes(x=ConversionPct, y=ConversionPct_Spin), size=3) +
  xlab("RACIPE Conversion %") +
  ylab("Spin Conversion %") +
  theme_sticcc() +
  theme(axis.line = element_line(color = "black"))

indices <- which(!is.na(sigEffs_1gene_racipe$ConversionPct_Spin))
cor(sigEffs_1gene_racipe$ConversionPct[indices], sigEffs_1gene_racipe$ConversionPct_Spin[indices], 
    use = "complete.obs", method = "spearman")



## Correlation of 2-gene signals
sigEffs_2gene_racipe <- resultSet[which(resultSet$NumGenes == 2 & resultSet$Noise == selectedNoise), c("Species 1", "Species 2", "ConversionPct")]
sigEffs_2gene_spin <- read.table(file.path(dataDir, "spin_2node_effs.dat"), sep = " ")
rownames(sigEffs_2gene_spin) <- genes_reordered
colnames(sigEffs_2gene_spin) <- genes_reordered

sigEffs_2gene_racipe$ConversionPct_Spin <- NA

find_index <- function(df, set) {
  for(row in rownames(df)) {
    rowGeneSet <- c(df[row,"Species 1"], df[row, "Species 2"])
    if(length(setdiff(rowGeneSet, set)) == 0) {
      return(row)
    }
  }
}

for(gene1 in genes_reordered) {
  for(gene2 in genes_reordered) {
    spinEff <- sigEffs_2gene_spin[gene1, gene2]
    idx <- find_index(sigEffs_2gene_racipe, c(gene1, gene2))
    
    sigEffs_2gene_racipe[idx, "ConversionPct_Spin"] <- spinEff
  }
}


diff_df <- data.frame(Signal=paste0(sigEffs_2gene_racipe$`Species 1`,"_",sigEffs_2gene_racipe$`Species 2`),
                      Rank=order(sigEffs_2gene_racipe$ConversionPct),
                      Diff=sigEffs_2gene_racipe$ConversionPct_Spin - sigEffs_2gene_racipe$ConversionPct)
ggplot(data=diff_df[which(!is.na(diff_df$Diff)),]) +
  geom_point(aes(x=Rank, y=Diff), size=3) +
  ylab("Conversion %, Spin - RACIPE") +
  xlab("Signal Rank") +
  theme_sticcc() +
  theme(axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle=90))


ggplot() +
  geom_point(data=sigEffs_2gene_racipe, aes(x=ConversionPct, y=ConversionPct_Spin), size=3) +
  xlab("RACIPE Conversion %") +
  ylab("Spin Conversion %") +
  theme_sticcc() +
  theme(axis.line = element_line(color = "black"))


## Now combine with info from 5B (betweenness centrality)
sigEffs_2gene_racipe$Signal <- paste0(sigEffs_2gene_racipe$`Species 1`, "_", sigEffs_2gene_racipe$`Species 2`)
signal_summary_compare <- merge(sigEffs_2gene_racipe, signal_summary_df, by="Signal")
image <- ggplot() +
  geom_point(data=signal_summary_compare, aes(x=ConversionPct, y=ConversionPct_Spin, color=GroupBetweenCentrality), size=3) +
  xlab("RACIPE Conversion %") +
  ylab("Spin Conversion %") +
  scale_color_gradient(name="Betweenness") +
  theme_sticcc() +
  theme(axis.line = element_line(color = "black"))
plot_fname <- file.path(plotDir,paste0("fig4b_spin_vs_racipe_betweenness_noise=",selectedNoise,".pdf"))
pdf(plot_fname, width = 10, height = 10)
print(image)
dev.off()



cor(signal_summary_compare$ConversionPct, signal_summary_compare$ConversionPct_Spin,
    use = "complete.obs", method = "spearman")
cor(signal_summary_compare$ConversionPct, signal_summary_compare$GroupBetweenCentrality,
    use = "complete.obs", method = "spearman")
cor(signal_summary_compare$ConversionPct_Spin, signal_summary_compare$GroupBetweenCentrality,
    use = "complete.obs", method = "spearman")



indices <- which(!is.na(sigEffs_2gene_racipe$ConversionPct_Spin))
cor(sigEffs_2gene_racipe$ConversionPct[indices], sigEffs_2gene_racipe$ConversionPct_Spin[indices], 
    use = "complete.obs", method = "spearman")
cor(sigEffs_2gene_racipe$ConversionPct[indices], sigEffs_2gene_racipe$ConversionPct_Spin[indices], 
    use = "complete.obs", method = "spearman")


# try removing signals where sigEff_spin is one
# sigEffs_2gene_racipe <- sigEffs_2gene_racipe[which(sigEffs_2gene_racipe$ConversionPct_Spin != 1),]
# cor(sigEffs_2gene_racipe$ConversionPct, sigEffs_2gene_racipe$ConversionPct_Spin, 
#     use = "complete.obs", method = "spearman")

## Bland-altman plot
ranks_racipe <- rank(sigEffs_2gene_racipe$ConversionPct)
ranks_spin <- rank(sigEffs_2gene_racipe$ConversionPct_Spin)

plot(ranks_racipe, ranks_spin)

# Calculate means and differences
means <- (ranks_racipe + ranks_spin) / 2
differences <- ranks_racipe - ranks_spin

spinEffOne <- as.factor(sigEffs_2gene_racipe$ConversionPct_Spin == 1)

# Create a Bland-Altman Plot
bland_altman_plot <- ggplot(data = NULL, aes(x = means, y = differences, color=spinEffOne)) +
  geom_point() +  # Add points
  geom_hline(yintercept = mean(differences), linetype = "dashed", color = "red") +  # Add mean line
  geom_hline(yintercept = mean(differences) + 1.96 * sd(differences), linetype = "dashed", color = "blue") +  # Add upper limit
  geom_hline(yintercept = mean(differences) - 1.96 * sd(differences), linetype = "dashed", color = "blue") +  # Add lower limit
  labs(x = "Average Rank", y = "Difference in Rank", title = "Bland-Altman Plot for Ranking Comparisons") +
  theme_minimal()

# Display the plot
print(bland_altman_plot)

#### Percentile plot has been cut, too confusing :(
# # Function to convert ranks to percentile ranks
# convert_to_percentile <- function(ranks) {
#   # Obtain the total number of items
#   total_items <- length(ranks)
#   
#   # Calculate percentile ranks
#   percentile_ranks <- (rank(ranks, ties.method = "min") - 1) / total_items * 100
#   
#   return(percentile_ranks)
# }
# ranks_racipe_pct <- convert_to_percentile(ranks_racipe)
# ranks_spin_pct <- convert_to_percentile(ranks_spin)
# plot(ranks_racipe_pct, ranks_spin_pct)
# cor(ranks_racipe_pct, ranks_spin_pct)
# 
# ggplot() +
#   geom_point(aes(x=ranks_racipe_pct, y=ranks_spin_pct), size=3) +
#   xlab("RACIPE Conversion %") +
#   ylab("Spin Conversion %") +
#   theme_sticcc() +
#   theme(axis.line = element_line(color = "black"))



######## ZEB1 BOOLEAN EFFICACY #####

sigEffs_2gene_spin_long <- sigEffs_2gene_spin
sigEffs_2gene_spin_long$gene1 <- rownames(sigEffs_2gene_spin_long)

sigEffs_2gene_spin_long <- melt(sigEffs_2gene_spin_long, id.vars = "gene1", 
                                variable.name = "gene2", value.name = "efficacy")
sigEffs_2gene_spin_long$gene2 <- as.character(sigEffs_2gene_spin_long$gene2)
unique_sets <- c()
keep_idx <- c()
for(row in 1:nrow(sigEffs_2gene_spin_long)) {
  gene1 <- sigEffs_2gene_spin_long[row, "gene1"]
  gene2 <- sigEffs_2gene_spin_long[row, "gene2"]
  geneset <- paste0(sort(c(gene1, gene2)), collapse = "_")
  
  if(geneset %in% unique_sets) {
    next
  } else {
    unique_sets <- c(unique_sets, geneset)
    keep_idx <- c(keep_idx, row)
  }
}

sigEffs_2gene_spin_long <- sigEffs_2gene_spin_long[keep_idx,]


## Stats for Zeb1-related signals
summary(sigEffs_2gene_spin_long[which(sigEffs_2gene_spin_long$gene1 == "Zeb1" |
                                        sigEffs_2gene_spin_long$gene2 == "Zeb1"), "efficacy"])
sigEffs_2gene_spin_long[which(sigEffs_2gene_spin_long$gene1 == "Zeb1" |
                                sigEffs_2gene_spin_long$gene2 == "Zeb1"),]

## Stats for Snai1-related signals
summary(sigEffs_2gene_spin_long[which(sigEffs_2gene_spin_long$gene1 == "Snai1" |
                                        sigEffs_2gene_spin_long$gene2 == "Snai1"), "efficacy"])
sigEffs_2gene_spin_long[which(sigEffs_2gene_spin_long$gene1 == "Snai1" |
                                        sigEffs_2gene_spin_long$gene2 == "Snai1"),]


######## TIMING SIMULATIONS IN PARALLEL #####

# Select signals to use
# top n by conversion pct,  and n more random
# or read old signals & reuse those
sigNames <- c(readRDS(file.path(dataDir,paste0("transitions_vs_time_signalNames.Rds"))), "Zeb1_noise=0.04")
detSigNames <- gsub("0.04", "0", sigNames)


sigNames <- c("Zeb1_noise=0.04")

setIDList <- c(which(resultSet$SetName %in% sigNames), which(resultSet$SetName %in% detSigNames))
setIDList <- c(which(resultSet$SetName %in% sigNames))
n <- 10
# setIDList <- c()
# 
# rankorder <- order(resultSet$ConversionPct, decreasing = T)
# topn <- rankorder[which(!rankorder %in% setIDList)][1:n]
# setIDList <- c(setIDList, topn)
# 
# randn <- sample(setdiff(which(resultSet$Simulated), setIDList), n)
# setIDList <- c(setIDList, randn)


nCores <- 4
requireNamespace("doFuture")
registerDoFuture()
plan(multisession, workers = nCores)

## CHECK THESE EVERY TIME
#noise <- 0.04
#resultSet$Noise <- rep(noise, nrow(resultSet))
resultSet$Tau <- rep(signal_tcorr, nrow(resultSet))
resultSet$ScaledNoise <- rep(T, nrow(resultSet))



trajSimTime <- 300
expName_new_list <- paste0("bhtopo_t=",trajSimTime,
                           "_relax_OUnoise=",0.04,
                           "_tau=",signal_tcorr,
                           "_genes=",paste0(signal_nGenes,collapse = "."),
                           "_SIG=",resultSet[setIDList, "SetName"])
times <- c(seq(2, 30, 2), seq(35, 100, 5), seq(120, trajSimTime, 20))
#times <- c(seq(10, 100, 10), seq(120, 500, 20), seq(550, 1000, 50))



for(setIDNo in 1:length(setIDList)) {
  setID = setIDList[setIDNo]
  setName = resultSet[setID, "SetName"]
  expName_new = expName_new_list[setIDNo]
  #debug(calcTransitionRate)
  sampleSet_times <- calcTransitionRate(paramSets = resultSet,
                                             setID = setID,
                                             racipe = racipe_bistable_raw, # bistable or bistable_raw?
                                             pca = pca,
                                             clust = clust,
                                             initialClust = 1,
                                             targetClust = 2,
                                             sigName = setName,
                                             outDir = file.path(topoDir),
                                             expName = expName_new,
                                             plot=F,
                                             noise = signal_noise,
                                             forceRerun = F,
                                             forceRecompute = F,
                                             anneal=F,
                                             relax=T,
                                             simTimes = times,
                                             simTimeRelax = 10,
                                             save=T,
                                             clamp=T) 
}


# x <- foreach(setIDNo = as.list(c(1:length(setIDList))),
#              setID = setIDList,
#              setName = resultSet[setIDList, "SetName"],
#              expName_new = expName_new_list,
#              .export = c("resultSet","racipe_bistable_raw","pca","clust","topoDir","signal_noise", 
#                          "times")) %dorng% {
#                            
#                            sampleSet_times <- calcTransitionRate(paramSets = resultSet,
#                                                                  setID = setID,
#                                                                  racipe = racipe_bistable_raw, # bistable or bistable_raw?
#                                                                  pca = pca,
#                                                                  clust = clust,
#                                                                  initialClust = 1,
#                                                                  targetClust = 2,
#                                                                  sigName = setName,
#                                                                  outDir = file.path(topoDir),
#                                                                  expName = expName_new,
#                                                                  plot=F,
#                                                                  noise = signal_noise,
#                                                                  forceRerun = F,
#                                                                  forceRecompute = F,
#                                                                  anneal=F,
#                                                                  relax=T,
#                                                                  simTimes = times,
#                                                                  simTimeRelax = 10,
#                                                                  save=T) 
#                            
#                          } # end parallel loop

multiSet_fname <- file.path(dataDir,
                            paste0("transitions_vs_time_20SignalsManual_DetVStoch",expName,".Rds"))
if(!file.exists(multiSet_fname)) {
  multiSet_times <- list()
  for(setIDNo in 1:length(setIDList)) {
    setID <- setIDList[setIDNo]
    expDir <- file.path(topoDir, expName_new_list[setIDNo])
    sampleSet_times <- readRDS(file.path(expDir, 
                                         paste0("transitionTimes_",resultSet[setID,"SetName"],".Rds")))
    multiSet_times[[setIDNo]] <- sampleSet_times
  }
  saveRDS(multiSet_times, multiSet_fname)
  sigNames <- resultSet[setIDList, "SetName"]
  saveRDS(sigNames, file.path(dataDir,paste0("transitions_vs_time_signalNames_DetVStoch.Rds")))
  saveRDS(setIDList, file.path(dataDir,paste0("transitions_vs_time_signalIDs_DetVStoch.Rds")))
} else {
  multiSet_times <- readRDS(multiSet_fname)
  sigNames <- readRDS(file.path(dataDir,paste0("transitions_vs_time_signalNames_DetVStoch.Rds")))
  setIDList <- readRDS(file.path(dataDir,paste0("transitions_vs_time_signalIDs_DetVStoch.Rds")))
}




## Plot conversions vs time
plot_df_list <- list()
idx <- 1
for(setIDNo in seq_along(setIDList)) {
  setID <- setIDList[setIDNo]
  setName <- resultSet[setID, "SetName"]
  sampleSet_times <- multiSet_times[[setIDNo]]
  
  for(timeID in seq_along(times)) {
    time <- times[timeID]
    newrow <- list(Signal=setName, Time=time, Conversions=sampleSet_times[[timeID]][[1]])
    plot_df_list[[idx]] <- newrow
    idx <- idx+1
  }
  
}

plot_df <- do.call(rbind, lapply(plot_df_list, function(x) as.data.frame(t(x))))
plot_df$Time <- as.numeric(plot_df$Time)
plot_df$Conversions <- as.numeric(plot_df$Conversions)
plot_df$Signal <- as.character(plot_df$Signal)

image <- ggplot(plot_df[which(plot_df$Signal %in% c("Zeb1_noise=0.04", "Zeb1_noise=0")),], aes(x=Time, y=Conversions, color=Signal)) +
  geom_point() +
  geom_path() +
  ggtitle(paste0(expName,"_20ManualSignals"))

image






######## TRANSIT TIME DISTRIBUTION ########
library(MASS)


sigSpeciesList <- unique(paste0(resultSet[setIDList,c("Species 1")], "_",
                         resultSet[setIDList,c("Species 2")]))

initial_e_models <- which(clust==1)
transit_time_df <- data.frame(Model=rep(initial_e_models, length(sigSpeciesList)), 
                              Signal=rep(sigSpeciesList, each=length(initial_e_models)), 
                              TransitType="None", 
                              TransitTime_Det=NA, TransitTime_Stoch_0.04=NA)
#det_times <- times
#stoch_times <- times

for(sigIdx in seq_along(setIDList)) {
  # get signal data
  sigNo <- setIDList[sigIdx]
  sigSpecies1 <- resultSet[sigNo,"Species 1"]
  sigSpecies2 <- resultSet[sigNo,"Species 2"]
  sigSpecies <- paste0(sigSpecies1, "_", sigSpecies2)
  
  sigName <- resultSet[sigNo,"SetName"]
  sigNoise <- resultSet[sigNo,"Noise"]
  
  timing_results_sig <- multiSet_times[[sigIdx]]
  
  # Loop over det results
  if(sigNoise == 0) {
    for(tID in seq_along(times)) {
      t <- times[tID]
      unsolved_models <- initial_e_models[which(is.na(
        transit_time_df[which(transit_time_df$Signal == sigSpecies), "TransitTime_Det"]))]
      for(model in unsolved_models) {
        if(timing_results_sig[[tID]]$NewClusters[model] == tgtClust) {
          transit_time_df[which(transit_time_df$Model == model &
                                  transit_time_df$Signal == sigSpecies),"TransitTime_Det"] <- t
          #transit_time_df[which(transit_time_df$Model == model &
          #                        transit_time_df$Signal == sigSpecies),"TransitType"] <- "Det"
          
        }
      }
    }
  } else if(sigNoise == 0.04) {
    # Loop over stoch results
    for(tID in seq_along(times)) {
      t <- times[tID]
      unsolved_models <- initial_e_models[which(is.na(
        transit_time_df[which(transit_time_df$Signal == sigSpecies), "TransitTime_Stoch_0.04"]))]
      for(model in unsolved_models) {
        if(timing_results_sig[[tID]][[2]][model] == tgtClust) {
          transit_time_df[which(transit_time_df$Model == model),"TransitTime_Stoch_0.04"] <- t
          #if(transit_time_df[which(transit_time_df$Model == model &
          #                         transit_time_df$Signal == sigSpecies),"TransitType"] == "None") {
            #transit_time_df[which(transit_time_df$Model == model &
            #                        transit_time_df$Signal == sigSpecies),"TransitType"] <- "Stoch"
          #}
          
        }
        
      }
    }
  }

  

  
  
}




transit_time_df_name <- file.path(dataDir, 
                                  "det_vs_stoch_transit_times_20SignalsManual.Rds")
#saveRDS(transit_time_df, transit_time_df_name)


#transit_time_df <- readRDS(transit_time_df_name)
hist(transit_time_df$TransitTime_Det)
hist(transit_time_df$TransitTime_Stoch_0.04[which(is.na(transit_time_df$TransitTime_Det))])


# estimate the parameters
nonNA_stoch_times <- transit_time_df$TransitTime_Stoch[which(!is.na(transit_time_df$TransitTime_Stoch_0.04) &
                                                               is.na(transit_time_df$TransitTime_Det) &
                                                               transit_time_df$TransitTime_Stoch_0.04 < 100)]
fit1 <- fitdistr(nonNA_stoch_times, "exponential") 

# goodness of fit test
ks.test(nonNA_stoch_times, "pexp", fit1$estimate) # p-value > 0.05 -> distribution not refused

# plot a graph
hist(nonNA_stoch_times, freq = FALSE, xlim = c(0, quantile(nonNA_stoch_times, 0.99)))
curve(dexp(x, rate = fit1$estimate), from = 0, col = "red", add = TRUE)
#curve(dpois(x, lambda = fit1$estimate), from = 0, col = "red", add = TRUE)

## histogram overlaying transit times
hist_df <- data.frame(TransitTime_Det = transit_time_df$TransitTime_Det,
                      TransitTime_Stoch = transit_time_df$TransitTime_Stoch_0.04,
                      Type=NA,
                      TransitTime = NA)
hist_df$TransitTime[which(!is.na(hist_df$TransitTime_Det))] <- 
  hist_df[which(!is.na(hist_df$TransitTime_Det)),"TransitTime_Det"]
hist_df[which(!is.na(hist_df$TransitTime_Det)),"Type"] <- "Deterministic"

hist_df$TransitTime[which(!is.na(hist_df$TransitTime_Stoch) & is.na(hist_df$TransitTime_Det))] <- 
  hist_df[which(!is.na(hist_df$TransitTime_Stoch) & is.na(hist_df$TransitTime_Det)),"TransitTime_Stoch"]
hist_df[which(!is.na(hist_df$TransitTime_Stoch) & is.na(hist_df$TransitTime_Det)),"Type"] <- "D=0.04"


#topn_short_names <- gsub(resultSet[topn,"SetName"])

image <- ggplot(hist_df[which(hist_df$TransitTime < 60 & transit_time_df$Signal %in% sigSpeciesList[1:10]),]) +
  #geom_histogram(aes(x=TransitTime, fill=Type), alpha=0.5, bins=20) +
  geom_boxplot(aes(x=Type, y=TransitTime, fill=Type)) +
  ylab("Transition Time") +
  xlab("Noise") +
  theme_sticcc() + 
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color="black"),
        axis.text.x = element_text(angle = 30, vjust = 0.6))
image

plot_fname <- file.path(plotDir,"det_vs_stoch_times_1024.pdf")
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()



## Plot just one signal
image <- ggplot(hist_df[which(hist_df$TransitTime < 60 & transit_time_df$Signal %in% c("miR200b_Twist2")),]) +
  #geom_histogram(aes(x=TransitTime, fill=Type), alpha=0.5, bins=20) +
  geom_boxplot(aes(x=Type, y=TransitTime, fill=Type)) +
  ylab("Transition Time") +
  xlab("Noise") +
  theme_sticcc() + 
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color="black"),
        axis.text.x = element_text(angle = 30, vjust = 0.6))
image



## Transit times only if deterministic transition also possible
hist_df_2 <- pivot_longer(transit_time_df[which(#!is.na(transit_time_df$TransitTime_Det) &
                                                  transit_time_df$Signal %in% sigSpeciesList[1:20]),
                                          c("Model","Signal",
                                             "TransitTime_Det","TransitTime_Stoch_0.04")], 
                          cols = all_of(c("TransitTime_Det","TransitTime_Stoch_0.04")), 
                          values_to = "TransitTime",
                          names_to = "Type")
image <- ggplot(hist_df_2[which(hist_df_2$TransitTime < 60),]) +
  #geom_histogram(aes(x=TransitTime, fill=Type), alpha=0.5, bins=20) +
  geom_violin(aes(x=Type, y=TransitTime, fill=Type)) +
  ylab("Transition Time") +
  xlab("Noise") +
  theme_sticcc() + 
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color="black"),
        axis.text.x = element_text(angle = 30, vjust = 0.6))
image





######## AVERAGE TIME TRAJECTORY ############
# Select a signal with a fairly high efficacy
# Collect average expression over all models at each timepoint
#avg_traj_sigID <- which(resultSet$ConversionPct == max(resultSet$ConversionPct))[1]
avg_traj_sigID <- which(resultSet$SetName == "Zeb1_noise=0.04")
avg_traj_setName <- resultSet[avg_traj_sigID,"SetName"]
avg_traj_expName <- paste0("bhtopo_t=",trajSimTime,
                           "_relax_OUnoise=",0.04,
                           "_tau=",signal_tcorr,
                           "_genes=",paste0(signal_nGenes,collapse = "."),
                           "_SIG=",resultSet[avg_traj_sigID, "SetName"])
avg_traj_expDir <- file.path(topoDir, avg_traj_expName)

avg_traj_useCells <- which(clust == 1)

tp_data_all <- data.frame(matrix(NA, nrow=0, ncol=length(genes)+2))
colnames(tp_data_all) <- c("Time", "Model", genes)


ics <- as.data.frame(t(assay(racipe_bistable_raw)))[avg_traj_useCells,]
ics$Time <- 0
ics$Model <- avg_traj_useCells


tp_data_all <- rbind(tp_data_all, ics[,c("Time", "Model", genes)])

for(tp in times) {
  
  fname_checkpoint <- file.path(avg_traj_expDir, paste0("racipe_",avg_traj_setName,"_t=",tp,".Rds"))
  racipe_checkpoint <- readRDS(fname_checkpoint)
  tp_data <- as.data.frame(t(assay(racipe_checkpoint)))[avg_traj_useCells,]
  tp_data$Time <- tp
  tp_data$Model <- avg_traj_useCells
  tp_data_all <- rbind(tp_data_all, tp_data[,c("Time", "Model", genes)])
}


tp_data_agg <- tp_data_all %>%
  group_by(Time) %>%
  summarise(across(all_of(c(genes)), mean, na.rm = TRUE))


tp_data_agg_long <- pivot_longer(tp_data_agg[which(tp_data_agg$Time < 25),], cols=c(genes), values_to = "Expression", names_to = "Gene")
ggplot(tp_data_agg_long) +
  geom_line(aes(x=Time,y=Expression,color=Gene))

saveRDS(tp_data_all, file.path(avg_traj_expDir,"Zeb1_fullTrajectory.Rds"))
saveRDS(tp_data_agg, file.path(avg_traj_expDir,"Zeb1_aggTrajectory.Rds"))
saveRDS(tp_data_agg_long, file.path(avg_traj_expDir,"Zeb1_longAggTrajectory.Rds"))



ggplot(tp_data_agg_long[which(tp_data_agg_long$Gene %in% c("Zeb1","Zeb2","Snai1","Snai2","Twist1","Twist2","Cdh1")),]) +
  geom_line(aes(x=Time,y=Expression,color=Gene))

######## TEST ONE MODEL AT A TIME ############

# Select signal
which(resultSet$SetName == "Zeb1_noise=0.04")#"miR34a_Zeb1_noise=0.04")
hd_traj_sigID <- 114#57#14#setID
hd_traj_sig <- resultSet[hd_traj_sigID, "SetName"]

# Identify transitions in ensemble case for selected signal
if(hd_traj_sig %in% sigNames) {
  times_idx <- which(sigNames == hd_traj_sig)
  sigTimes <- multiSet_times[[times_idx]]
  transition_time_df <- as.data.frame(matrix(data = NA, 
                                                        nrow = dim(racipe_bistable_raw)[2], 
                                                        ncol = length(times)+1))
  colnames(transition_time_df) <- c("Model",times)
  transition_time_df$Model <- 1:nrow(transition_time_df)
  for(time in seq_along(times)) {
    transition_time_df[,as.character(times[time])] <- sigTimes[[time]]$NewClusters
  }
  
}

transition_time_summary <- data.frame(Model=transition_time_df$Model,
                                      InitClust=clust,
                                      FinalClust=transition_time_df$`300`,
                                      EMTTime=apply(transition_time_df, 1, 
                                                    function(x) colnames(transition_time_df)[(which(x == 2)[1])]))


# Select model
#multiSet_times
emt_models <- transition_time_df[which(transition_time_df$`2`==1 & transition_time_df$`300`==2),"Model"]
hd_traj_model <- 4
#hd_traj_model <- sample(which(clust==1), 1)

# Sim parameters
hd_traj_startTime <- 2
hd_traj_simTime <- 100
hd_traj_resolution <- 0.2
set.seed(1234)

# Load checkpoint
#hd_traj_checkpoint <- readRDS(file.path(topoDir,expName_hd_traj,paste0("racipe_",hd_traj_sig,"_t=",hd_traj_startTime,".Rds")))
hd_traj_checkpoint <- racipe_bistable_raw # actually just use WT object
dim(hd_traj_checkpoint)

# Filter selected model & set ICs (to WT steady state)
hd_traj_p2 <- hd_traj_checkpoint[,hd_traj_model]
hd_traj_p2@metadata$config$simParams["numModels"] <- 1
sracipeIC(hd_traj_p2) <- assay(racipe_bistable_raw)[,as.numeric(hd_traj_model)] 



## TODO
# QC simulation (ensure steady state is the same, simulate for 10 @ 0.1)
hd_traj_p2@metadata$config$simParams["nIC"] <- 1
hd_traj_p2 <- sracipeSimulate(hd_traj_p2, timeSeries = T, printStart = 0, 
                              printInterval = 0.1, simulationTime = 10,
                              genParams = F, genIC = F, integrateStepSize = 0.1, initialNoise=0.0, scaledNoise = T, 
                              stepper = "EM", ouNoise_t = signal_tcorr)


qc_ts <- unlist(as.data.frame(hd_traj_p2@metadata$timeSeriesDet)[,1])
all.equal(sracipeIC(hd_traj_p2), qc_ts)


# Implement signal
sig_clamp_genes <- getClampGenes(resultSet, hd_traj_sigID)
sig_clamp_gene_ids <- unlist(lapply(sig_clamp_genes, function(x) which(rownames(racipe) == x)))
sig_clamp_df <- getClampDF(clamp_df, sig_clamp_genes, 2)
colnames(sig_clamp_df) <- as.numeric(sig_clamp_gene_ids)
sig_clamp_df <- as.matrix(sig_clamp_df)

sig_1step_tgts <- topo[which(topo$Source %in% sig_clamp_genes), "Target"]
sig_2step_tgts <- unique(topo[which(topo$Source %in% sig_1step_tgts), "Target"])
sig_3step_tgts <- unique(topo[which(topo$Source %in% sig_2step_tgts), "Target"])

## Simulate with signal
# (OLD) TEST: remove signal now. Simulate for 10 time units (relaxation) and see if it ends up in M, as transit_time_df indicates it should
all.equal(sracipeParams(hd_traj_p2), sracipeParams(racipeNorm)[hd_traj_model,]) 

hd_traj_p2@metadata$config$simParams["nIC"] <- 1
hd_traj_p2 <- sracipeSimulate(hd_traj_p2, timeSeries = T, printStart = 0, 
                              printInterval = hd_traj_resolution, simulationTime = hd_traj_simTime,
                              genParams = F, genIC = F, integrateStepSize = 0.1, initialNoise=0.04, scaledNoise = T, 
                              stepper = "EM_Clamp", ouNoise_t = signal_tcorr, clampGenes=sig_clamp_gene_ids,
                              clampValues=sig_clamp_df)


hd_traj_df <- as.data.frame(t(hd_traj_p2@metadata$timeSeries)) # for stoch
#hd_traj_df <- as.data.frame(t(hd_traj_p2@metadata$timeSeriesDet)) # for det


## TODO
# relaxation simulation (check stability of final state, simulate for 40 @ 0.1)
hd_traj_relaxation <- hd_traj_p2
sracipeIC(hd_traj_relaxation) <-assays(hd_traj_relaxation)
hd_traj_relaxation <- sracipeSimulate(hd_traj_relaxation, timeSeries = T, printStart = 0, 
                              printInterval = hd_traj_resolution, simulationTime = hd_traj_simTime,
                              genParams = F, genIC = F, integrateStepSize = 0.1, initialNoise=0.0, scaledNoise = T, 
                              stepper = "EM_Clamp", ouNoise_t = signal_tcorr)






hd_traj_df <- hd_traj_df[,c(names(tmpMeans))]
for(gene in genes) {
  hd_traj_df[,gene] <- as.numeric(hd_traj_df[,gene])
}
hd_traj_df[,1:length(genes)] <- log2(1+hd_traj_df[,1:length(genes)]) # Log transform
hd_traj_df[,1:length(genes)] <- sweep(hd_traj_df[,1:length(genes)], 2, tmpMeans, FUN = "-") # scale
hd_traj_df[,1:length(genes)] <- sweep(hd_traj_df[,1:length(genes)], 2, tmpSds, FUN = "/") # scale

hd_traj_df$Time <- as.numeric(gsub("X", "", rownames(hd_traj_df)))
rownames(hd_traj_df) <- hd_traj_df$Time
hd_traj_df <- hd_traj_df[,c("Time",genes_reordered)]


# Add ICs
#modelICs <- as.data.frame(t(sracipeIC(hd_traj_checkpoint[,hd_traj_model])))[,genes_reordered] # use checkpoint ICs
modelICs <- as.data.frame(t(sracipeIC(hd_traj_p2)))[,genes_reordered] # use checkpoint SS
modelICs[,names(tmpMeans)] <- log2(1+modelICs[,names(tmpMeans)]) # Log transform
modelICs[,names(tmpMeans)] <- sweep(modelICs[,names(tmpMeans)], 2, tmpMeans, FUN = "-") # scale
modelICs[,names(tmpMeans)] <- sweep(modelICs[,names(tmpMeans)], 2, tmpSds, FUN = "/") # scale

modelICs$Time <- 0
modelICs <- modelICs[,c("Time",genes_reordered)]
hd_traj_df <- as.data.frame(rbind(modelICs, hd_traj_df))
rownames(hd_traj_df)[1] <- "0"
rownames(hd_traj_df) <- hd_traj_df$Time
#heatmap(as.matrix(hd_traj_df[rev(rownames(hd_traj_df)),genes_reordered]), Rowv = NA)#, Colv = NA)


# Plot zeb trajectory
#ggplot(hd_traj_df, aes(x=Time,y=Zeb1)) +
#  geom_line()

hd_traj_df_long <- pivot_longer(hd_traj_df, cols=all_of(genes), 
                                names_to = "Gene", values_to = "Expression")

ggplot(hd_traj_df_long[which(hd_traj_df_long$Gene %in% sig_1step_tgts),], 
       aes(x=Time,y=Expression, group=Gene, color=Gene)) +
  geom_line()


# Plot trajectory
hd_traj_pca <- as.data.frame(predict(pca, hd_traj_df[,names(tmpMeans)]))
hd_traj_pca$Time <- hd_traj_df$Time
hd_traj_pca$Model <- hd_traj_model

ggplot(pca_df) +
  geom_point(aes(x=PC1,y=PC2)) +
  geom_point(data = hd_traj_pca, aes(x=PC1,y=PC2,color=Time)) +
  geom_path(data = hd_traj_pca, aes(x=PC1,y=PC2,color=Time)) +
  scale_color_gradient()















ggplot(pca_df) +
  geom_point(aes(x=PC1,y=PC2, color=clust)) +
  #geom_point(data = hd_traj_pca, aes(x=PC1,y=PC2,color=Time)) +
  #geom_path(data = hd_traj_pca, aes(x=PC1,y=PC2,color=Time)) +
  scale_color_gradient()






plotCondition(racipe = racipe_bistable_raw,
              pca = pca,
              clust = clust,
              initialClust = 1,
              targetClust = 2,
              expName = expName,
              setName = hd_traj_sig,
              suffix="_09112024_example",
              outDir = file.path(topoDir))







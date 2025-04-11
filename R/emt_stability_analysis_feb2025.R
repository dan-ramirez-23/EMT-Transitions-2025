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
library(ComplexHeatmap)
library(extrafont)
loadfonts(device = "postscript")
source("R/utils.R")
source("R/utils_clamping.R")
source("R/scratch.R")

library(Ryacas0)  # For symbolic algebra
library(pracma) # numeric methods
library(readr)     # For saving progress

set.seed(1234)
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(3,2,4:8)]

# set up directories
topoName <- "emt_bhtopo_26node_CLAMP"
topoDir <- file.path(getwd(),topoName)
plotDir <- file.path(topoDir,"plots")
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
summary_df_fname <- file.path(dataDir,"state_summary_df_jan25test.Rds")
ss_unique_fname <- file.path(dataDir,"ss_unique_df_jan25test.Rds")
if(!file.exists(summary_df_fname)) {
  # Find unique steady states per model
  ss_rounded <- round(as.data.frame(unnormData), 1)
  ss_rounded$Model <- rep(1:numModels, each=numICs)
  ss_unique <- ss_rounded %>%
    #group_by(Model) %>%
    distinct()
  
  
  clust_fname <- file.path(dataDir,"clust_jan25test.Rds") 
  if(!file.exists(clust_fname)) {
    
    gmm = GMM(ss_unique[,1:nGenes], 2, dist_mode = "eucl_dist", seed_mode = "random_subset", km_iter = 10,
              em_iter = 10, verbose = F)          
    
    # predict centroids, covariance matrix and weights
    clust = predict(gmm, newdata = ss_unique[,1:nGenes])
    
    # check that clust 1 has higher Cdh
    summary(ss_unique[which(clust==1),"Cdh1"])
    summary(ss_unique[which(clust==2),"Cdh1"])
    
    # Im reversing it so the initial state is cluster 1
    revClust <- clust
    revClust <- ifelse(clust == 1, 2, ifelse(clust == 2, 1, clust))
    
    #mean(ss_unique$Cdh1[which(clust == 2)])
    #ggplot(pca_df, aes(x=PC1, y=PC2, color=as.factor(revClust))) + geom_point()
    
    saveRDS(revClust, file = clust_fname)
    clust <- revClust
  } else {
    clust <- readRDS(clust_fname)
  }
  ss_unique$Cluster <- clust
  
  
  summary_df <- ss_unique %>%
    dplyr::select(Model, Cluster) %>%
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


racipe_bistable_fname <- file.path(dataDir,"racipe_bistable_jan25test.Rds")
models_selected_fname <- file.path(dataDir,"racipe_bistable_indexMap_jan25test.Rds")
if(!file.exists(racipe_bistable_fname)) {
  # filter for models with 2 states, only bistable
  # save new racipe object w/ only bistable
  models_selected <- unlist(summary_df[which(summary_df$NumStates == 2 & summary_df$Stability == "bistable"),"Model"])
  
  keepIdx <- c()
  for(model in models_selected) {
    addIdx <- sample((numICs*(model-1)+1):(numICs*model), 1) # select random state from each selected model
    keepIdx <- c(keepIdx, addIdx)
  }
  
  racipe_bistable <- racipe[,keepIdx]
  sracipeParams(racipe_bistable) <- sracipeParams(racipe)[as.numeric(unname(models_selected)),]
  racipe_bistable@metadata$config$simParams[["numModels"]] <- length(keepIdx)
  
  saveRDS(racipe_bistable, racipe_bistable_fname)
  saveRDS(models_selected, models_selected_fname)
} else {
  racipe_bistable <- readRDS(racipe_bistable_fname)
  models_selected <- readRDS(models_selected_fname)
}

racipe_bistable_raw <- racipe_bistable
racipe_bistable_raw@metadata$config$simParams["nIC"] <- 1
unnormData <- t(assay(racipe_bistable_raw))
racipe_bistable <- sracipeNormalize(racipe_bistable)
exprMat <- as.data.frame(t(assay(racipe_bistable)))



########## PCA & CLUSTERING ############

# PCA-based "cluster" assignment
pca_fname <- file.path(dataDir,"pca_jan25test.Rds")
if(!file.exists(pca_fname)) {
  pca <- prcomp(exprMat[,1:nGenes])
  pca$rotation[,"PC1"] <- -1*pca$rotation[,"PC1"]
  pca$x[,1] <- -1*pca$x[,1]
  
  pca_df <- as.data.frame(pca$x)
  
  ggplot(pca_df, aes(x=PC1, y=PC2, color=exprMat$Cdh1)) + geom_point()
  
  saveRDS(pca, pca_fname)
} else {
  pca <- readRDS(pca_fname)
  pca_df <- as.data.frame(pca$x)
}

# Clustering bistable models
clust_fname <- file.path(dataDir,"clust_bistable_jan25test.Rds") 
if(!file.exists(clust_fname)) {
  
  gmm = GMM(exprMat[,1:nGenes], 2, dist_mode = "eucl_dist", seed_mode = "random_subset", km_iter = 10,
            em_iter = 10, verbose = F)          
  
  # predict centroids, covariance matrix and weights
  clust = predict(gmm, newdata = exprMat[,1:nGenes])
  
  # Im reversing it so the initial state is cluster 1
  #revClust <- clust
  #revClust <- ifelse(clust == 1, 2, ifelse(clust == 2, 1, clust))
  
  ggplot(pca_df, aes(x=PC1, y=PC2, color=as.factor(clust))) + geom_point()
  
  saveRDS(clust, file = clust_fname)
} else {
  clust <- readRDS(clust_fname)
}
exprMat$Cluster <- clust


# Plot WT as baseline
pca_df <- as.data.frame(pca$x)
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
clamp_df_fname <- file.path(dataDir,"clamp_values_jan25test.Rds")
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



ggplot(data = clamp_df, aes(x = Gene, y = Expression, fill = as.factor(Cluster))) + 
  geom_boxplot() +
  labs(title = "Gene Expression by Cluster", 
       x = "Gene", 
       y = "Expression",
       fill = "Cluster") +
  theme_minimal() +
  #scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust text angle for better visibility


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
# pca_plot_fname <- file.path(plotDir,"fig1b_pca_wt_dist.pdf")
# pdf(pca_plot_fname, height = 10, width = 10)
# print(image)
# dev.off()

# sanity check - Cdh1 should be high on the left side, cluster 1
#ggplot(pca_df, aes(x=PC1, y=PC2, color=racipeData$Cdh1)) +
#  geom_point(size=3) 


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

# wt_hmap_fname <- file.path(plotDir,"fig1c_wt_hmap.pdf")
# pdf(wt_hmap_fname, height = 10, width = 10)
# print(image)
# dev.off()


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
  geom_point(size=3, alpha=0.5) +
  geom_point(data=pca_plot_df, aes(x=PC1, y=PC2, fill=Cluster), color="black", pch=21, size=3) +
  guides(color=guide_legend(title = "Cluster")) +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(values=cbPalette) +
  xlab(plot_xlab) +
  ylab(plot_ylab) +
  theme_sticcc() +
  guides(color="none") +
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color="black"))
image
pca_plot_fname <- file.path(plotDir,"figxx_pca_wt_dist_ALL_hilites.pdf")
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
                 row_names_gp=gpar(fontsize=12))

wt_hmap_fname <- file.path(plotDir,"figxx_wt_hmap_ALL.pdf")
pdf(wt_hmap_fname, height = 10, width = 10)
print(image)
dev.off()



## Multistability histogram
# Define colors


summary_hist_df_in <- read.table(file.path("/Users/danramirez/Desktop/NEU/projects/emt_state_transitions",
                                           "pyracipe","emt_bhtopo_26node_5k200ic_rep2",
                                           "emt_bhtopo_26node_5k200ic_rep2_pyracipe_summary.csv"),
                                 header = T, row.names = 1, sep = ",")


# Ensure summary has the expected structure
count_df <- summary_hist_df_in %>%
  count(NO_STATES, StateIdentity, name = "Count") %>%  # Count occurrences safely
  arrange(NO_STATES)
count_long <- count_df[which(count_df$NO_STATES > 0),]
# Convert to long format (no need for pivot_wider, just use the counted data)
#count_long <- count_df %>%
#  mutate(NO_STATES = as.factor(NO_STATES))  # Ensure NO_STATES is treated as categorical

# Define colors
count_long$StateIdentity <- factor(count_long$StateIdentity, levels=c("E", "Bistable", "M"))
colors <- c("E" = "blue", "Bistable" = "purple", "M" = "red")

max_bin <- 7
count_long <- count_long %>%
  mutate(NO_STATES = ifelse(NO_STATES >= max_bin, paste0(max_bin, "+"), as.character(NO_STATES)))

# Convert NumStates to factor for proper ordering
count_long$NO_STATES <- factor(count_long$NO_STATES, levels = c(as.character(1:max_bin), paste0(max_bin, "+")))


# Create stacked bar plot
ggplot(count_long, aes(x = factor(NO_STATES), y = Count, fill = StateIdentity)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  labs(title = "RACIPE Models by Steady States & Cluster",
       x = "States",
       y = "Count",
       fill = "Cluster Identity") +
  theme_sticcc() +
  theme(title = element_text(size=10, inherit.blank = FALSE))

# summary_hist_df <- summary_df
# summary_hist_df[which(summary_hist_df$Stability == 1),"Stability"] <- "E"
# summary_hist_df[which(summary_hist_df$Stability == 2),"Stability"] <- "M"
# summary_hist_df[which(summary_hist_df$Stability == "bistable"),"Stability"] <- "Bistable"
# summary_hist_df$Stability <- factor(summary_hist_df$Stability, levels = c("E", "Bistable", "M"))
# colors <- c("E" = "blue", "Bistable" = "purple", "M" = "red")
# 
# # Recode NumStates to group values above max_bin
# max_bin <- 7
# summary_hist_df <- summary_hist_df %>%
#   mutate(NumStates = ifelse(NumStates >= max_bin, paste0(max_bin, "+"), as.character(NumStates)))
# 
# # Convert NumStates to factor for proper ordering
# summary_hist_df$NumStates <- factor(summary_hist_df$NumStates, levels = c(as.character(1:max_bin), paste0(max_bin, "+")))
# 
# 
# # Create the plot
# ggplot(summary_hist_df, aes(x = factor(NumStates), fill = Stability)) +
#   geom_bar() +
#   scale_fill_manual(values = colors) +
#   labs(title = "RACIPE Models by Steady States & Cluster", 
#        x = "States", 
#        y = "Counts",
#        fill = "Cluster Identity") +
#   theme_sticcc()


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
                  "_tau=",signal_tcorr,"_genes=",paste0(signal_nGenes,collapse = "."),"_CLAMPS_jan2025")


# undebug(optimizeST_parallel)
# paramSets <- optimizeST_parallel(racipe_bistable_raw,
#                                  pca,
#                                  clust,
#                                  totalPerturbations = 500,
#                                  initialClust=1,
#                                  targetClust=2,
#                                  nSigGenes=signal_nGenes,
#                                  outDir = file.path(topoDir),
#                                  expName = expName,
#                                  paramType="G",
#                                  nPerturbations = 351,
#                                  plot=F,
#                                  fcVals = 1,
#                                  noise = signal_noise,
#                                  checkpointSize = 5,
#                                  forceRerun = F,
#                                  forceRecompute=F,
#                                  anneal = F,
#                                  randomParams = F,
#                                  relax=T,
#                                  simTime = signal_simTime,
#                                  simTimeRelax = signal_relaxTime,
#                                  nCores = 6,
#                                  saveTrajectory = F,
#                                  printInterval = 10,
#                                  noise_tcorr = signal_tcorr,
#                                  clamp=T,
#                                  clamp_df=clamp_df)



######## AGGREGATE ANALYSIS #####
resultSet_fname <- file.path(topoDir,expName,"result_summary.Rds")
resultSet_full <- readRDS(resultSet_fname)
resultSet <- resultSet_full
rs_det <- resultSet_full[which(resultSet_full$Noise == 0),]
## this breaks sequence, must run after rs_list - only here as a reminder
#resultSet_full <- rs_list_full[[2]]
#resultSet_full$ConversionPct <- resultSet_full$ConvertingModels / resultSet_full$startingInitPopulation
#saveRDS(resultSet_full, file.path(topoDir,expName,"result_summary.Rds"))

ggplot(resultSet_full[,], aes(x=factor(Noise), y=ConversionPct, fill=factor(Noise))) +
  #geom_violin() +
  geom_boxplot() +
  xlab("Noise") +
  ylab("Models Undergoing EMT (%)") +
  scale_fill_discrete(name="Noise") +
  theme_sticcc() +
  theme(axis.line = element_line(color="black"), axis.title = element_text(size=20))



#selectedNoise <- c(0, 0.04) 
selectedNoise <- c(0.04)
#selectedNoise <- c(0.04)
resultSet <- resultSet_full[which(resultSet_full$Noise %in% selectedNoise),]



# 
# rs_list <- genes_x_transitions(resultSet, # dataframe w/ columns: ModelNo, SetName
#                                topoName = topoName, # string
#                                collectionName = expName, # string
#                                initialClust = 1, # int
#                                targetClust = 2, # int
#                                wt_data = exprMat[,genes], # matrix/dataframe of size numSamples x numFeatures
#                                clust = clust, # vector of length numSamples containing integers
#                                tmpMeans = tmpMeans,
#                                tmpSds = tmpSds
# )
# 
# 
# rsMatrix <- rs_list[[1]]
# # resultSet <- rs_list[[2]]
# # resultSet$ConversionPct <- resultSet$ConvertingModels / resultSet$startingInitPopulation
# rownames(rsMatrix) <- resultSet$SetName[which(!is.na(resultSet$ConversionPct))]
# # summary(resultSet$ConversionPct)



# Plot 2D heatmap at different noise levels
eff_hmap_0.04 <- resultSet_full[which(resultSet_full$Noise == 0.2),
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


#image(eff_hmap_0.04_matrix)
eff_hmap_0.04_matrix[lower.tri(eff_hmap_0.04_matrix)] <- NA

# Define the color mapping
library(circlize)
library(viridis)
# Define the color mapping using viridis
col_fun <- colorRamp2(c(min(eff_hmap_0.04_matrix, na.rm = TRUE), 
                        mean(eff_hmap_0.04_matrix, na.rm = TRUE), 
                        max(eff_hmap_0.04_matrix, na.rm = TRUE)), 
                      viridis(3))  # Generates 3 colors from viridis

# Generate the heatmap
Heatmap(eff_hmap_0.04_matrix, 
        name = "ConversionPct", 
        col = col_fun, 
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        show_row_names = TRUE, 
        show_column_names = TRUE)


## 1-gene efficacy barchart mirroring boolean model

barchart_df <- resultSet_full[which(resultSet_full$NumGenes == 1 & resultSet_full$Noise == 0.04),
                              c("Species 1", "ConversionPct")]
barchart_df$`Species 1` <- factor(barchart_df$`Species 1`, levels = genes_reordered)

ggplot(barchart_df, aes(x=`Species 1`, y=ConversionPct)) +
  geom_bar(stat="identity") +
  labs(y="Models Undergoing EMT (%)") +
  theme_sticcc() +
  theme(axis.text.x = element_text(size=18, angle = 90),
        axis.line = element_line(),
        axis.ticks = element_line())



######## TIMING SIMULATIONS IN PARALLEL #####

# Select signals to use
# top n by conversion pct,  and n more random
# or read old signals & reuse those
# resultSet <- resultSet_full[which(resultSet_full$NumGenes == 1 & resultSet_full$Noise == 0.04),] 
# n <- 1
# nPer <- 10
# 
# setIDList <- c()
# rankorder <- order(resultSet$ConversionPct, decreasing = T)
# topn <- rankorder[which(!rankorder %in% setIDList)][1:n]
# setIDList <- rep(c(setIDList, topn), 10)
# sigNames <- resultSet[setIDList,"SetName"]
# names(sigNames) <- setIDList
# 
# 
# 
# nCores <- 4
# requireNamespace("doFuture")
# registerDoFuture()
# plan(multisession, workers = nCores)
# 
# ## CHECK THESE EVERY TIME
# 
# resultSet$Tau <- rep(signal_tcorr, nrow(resultSet))
# resultSet$ScaledNoise <- rep(T, nrow(resultSet))
# resultSet$Noise <- 0.04
# resultSet$SetName <- paste0(resultSet$`Species 1`,"_noise=",resultSet$Noise)
# 
# trajSimTime <- 300
# expName_new_list <- paste0("bhtopo_t=",trajSimTime,
#                            "_relax_OUnoise=",0.04,
#                            "_tau=",signal_tcorr,
#                            "_genes=",paste0(signal_nGenes,collapse = "."),
#                            "_SIG=",resultSet[setIDList, "SetName"],
#                            "_runNo=",rep(seq(nPer),each=n))
# times <- c(seq(2, 30, 2), seq(35, 100, 5), seq(120, trajSimTime, 20))
# #times <- c(seq(10, 100, 10), seq(120, 500, 20), seq(550, 1000, 50))

######## TIME SERIES, ONE MODEL ############

trajSimTime <- 10
step <- 2
hd_traj_sig <- "Zeb1_noise=0.04"#resultSet[setID, "SetName"]
hd_traj_sigID <- which(resultSet$SetName == hd_traj_sig)#setID
setName <- hd_traj_sig

zeb_1step_tgts <- unique(topo[which(topo$Source == "Zeb1"), "Target"])
zeb_2step_tgts <- unique(topo[which(topo$Source %in% zeb_1step_tgts), "Target"])
zeb_3step_tgts <- unique(topo[which(topo$Source %in% zeb_2step_tgts), "Target"])


genes_reordered <- c("Foxc2","Zeb1","Klf8","Cdh1","miR101", "Zeb2", "Snai1", "miR141",
                     "Tgfbeta","miR200a","miR200b","miR200c","miR205","miR30c","Snai2",
                     "miR34a","Twist2","miR9","Vim","Twist1","Tcf3","Gsc", "Ovol2", "Grhl2",  "Np63a", "Cldn7")


traj_model_list <- c(423)#c(423, 461, 469, 475)#c(sample(noise_models_list, 3))#c(7,8)#c(14, 26, 41, 156)
# jan 2025: 16, 26 transition deterministically. 

# 73, 81, 92 are stochastic-10/10 and stay in 2.
# t=1000: 423, 461 same as above

# 28, 32, 82 are stochastic-10/10 but don't stay in cluster 2
# t=1000: 469, 475


noAttemptsList <- c(20)
# Sim parameters
hd_traj_startTime <- 2
hd_traj_simTime <- 1000
hd_traj_relaxTime <- 50
hd_traj_resolution <- 0.2


######## DETAILED TRAJECTORY ANALYSIS ############
hd_traj_model <- 423
attemptNo <- 1

# set up data storage
modelDataDir <- file.path(dataDir,paste0("trajectories_model",hd_traj_model))
modelPlotDir <- file.path(plotDir,paste0("trajectories_model",hd_traj_model))
if(!dir.exists(modelDataDir)) {
  dir.create(modelDataDir)
}

# Setup
hd_traj_model_idx <- models_selected[paste0("Model",hd_traj_model)] # for retrieving relevant state
model_states <- ss_unique[which(ss_unique$Model == hd_traj_model_idx), ]
#hd_traj_model <- sample(which(clust==1), 1)
# Sim parameters


traj_fname <- file.path(modelDataDir,paste0("trajectory_model",hd_traj_model,
                                            "_simTime=",hd_traj_simTime,
                                            "_SIG=",hd_traj_sig,"_2025-01-28_v",attemptNo,".Rds"))

traj_pca_fname <- file.path(modelDataDir,paste0("trajectoryPCA_model",hd_traj_model,
                                                "_simTime=",hd_traj_simTime,
                                                "_SIG=",hd_traj_sig,"_2025-01-28_v",attemptNo,".Rds"))

hd_traj_pca <- readRDS(traj_pca_fname)
hd_traj_df_long <- readRDS(traj_fname)
hd_traj_df_full <- pivot_wider(
  hd_traj_df_long,
  names_from = "Gene",      # Column with names to become new columns
  values_from = "Expression" # Column with values to fill the new columns
)


# reconstruct raw model trajectory
hd_traj_df_raw <- hd_traj_df_full[,genes]
hd_traj_df_raw[,1:length(genes)] <- sweep(hd_traj_df_raw[,1:length(genes)], 2, tmpSds, FUN = "*") # scale
hd_traj_df_raw[,1:length(genes)] <- sweep(hd_traj_df_raw[,1:length(genes)], 2, tmpMeans, FUN = "+") # scale
hd_traj_df_raw[,1:length(genes)] <- 2^hd_traj_df_raw[,1:length(genes)] - 1
hd_traj_df_raw$Time <- hd_traj_df_full$Time  



######## CLAMPED STABILITY ANALYSIS ############
## Here, we'll take one model of interest, and simulate 2000 random ICs with & without clamps
## Goal is to see if a third intermediate state is present in WT (explaining jumps)

modelDataDir <- file.path(dataDir,paste0("trajectories_model",hd_traj_model))
if(!dir.exists(modelDataDir)) {
  dir.create(modelDataDir)
}

# Setup
hd_traj_model_idx <- models_selected[paste0("Model",hd_traj_model)] # for retrieving relevant state
model_states <- ss_unique[which(ss_unique$Model == hd_traj_model_idx), ]
num_ics <- 500
num_models <- 1

hd_traj_checkpoint <- racipe_bistable_raw # actually just use WT object
dim(hd_traj_checkpoint)

# Filter selected model & set ICs
hd_traj_p2 <- hd_traj_checkpoint[,hd_traj_model]
hd_traj_p2@metadata$config$simParams["numModels"] <- num_models
hd_traj_p2@metadata$config$simParams["nIC"] <- num_ics

# confirm parameters match
orig_params <- as.numeric(sracipeParams(racipe_bistable_raw)[as.numeric(hd_traj_model),])
new_params <- as.numeric(sracipeParams(hd_traj_p2))
all.equal(new_params, orig_params)


# Generate ICs
hd_traj_p2 <- sracipeSimulate(hd_traj_p2, timeSeries = F, simulationTime = 0.1, integrate = F,
                              genParams = F, genIC = T, integrateStepSize = 0.01, initialNoise=0.0, scaledNoise = T, 
                              stepper = "EM", ouNoise_t = signal_tcorr)


# Implement signal
sig_clamp_genes <- getClampGenes(resultSet, hd_traj_sigID)
sig_clamp_gene_ids <- unlist(lapply(sig_clamp_genes, function(x) which(rownames(racipe) == x)))
sig_clamp_df <- getClampDF(clamp_df, sig_clamp_genes, 2)
colnames(sig_clamp_df) <- as.numeric(sig_clamp_gene_ids)
sig_clamp_df <- as.matrix(sig_clamp_df)


## Simulate with signal
hd_traj_p2@metadata$config$simParams["numModels"] <- num_models
hd_traj_p2@metadata$config$simParams["nIC"] <- num_ics
hd_traj_p2 <- sracipeSimulate(hd_traj_p2, simulationTime = 200,
                              genParams = F, genIC = F, integrateStepSize = 0.2, initialNoise=0.0, scaledNoise = T, 
                              integrate = T,
                              stepper = "EM_Clamp", 
                              ouNoise_t = signal_tcorr, 
                              clampGenes=sig_clamp_gene_ids,
                              clampValues=sig_clamp_df[hd_traj_model,]
)
dim(assay(hd_traj_p2))
as.data.frame(sracipeParams(hd_traj_p2))[40,1:5]


# normalize results
model_ss_df <- as.data.frame(t(assay(hd_traj_p2)))
clamped_states_unique_raw <- distinct(round(model_ss_df, 3))
model_ss_df_norm <- model_ss_df
for(gene in genes) {
  model_ss_df_norm[,gene] <- as.numeric(model_ss_df_norm[,gene])
}
model_ss_df_norm[,1:length(genes)] <- log2(1+model_ss_df_norm[,1:length(genes)]) # Log transform
model_ss_df_norm[,1:length(genes)] <- sweep(model_ss_df_norm[,1:length(genes)], 2, tmpMeans, FUN = "-") # scale
model_ss_df_norm[,1:length(genes)] <- sweep(model_ss_df_norm[,1:length(genes)], 2, tmpSds, FUN = "/") # scale
model_ss_pca <- as.data.frame(predict(pca, model_ss_df_norm[,names(tmpMeans)]))
model_ss_pca_unique <- unique(round(model_ss_pca, 2))
model_ss_clamped_unique <- unique(round(model_ss_df_norm, 2))

# Plot states, with known model states in red
model_states_norm <- model_states
model_states_norm[,1:length(genes)] <- log2(1+model_states_norm[,1:length(genes)]) # Log transform
model_states_norm[,1:length(genes)] <- sweep(model_states_norm[,1:length(genes)], 2, tmpMeans, FUN = "-") # scale
model_states_norm[,1:length(genes)] <- sweep(model_states_norm[,1:length(genes)], 2, tmpSds, FUN = "/") # scale
model_states_pca <- as.data.frame(predict(pca, model_states_norm[,names(tmpMeans)]))

# Plot available states under clamping
ggplot(model_ss_pca) +
  geom_point(aes(x=PC1, y=PC2), size=4) +
  geom_point(data=pca_df, aes(x=PC1, y=PC2), color="lightblue") +
  geom_point(data=model_states_pca, aes(x=PC1, y=PC2), color="red", size=3)


########## ODE EXTRACTION FUNCTIONS ############

# Function to define regulation using Hill functions in symbolic form
hill_function <- function(gene_value, threshold, n, lambda, mode) {
  gene_value <- Sym(gene_value)  # Convert gene symbolically
  threshold <- Sym(threshold)    # Convert threshold to symbol
  lambda <- Sym(lambda)          # Convert lambda to symbol
  n <- Sym(n)                    # Convert n to symbol
  
  if (mode == "activation") {
    return( (lambda + (1 - lambda) / (1 + (gene_value / threshold)^n)) / lambda )
  } else if (mode == "repression") {
    return((lambda + (1 - lambda) / (1 + (gene_value / threshold)^n)))
  } else {
    return(Sym("1"))  # No interaction
  }
}

# Function to generate ODEs from the topology matrix
generate_odes <- function(topology, gGene, kGene, threshold_matrix, n_matrix, lambda_matrix) {
  num_genes <- nrow(topology)
  equations <- list()
  gene_names <- rownames(topology)
  
  for (i in 1:num_genes) {
    gene_name <- gene_names[i]  # Current gene being processed
    interaction_terms <- list()
    
    for (j in 1:num_genes) {
      interaction_type <- topology[j, i]  # Incoming edge: j → i
      
      if (interaction_type == 1) {  # Activation
        term <- hill_function(paste0("X", gene_names[j]), threshold_matrix[j, i], n_matrix[j, i], lambda_matrix[j, i], "activation")
        interaction_terms <- append(interaction_terms, term)
      } else if (interaction_type == 2) {  # Repression
        term <- hill_function(paste0("X", gene_names[j]), threshold_matrix[j, i], n_matrix[j, i], lambda_matrix[j, i], "repression")
        interaction_terms <- append(interaction_terms, term)
      }
    }
    
    # Multiply all interaction terms
    if (length(interaction_terms) > 0) {
      final_multiplier <- Reduce(function(a, b) Sym(paste0("(", a, ") * (", b, ")")), interaction_terms)
    } else {
      final_multiplier <- Sym("1")
    }
    
    # Construct full ODE
    expr <- Sym(paste0(gGene[i], " * (", final_multiplier, ") - ", kGene[i], " * X", gene_name))
    equations[[i]] <- expr
  }
  
  return(equations)
}

# Function to convert RACIPE output into matrices for ODE generation
convert_racipe_to_ode_input <- function(topology_df, param_vector) {
  # Get unique gene names
  genes <- unique(c(topology_df$Source, topology_df$Target))
  num_genes <- length(genes)
  
  # Create empty matrices for parameters
  topology <- matrix(0, nrow = num_genes, ncol = num_genes, dimnames = list(genes, genes))
  threshold_matrix <- matrix(NA, nrow = num_genes, ncol = num_genes, dimnames = list(genes, genes))
  n_matrix <- matrix(NA, nrow = num_genes, ncol = num_genes, dimnames = list(genes, genes))
  lambda_matrix <- matrix(NA, nrow = num_genes, ncol = num_genes, dimnames = list(genes, genes))
  
  # Initialize vectors for gene-level parameters
  gGene <- setNames(rep(NA, num_genes), genes)
  kGene <- setNames(rep(NA, num_genes), genes)
  
  # Fill the topology matrix based on the input table
  for (i in 1:nrow(topology_df)) {
    source <- topology_df$Source[i]
    target <- topology_df$Target[i]
    interaction_type <- topology_df$Type[i]
    
    # Store interaction type (1 = activation, 2 = repression)
    topology[source, target] <- interaction_type
  }
  
  # Extract gene-level parameters
  for (gene in genes) {
    gGene[gene] <- param_vector[paste0("G_", gene)]  # Max production
    kGene[gene] <- param_vector[paste0("K_", gene)]  # Degradation
  }
  
  # Extract edge-level parameters
  for (i in 1:nrow(topology_df)) {
    source <- topology_df$Source[i]
    target <- topology_df$Target[i]
    
    threshold_matrix[source, target] <- param_vector[paste0("TH_", source, "_", target)]
    n_matrix[source, target] <- param_vector[paste0("N_", source, "_", target)]
    
    # Fold-change is used to compute lambda (lambda = 1 / FC for repression)
    fold_change <- param_vector[paste0("FC_", source, "_", target)]
    
    if (topology[source, target] == 1) {  # Activation
      lambda_matrix[source, target] <- fold_change
    } else if (topology[source, target] == 2) {  # Repression
      lambda_matrix[source, target] <- 1 / fold_change
    }
  }
  
  return(list(
    topology = topology,
    gGene = gGene,
    kGene = kGene,
    threshold_matrix = threshold_matrix,
    n_matrix = n_matrix,
    lambda_matrix = lambda_matrix
  ))
}

########## ODE EXTRACTION SETUP ############

# # Example topology matrix (3 genes)
# topology <- matrix(c(0, 1, 2,
#                      2, 0, 1,
#                      1, 2, 0), nrow=3, byrow=TRUE)
# 
# # Example parameters (these would come from the RACIPE model)
# gGene <- c(1.5, 1.2, 1.8)  # Production rates
# kGene <- c(0.5, 0.6, 0.4)  # Degradation rates
# threshold_matrix <- matrix(c(1, 1, 1, 
#                              1, 1, 1,
#                              1, 1, 1), nrow=3, byrow=TRUE)
# n_matrix <- matrix(c(2, 2, 2, 
#                      2, 2, 2, 
#                      2, 2, 2), nrow=3, byrow=TRUE)
# lambda_matrix <- matrix(c(0.5, 0.5, 0.5,
#                           0.5, 0.5, 0.5,
#                           0.5, 0.5, 0.5), nrow=3, byrow=TRUE)
# 
# # Generate ODEs
# odes <- generate_odes(topology, gGene, kGene, threshold_matrix, n_matrix, lambda_matrix)
# 
# # Print ODEs
# print(odes)


ode_model_params <- unlist(sracipeParams(racipe_bistable_raw)[hd_traj_model,])
ode_model_inputs <- convert_racipe_to_ode_input(topo, ode_model_params)
odes <- generate_odes(ode_model_inputs$topology, 
                      ode_model_inputs$gGene, 
                      ode_model_inputs$kGene, 
                      ode_model_inputs$threshold_matrix, 
                      ode_model_inputs$n_matrix, 
                      ode_model_inputs$lambda_matrix)
#print(odes)



# # spot check, see if interactions imported correctly
# 
# #miR9 is an input, gene 4
# odes[[4]]
# 
# # Vim is an output, gene 26
# topo[which(topo$Target == "Vim"),]
# ode_model_params[grep("Vim", names(ode_model_params))]
# odes[[26]]
# 
# # Twist2
# topo[which(topo$Target == "Twist2"),]
# which(genes == "Twist2")
# ode_model_params[grep("Twist2", names(ode_model_params))]
# odes[[18]]


########## NEWTON'S METHOD SETUP ############

# Define function F(X) for steady-state equations with clamping
make_numeric_ode_system <- function(equations, vars, X_clamp = list()) {
  function(x) {
    subs_list <- setNames(as.list(x), vars)
    
    # Evaluate ODE system normally
    f_val <- sapply(equations, function(eq) Eval(eq, subs_list))
    
    # Override with clamped values
    for (clamp_var in names(X_clamp)) {
      idx <- which(vars == clamp_var)
      if (length(idx) > 0) {
        f_val[idx] <- X_clamp[[clamp_var]] - x[idx]  # Enforce clamp condition
      }
    }
    
    return(f_val)
  }
}

# Compute Jacobian numerically via finite differences with clamping
finite_difference_jacobian <- function(f_numeric, x, vars, X_clamp = list(), epsilon = 1e-6) {
  n <- length(x)
  J <- matrix(0, n, n)

  for (j in 1:n) {
    if (vars[j] %in% names(X_clamp)) {
      J[, j] <- 0  # Clamped variables have zero influence on the system
      J[j, j] <- -1  # Keep diagonal term to enforce clamp condition
    } else {
      x_plus <- x
      x_plus[j] <- x_plus[j] + epsilon  # Perturb x_j

      f_x <- f_numeric(x)               # F(X)
      f_x_plus <- f_numeric(x_plus)      # F(X + εe_j)

      J[, j] <- (f_x_plus - f_x) / epsilon  # Compute finite difference
    }
  }

  return(J)
}


# finite_difference_jacobian_num <- function(f_numeric, x, vars, X_clamp = list()) {
#   n <- length(x)
#   J <- matrix(0, n, n)  # Initialize Jacobian
#   
#   # Identify clamped variable indices
#   clamped_indices <- which(vars %in% names(X_clamp))
#   unclamped_indices <- setdiff(1:n, clamped_indices)
#   
#   # Define wrapper function that evaluates only unclamped variables
#   f_unclamped <- function(x_unclamped) {
#     x_full <- x
#     x_full[unclamped_indices] <- x_unclamped  # Update only unclamped values
#     f_numeric(x_full)
#   }
#   
#   # Compute Jacobian only for unclamped variables
#   if (length(unclamped_indices) > 0) {
#     J_unclamped <- jacobian(f_unclamped, x[unclamped_indices])
#     
#     # Assign computed values to the correct locations in J
#     J[unclamped_indices, unclamped_indices] <- J_unclamped
#   }
#   
#   # Ensure clamped variable columns are zero
#   J[, clamped_indices] <- 0
#   
#   return(J)
# }

# Newton's Method with finite-difference Jacobian, clamping, damping, and bisection
newtons_method_fd <- function(f_numeric, x0, vars, X_clamp = list(), tol = 1e-6, max_iter = 100) {
  x <- x0
  alpha <- 1  # Initial step size for damping
  prev_f_val <- NULL  # Track previous function values for oscillation detection
  prev2_f_val <- NULL
  damping_count <- 0 # Number of consecutive damped steps
  
  for (i in 1:max_iter) {
    f_val <- f_numeric(x)  # Evaluate function
    J_val <- finite_difference_jacobian(f_numeric, x, vars, X_clamp)  # Compute Jacobian
    
    print(paste("Iteration", i))
    print(paste("Max f_x val:", max(abs(f_val))))
    
    # Check for convergence
    if (i == 30) { # Relax constraints if it's converging slowly
      tol <- tol * 10
    }
    if (i == 50) {
      tol <- tol * 10
    }
    if (max(abs(f_val)) < tol) {
      return(list(solution = x, converged = TRUE, iterations = i))
    }
    
    # Solve for delta_x
    delta_x <- tryCatch(
      qr.solve(J_val, -f_val),
      error = function(e) {
        warning("Singular Jacobian encountered. Using pseudo-inverse.")
        return(MASS::ginv(J_val) %*% -f_val)  # Generalized inverse as fallback
      }
    )
    
    if (any(is.na(delta_x))) {
      stop("Newton update resulted in NaN values. Stopping iteration.")
    }
    
    # Detect oscillation: Check if f_val is alternating between two values
    if (!is.null(prev_f_val) && !is.null(prev2_f_val)) {
      if (all(abs(f_val - prev2_f_val) < tol) && all(abs(prev_f_val - prev3_f_val) < tol)) {
        print("Oscillation detected. Reducing step size.")
        alpha <- alpha / 2  # Reduce step size if oscillation is detected
      }
    }
    
    # Store function values for oscillation detection
    prev2_f_val <- prev_f_val
    prev_f_val <- f_val
    
    
    # Damped update step
    cond_num <- kappa(J_val)
    if (cond_num > 1e6 && damping_count < 5) {  # Allow only 5 consecutive damped steps
      alpha <- 0.5
      damping_count <- damping_count + 1
      print(paste0("Ill-conditioned Jacobian. Damping for the ",damping_count," consecutive step"))
    } else {
      alpha <- 1  # Restore full Newton step occasionally
      damping_count <- 0
    }
    x_new <- x + alpha * delta_x
    

    
    # Bisection step if f_x is not decreasing significantly
    # if (max(abs(f_numeric(x_new))) > max(abs(f_val))) {
    #   print("Overshooting detected. Using bisection step.")
    #   x_new <- (x + x_new) / 2  # Bisection step to reduce overshoot
    # }
    
    x <- x_new  # Update estimate
  }
  
  return(list(solution = x, converged = FALSE, iterations = max_iter))
}

########## NEWTON'S METHOD WT ############
# Compute Jacobian symbolically
vars = paste0("X",genes)

# # Initial guess for steady state
x0 <- as.vector(unlist(model_states[1,genes])) * rnorm(length(genes), mean=1, sd=0.1)  # Start from equal expression levels

# Create numeric function F(X)
f_numeric <- make_numeric_ode_system(odes, vars)


# Run Newton's method with finite-difference Jacobian
result <- newtons_method_fd(f_numeric, x0, vars)

X_clamp <- list("XZeb1" = 7.0)
x0_clamp <- as.numeric(clamped_states_unique_raw[1,]) * rnorm(length(genes), mean=1, sd=0.01) 
result_clamp <- newtons_method_fd(f_numeric, x0_clamp, vars, X_clamp)
x0_clamp_h <- as.numeric(clamped_states_unique_raw[3,]) * rnorm(length(genes), mean=1, sd=0.01) 
result_clamp_h <- newtons_method_fd(f_numeric, x0_clamp_h, vars, X_clamp)



# Benchmark shows numeric differentiation doesn't give any speedup
# J_val_fd <- finite_difference_jacobian(f_numeric, x0_clamp, vars, X_clamp)
# 
# #J_val_fd_num <- finite_difference_jacobian_num(f_numeric, x0_clamp, vars, X_clamp)
# library(microbenchmark)
# microbenchmark(
#   manual = finite_difference_jacobian(f_numeric, x0_clamp, vars, X_clamp),
#   numDeriv =  finite_difference_jacobian_num(f_numeric, x0_clamp, vars, X_clamp),
#   times = 5
# )
# Unit: seconds
# expr      min       lq     mean   median       uq      max neval
# manual 36.85895 36.90693 37.30315 36.92896 37.06220 38.75871     5
# numDeriv 37.51600 37.57494 38.26062 37.67757 39.21903 39.31557     5



# Print results
if (result$converged) {
  print("Steady state found:")
  print(result$solution)
} else {
  print("Newton's method did not converge")
}



########## NEWTON BIFURCATION SETUP ############
library(doFuture)
library(doRNG)
bifurcation_analysis <- function(f_numeric, vars, clamp_var, clamp_values, 
                                 initial_steady_states, num_initial_guesses = 8, tol = 1e-6, 
                                 max_iter = 100, perturb_scale = 0.1, eigen_threshold = -0.5,
                                 save_results = TRUE, save_file = "bifurcation_results.csv") {
  registerDoFuture()
  options(future.availableCores = 4)
  plan(multisession, workers = 4)  # Use multiple cores
  
  # Load existing results if the file exists
  if (file.exists(save_file)) {
    existing_results <- read_csv(save_file)
    processed_clamp_values <- unique(existing_results$clamp_value)
  } else {
    existing_results <- NULL
    processed_clamp_values <- c()
  }
  
  # Function to generate diverse initial guesses
  generate_initial_guesses <- function(previous_states, clamp_val, vars) {
    
    previous_states <- Filter(Negate(is.null), previous_states)  # Remove NULL values
    
    if (length(previous_states) == 0) {
      return(replicate(num_initial_guesses, rep(1, length(vars)), simplify = FALSE))  # Default guesses
    }
    
    num_states <- length(previous_states)
    if (num_states >= num_initial_guesses) {
      return(previous_states[1:num_initial_guesses])
    }
    
    interpolated_states <- previous_states
    while (length(interpolated_states) < num_initial_guesses) {
      if (num_states > 1) {
        weights <- runif(2)
        weights <- weights / sum(weights)  # Random convex combination
        i <- sample(1:num_states, 2, replace = TRUE)
        new_guess <- weights[1] * previous_states[[i[1]]] + weights[2] * previous_states[[i[2]]]
      } else {
        new_guess <- previous_states[[1]]
      }
      names(new_guess) <- vars
      
      # Compute Jacobian at existing steady states and perturb along eigenvectors
      jacobian <- finite_difference_jacobian(f_numeric, new_guess, vars, list(clamp_var = clamp_val))
      eigvals <- eigen(jacobian)
      if (!any(is.na(eigvals$vectors))) {
        unstable_indices <- which(Re(eigvals$values) > eigen_threshold)  # Select eigenvectors above threshold
        if (length(unstable_indices) > 0) {
          chosen_vector <- eigvals$vectors[, sample(unstable_indices, 1)]
          perturbation <- chosen_vector * perturb_scale # Perturb along an unstable eigenvector
          new_guess <- pmax(new_guess + perturbation, 0)  # Ensure non-negative values
        }
      }
      
      new_guess[clamp_var] <- clamp_val  # Fix clamp variable
      interpolated_states <- append(interpolated_states, list(new_guess))
    }
    return(interpolated_states)
  }
  
  # Function to check stability using Jacobian eigenvalues
  is_stable <- function(jacobian) {
    eigvals <- eigen(jacobian)$values
    return(all(Re(eigvals) < 0))  # Stable if all real parts are negative
  }
  
  previous_states <- initial_steady_states
  results <- list()
  
  #adaptive_step <- diff(clamp_values)[1]  # Initial step size
  for (i in seq_along(clamp_values)) {
    
    clamp_val <- clamp_values[i]
    if (clamp_val %in% processed_clamp_values) {
      next  # Skip values already processed
    }
    
    #X_clamp <- list(clamp_var = clamp_val)
    X_clamp <- setNames(list(clamp_val), clamp_var)
    f_numeric_clamped <- make_numeric_ode_system(odes, vars, X_clamp)
    
    initial_guesses <- generate_initial_guesses(previous_states, clamp_val, vars)
    for(guess in seq_along(initial_guesses)) {
      initial_guesses[[guess]][clamp_var] <- clamp_val
    }
    
    steady_states <- foreach(x0 = initial_guesses, 
                             .combine = 'list', 
                             .packages = c("pracma", "MASS"),
                             .options.future = list(scheduling = 5)) %dorng% {
      result <- newtons_method_fd(f_numeric_clamped, x0, vars, X_clamp, tol, max_iter)
      if (result$converged) {
        jacobian <- finite_difference_jacobian(f_numeric_clamped, result$solution, vars, X_clamp)
        stable <- is_stable(jacobian)
        return(list(solution = result$solution, stability = stable, initial_guess=x0))
      }
      return(NULL)
    }
    
    steady_states <- Filter(Negate(is.null), steady_states)  # Remove failed solutions
    
    if (length(steady_states) > 0) {
      previous_states <- lapply(steady_states, function(ss) ss$solution)  # Use solutions for next step
    } else {
      print(paste0("No steady states for clamp value ",clamp_val," - proceeding to next clamp value"))
      #adaptive_step <- adaptive_step / 2  # Reduce step size if solutions disappear
      next
    }
    
    new_results <- data.frame(
      clamp_value = clamp_val, 
      steady_states = I(lapply(steady_states, `[[`, "solution")),
      stability = I(lapply(steady_states, `[[`, "stability")),
      initial_guess = I(lapply(steady_states, `[[`, "initial_guess"))
    )
    
    results <- append(results, list(new_results))
    
    if (save_results) {
      write_csv(bind_rows(existing_results, do.call(rbind, results)), save_file)
    }
  }
  
  return(do.call(rbind, results))
}





########## NEWTON'S METHOD BIFURCATION ############


# Define range of clamp values
clamp_var <- "XZeb1"  # Example

cluster1_expr <- unlist(clamp_df[which(clamp_df$Model == hd_traj_model_idx &
                                         clamp_df$Gene == sig_clamp_genes &
                                         clamp_df$Cluster == 1), "Expression"])
cluster2_expr <- unlist(clamp_df[which(clamp_df$Model == hd_traj_model_idx &
                                         clamp_df$Gene == sig_clamp_genes &
                                         clamp_df$Cluster == 2), "Expression"])
#clamp_values <- seq(0, 1.5e-05, length.out = 20)
clamp_values <- seq(1.5e-05, 0, length.out = 20)


# Known steady states from unclamped model
known_steady_states <- list(
  #setNames(as.numeric(model_states[1,genes]), vars),
  #setNames(as.numeric(model_states[2,genes]), vars),
  setNames(as.numeric(clamped_states_unique_raw[1,]), vars),
  setNames(as.numeric(clamped_states_unique_raw[2,]), vars),
  setNames(as.numeric(clamped_states_unique_raw[3,]), vars)
)

# Run bifurcation analysis
debug(bifurcation_analysis)
bifurcation_results_2vals <- bifurcation_analysis(f_numeric, vars, clamp_var, 
                                                  clamp_values, known_steady_states,
                                                  num_initial_guesses = 8,
                                                  save_file = file.path(modelDataDir,
                                                                        "bifurcation_results_newton_feb152025.csv"))



bifurcation_results <- read.table(file.path(modelDataDir,
                                            "bifurcation_results_newton_feb152025.csv"), sep = ",")
colnames(bifurcation_results) <- bifurcation_results[1,]
bifurcation_results <- bifurcation_results[-1,]

parse_vector_column <- function(df, col_name, suffix = "") {
  parsed_list <- lapply(df[[col_name]], function(x) eval(parse(text = x)))  # Convert string to numeric
  parsed_matrix <- do.call(rbind, parsed_list)  # Convert list to matrix
  colnames(parsed_matrix) <- paste0(colnames(parsed_matrix), suffix)  # Add suffix if needed
  return(as.data.frame(parsed_matrix))
}

steady_state_df <- parse_vector_column(bifurcation_results, "steady_states")
initial_guess_df <- parse_vector_column(bifurcation_results, "initial_guess", suffix = "_init")

# Combine all results into a final dataframe
bifurcation_results_final <- cbind(steady_state_df, initial_guess_df)
bifurcation_results_final[, c("clamp_value", "stability")] <- bifurcation_results[which(bifurcation_results$stability != "NULL"), c("clamp_value", "stability")]

bifurcation_results_final_norm <- bifurcation_results_final
colnames(bifurcation_results_final_norm)[1:length(genes)] <- genes
bifurcation_results_final_norm[,1:length(genes)] <- log2(1+bifurcation_results_final_norm[,1:length(genes)]) # Log transform
bifurcation_results_final_norm[,1:length(genes)] <- sweep(bifurcation_results_final_norm[,1:length(genes)], 2, tmpMeans, FUN = "-") # scale
bifurcation_results_final_norm[,1:length(genes)] <- sweep(bifurcation_results_final_norm[,1:length(genes)], 2, tmpSds, FUN = "/") # scale
bifurcation_results_final_pca <- as.data.frame(predict(pca, bifurcation_results_final_norm[,names(tmpMeans)]))
bifurcation_results_final_pca$clamp_value <- bifurcation_results_final$clamp_value
bifurcation_results_final_pca$stability <- bifurcation_results_final$stability


bifurcation_plot <- bifurcation_results_final_pca
bifurcation_plot$color <- ifelse(bifurcation_plot$stability, "black", "red")

ggplot(bifurcation_plot, aes(x = as.numeric(clamp_value), y = PC1, color = color)) +
  geom_point() +
  scale_color_identity() +
  xlim(0, 2e-05) +
  ylim(-4, 4) +
  labs(title = "Bifurcation Diagram", x = clamp_var, y = "Steady-State Values")

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
source("R/utils.R")
source("R/utils_clamping.R")
source("R/scratch.R")



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



# small test of clamping code
#clamp_test_df <- as.matrix(data.frame('6'=rnorm(numModels*0.1, mean = 10), '19'=rnorm(numModels*0.1, mean = 100)))
# clamp_test_df <- as.matrix(data.frame('6'=rnorm(numModels*0.1, mean = 10), '19'=rnorm(numModels*0.1, mean = 100)))
# 
# racipetest <- racipe[,seq(1, numModels*numICs*0.1, numICs)]
# racipetest@metadata$config$simParams[["numModels"]] <- numModels*0.1
# racipetest <- sracipeSimulate(racipetest, simulationTime = 10, nIC = 1, numModels = numModels*0.1,
#                               integrateStepSize = 0.1, initialNoise=0.04, scaledNoise = T, 
#                               stepper = "EM_Clamp", ouNoise_t = 10, clampGenes=c(6, 19), 
#                               clampValues=clamp_test_df)#c(0, 100))
# 

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
  # filter for models with <10 states, only bistable
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


########## WT VISUALIZATIONS ############

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

pca_plot_fname <- file.path(plotDir,"fig1b_pca_wt_dist.pdf")
pdf(pca_plot_fname, height = 10, width = 10)
print(image)
dev.off()

# sanity check - Cdh1 should be high on the left side, cluster 1
#ggplot(pca_df, aes(x=PC1, y=PC2, color=racipeData$Cdh1)) +
#  geom_point(size=3) 



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




rs_list <- genes_x_transitions(resultSet, # dataframe w/ columns: ModelNo, SetName
                               topoName = topoName, # string
                               collectionName = expName, # string
                               initialClust = 1, # int
                               targetClust = 2, # int
                               wt_data = exprMat[,genes], # matrix/dataframe of size numSamples x numFeatures
                               clust = clust, # vector of length numSamples containing integers
                               tmpMeans = tmpMeans,
                               tmpSds = tmpSds
)


rsMatrix <- rs_list[[1]]
# resultSet <- rs_list[[2]]
# resultSet$ConversionPct <- resultSet$ConvertingModels / resultSet$startingInitPopulation
rownames(rsMatrix) <- resultSet$SetName[which(!is.na(resultSet$ConversionPct))]
# summary(resultSet$ConversionPct)


######## TIMING SIMULATIONS IN PARALLEL #####

# Select signals to use
# top n by conversion pct,  and n more random
# or read old signals & reuse those
resultSet <- resultSet_full[which(resultSet_full$NumGenes == 1 & resultSet_full$Noise == 0.04),] 
n <- 1
nPer <- 10

setIDList <- c()
rankorder <- order(resultSet$ConversionPct, decreasing = T)
topn <- rankorder[which(!rankorder %in% setIDList)][1:n]
setIDList <- rep(c(setIDList, topn), 10)
sigNames <- resultSet[setIDList,"SetName"]
names(sigNames) <- setIDList



nCores <- 4
requireNamespace("doFuture")
registerDoFuture()
plan(multisession, workers = nCores)

## CHECK THESE EVERY TIME

resultSet$Tau <- rep(signal_tcorr, nrow(resultSet))
resultSet$ScaledNoise <- rep(T, nrow(resultSet))
resultSet$Noise <- 0.04
resultSet$SetName <- paste0(resultSet$`Species 1`,"_noise=",resultSet$Noise)

trajSimTime <- 300
expName_new_list <- paste0("bhtopo_t=",trajSimTime,
                           "_relax_OUnoise=",0.08,
                           "_tau=",signal_tcorr,
                           "_genes=",paste0(signal_nGenes,collapse = "."),
                           "_SIG=",resultSet[setIDList, "SetName"],
                           "_runNo=",rep(seq(nPer),each=n))
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
                                        noise = 0.08,
                                        forceRerun = F,
                                        forceRecompute = F,
                                        anneal=F,
                                        relax=T,
                                        simTimes = times,
                                        simTimeRelax = 40,
                                        save=T,
                                        clamp=T) 
}


multiSet_fname <- file.path(dataDir,
                            paste0("transitions_vs_time_Zeb1_10replicates_noise=0.04_",expName,".Rds"))
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
  #saveRDS(sigNames, file.path(dataDir,paste0("transitions_vs_time_signalNames_DetVStoch.Rds")))
  #saveRDS(setIDList, file.path(dataDir,paste0("transitions_vs_time_signalIDs_DetVStoch.Rds")))
} else {
  multiSet_times <- readRDS(multiSet_fname)
  #sigNames <- readRDS(file.path(dataDir,paste0("transitions_vs_time_signalNames_DetVStoch.Rds")))
  #setIDList <- readRDS(file.path(dataDir,paste0("transitions_vs_time_signalIDs_DetVStoch.Rds")))
}




## Plot conversions vs time
plot_df_list <- list()
idx <- 1
attemptIdxList <- rep(1:nPer, each=n)
for(setIDNo in seq_along(setIDList)) {
  setID <- setIDList[setIDNo]
  setName <- resultSet[setID, "SetName"]
  sampleSet_times <- multiSet_times[[setIDNo]]
  attemptNo <- attemptIdxList[setIDNo]
  
  for(timeID in seq_along(times)) {
    time <- times[timeID]
    newrow <- list(Signal=setName, Time=time, Conversions=sampleSet_times[[timeID]][[1]], Attempt=attemptNo)
    plot_df_list[[idx]] <- newrow
    idx <- idx+1
  }
  
}

plot_df <- do.call(rbind, lapply(plot_df_list, function(x) as.data.frame(t(x))))
plot_df$Time <- as.numeric(plot_df$Time)
plot_df$Conversions <- as.numeric(plot_df$Conversions)
plot_df$Signal <- as.character(plot_df$Signal)
plot_df$Attempt <- as.numeric(plot_df$Attempt)
plot_df$Signal_Full <- paste0(plot_df$Signal,"_",plot_df$Attempt)

image <- ggplot(plot_df[,], aes(x=Time, y=Conversions, color=Signal_Full)) +
  geom_point() +
  geom_path() +
  ggtitle(paste0(expName,"_10TopSignals"))
image

plot_fname <- file.path(plotDir,"figxx_zeb1_noise=0.08_transitionsVsTime_lineplot.pdf")
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()

# Summarize plot_df to aggregate data by Signal and Time
summary_plot_df <- plot_df %>%
  group_by(Signal, Time) %>%
  summarize(
    MeanConversions = mean(Conversions),
    SDConversions = sd(Conversions),
    SEConversions = SDConversions / sqrt(n()),
    .groups = "drop"
  )

# Create the line plot with error bars
image <- ggplot(summary_plot_df, aes(x = Time, y = MeanConversions, color = Signal, group = Signal)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = MeanConversions - SEConversions, ymax = MeanConversions + SEConversions), width = 0.5) +
  labs(
    title = paste0(expName, "_10TopSignals"),
    x = "Time",
    y = "Mean Conversions ± SE"
  ) +
  theme_minimal()
image

plot_fname <- file.path(plotDir,"figxx_zeb1_10reps_noise=0.08_transitionsVsTime_lineplot_AGG.pdf")
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()






######## MODEL-WISE TRANSITION PROB #####

# Identify unique signals
unique_signals <- unique(sigNames)

# Reorganize multiSet_times to group by signals and trials
signal_trials <- lapply(unique_signals, function(signal) {
  indices <- which(sigNames == signal)
  trials <- multiSet_times[indices]
  trials
})
names(signal_trials) <- unique_signals

# Initialize storage for results
variance_results <- list()
transition_results <- list()
outcome_counts <- list()
model_outcomes <- list()

# Loop through each signal
for (signal_name in unique_signals) {
  # Extract trials for the current signal
  trials <- signal_trials[[signal_name]]
  
  # Variance in NumTransitions across trials
  num_transitions <- sapply(trials, function(trial) {
    sapply(trial, function(timepoint) timepoint$NumTransitions)
  })
  variance_results[[signal_name]] <- apply(num_transitions, 1, var) # Variance across trials for each timepoint
  
  # Analyze transitions for each trial
  transition_results[[signal_name]] <- lapply(trials, function(trial) {
    sapply(trial, function(timepoint) {
      # Identify transitions for model clusters
      initial_clusters <- clust
      final_clusters <- timepoint$NewClusters
      data.frame(
        Initial = initial_clusters,
        Final = final_clusters
      )
    })
  })
  
  # Count outcomes at the final timepoint (1 -> 1, 1 -> 2, etc.)
  outcome_counts[[signal_name]] <- lapply(trials, function(trial) {
    final_timepoint <- trial[[length(trial)]]
    initial_clusters <- clust
    final_clusters <- final_timepoint$NewClusters
    
    table(factor(initial_clusters, levels = c(1, 2)),
          factor(final_clusters, levels = c(1, 2)))
  })
  
  # Analyze model outcomes across trials, focusing on 1->2 transitions
  model_outcomes[[signal_name]] <- do.call(rbind, lapply(1:length(clust), function(model_idx) {
    model_transitions <- sapply(trials, function(trial) {
      final_timepoint <- trial[[length(trial)]]
      initial <- clust[model_idx]
      final <- final_timepoint$NewClusters[model_idx]
      paste0(initial, "->", final)
    })
    
    transition_times <- lapply(trials, function(trial) {
      time_reached <- which(sapply(trial, function(tp) tp$NewClusters[model_idx]) == 2)
      if (length(time_reached) > 0) {
        list(FirstTime = times[min(time_reached)],
             Variance = var(times[time_reached]),
             EarliestTime = times[min(time_reached)],
             StaysIn2 = all(time_reached == seq(min(time_reached), length(trial))))
      } else {
        list(FirstTime = NA, Variance = NA, EarliestTime = NA, StaysIn2 = FALSE)
      }
    })
    
    # Focus on 1->2 transitions only
    count_1_to_2 <- sum(model_transitions == "1->2")
    
    data.frame(
      Model = model_idx,
      Count1to2 = count_1_to_2,
      AvgFirstTime = mean(unlist(lapply(transition_times, "[[", "FirstTime")), na.rm = TRUE),
      VarianceFirstTime = mean(unlist(lapply(transition_times, "[[", "Variance")), na.rm = TRUE),
      EarliestFirstTime = min(unlist(lapply(transition_times, "[[", "EarliestTime")), na.rm = TRUE),
      StaysIn2 = all(unlist(lapply(transition_times, "[[", "StaysIn2")))
    )
  }))
}

# Combine model outcomes into a single dataframe
model_outcome_df <- do.call(rbind, lapply(names(model_outcomes), function(signal) {
  data.frame(Signal = signal, model_outcomes[[signal]])
}))

# Filter for models with initial cluster 1
model_outcome_df <- model_outcome_df %>%
  group_by(Signal, Model) %>%
  mutate(Probability1to2 = Count1to2 / 10) # Calculate probability of 1->2 transitions

# Plot boxplot of 1->2 transition probabilities for models with initial cluster 1
image <- ggplot(model_outcome_df, aes(x = Signal, y = Probability1to2, fill = Signal)) +
  geom_violin() +
  labs(title = "Distribution of EMT Probabilities",
       x = "Signal",
       y = "Model EMT Probability") +
  theme_sticcc() +
  theme(axis.text.x = element_text(angle = 90, size=12))
image

plot_fname <- file.path(plotDir,"figxx_zeb1_noise=0.08_transitionProb_violin.pdf")
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()

# Visualize variance in NumTransitions across signals
variance_df <- do.call(rbind, lapply(names(variance_results), function(signal) {
  data.frame(
    Signal = signal,
    Timepoint = times,
    Variance = variance_results[[signal]]
  )
}))

image <- ggplot(variance_df, aes(x = Timepoint, y = Variance, color = Signal)) +
  geom_line() +
  labs(title = "EMT Variance Across Signals",
       x = "Timepoint",
       y = "Variance") +
  theme_sticcc()
image

plot_fname <- file.path(plotDir,"figxx_zeb1_noise=0.08_variance_linechart.pdf")
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()

# Visualize transition outcomes
outcome_df <- do.call(rbind, lapply(names(outcome_counts), function(signal) {
  outcome <- outcome_counts[[signal]]
  do.call(rbind, lapply(seq_along(outcome), function(trial_idx) {
    trial_outcome <- outcome[[trial_idx]]
    data.frame(
      Signal = signal,
      Trial = trial_idx,
      Initial = rep(c(1, 2), 2),
      Final = rep(c(1, 2), each = 2),
      Count = as.vector(trial_outcome)
    )
  }))
}))

image <- ggplot(outcome_df, aes(x = interaction(Initial, Final), y = Count, fill = Signal)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Transition Outcome Counts",
       x = "Transition (Initial -> Final)",
       y = "Count") +
  theme_sticcc() 
image

plot_fname <- file.path(plotDir,"figxx_zeb1_noise=0.08_outcomes_barchart.pdf")
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()


######## WHAT MODELS ARE NOISE-DRIVEN? ############

# How many models transit somewhere between 0 and 10 times?
length(which(model_outcome_df$Signal == "Zeb1_noise=0.08" & 
               model_outcome_df$Probability1to2 > 0 & 
               model_outcome_df$Probability1to2 < 1))
noise_models_list <- which(model_outcome_df$Signal == "Zeb1_noise=0.08" & 
                             model_outcome_df$Probability1to2 > 0 & 
                             model_outcome_df$Probability1to2 < 1)
summary(model_outcome_df[noise_models_list,"Probability1to2"])


######## ZEB1: EFFECT OF NOISE ############

det_compare_signal <- "Zeb1_noise=0"
## Read in deterministic time-series results
zeb_det_transitionTimes <- readRDS(file.path(topoDir,
                                             paste0("bhtopo_t=300_relax_OUnoise=0_tau=10_genes=1.2_SIG=",det_compare_signal,"_runNo=1"),
                                             paste0("transitionTimes_",det_compare_signal,".Rds")))

#### WHY WHY WHY is the number of transitions different in this new run? Is resultSet accurate?
# It matches for 3 of the signals (Snai1, Snai2, Zeb2) but not Zeb1 or Tgfbeta... whatthefuck
resultSet_full[which(resultSet_full$SetName == det_compare_signal),"ConvertingModels"]
resultSet_full[which(resultSet_full$SetName == gsub("0","0.04",det_compare_signal)),"ConvertingModels"]
zeb_det_transitionModels <- which(clust == 1 & zeb_det_transitionTimes[[length(times)]]$NewClusters == 2)


length(zeb_det_transitionModels) == resultSet_full[which(resultSet_full$SetName == det_compare_signal),"ConvertingModels"]


# Plot boxplot of 1->2 transition probabilities for models with initial cluster 1
image <- ggplot(model_outcome_df[which(model_outcome_df$Signal == gsub("0","0.08",det_compare_signal) & 
                                         !model_outcome_df$Model %in% zeb_det_transitionModels),], 
                aes(x = Signal, y = Probability1to2, fill = Signal)) +
  geom_violin() +
  labs(title = "Distribution of EMT Probabilities",
       x = "Signal",
       y = "Model EMT Probability") +
  theme_sticcc() +
  theme(axis.text.x = element_text(angle = 90, size=12))
image

plot_fname <- file.path(plotDir,"figxx_zeb1_noise=0.08_transitionProb_violin_stochOnly.pdf")
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()



# Deterministic
# interestingly, two models which DO transition deterministically, DO NOT always with noise
table(model_outcome_df[which(clust == 1 &
                              model_outcome_df$Signal == "Zeb1_noise=0.08" & 
                              model_outcome_df$Model %in% zeb_det_transitionModels),"Probability1to2"])
# stochastic
# also interesting: 76 models require noise, but at this noise level STRONGLY favor EMT over remaining
# noise...deterministically..amplifies the signal for these??
table(model_outcome_df[which(clust == 1 &
                               model_outcome_df$Signal == "Zeb1_noise=0.08" & 
                              !model_outcome_df$Model %in% zeb_det_transitionModels),"Probability1to2"])




# comparing to noise=0.04, I'm curious if models with transition probs 0 < p < 1 go to 1 at higher noise
noise0.04_uncertain_models <- c(31,  38,  91, 121, 187, 198, 205, 206, 224, 268,
                                281, 282, 285, 300, 303, 320, 335, 355, 361, 362,
                                384, 389, 402, 412, 430, 436, 444, 452, 457, 465, 476, 480)
summary(model_outcome_df[noise0.04_uncertain_models, "Probability1to2"])


# hypothesis: time-series simulations have a short relaxation phase (10), and some models may not fully relax
# this would lead to inflated transition counts in the time-series results
# and should be visible in the last timepoint

# plot final timepoint
zeb_det_lastT <- readRDS(file.path(topoDir,
                                   paste0("bhtopo_t=300_relax_OUnoise=0_tau=10_genes=1.2_SIG=",det_compare_signal,"_runNo=1"),
                                   paste0("racipe_",det_compare_signal,"_t=300_relaxed.Rds")))
zeb_det_lastT_norm <- assay(zeb_det_lastT)





zeb_det_lastT_norm <- as.data.frame(t(assay(zeb_det_lastT)))
zeb_det_lastT_norm[,1:length(genes)] <- log2(1+zeb_det_lastT_norm[,1:length(genes)]) # Log transform
zeb_det_lastT_norm[,1:length(genes)] <- sweep(zeb_det_lastT_norm[,1:length(genes)], 2, tmpMeans, FUN = "-") # scale
zeb_det_lastT_norm[,1:length(genes)] <- sweep(zeb_det_lastT_norm[,1:length(genes)], 2, tmpSds, FUN = "/") # scale
zeb_det_lastT_pca <- as.data.frame(predict(pca, zeb_det_lastT_norm[,names(tmpMeans)]))
zeb_det_lastT_pca


image <- ggplot(pca_df) +
  geom_point(aes(x=PC1,y=PC2)) + ## does PC2 need to flip?
  geom_point(data = zeb_det_lastT_pca, aes(x=PC1,y=PC2),color="red") +
  #geom_path(data = zeb_det_lastT_pca, aes(x=PC1,y=PC2)) +
  #geom_point(data = model_states_pca, aes(x=PC1,y=PC2),color="red", size=3) +
  theme_sticcc() +
  theme(axis.line = element_line(linewidth = 1, color = "black"), 
        axis.ticks = element_line(linewidth = 1, color="black"),
        axis.title = element_text(size=22),
        axis.text = element_text(size=18)
  ) +
  scale_color_gradient()
image


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

traj_model_list <- c(sample(noise_models_list, 3))#c(7,8)#c(14, 26, 41, 156)
noAttempts <- 10

#hd_traj_model <- 14
#attemptNo <- 4
trajPlots <- T

for(hd_traj_model in traj_model_list) {
  for(attemptNo in seq(1, noAttempts)) {
    
    # set up data storage
    modelDataDir <- file.path(dataDir,paste0("trajectories_model",hd_traj_model))
    if(!dir.exists(modelDataDir)) {
      dir.create(modelDataDir)
    }
    
    # Setup
    hd_traj_model_idx <- models_selected[paste0("Model",hd_traj_model)] # for retrieving relevant state
    model_states <- ss_unique[which(ss_unique$Model == hd_traj_model_idx), ]
    #hd_traj_model <- sample(which(clust==1), 1)
    # Sim parameters
    hd_traj_startTime <- 2
    hd_traj_simTime <- 100
    hd_traj_resolution <- 0.2
    
    traj_fname <- file.path(modelDataDir,paste0("trajectory_model",hd_traj_model,
                                                "_simTime=",hd_traj_simTime,
                                                "_SIG=",hd_traj_sig,"_2025-01-23_v",attemptNo,".Rds"))
    
    traj_pca_fname <- file.path(modelDataDir,paste0("trajectoryPCA_model",hd_traj_model,
                                                    "_simTime=",hd_traj_simTime,
                                                    "_SIG=",hd_traj_sig,"_2025-01-23_v",attemptNo,".Rds"))
    
    
    
    # run simulations if not already done
    if(!file.exists(traj_pca_fname) | !file.exists(traj_fname)) {
      # Load checkpoint
      #hd_traj_checkpoint <- readRDS(file.path(topoDir,expName_hd_traj,paste0("racipe_",hd_traj_sig,"_t=",hd_traj_startTime,".Rds")))
      hd_traj_checkpoint <- racipe_bistable_raw # actually just use WT object
      dim(hd_traj_checkpoint)
      
      # Filter selected model & set ICs (to WT steady state)
      hd_traj_p2 <- hd_traj_checkpoint[,hd_traj_model]
      hd_traj_p2@metadata$config$simParams["numModels"] <- 1
      sracipeIC(hd_traj_p2) <- assay(racipe_bistable_raw)[,as.numeric(hd_traj_model)] 
      
      # confirm parameters match
      orig_params <- as.numeric(sracipeParams(racipe_bistable_raw)[as.numeric(hd_traj_model),])
      new_params <- as.numeric(sracipeParams(hd_traj_p2))
      all.equal(new_params, orig_params)
      
      # confirm ICs match previous SS
      orig_SS <- as.numeric(unname(assay(racipe_bistable_raw)[,as.numeric(hd_traj_model)]))
      new_IC <- as.numeric(sracipeIC(hd_traj_p2))
      all.equal(orig_SS, new_IC)
      
      ## TODO
      # QC simulation (ensure steady state is the same, simulate for 10 @ 0.1)
      hd_traj_p2@metadata$config$simParams["nIC"] <- 1
      hd_traj_p2 <- sracipeSimulate(hd_traj_p2, timeSeries = T, printStart = 0, 
                                    printInterval = 0.1, simulationTime = 10,
                                    genParams = F, genIC = F, integrateStepSize = 0.01, initialNoise=0.0, scaledNoise = T, 
                                    stepper = "EM", ouNoise_t = signal_tcorr)
      
      
      
      # confirm parameters still match
      all.equal(as.numeric(sracipeParams(hd_traj_p2)), 
                orig_params)
      
      # confirm ICs still match previous SS
      all.equal(as.numeric(sracipeIC(hd_traj_p2)), 
                orig_SS )
      
      # check if new SS matches ICs and previous SS
      all.equal(rownames(sracipeIC(hd_traj_p2)), rownames(hd_traj_p2@metadata$timeSeriesDet))
      
      # early
      qc_ts <- as.numeric(unlist(as.data.frame(hd_traj_p2@metadata$timeSeriesDet)[,1]))
      all.equal(orig_SS, as.numeric(qc_ts))
      # mid
      qc_ts <- as.numeric(unlist(as.data.frame(hd_traj_p2@metadata$timeSeriesDet)[,10]))
      all.equal(orig_SS, as.numeric(qc_ts))
      # late
      qc_ts <- as.numeric(unlist(as.data.frame(hd_traj_p2@metadata$timeSeriesDet)[,101]))
      all.equal(orig_SS, as.numeric(qc_ts))
      
      # does it approach its other state?
      model_states <- ss_unique[which(ss_unique$Model == hd_traj_model_idx), ]
      
      model_states
      round(orig_SS, 1)
      round(qc_ts, 1)
      
      
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
      #all.equal(sracipeParams(hd_traj_p2), sracipeParams(racipeNorm)[hd_traj_model,]) 
      
      hd_traj_p2@metadata$config$simParams["nIC"] <- 1
      hd_traj_p2 <- sracipeSimulate(hd_traj_p2, timeSeries = T, printStart = 0, 
                                    printInterval = hd_traj_resolution, simulationTime = hd_traj_simTime,
                                    genParams = F, genIC = F, integrateStepSize = 0.1, initialNoise=0.04, scaledNoise = T, 
                                    stepper = "EM_Clamp", ouNoise_t = signal_tcorr, clampGenes=sig_clamp_gene_ids,
                                    clampValues=sig_clamp_df[hd_traj_model,])
      
      
      hd_traj_df <- as.data.frame(t(hd_traj_p2@metadata$timeSeries)) # for stoch
      #hd_traj_df <- as.data.frame(t(hd_traj_p2@metadata$timeSeriesDet)) # for det
      
      
      ## TODO
      # relaxation simulation (check stability of final state, simulate for 40 @ 0.1)
      hd_traj_relaxation <- hd_traj_p2
      sracipeIC(hd_traj_relaxation) <- as.numeric(t(as.data.frame(t(hd_traj_p2@metadata$timeSeries))[501,]))
      hd_traj_relaxation <- sracipeSimulate(hd_traj_relaxation, timeSeries = T, printStart = 0, 
                                            printInterval = hd_traj_resolution, simulationTime = hd_traj_simTime,
                                            genParams = F, genIC = F, integrateStepSize = 0.1, initialNoise=0.0, scaledNoise = T, 
                                            stepper = "EM_Clamp", ouNoise_t = signal_tcorr)
      
      #model_states
      #round(as.numeric(unlist(as.data.frame(hd_traj_relaxation@metadata$timeSeriesDet)[,101])), 1)
      
      hd_traj_df_relaxation <- as.data.frame(t(hd_traj_relaxation@metadata$timeSeries))
      hd_traj_df_relaxation$Phase <- "Relax"
      hd_traj_df_relaxation$Time <- as.numeric(gsub("X", "", rownames(hd_traj_df_relaxation)))
      
      hd_traj_df$Phase <- "Signal"
      hd_traj_df$Time <- as.numeric(gsub("X", "", rownames(hd_traj_df)))
      hd_traj_df_full <- rbind(hd_traj_df, hd_traj_df_relaxation)
      hd_traj_df_full <- hd_traj_df_full[,c(names(tmpMeans), "Phase","Time")]
      for(gene in genes) {
        hd_traj_df_full[,gene] <- as.numeric(hd_traj_df_full[,gene])
      }
      hd_traj_df_full[,1:length(genes)] <- log2(1+hd_traj_df_full[,1:length(genes)]) # Log transform
      hd_traj_df_full[,1:length(genes)] <- sweep(hd_traj_df_full[,1:length(genes)], 2, tmpMeans, FUN = "-") # scale
      hd_traj_df_full[,1:length(genes)] <- sweep(hd_traj_df_full[,1:length(genes)], 2, tmpSds, FUN = "/") # scale
      
      
      hd_traj_df_full[which(hd_traj_df_full$Phase == "Relax"),"Time"] <- hd_traj_df_full[which(hd_traj_df_full$Phase == "Relax"),"Time"] + max(hd_traj_df_full[which(hd_traj_df_full$Phase == "Signal"),"Time"])
      rownames(hd_traj_df_full) <- hd_traj_df_full$Time
      hd_traj_df_full <- hd_traj_df_full[,c("Time","Phase",genes_reordered)]
      
      
      # Add ICs
      #modelICs <- as.data.frame(t(sracipeIC(hd_traj_checkpoint[,hd_traj_model])))[,genes_reordered] # use checkpoint ICs
      modelICs <- as.data.frame(t(sracipeIC(hd_traj_p2)))[,genes_reordered] # use checkpoint SS
      modelICs[,names(tmpMeans)] <- log2(1+modelICs[,names(tmpMeans)]) # Log transform
      modelICs[,names(tmpMeans)] <- sweep(modelICs[,names(tmpMeans)], 2, tmpMeans, FUN = "-") # scale
      modelICs[,names(tmpMeans)] <- sweep(modelICs[,names(tmpMeans)], 2, tmpSds, FUN = "/") # scale
      
      modelICs$Time <- 0
      modelICs$Phase <- "IC"
      modelICs <- modelICs[,c("Time","Phase",genes_reordered)]
      
      hd_traj_df_full <- as.data.frame(rbind(modelICs, hd_traj_df_full))
      rownames(hd_traj_df_full)[1] <- "0"
      rownames(hd_traj_df_full) <- hd_traj_df_full$Time
      #heatmap(as.matrix(hd_traj_df[rev(rownames(hd_traj_df)),genes_reordered]), Rowv = NA)#, Colv = NA)
      
      
      # Plot zeb trajectory
      #ggplot(hd_traj_df, aes(x=Time,y=Zeb1)) +
      #  geom_line()
      
      hd_traj_df_long <- pivot_longer(hd_traj_df_full, cols=all_of(genes), 
                                      names_to = "Gene", values_to = "Expression")
      
    } else {
      hd_traj_pca <- readRDS(traj_pca_fname)
      hd_traj_df_long <- readRDS(traj_fname)
    }
    
    
    image <- ggplot(hd_traj_df_long[which(hd_traj_df_long$Gene %in% c("Zeb1","Cdh1","Vim","Cldn7","Np63a")),], 
                    aes(x=Time,y=Expression, group=Gene, color=Gene)) +
      geom_line() +
      theme_sticcc() +
      theme(axis.line = element_line(linewidth = 1, color = "black"), 
            axis.ticks = element_line(linewidth = 1, color="black"),
            axis.title = element_text(size=22),
            axis.text = element_text(size=18)
      )
    image
    
    if(trajPlots) {
      modelPlotDir <- file.path(plotDir,paste0("trajectories_model",hd_traj_model))
      if(!dir.exists(modelPlotDir)) {
        dir.create(modelPlotDir)
      }
      
      traj_plot_fname <- file.path(modelPlotDir, paste0("trajectory_model",hd_traj_model,
                                                        "_simTime=",hd_traj_simTime,
                                                        "_SIG=",hd_traj_sig,"_2025-01-23_v",attemptNo,".pdf"))
      pdf(traj_plot_fname, height = 10, width = 10)
      print(image)
      dev.off()
    }
    
    
    
    # Plot trajectory, with model states in red
    model_states_norm <- model_states
    model_states_norm[,1:length(genes)] <- log2(1+model_states_norm[,1:length(genes)]) # Log transform
    model_states_norm[,1:length(genes)] <- sweep(model_states_norm[,1:length(genes)], 2, tmpMeans, FUN = "-") # scale
    model_states_norm[,1:length(genes)] <- sweep(model_states_norm[,1:length(genes)], 2, tmpSds, FUN = "/") # scale
    model_states_pca <- as.data.frame(predict(pca, model_states_norm[,names(tmpMeans)]))
    model_states_pca
    
    ## NB LOGIC NEEDS FIXING HERE - these should be uncommented probably
    hd_traj_pca <- as.data.frame(predict(pca, hd_traj_df_full[,names(tmpMeans)]))
    hd_traj_pca$Time <- hd_traj_df_full$Time
    hd_traj_pca$Model <- hd_traj_model
    
    image <- ggplot(pca_df) +
      geom_point(aes(x=PC1,y=PC2)) + ## does PC2 need to flip?
      geom_point(data = hd_traj_pca, aes(x=PC1,y=PC2,color=Time)) +
      geom_path(data = hd_traj_pca, aes(x=PC1,y=PC2,color=Time)) +
      geom_point(data = model_states_pca, aes(x=PC1,y=PC2),color="red", size=3) +
      theme_sticcc() +
      theme(axis.line = element_line(linewidth = 1, color = "black"), 
            axis.ticks = element_line(linewidth = 1, color="black"),
            axis.title = element_text(size=22),
            axis.text = element_text(size=18)
      ) +
      scale_color_gradient()
    
    if(trajPlots) {
      traj_plot_pca_fname <- file.path(modelPlotDir, paste0("trajectoryPCA_model",hd_traj_model,
                                                            "_simTime=",hd_traj_simTime,
                                                            "_SIG=",hd_traj_sig,"_2025-01-23_v",attemptNo,".pdf"))
      pdf(traj_plot_pca_fname, height = 10, width = 10)
      print(image)
      dev.off()
    }
    
    
    
    # save trajectories
    if(trajPlots) {
      modelDataDir <- file.path(dataDir,paste0("trajectories_model",hd_traj_model))
      if(!dir.exists(modelDataDir)) {
        dir.create(modelDataDir)
      }
      
      traj_fname <- file.path(modelDataDir,paste0("trajectory_model",hd_traj_model,
                                                  "_simTime=",hd_traj_simTime,
                                                  "_SIG=",hd_traj_sig,"_2025-01-23_v",attemptNo,".Rds"))
      saveRDS(hd_traj_df_long, traj_fname)
      
      traj_pca_fname <- file.path(modelDataDir,paste0("trajectoryPCA_model",hd_traj_model,
                                                      "_simTime=",hd_traj_simTime,
                                                      "_SIG=",hd_traj_sig,"_2025-01-23_v",attemptNo,".Rds"))
      saveRDS(hd_traj_pca, traj_pca_fname)
      
    }
    
  }
}



######## DETAILED TRAJECTORY PLOTS ############
hd_traj_model <- 8
attemptNo <- 8

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
hd_traj_startTime <- 2
hd_traj_simTime <- 100
hd_traj_resolution <- 0.2

traj_fname <- file.path(modelDataDir,paste0("trajectory_model",hd_traj_model,
                                            "_simTime=",hd_traj_simTime,
                                            "_SIG=",hd_traj_sig,"_2025-01-05_v",attemptNo,".Rds"))

traj_pca_fname <- file.path(modelDataDir,paste0("trajectoryPCA_model",hd_traj_model,
                                                "_simTime=",hd_traj_simTime,
                                                "_SIG=",hd_traj_sig,"_2025-01-05_v",attemptNo,".Rds"))

hd_traj_pca <- readRDS(traj_pca_fname)
hd_traj_df_long <- readRDS(traj_fname)



# Define named subsets of genes
my_e_markers <- c("Cdh1","miR141","miR34a","miR101",
                  "miR200a","miR200b","miR200c","Ovol2","Grhl2","Cldn7")
my_m_markers <-  c("Vim","Snai1","Snai2",
                   "Twist1","Twist2","Foxc2","Gsc",
                   "Klf8","Tcf3","Tgfbeta")

e_markers_step1tgts <- my_e_markers[which(my_e_markers %in% zeb_1step_tgts)]
e_markers_step2tgts <- my_e_markers[which(my_e_markers %in% zeb_2step_tgts & 
                                            !my_e_markers %in% zeb_1step_tgts)]
e_markers_step3tgts <- my_e_markers[which(my_e_markers %in% zeb_3step_tgts & 
                                            !my_e_markers %in% zeb_1step_tgts &
                                            !my_e_markers %in% zeb_2step_tgts)]

m_markers_step1tgts <- my_m_markers[which(my_m_markers %in% zeb_1step_tgts)]
m_markers_step2tgts <- my_m_markers[which(my_m_markers %in% zeb_2step_tgts & 
                                            !my_m_markers %in% zeb_1step_tgts)]
m_markers_step3tgts <- my_m_markers[which(my_m_markers %in% zeb_3step_tgts & 
                                            !my_m_markers %in% zeb_1step_tgts &
                                            !my_m_markers %in% zeb_2step_tgts)]

gene_subsets <- list(
  'E Markers 1' = e_markers_step1tgts,
  #'E Markers 2'= e_markers_step2tgts,
  'E Markers 3' = e_markers_step3tgts,
  'M Markers 1' = m_markers_step1tgts,
  #'M Markers 2' = m_markers_step2tgts,
  'M Markers 3' = m_markers_step3tgts
)

gene_subsets <- list(
  'E Markers' = c("Cdh1","miR141","miR34a","miR101",
                  "miR200a","miR200b","miR200c","Ovol2","Grhl2","Cldn7"),
  'M Markers' = c("Vim","Snai1","Snai2",
                  "Twist1","Twist2","Foxc2","Gsc",
                  "Klf8","Tcf3","Tgfbeta"),
  'Clamp Genes' = c("Zeb1"),
  'Zeb2' = c('Zeb2')#,
  #'Disconnected' = c("Np63a","miR9","miR30c","miR205")
)

# Assign each gene to a subset
hd_traj_df_long$Subset <- sapply(hd_traj_df_long$Gene, function(gene) {
  subset_name <- names(gene_subsets)[sapply(gene_subsets, function(subset) gene %in% subset)]
  if (length(subset_name) > 0) subset_name else NA
})

# Filter out genes not in any subset
filtered_df <- hd_traj_df_long[!is.na(hd_traj_df_long$Subset),]

# Calculate average and standard error for each subset and time point
library(dplyr)
avg_traj_df <- filtered_df %>%
  group_by(Subset, Time) %>%
  summarize(
    Mean_Expression = mean(Expression),
    SE_Expression = sd(Expression) / sqrt(n()),
    .groups = "drop"
  )

# Plot the average trajectory with error bars
image <- ggplot(avg_traj_df, aes(x = Time, y = Mean_Expression, color = Subset)) +
  geom_line() +
  geom_ribbon(aes(ymin = Mean_Expression - SE_Expression, ymax = Mean_Expression + SE_Expression, fill = Subset), alpha = 0.2, color = NA) +
  theme_sticcc() +
  theme(
    axis.line = element_line(linewidth = 1, color = "black"), 
    axis.ticks = element_line(linewidth = 1, color = "black"),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18)
  ) +
  labs(y = "Average Expression", x = "Time", color = "Gene Subset", fill = "Gene Subset")

image

plot_fname <- file.path(modelPlotDir, paste0("trajectory_model",hd_traj_model,
                                             "_simTime=",hd_traj_simTime,
                                             "_SIG=",hd_traj_sig,"_2025-01-04_v",attemptNo,"_aggregated_zeb2.pdf"))
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()


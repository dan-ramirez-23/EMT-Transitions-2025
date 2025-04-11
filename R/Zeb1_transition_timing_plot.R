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


signal_simTime <- 500
signal_relaxTime <- 50
signal_nGenes <- c(1,2)
signal_noise <- c(0, 0.04, 0.2) #0.04
signal_tcorr <- 10
initClust <- 1
tgtClust <- 2


expName <- paste0("bhtopo_t=",signal_simTime,"_relax_OUnoise=",paste0(signal_noise, collapse = "."),
                  "_tau=",signal_tcorr,"_genes=",paste0(signal_nGenes,collapse = "."),"_CLAMPS_jan2025")

resultSet_fname <- file.path(topoDir,expName,"result_summary.Rds")
resultSet_full <- readRDS(resultSet_fname)
resultSet <- resultSet_full
rs_det <- resultSet_full[which(resultSet_full$Noise == 0),]

#selectedNoise <- c(0, 0.04) 
selectedNoise <- c(0.04)
#selectedNoise <- c(0.04)
resultSet <- resultSet_full[which(resultSet_full$Noise %in% selectedNoise),]

resultSet <- resultSet_full[which(resultSet_full$NumGenes == 1 & resultSet_full$Noise == 0.04),] 
n <- 1
nPer <- 10

setIDList <- c()
rankorder <- order(resultSet$ConversionPct, decreasing = T)
topn <- rankorder[which(!rankorder %in% setIDList)][1:n]
setIDList <- rep(c(setIDList, topn), 10)
sigNames <- resultSet[setIDList,"SetName"]
names(sigNames) <- setIDList

trajSimTime <- 300
times <- c(seq(2, 30, 2), seq(35, 100, 5), seq(120, trajSimTime, 20))


plotNoise_list <- c(0.04, 0.08)
plot_df_comb_list <- list()
for(plotNoise in plotNoise_list) {
  multiSet_fname <- file.path(dataDir,
                              paste0("transitions_vs_time_Zeb1_10replicates_noise=",plotNoise,"_",expName,".Rds"))
  multiSet_times <- readRDS(multiSet_fname)
  
  
  ## Plot conversions vs time
  plot_df_list <- list()
  idx <- 1
  attemptIdxList <- rep(1:nPer, each=n)
  for(setIDNo in seq_along(setIDList)) {
    setID <- setIDList[setIDNo]
    setName <- gsub("0.04",plotNoise,resultSet[setID, "SetName"])
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
  plot_df$Transitions <- as.numeric(plot_df$Conversions)
  plot_df$Signal <- as.character(plot_df$Signal)
  plot_df$Trial <- factor(as.numeric(plot_df$Attempt))
  plot_df$Signal_Full <- paste0(plot_df$Signal,"_",plot_df$Attempt)
  plot_df$Noise <- plotNoise
  
  
  plot_df_comb_list[[as.character(plotNoise)]] <- plot_df
  
}
plot_df_comb <- do.call(rbind, lapply(plot_df_comb_list, function(x) as.data.frame(x)))
plot_df_comb$Noise <- factor(plot_df_comb$Noise)
plot_df_comb$Trial <- factor(plot_df_comb$Trial)


image <- ggplot(plot_df_comb[,], aes(x=Time, y=Transitions, color=Noise, group=interaction(Noise, Trial))) +
  geom_point() +
  geom_line() +
  #ggtitle(paste0(expName,"_Zeb1"))
  ggtitle("Zeb1 EMT Induction Time") +
  theme_sticcc() +
  theme(panel.grid = element_line(color = "black",
                                  size = 0.5,
                                  linetype = 1),
        axis.line = element_line(color="black", size=1),
        axis.ticks = element_line(color="black", size=1),
        title = element_text(size=22))
image

plot_fname <- file.path(plotDir,"figS4_zeb1_10reps_1gene_noise=0.04.0.08_transitionsVsTime_lineplot.pdf")
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()

# Summarize plot_df to aggregate data by Signal and Time
summary_plot_df <- plot_df_comb %>%
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

plot_fname <- file.path(plotDir,"figxx_zeb1_10reps_noise=0.04.0.08_transitionsVsTime_lineplot_AGG.pdf")
pdf(plot_fname, height = 10, width = 10)
print(image)
dev.off()


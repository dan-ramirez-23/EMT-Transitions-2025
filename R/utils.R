
library(dplyr)
library(sRACIPE)#, lib.loc = "/Users/danramirez/localR/4.2.2-arm")
library(doFuture)
library(doRNG)
library(igraph)
source("R/utils_clamping.R")

# class used for KNN classification
library(class)

# Load topo in as dataframe
loadTopo <- function(topoName,
                     topoDir=NA) {
  ## Import topology
  if(is.na(topoDir)) {
    topoDir <- file.path(getwd(),"inputs")
  }
  
  if(!dir.exists(topoDir)) {
    dir.create(topoDir)
  }
  topo <- read.table(file.path(topoDir,paste0(topoName,".tpo")), header = T)
  topo$Source <- gsub("-","",topo$Source)
  topo$Target <- gsub("-","",topo$Target)
  return(topo)
}



plotNetwork <- function(topo, topoName, outputDir=NA, plot_suffix=NA) {
  if(is.na(outputDir)) {
    outputDir <- file.path(getwd(),topoName)
  }
  if(!dir.exists(outputDir)) {
    dir.create(outputDir)
  }
  
  if(is.na(plot_suffix)) {
    plot_suffix <- Sys.time()
  }
  
  
  net_file <- file.path(outputDir,paste0("networkVis",plot_suffix,".html"))
  
  topo[which(as.numeric(topo$Type) %% 2 == 0),"Type"] <- 2
  topo[which(as.numeric(topo$Type) %% 2 == 1),"Type"] <- 1
  
  
  node_list <-
    unique(c(topo[, 1], topo[, 2]))
  
  nodes <-
    data.frame(
      id = node_list,
      label = node_list,
      font.size = 50,
      value = c(rep(1, length(node_list)))
    )
  edge_col <- data.frame(c(1, 2), c("blue", "darkred"))
  arrow_type <- data.frame(c(1, 2), c("arrow", "circle"))
  colnames(arrow_type) <- c("type", "color")
  colnames(edge_col) <- c("type", "color")
  edges <-
    data.frame(
      from = c(topo[, 1]),
      to = c(topo[, 2]),
      arrows.to.type	= arrow_type$color[c(as.numeric(topo[, 3]))],
      width = 3,
      color = edge_col$color[c(as.numeric(topo[, 3]))]
    )
  
  networkPlot <-
    visNetwork::visNetwork(nodes, edges, height = "1000px", width = "100%") %>%
    visEdges(arrows = "to") %>%
    visOptions(manipulation = FALSE) %>%
    visLayout(randomSeed = 123) %>%
    #visNodes(scaling = list(label = list(enabled = T))) %>%
    visPhysics(solver = "forceAtlas2Based", stabilization = FALSE)
  visNetwork::visSave(networkPlot, file = net_file, selfcontained = FALSE)
  
}




gmmClust <- function(data, k) {
  
}

# take a matrix newdata and assign clusters to each row based on the supplied clustering
# return a vector of cluster IDs
assignClusters <- function(newdata, gmm, tmpMeans, tmpSds) {
  
  nCols <- ncol(newdata)
  
  newdata[,1:nCols] <- log2(1+newdata[,1:nCols]) # Log transform
  newdata[,1:nCols] <- sweep(newdata[,1:nCols], 2, tmpMeans, FUN = "-") # scale
  newdata[,1:nCols] <- sweep(newdata[,1:nCols], 2, tmpSds, FUN = "/") # scale
  
  
  out_clust <- predict_GMM(newdata, gmm$centroids, gmm$covariance_matrices, gmm$weights)
  
  ### OLD Method used KNN, now using GMM
  # out_clust <- rep(NA, nrow(newdata))
  # # loop over each model
  # for(model in 1:nrow(newdata)) {
  #   query <- t(as.data.frame(newdata[model,]))
  #   neighbors <- get.knnx(expr_norm, query, k = myK)$nn.index
  #   neighbors <- neighbors[which(neighbors != model)]
  #   neighborsTopClust <- as.numeric(names(sort(table(clust$cluster[neighbors]), decreasing = T))[1])
  #   out_clust[model] <- neighborsTopClust
  # }
  
  return(out_clust)
  
}




assignClustersMultinom <- function(newdata, model, tmpMeans, tmpSds) {
  
  nCols <- ncol(newdata)
  
  newdata[,1:nCols] <- log2(1+newdata[,1:nCols]) # Log transform
  newdata[,1:nCols] <- sweep(newdata[,1:nCols], 2, tmpMeans, FUN = "-") # scale
  newdata[,1:nCols] <- sweep(newdata[,1:nCols], 2, tmpSds, FUN = "/") # scale
  
  out_clust <- predict(model, data = newdata)
  #out_clust <- predict_GMM(newdata, gmm$centroids, gmm$covariance_matrices, gmm$weights)

  return(out_clust)
  
}



assignClustersPC1 <- function(newpca) {
  
  
  x <- newpca[,1]
  out_clust <- unlist(lapply(x, PC1Clustering))
  #out_clust <- predict_GMM(newdata, gmm$centroids, gmm$covariance_matrices, gmm$weights)
  
  return(out_clust)
  
}

PC1Clustering <- function(pc) {
  if(pc <= 0) {
    return(1)
  } else {
    return(2)
  }
}


# Take as input a clustered gene expression dataset and return a list of its genes, 
#         associating each gene with a cluster
idGenes <-  function(clust, # vector of cluster ids
                     initialClust,
                     targetClust, 
                     expr # expr should have samples/cells in rows, genes in columns
                     ) { 
  
  rs <- data.frame(gene=colnames(expr),dir=NA)
  
  for(gene in colnames(expr)) {
    
    avg_exp_c1 <- mean(expr[which(clust == initialClust),gene], na.rm=T)
    avg_exp_c2 <- mean(expr[which(clust == targetClust),gene], na.rm=T)
    diff <- avg_exp_c2 - avg_exp_c1
    
    
    if(is.na(diff)) {
      rs[which(rs$gene==gene),"dir"] <- NA
    } else if(diff >= 0) {
      rs[which(rs$gene==gene),"dir"] <- "Up"
    } else {
      rs[which(rs$gene==gene),"dir"] <- "Down"
    }
  }
  # Manual correction, #TODO: find better solution for this
  rs[which(rs$gene == "Goosecoid"),"dir"] <- "Up"
  rs$GeneID <- 1:nrow(rs)
  
  return(rs)
  
}


# Compute out-degree for every gene in a topology
topoOutDegrees <- function(geneinfo, 
                           topo) {
  
  outdeg <- rep(NA, dim(geneinfo)[1])
  for(i in 1:length(outdeg)) {
    outdeg[i] <- length(which(topo$Source == geneinfo[i,"gene"]))
  }
  
  geneinfo$outdegree <- outdeg
  
  
  # convert to graph and compute pagerank
  g <- graph_from_edgelist(as.matrix(topo[,c(1:2)]))
  pr <- page.rank(g)$vector
  
  # add pagerank to geneinfo
  # this is inefficient, but I can't remember the better way right now
  geneinfo$PageRank <- NA
  for(gene in 1:length(pr)) {
    geneinfo[which(geneinfo$gene == names(pr)[gene]),"PageRank"] <- unname(pr[gene])
  }
  
  
  return(geneinfo)
}

# use hyperparameters and output from idGenes to create parameter sets
# return a nested list where each top-level element is a paramset, containing:
  # sigGenes = list() of length nSigGenes
  # fcVals = list() of length nSigGenes
createParamSets <- function(fcVals = c(2,5,15),
                            targetClust=2,
                            nSigGenes = 2,
                            geneInfo,
                            nPerturbations=NA,
                            inverted=FALSE,
                            constantFC=T, # Disabling this vastly increase parameter search space, and is probably unrealistic
                            order = T
                            ) {
  

  ## Combination generation code from chatGPT
  ## prompt:
        # I am searching a parameter space of possible perturbations to a system. 
        # Given n total species, I would like to select k species to perturb, with a 
        # set of S different perturbed conditions for each species. Write R code to 
        # enumerate all combinations of k species at each of S signaling values.
  k = nSigGenes
  nr = dim(geneInfo)[1]
  S = fcVals
  
  # Generate all combinations of k species from n total species
  species_combinations <- combn(nr, k)
  
  if(!is.na(nPerturbations)) {
    # Select top ranked combinations by total out-degree
    comb_df <- as.data.frame(t(species_combinations))
    comb_df$OutDegree <- 0
    comb_df$PageRank <- 0
    for(c in 1:nSigGenes) {
      # complicated, but here we are adding the out-degree for each constituent gene in a signal
      comb_df$OutDegree <- comb_df$OutDegree + geneInfo[comb_df[,paste0("V",c)],"outdegree"] 
      comb_df$PageRank <- comb_df$PageRank + geneInfo[comb_df[,paste0("V",c)],"PageRank"] 
    }
    
    if(order) {
      keep_sets <- order(comb_df$OutDegree, decreasing = T)[1:nPerturbations]  
    } else {
      keep_sets <- order(comb_df$OutDegree, decreasing = T)[sample(1:nrow(comb_df), nPerturbations),]  
    }
    
    keep_sets <- keep_sets[which(!is.na(keep_sets))]
    out_deg_list <- comb_df[keep_sets,"OutDegree"]
    pr_list <- comb_df[keep_sets,"PageRank"]
    
    species_combinations <- as.data.frame(species_combinations[,keep_sets])
    if(c==1) {
      species_combinations <- t(species_combinations)
    }
    
  }
  
  # Initialize a list to hold the results
  all_combinations <- list()
  
  # Loop through each combination of species to perturb
  for (i in seq_len(ncol(species_combinations))) {
    current_combination <- species_combinations[, i]
    
    # Create all combinations of signaling values for the current combination of species
    signaling_combinations <- expand.grid(rep(list(S), k)) 
    if(constantFC) {
      signaling_combinations <- as.data.frame(signaling_combinations[apply(signaling_combinations, 1, function(row) length(unique(row)) == 1), ])
    } 
    
    # Attach the species combination to the signaling combinations
    result <- cbind(matrix(current_combination, nrow=nrow(signaling_combinations), ncol=k, byrow=TRUE),
                    signaling_combinations)
    
    colnames(result) <- c(paste("Species", seq_len(k)), paste("Signal", seq_len(k)))
    
    # Change signal sign based on direction
    for(g in current_combination) {
      if(geneInfo[g,"dir"] == "Down") {
        sigID = which(current_combination == g)
        result[,paste0("Signal ",sigID)] <-  1/result[,paste0("Signal ",sigID)]
      }
    }
    
    
    # Change signal based on parameter overried
    if(inverted) {
      for(idx in 1:nSigGenes) {
        result[,paste0("Signal ",idx)] <-  1/result[,paste0("Signal ",idx)]
      }
    }
    
    # Store the result in the list
    all_combinations[[i]] <- result
  }
  
  # If you want to merge all combinations into one data frame, you can do the following
  final_result <- do.call(rbind, all_combinations)
  
  
  # Rename signal genes
  for(r in 1:nSigGenes) {
    # Use match to find the index of each value in df1$int_column in df2$int_value
    index <- match(final_result[,paste0("Species ",r)], geneInfo$GeneID)
    
    # Replace df1$int_column with corresponding string values
    final_result[,paste0("Species ",r)] <- geneInfo$gene[index]
  }
  
  if(!is.na(nPerturbations)) {
    combinations_per_geneCombo = length(S)^nSigGenes
    if(constantFC) {
      combinations_per_geneCombo = 1
    }
    final_result$TotalOutDegree <- rep(out_deg_list, each=combinations_per_geneCombo)
    final_result$TotalPageRank <- rep(pr_list, each=combinations_per_geneCombo)
    final_result$FoldChange <- rep(fcVals, length(all_combinations))
  }

  
  return(final_result)
  
  
}



# Simulate from previous steady state under new conditions:
# Remove an edge by setting its threshold value to 1(?)

simulateST_EKD <- function(racipe,
                               pca,
                               clust,
                               initialClust,
                               targetClust,
                               outDir,
                               kd_src,
                               kd_tgt = NA,
                               expName = NA,
                               plot=F,
                               noise = 0.2,
                               ...
) {
  
  # Directory setup, variable QC
  if(!dir.exists(outDir)) {
    print("Error: outdir could not be found")
    return(NULL)
  } 
  
  if(is.na(expName)) {
    expName <- paste0("Param_Analysis_",Sys.Date())
    print(paste("Using experiment name",expName,"- please specify this to 
                consolidate future results for this analysis"))
  } 
  expDir <- file.path(outDir,expName)
  if(!dir.exists(expDir)) {
    dir.create(expDir, recursive = T)
  }
  
  plotDir <- file.path(expDir,"plots")
  if(!dir.exists(plotDir)) {
    dir.create(plotDir, recursive = T)
  }
  
  # Pre-treat racipe results, save normalization moments
  simExp <- assay(racipe, 1)
  simExp <- log2(1+simExp)
  tmpMeans <- rowMeans(simExp)
  tmpSds <- apply(simExp,1,sd)
  
  racipeNorm <- sracipeNormalize(racipe)
  racipeData <- as.data.frame(t(assay(racipeNorm)))
  
  
  # get gene info, out-degree
  geneinfo <- idGenes(clust, initialClust, targetClust, racipeData)
  geneinfo <- topoOutDegrees(geneinfo, 
                             topo = sracipeCircuit(racipe))
  
  
  # Get parameters across different clusters
  params_full <- sracipeParams(racipe)
  params_init <- params_full[which(clust == initialClust),]
  params_tgt <- params_full[which(clust == targetClust),]
  
  # Find determinative parameters by mean
  colmeans_init <- colMeans(params_init)
  colmeans_tgt <- colMeans(params_tgt)
  param_diffs <- colmeans_tgt - colmeans_init
  names(param_diffs) <- colnames(params_full)
  
  # Set new simulation parameters
  if(type == "exact") {
    # Check if there are enough existing models. If not, we will duplicate some (for now)
    nTgtModels <- length(which(clust == targetClust))
    nInitModels <- length(which(clust == initialClust))
    shortfall = nTgtModels - nInitModels
    
    if(shortfall >= 0) {
      # Copy over target model parameters to initial models
      keep_idx <- sample(1:nrow(params_tgt), nrow(params_init))
      newparams_init <- params_tgt[keep_idx,]
    } else {
      # Duplicate enough to make up the shortfall
      dupes <- params_tgt[sample(1:nrow(params_tgt), shortfall),]
      newparams_init <- rbind(params_init, dupes)
    }
    
    # Update parameter table in a duplicate racipe object
    newparams <- as.matrix(rbind(params_tgt, newparams_init))
    
    
  } else if(type == "mean") {
    
  }
  
  
  # Prepare new simulations
  racipe2 <- racipe 
  #sracipeParams(racipe2) <- params[,]
  sracipeParams(racipe2) <- newparams
  
  # Set ICs to WT steady states
  sracipeIC(racipe2) <- assay(racipe)[,]
  
  # Simulate if not done already
  fname_racipe2 <- file.path(expDir, paste0("racipe_",expName,"_","_noise=",noise,".Rds"))
  if(!file.exists(fname_racipe2)) {
    racipe2 <- sracipeSimulate(racipe2, genIC = F, genParams = F, simulationTime = 10, initialNoise=noise, nNoise=1, scaledNoise = T)
    saveRDS(racipe2, fname_racipe2)
  } else {
    racipe2 <- readRDS(fname_racipe2)
  }
  
  
  # Rescale data & assign clusters to perturbed data
  racipe2Norm <- as.data.frame(t(assay(racipe2,2)))
  nCols <- ncol(racipe2Norm)
  racipe2Norm[,1:nCols] <- log2(1+racipe2Norm[,1:nCols]) # Log transform
  racipe2Norm[,1:nCols] <- sweep(racipe2Norm[,1:nCols], 2, tmpMeans, FUN = "-") # scale
  racipe2Norm[,1:nCols] <- sweep(racipe2Norm[,1:nCols], 2, tmpSds, FUN = "/") # scale
  newpca <- scale(racipe2Norm, pca$center, pca$scale) %*% pca$rotation 
  
  newlabels <- assignClustersPC1(newpca)
  
  
  # Compute number of models in target cluster
  tgt_count <- length(which(newlabels == targetClust))
  
  # Compute effectiveness
  INIT_VIABLE <- length(which(clust == initialClust))
  INIT_TARGET <- length(which(clust == targetClust))
  effectiveness <- (tgt_count - INIT_TARGET) / INIT_VIABLE
  print(paste("Perfect ",type," perturbation is ",round(effectiveness, 2),
              " effective in ",expName))  
  
  
  
  # Plot if specified
  if(plot) {
    label_df <- data.frame(x=c(-3,3),y=c(3,3),text=NA)
    label_df$text <- c(unname(table(newlabels)[1]), unname(table(newlabels)[2]))
    
    
    image <- ggplot() +
      geom_point(data=pca_df, aes(x=PC1,y=PC2),color="gray", alpha=0.6) +
      geom_point(data=as.data.frame(newpca), aes(x=PC1, y=PC2, color=as.factor(clust))) +
      ggtitle(paste0("Parameter set: ",expName," with noise=",noise, 
                     "(Eff=",round(effectiveness, 2),")")) +
      guides(color=guide_legend("Cluster")) +
      geom_text(data=label_df, aes(x=x,y=y,label=text))
    
    plot_fname <- file.path(plotDir, paste0("pca_",expName,"_noise=",noise,"_byOldClusters.pdf"))
    pdf(plot_fname, 10, 10)
    print(image)
    dev.off()
  }
  
  # Save results (cluster assignment) to output list
  #results[[i]] <- newlabels
  
  
}



simulateST_Perfect <- function(racipe,
                             pca,
                             clust,
                             initialClust,
                             targetClust,
                             outDir,
                             type = "exact",
                             expName = NA,
                             plot=F,
                             noise = 0.2,
                             ...
) {
  
  # Directory setup, variable QC
  if(!dir.exists(outDir)) {
    print("Error: outdir could not be found")
    return(NULL)
  } 
  
  if(is.na(expName)) {
    expName <- paste0("Param_Analysis_",Sys.Date())
    print(paste("Using experiment name",expName,"- please specify this to 
                consolidate future results for this analysis"))
  } 
  expDir <- file.path(outDir,expName)
  if(!dir.exists(expDir)) {
    dir.create(expDir, recursive = T)
  }
  
  plotDir <- file.path(expDir,"plots")
  if(!dir.exists(plotDir)) {
    dir.create(plotDir, recursive = T)
  }
  
  # Pre-treat racipe results, save normalization moments
  simExp <- assay(racipe, 1)
  simExp <- log2(1+simExp)
  tmpMeans <- rowMeans(simExp)
  tmpSds <- apply(simExp,1,sd)
  
  racipeNorm <- sracipeNormalize(racipe)
  racipeData <- as.data.frame(t(assay(racipeNorm)))
  
  
  # get gene info, out-degree
  geneinfo <- idGenes(clust, initialClust, targetClust, racipeData)
  geneinfo <- topoOutDegrees(geneinfo, 
                             topo = sracipeCircuit(racipe))
  
  
  # Get parameters across different clusters
  params_full <- sracipeParams(racipe)
  params_init <- params_full[which(clust == initialClust),]
  params_tgt <- params_full[which(clust == targetClust),]
  
  # Find determinative parameters by mean
  colmeans_init <- colMeans(params_init)
  colmeans_tgt <- colMeans(params_tgt)
  param_diffs <- colmeans_tgt - colmeans_init
  names(param_diffs) <- colnames(params_full)
  
  # Set new simulation parameters
  if(type == "exact") {
    # Check if there are enough existing models. If not, we will duplicate some (for now)
    nTgtModels <- length(which(clust == targetClust))
    nInitModels <- length(which(clust == initialClust))
    shortfall = nTgtModels - nInitModels
    
    if(shortfall >= 0) {
      # Copy over target model parameters to initial models
      keep_idx <- sample(1:nrow(params_tgt), nrow(params_init))
      newparams_init <- params_tgt[keep_idx,]
    } else {
      # Duplicate enough to make up the shortfall
      dupes <- params_tgt[sample(1:nrow(params_tgt), shortfall),]
      newparams_init <- rbind(params_init, dupes)
    }
    
    # Update parameter table in a duplicate racipe object
    newparams <- as.matrix(rbind(params_tgt, newparams_init))
    
    
  } else if(type == "mean") {
    
  }
  
  
  # Prepare new simulations
  racipe2 <- racipe 
  #sracipeParams(racipe2) <- params[,]
  sracipeParams(racipe2) <- newparams
  
  # Set ICs to WT steady states
  sracipeIC(racipe2) <- assay(racipe)[,]
  
  # Simulate if not done already
  fname_racipe2 <- file.path(expDir, paste0("racipe_",expName,"_","_noise=",noise,".Rds"))
  if(!file.exists(fname_racipe2)) {
    racipe2 <- sracipeSimulate(racipe2, genIC = F, genParams = F, simulationTime = 10, initialNoise=noise, nNoise=1, scaledNoise = T)
    saveRDS(racipe2, fname_racipe2)
  } else {
    racipe2 <- readRDS(fname_racipe2)
  }
  
  
  # Rescale data & assign clusters to perturbed data
  racipe2Norm <- as.data.frame(t(assay(racipe2,2)))
  nCols <- ncol(racipe2Norm)
  racipe2Norm[,1:nCols] <- log2(1+racipe2Norm[,1:nCols]) # Log transform
  racipe2Norm[,1:nCols] <- sweep(racipe2Norm[,1:nCols], 2, tmpMeans, FUN = "-") # scale
  racipe2Norm[,1:nCols] <- sweep(racipe2Norm[,1:nCols], 2, tmpSds, FUN = "/") # scale
  newpca <- scale(racipe2Norm, pca$center, pca$scale) %*% pca$rotation 
  
  newlabels <- assignClustersPC1(newpca)
  
  
  # Compute number of models in target cluster
  tgt_count <- length(which(newlabels == targetClust))
  
  # Compute effectiveness
  INIT_VIABLE <- length(which(clust == initialClust))
  INIT_TARGET <- length(which(clust == targetClust))
  effectiveness <- (tgt_count - INIT_TARGET) / INIT_VIABLE
  print(paste("Perfect ",type," perturbation is ",round(effectiveness, 2),
              " effective in ",expName))  
  

  
  # Plot if specified
  if(plot) {
    label_df <- data.frame(x=c(-3,3),y=c(3,3),text=NA)
    label_df$text <- c(unname(table(newlabels)[1]), unname(table(newlabels)[2]))
    
    
    image <- ggplot() +
      geom_point(data=pca_df, aes(x=PC1,y=PC2),color="gray", alpha=0.6) +
      geom_point(data=as.data.frame(newpca), aes(x=PC1, y=PC2, color=as.factor(clust))) +
      ggtitle(paste0("Parameter set: ",expName," with noise=",noise, 
                     "(Eff=",round(effectiveness, 2),")")) +
      guides(color=guide_legend("Cluster")) +
      geom_text(data=label_df, aes(x=x,y=y,label=text))
    
    plot_fname <- file.path(plotDir, paste0("pca_",expName,"_noise=",noise,"_byOldClusters.pdf"))
    pdf(plot_fname, 10, 10)
    print(image)
    dev.off()
  }
  
  # Save results (cluster assignment) to output list
  #results[[i]] <- newlabels
  
  
  
  
  
}


simulateStateTransitions <- function(racipe,
                                     pca,
                                     clust,
                                     initialClust,
                                     targetClust,
                                     nSigGenes,
                                     outDir,
                                     expName = NA,
                                     fcVals = NA,
                                     paramType="G",
                                     plot=F,
                                     noise = 0.2,
                                     ...) {
  
  if(!dir.exists(outDir)) {
    print("Error: outdir could not be found")
    return(NULL)
  } 
  # dataDir <- file.path(outDir,"data","Parameter Perturbation")
  # if(!dir.exists(dataDir)) {
  #   dir.create(dataDir, recursive = T)
  # }
  
  if(is.na(expName)) {
    expName <- paste0("Param_Analysis_",Sys.Date())
    print(paste("Using experiment name",expName,"- please specify this to 
                consolidate future results for this analysis"))
  } 
  expDir <- file.path(outDir,expName)
  if(!dir.exists(expDir)) {
    dir.create(expDir, recursive = T)
  }
  
  plotDir <- file.path(expDir,"plots")
  if(!dir.exists(plotDir)) {
    dir.create(plotDir, recursive = T)
  }
  
  if(all(is.na(fcVals))) {
    fcVals = c(5)
  }
  
  # Pre-treat racipe results, save normalization moments
  simExp <- assay(racipe, 1)
  simExp <- log2(1+simExp)
  tmpMeans <- rowMeans(simExp)
  tmpSds <- apply(simExp,1,sd)
  
  racipeNorm <- sracipeNormalize(racipe)
  racipeData <- as.data.frame(t(assay(racipeNorm)))
  
  
  # get gene info, out-degree
  geneinfo <- idGenes(clust, initialClust, targetClust, racipeData)
  geneinfo <- topoOutDegrees(geneinfo, 
                             topo = sracipeCircuit(racipe))
  
  # prepare parameter sets
  inverted = FALSE
  if(paramType == "K") {
    inverted = TRUE
  }
  paramSets <- createParamSets(fcVals = fcVals,
                               targetClust=targetClust,
                               nSigGenes = nSigGenes,
                               geneInfo = geneinfo,
                               inverted = inverted,
                               ...)
  params <- sracipeParams(racipe)
  
  # simulate parameter sets
  results = list()
  for(i in 1:nrow(paramSets)) {
    print(paste0("Beginning parameter set ",i))
    
    racipe2 <- racipe 
    sracipeParams(racipe2) <- params[,]
    
    for(sigID in 1:nSigGenes) {
      sracipeParams(racipe2)[,paste0(paramType,"_",paramSets[i,paste0("Species ",sigID)])] <- 
        params[,paste0(paramType,"_",paramSets[i,paste0("Species ",sigID)])] * 
        paramSets[i,paste0("Signal ",sigID)]
    }
    # Set ICs to WT steady states
    sracipeIC(racipe2) <- assay(racipe)[,]
    
    # Simulate if not done already
    keep_idx <- which(!colnames(paramSets) %in% c("TotalOutDegree","TotalPageRank"))
    fname_racipe2 <- file.path(expDir, paste0("racipe_",paste(paramSets[i,keep_idx], collapse = "_"),"_","_noise=",noise,".Rds"))
    if(nSigGenes > 5) {
      fname_racipe2 <- file.path(expDir, paste0("racipe_",i,"_","_noise=",noise,".Rds"))
    }
    if(!file.exists(fname_racipe2)) {
      racipe2 <- sracipeSimulate(racipe2, genIC = F, genParams = F, simulationTime = 10, initialNoise=noise, nNoise=1, scaledNoise = T)
      saveRDS(racipe2, fname_racipe2)
    } else {
      racipe2 <- readRDS(fname_racipe2)
    }
    
    
    # Rescale data & assign clusters to perturbed data
    racipe2Norm <- as.data.frame(t(assay(racipe2,2)))
    nCols <- ncol(racipe2Norm)
    racipe2Norm[,1:nCols] <- log2(1+racipe2Norm[,1:nCols]) # Log transform
    racipe2Norm[,1:nCols] <- sweep(racipe2Norm[,1:nCols], 2, tmpMeans, FUN = "-") # scale
    racipe2Norm[,1:nCols] <- sweep(racipe2Norm[,1:nCols], 2, tmpSds, FUN = "/") # scale
    newpca <- scale(racipe2Norm, pca$center, pca$scale) %*% pca$rotation 
    
    newlabels <- assignClustersPC1(newpca)
    
    # Plot if specified
    if(plot) {
      label_df <- data.frame(x=c(-3,3),y=c(3,3),text=NA)
      label_df$text <- c(unname(table(newlabels)[1]), unname(table(newlabels)[2]))
      
      
      image <- ggplot() +
        geom_point(data=pca_df, aes(x=PC1,y=PC2),color="gray", alpha=0.6) +
        geom_point(data=as.data.frame(newpca), aes(x=PC1, y=PC2, color=as.factor(clust))) +
        ggtitle(paste0("Parameter set: ",paste(paramSets[i,], collapse = "_")," with noise=",noise)) +
        guides(color=guide_legend("Cluster")) +
        geom_text(data=label_df, aes(x=x,y=y,label=text))
      
      
      plot_fname <- file.path(plotDir, paste0("pca_",paste(paramSets[i,], collapse = "_"),"_noise=",noise,"_byOldClusters.pdf"))
      if(nSigGenes > 5) {
        plot_fname <- file.path(plotDir, paste0("pca_paramSet",i,"_noise=",noise,"_byOldClusters.pdf"))  
      }
      
      pdf(plot_fname, 10, 10)
      print(image)
      dev.off()
    }
    
    # Save results (cluster assignment) to output list
    results[[i]] <- newlabels
    
    
  }
  
  # Save all data
  saveRDS(results, file.path(expDir,"cluster_assignments_all.Rds"))
  
  
  
  # Compute efficiency and save paramSets
  clust_all_df <- as.matrix(do.call(rbind, results))
  
  # Compute number of models in target cluster
  tgt_counts <- countState(clust_all_df, targetClust)
  
  paramSets$EndTargetStates <- tgt_counts
  
  # Compute effectiveness
  INIT_VIABLE <- length(which(clust == initialClust))
  INIT_TARGET <- length(which(clust == targetClust))
  
  
  paramSets$Effectiveness <- (paramSets$EndTargetStates - INIT_TARGET) / INIT_VIABLE
  
  
  saveRDS(paramSets, file.path(topoDir,expName,"result_summary.Rds"))
  
}




countState <- function(df, target) {
  counts <- apply(df, 1, countOccurrence, target)
  return(counts)
}

countOccurrence <- function(col, target) {
  return(length(which(col==target)))
}




# New methods to compute signal efficacy
# 1. Function to calculate the Euclidean distance between two vectors.
euclidean_distance <- function(v1, v2) {
  return(sqrt(sum((v1 - v2)^2)))
}

# 2. Function to compute the nearest neighbor distance for one sample from matrix A to matrix B.
nearest_neighbor_distance <- function(sample, matrixB) {
  distances <- apply(matrixB, 2, function(col) euclidean_distance(sample, col))
  return(min(distances))
}

# 3. Compute the Chamfer distance between two matrices.
chamfer_distance <- function(matrixA, matrixB) {
  distances_AtoB <- apply(matrixA, 2, function(col) nearest_neighbor_distance(col, matrixB))
  distances_BtoA <- apply(matrixB, 2, function(col) nearest_neighbor_distance(col, matrixA))
  return(sum(distances_AtoB) + sum(distances_BtoA))
}





# @racipe - perturbed simulation
# @clust - cluster assignment vector
# @initialClust, targetClust - integers
# @pca - pca object - I wonder if we should compute cluster overlap in high-dim space or first PCs?
# 
# First, compare the final position of each cluster with the target cluster
  # Use indices to select them from racipe
  # Add pairwise distances to target to the loss function
# Second, compute the overlap between both clusters after signaling
  # Add distance b/w medians to loss
transitionLoss <- function(racipe,
                        clust,
                        initialClust,
                        targetClust,
                        noise,
                        anneal=F,
                        relax=F) {
  
  if(relax) {
    finalState <- assay(racipe)
  } else {
    finalState <- assay(racipe, as.character(noise))  
  }
  
  
  targetWT <- sracipeIC(racipe)[,which(clust == targetClust)]
  targetEndpt <- finalState[,which(clust == targetClust)]
  initialEndpt <- finalState[,which(clust == initialClust)]
  if(anneal) {
    targetEndpt <- assay(racipe, as.character(noise/2^7))[,which(clust == targetClust)]
    initialEndpt <- assay(racipe, as.character(noise/2^7))[,which(clust == initialClust)]
  }
  
  # Now, compute the loss values for the three components:
  loss_1 <- chamfer_distance(initialEndpt, targetWT)
  loss_2 <- chamfer_distance(targetEndpt, targetWT)
  loss_3 <- chamfer_distance(initialEndpt, targetEndpt)
  
  # 4. Total loss value
  overlap_loss <- loss_3
  target_loss <- max(loss_1, loss_2)
  total_loss <- overlap_loss + target_loss
  
  return(c(overlap_loss, target_loss, total_loss))
  
  
}
  



optimizeST <- function(racipe,
                       pca,
                       clust,
                       initialClust,
                       targetClust,
                       nSigGenes,
                       outDir,
                       expName = NA,
                       fcVals = NA,
                       paramType="G",
                       plot=F,
                       noise = NA,
                       forceRerun = F,
                       forceRecompute = F,
                       checkpointSize=25,
                       anneal=F,
                       randomParams=F, # If true, fcVals, nSigGenes, and noise should all be vectors of length 2 containing min and max values
                       totalPerturbations = 500, # Only considered if randomParams is TRUE
                       relax=F,
                       simTime = 10,
                       simTimeRelax = 10,
                       onlyParams = F,
                       ...) {
  
  if(!dir.exists(outDir)) {
    print("Error: outdir could not be found")
    return(NULL)
  } 
  
  if(is.na(expName)) {
    expName <- paste0("Param_Analysis_",Sys.Date())
    print(paste("Using experiment name",expName,"- please specify this to 
                consolidate future results for this analysis"))
  } 
  expDir <- file.path(outDir,expName)
  if(!dir.exists(expDir)) {
    dir.create(expDir, recursive = T)
  }
  
  plotDir <- file.path(expDir,"plots")
  if(!dir.exists(plotDir)) {
    dir.create(plotDir, recursive = T)
  }
  
  # Defaults
  if(all(is.na(fcVals))) {
    fcVals = c(5,15,50)
  }
  if(all(is.na(nSigGenes))) {
    nSigGenes = c(1,2,3)
  }
  if(all(is.na(noise))) {
    noise = c(0.2,0.5,2,5)
  }
  
  cbPalette <- palette.colors(palette = "Okabe-Ito")[2:9]
  names(cbPalette) <- c(1:8)
  
  # Pre-treat racipe results, save normalization moments
  simExp <- assay(racipe, 1)
  simExp <- log2(1+simExp)
  tmpMeans <- rowMeans(simExp)
  tmpSds <- apply(simExp,1,sd)
  
  racipeNorm <- sracipeNormalize(racipe)
  racipeData <- as.data.frame(t(assay(racipeNorm)))
  
  
  # get gene info, out-degree
  geneinfo <- idGenes(clust, initialClust, targetClust, racipeData)
  geneinfo <- topoOutDegrees(geneinfo, 
                             topo = sracipeCircuit(racipe))
  
  # Get WT parameters
  params <- sracipeParams(racipe)  

  # iterate over nSigGenes to make a complete paramSet
  inverted = FALSE
  if(paramType == "K") {
    inverted = TRUE
  }
  
  # Generate parameter sets
  if(randomParams) {
    paramSets <- createRandomParamSets(numConditions = totalPerturbations,
                                       fcVals = fcVals,
                                       targetClust=2,
                                       nSigGenes = nSigGenes,
                                       noiseInterval = noise,
                                       geneInfo = geneinfo,
                                       ...
    ) 
  } else {
    paramSets <- optimizeSTParams(fcVals,
                                  targetClust,
                                  nSigGenes,
                                  geneinfo,
                                  inverted,
                                  noise,
                                  ...)   
  }
  
  
  
  # Check for an existing paramSet file under this expName
  # If one is found, merge and remove duplicates
  old_pset_fname <- file.path(topoDir,expName,"result_summary.Rds")
  if(file.exists(old_pset_fname)) {
    paramSet_old <- readRDS(old_pset_fname)
    # If random, discard new params to keep the total size true
    if(randomParams) {
      paramSets <- paramSet_old
    } 
    
    # By default, keep old results. If forcing parameter is true, discard these and recompute
    if(!forceRerun) {
      paramSets <- dplyr::bind_rows(paramSet_old, paramSets)
      paramSets <- paramSets[!duplicated(paramSets$SetName),]  
    } else {
      paramSets <- dplyr::bind_rows(paramSets, paramSet_old)
      paramSets <- paramSets[!duplicated(paramSets$SetName),] 
    }
    
  }
  
  paramSets <- paramSets[order(paramSets$SignalPower, decreasing = T),]
  if(all(noise) == 0) {
    paramSets <- paramSets[order(paramSets$TotalOutDegree, decreasing = T),]
  }
  
  if(onlyParams) {
    return(paramSets)
  }
  
  
  
  # simulate parameter sets
  currentMinLoss <- min(c(paramSets$TotalLoss, Inf), na.rm=T)
  newMinLoss <- F
  results = list()
  print(paste0("Beginning parameter set list of size ",nrow(paramSets)))
  for(i in 1:nrow(paramSets)) {
    
    if(paramSets[i,"Simulated"] & !forceRecompute){
      print(paste0("Skipping set ", i))
      next
    }
    
    
    racipe2 <- racipe 
    sracipeParams(racipe2) <- params[,]
    num <- paramSets[i,"NumGenes"]
    noisei <- paramSets[i,"Noise"]
    
    for(sigID in 1:num) {
      sracipeParams(racipe2)[,paste0(paramType,"_",paramSets[i,paste0("Species ",sigID)])] <- 
        params[,paste0(paramType,"_",paramSets[i,paste0("Species ",sigID)])] * 
        paramSets[i,paste0("Signal ",sigID)]
    }
    # Set ICs to WT steady states
    sracipeIC(racipe2) <- assay(racipe)[,]
    
    # Simulate if not done already
    #keep_idx <- which(!colnames(paramSets) %in% c("TotalOutDegree","TotalPageRank"))
    #fname_racipe2 <- file.path(expDir, paste0("racipe_",paste(paramSets[i,keep_idx], collapse = "_"),"_","_noise=",noise,".Rds"))
    #if(nSigGenes > 5) {
    fname_racipe2 <- file.path(expDir, paste0("racipe_",paramSets[i,"SetName"],".Rds"))
    #}
    if(!file.exists(fname_racipe2) | forceRerun) {
      if(anneal) {
        racipe2 <- sracipeSimulate(racipe2, genIC = F, genParams = F, simulationTime = simTime, 
                                   initialNoise=noisei, nNoise=8, scaledNoise = T, anneal = T, 
                                   integrateStepSize = 0.2)  
      } else {
        racipe2 <- sracipeSimulate(racipe2, genIC = F, genParams = F, simulationTime = simTime,
                                   initialNoise=noisei, nNoise=1, scaledNoise = T, 
                                   integrateStepSize = 0.2, simDet=F)
      }
      
      saveRDS(racipe2, fname_racipe2)
    } else {
      racipe2 <- readRDS(fname_racipe2)
    }
    
    
    if(relax) {
      racipe_relax_fname <- file.path(expDir, paste0("racipe_",paramSets[i,"SetName"],"_relaxed.Rds"))
      if(!file.exists(racipe_relax_fname) | forceRerun) {
        
        # Remove signal
        racipe_relax <- racipe2 
        sracipeParams(racipe_relax) <- params[,]
        # Set ICs to perturbed steady states
        sracipeIC(racipe_relax) <- assay(racipe2, as.character(noisei))[,]
        
        
        racipe_relax <- sracipeSimulate(racipe_relax, genIC = F, genParams = F, simulationTime = simTimeRelax,
                                   integrateStepSize = 0.2, nNoise = 0, iNoise = 0, nCores = 1, simDet = T)
        saveRDS(racipe_relax, racipe_relax_fname)
      } else {
        racipe_relax <- readRDS(racipe_relax_fname)
      }
      
      racipe2 <- racipe_relax
      racipe2Norm <- as.data.frame(t(assay(racipe2)))
      
    } else {
      racipe2Norm <- as.data.frame(t(assay(racipe2,as.character(noisei))))
    }
    
    
    
    # Rescale data & assign clusters to perturbed data
    
    nCols <- ncol(racipe2Norm)
    racipe2Norm[,1:nCols] <- log2(1+racipe2Norm[,1:nCols]) # Log transform
    racipe2Norm[,1:nCols] <- sweep(racipe2Norm[,1:nCols], 2, tmpMeans, FUN = "-") # scale
    racipe2Norm[,1:nCols] <- sweep(racipe2Norm[,1:nCols], 2, tmpSds, FUN = "/") # scale
    newpca <- scale(racipe2Norm, pca$center, pca$scale) %*% pca$rotation 
    
    # newlabels <- assignClustersPC1(newpca) # PC1-based clustering replaced with KNN
    newlabels <- knn_classifier(racipe2Norm, racipeData, clust, k=25)
    
    # Save results (cluster assignment) to output list
    results[[i]] <- newlabels
    
    
    # Compute loss - temporarily removed 11/20/23 for efficiency
    # Also, we can get similar information from # transitions
    # lossVec <- transitionLoss(racipe2,
    #                        clust,
    #                        initialClust,
    #                        targetClust,
    #                        noise=noisei,
    #                        anneal=anneal,
    #                        relax=relax)
    numTransitions <- length(which(newlabels == targetClust)) - length(which(clust == targetClust))
    
    #paramSets[i,"OverlapLoss"] <- lossVec[1]
    #paramSets[i,"TargetLoss"] <- lossVec[2]
    #paramSets[i,"TotalLoss"] <- lossVec[3]
    paramSets[i,"Simulated"] <- TRUE  
    paramSets[i,"EndTargetStates"] <- length(which(newlabels == targetClust))
    
    
    if(i %% checkpointSize == 0) {
      print(paste0("Checkpoint at parameter set ",i))
      saveRDS(paramSets, file.path(topoDir,expName,"result_summary.Rds"))
    }
    
    
    
    
    
    # Plot if specified, or if there's a new best
    if(i == 1) {
      #currentMinLoss <- lossVec[3] + 1
      currentMinLoss <- numTransitions - 1
      newMinLoss <- TRUE
    }
    if(numTransitions > currentMinLoss) {
      newMinLoss <- TRUE
      currentMinLoss <- numTransitions
    }
    
    if(plot | newMinLoss) {
      label_df <- data.frame(x=c(-3,3),y=c(3,3),text=NA)
      label_df$text <- c(unname(table(newlabels)[1]), unname(table(newlabels)[2]))
      
      
      image <- ggplot() +
        geom_point(data=as.data.frame(pca$x), aes(x=PC1,y=PC2),color="gray", alpha=0.6) +
        geom_point(data=as.data.frame(newpca), aes(x=PC1, y=PC2, color=as.factor(clust))) +
        ggtitle(paste0("Parameter set: ",paramSets[i,"SetName"])) +
        guides(color=guide_legend("Cluster")) +
        scale_color_manual(values=cbPalette) +
        geom_text(data=label_df, aes(x=x,y=y,label=text))
      
      
      plot_fname <- file.path(plotDir, paste0("pca_",paramSets[i,"SetName"],"_byOldClusters.pdf"))
      
      pdf(plot_fname, 10, 10)
      print(image)
      dev.off()
      
      newMinLoss <- FALSE
    }
    
    
  }
  
  saveRDS(paramSets, file.path(topoDir,expName,"result_summary.Rds"))
  
  

}







optimizeSTParams <- function(fcVals,
                             targetClust,
                             nSigGenes,
                             geneinfo,
                             inverted,
                             noise,
                             ...) {
  paramSets = data.frame()
  for(n in nSigGenes) {
    paramSets_new <- createParamSets(fcVals = fcVals,
                                     targetClust=targetClust,
                                     nSigGenes = n,
                                     geneInfo = geneinfo,
                                     inverted = inverted,
                                     ...)
    
    paramSets_new$NumGenes <- n
    paramSets = dplyr::bind_rows(paramSets, paramSets_new)
  }
  
  # Create a duplicate set for each noise level
  paramSets$Noise <- 0
  tmpList = list()
  for(n in seq_along(noise)) {
    ni <- noise[n]
    paramSets$Noise <- ni
    tmpList[[n]] <- paramSets
  }
  
  paramSets <- do.call(rbind, tmpList)
  
  # Create names for later
  relevant_columns <- grep("^Species|^Signal", colnames(paramSets))
  sub_mat <- paramSets[, relevant_columns]
  
  # Apply the function to each row
  names_vector <- apply(sub_mat, 1, create_name, colnames=colnames(sub_mat))
  paramSets$SetName <- paste0(names_vector,"_noise=",paramSets$Noise)
  
  
  # Leave empty columns to fill in
  paramSets$OverlapLoss <- NA
  paramSets$TargetLoss <- NA
  paramSets$TotalLoss <- NA
  paramSets$Simulated <- FALSE
  paramSets$EndTargetStates <- NA
  
  # experimental
  paramSets$SignalPower <- paramSets$TotalPageRank * paramSets$FoldChange * paramSets$Noise
  
  # reorder
  paramSets <- paramSets[order(paramSets$SignalPower, decreasing = T),]
  
  return(paramSets)
  
  
}


# Create the unified name function
create_name <- function(row, colnames) {
  # Modify values for 'Signal' columns
  for (i in seq_along(row)) {
    col_name <- colnames[i]
    
    # If the column name starts with "Signal" and the value is not NA, round to 2 significant figures
    if (startsWith(col_name, "Signal") && !is.na(row[i])) {
      row[i] <- signif(as.numeric(row[i]), 2)
    }
  }
  
  # Extract only non-NA values
  values <- row[!is.na(row)]
  
  # Paste the values together with "_" separator
  paste(values, collapse = "_")
}





createRandomParamSets <- function(numConditions,
                                  fcVals = c(1,100),
                                  targetClust=2,
                                  nSigGenes = c(1,10),
                                  noiseInterval = c(0,10),
                                  geneInfo,
                                  nPerturbations=NA,
                                  inverted=FALSE,
                                  constantFC=T # Disabling this vastly increases parameter search space, and is probably unrealistic
) {
  
  
  ## Find out how many genes (possible signal constituents) are present
  nr = dim(geneInfo)[1]
  
  ## First, generate numConditions uniform samples on fcVals, nSigGenes, and noiseInterval
  fcList <- runif(numConditions, min=min(fcVals), max=max(fcVals))
  ngList <- sample(nSigGenes[1]:nSigGenes[2], numConditions, replace=T)
  noiseList <- runif(numConditions, min=min(noiseInterval), max=max(noiseInterval))
  
  
  ## Combine these into a dataframe paramSets with column names FoldChange, NumGenes, Noise.
  ## Add metadata columns OutDegree, PageRank
  paramSets <- data.frame(FoldChange=fcList, NumGenes=ngList, Noise=noiseList)
  paramSets$OutDegree <- 0
  paramSets$PageRank <- 0
  
  ## Add columns 2 times the max NumGenes (Species 1, Signal 1, Species 2, Signal 2, etc), initialize with NA
  for(k in 1:max(nSigGenes)) {
    paramSets[,paste0("Species ",k)] <- rep(NA, numConditions)
    paramSets[,paste0("Signal ",k)] <- rep(NA, numConditions)
  }
  
  ## Iterate over each row
  for(i in 1:nrow(paramSets)) {
    fc <- paramSets[i,"FoldChange"]
    nSigGenes <- paramSets[i,"NumGenes"]
    
    # Sample gene indices from 1:nr
    sigGenes <- sample(1:nr, nSigGenes)
    
    # Order genes to make the sampling order-invariant
    sigGenes <- sigGenes[order(sigGenes)]
    
    # Iterate over signal genes and populate info
    for(c in seq_along(1:length(sigGenes))) {
      # Get gene index and name
      g <- sigGenes[c]
      gName <- geneInfo[g,"gene"]
      
      # Fill in gene indices
      paramSets[i,paste0("Species ",c)] <- gName
      paramSets[i,paste0("Signal ",c)] <- fc
      
      # Compute total outdegree and pagerank by adding info from geneinfo
      # here we are adding the out-degree for each constituent gene in a signal
      paramSets[i,"OutDegree"] <- paramSets[i,"OutDegree"] + geneInfo[g,"outdegree"] 
      paramSets[i,"PageRank"] <- paramSets[i,"PageRank"] + geneInfo[g,"PageRank"] 
      
      # Change signal sign based on direction
      if(geneInfo[g,"dir"] == "Down") {
        paramSets[i,paste0("Signal ",c)] <-  1/paramSets[i,paste0("Signal ",c)]
      }

      # Change signal based on parameter overried
      if(inverted) {
        paramSets[i,paste0("Signal ",c)] <-  1/paramSets[i,paste0("Signal ",c)]
      }
    }
  }
  
  
  
  
  
  
  # Create names for later
  relevant_columns <- grep("^Species|^Signal", colnames(paramSets))
  sub_mat <- paramSets[, relevant_columns]
  
  # Apply the function to each row
  names_vector <- apply(sub_mat, 1, create_name, colnames=colnames(sub_mat))
  paramSets$SetName <- paste0(names_vector,"_noise=",round(paramSets$Noise),2)
  
  
  # Leave empty columns to fill in
  paramSets$OverlapLoss <- NA
  paramSets$TargetLoss <- NA
  paramSets$TotalLoss <- NA
  paramSets$Simulated <- FALSE
  paramSets$EndTargetStates <- NA
  
  # experimental
  paramSets$SignalPower <- paramSets$PageRank * paramSets$FoldChange * paramSets$Noise
  
  return(paramSets)
  
}


# From phind.com, 11/2/2023
# Prompt: I have a topology file with 3 columns, Source, Target and Type. 
#         Write a function in R that takes two inputs: the topology dataframe 
#         and a vector list of nodes. Return a named list of nodes with the names being 
#         the same as the second input (the order must be preserved) and the vaues being 
#         "Source" for nodes with only outgoing edges, "Sink" for those with only incoming 
#         edges, and "Node" for those with both outgoing and incoming edges.
## NB: Per Phind's terms of service (https://www.phind.com/terms), Phind answers are not considered Phind content 
classify_nodes <- function(df, nodes) {
  # Initialize an empty list to store the results
  result <- list()
  
  # Iterate over each node
  for (node in nodes) {
    # Get the rows where the current node is the source
    source_rows <- df[df$Source == node,]
    
    # Get the rows where the current node is the target
    target_rows <- df[df$Target == node,]
    
    # Classify the node based on the type of edges
    if (nrow(source_rows) > 0 & nrow(target_rows) > 0) {
      result[[node]] <- "Node"
    } else if (nrow(source_rows) > 0) {
      result[[node]] <- "Source"
    } else if (nrow(target_rows) > 0) {
      result[[node]] <- "Sink"
    } else {
      result[[node]] <- "Isolated"
    }
  }
  
  # Return the result list
  return(result)
}



# utility function
knn_classifier <- function(newdata, wt_data, clust, k) {
  # Ensure that clust is a factor, as required by knn function
  clust <- as.factor(clust)
  
  # Perform KNN classification
  predicted_clusters <- knn(train = wt_data, test = newdata, cl = clust, k = k)
  
  # Return the list of predicted cluster IDs
  return(predicted_clusters)
}





# Compute transition-correlated genes

# Desired output is a matrix where each row is a paramSet, each column is a gene
# Workflow:
# For each paramset/signal (iterate over indices, but extract expName as paramSets[i,"SetName"]:
#   read simulation data from this condition (condition_data) which will be stored in file.path(getwd(),topoName,collectionName,paste0("racipe_))
#   re-assign clusters with function knn_classifier(condition_data, wt_data, clust, k=25)
#   generate a list of the model IDs which moved from initialClust to targetClust
#       targetClust (convertingModels), as well as those that did not (remainingModels)
#   for each gene/parameter (using lapply):
#     compute average expression on convertingModels and remainingModels
#     compute fold change between convertingModels/remainingModels
#   add FC values to the row for this paramSet
# return matrix of FC values
genes_x_transitions <- function(paramSets, # dataframe w/ columns: ModelNo, SetName
                                topoName, # string
                                collectionName, # string
                                initialClust, # int
                                targetClust, # int
                                wt_data, # matrix/dataframe of size numSamples x numFeatures
                                clust, # vector of length numSamples containing integers
                                clust_all,
                                tmpMeans,
                                tmpSds,
                                control = "non-target"
) {
  # Initialize a matrix to store the fold change values
  fc_matrix <- list()
  # Add results columns
  paramSets$InitialTargetStates <- NA
  paramSets$ConvertingModels <- NA
  paramSets$RemainingModels <- NA
  paramSets$RebelliousModels <- NA
  
  # Iterate over each parameter set
  for (i in 1:nrow(paramSets)) {
    # Extract the experiment name
    expName <- paramSets$SetName[i]
    
    # Extract noise level
    noise <- paramSets$Noise[i]
    
    # Read simulation data from file
    filepath <- file.path(getwd(), topoName, collectionName, paste0("racipe_", expName,".Rds"))
    filepath_final <- file.path(getwd(), topoName, collectionName, paste0("racipe_", expName,"_relaxed.Rds"))
    if(!file.exists(filepath) | !file.exists(filepath_final)) {
      next
    }
    condition_data <- readRDS(filepath)
    condition_data_final <- readRDS(filepath_final)
    #condition_data <- sracipeNormalize(condition_data)
    #condition_data <- assay(racipe, i=as.numeric(which(names(condition_data@assays) == noise)))
    
    condition_data <- as.data.frame(t(sracipeIC(condition_data)))
    condition_data_final <- as.data.frame(t(assay(condition_data_final)))
    
    
    #condition_data <- as.data.frame(t(assay(condition_data,as.character(noise))))

    # Rescale data & assign clusters to perturbed data
    nCols <- ncol(condition_data)
    condition_data[,1:nCols] <- log2(1+condition_data[,1:nCols]) # Log transform
    condition_data[,1:nCols] <- sweep(condition_data[,1:nCols], 2, tmpMeans, FUN = "-") # scale
    condition_data[,1:nCols] <- sweep(condition_data[,1:nCols], 2, tmpSds, FUN = "/") # scale
    
    condition_data_final[,1:nCols] <- log2(1+condition_data_final[,1:nCols]) # Log transform
    condition_data_final[,1:nCols] <- sweep(condition_data_final[,1:nCols], 2, tmpMeans, FUN = "-") # scale
    condition_data_final[,1:nCols] <- sweep(condition_data_final[,1:nCols], 2, tmpSds, FUN = "/") # scale
    #newpca <- scale(condition_data, pca$center, pca$scale) %*% pca$rotation 

    # Re-assign clusters
    new_clust <- knn_classifier(condition_data_final, wt_data, clust_all, k = 25)
    
    # Identify converting and remaining models
    paramSets$startingTargetPopulation[i] <- length(which(clust == targetClust))
    paramSets$startingInitPopulation[i] <- length(which(clust == initialClust))
    
    convertingModels <- condition_data[new_clust == targetClust & clust != targetClust, ]
    paramSets$ConvertingModels[i] <- nrow(convertingModels)
    
    if(control == "initial") {
      remainingModels <- condition_data[new_clust == initialClust & clust == initialClust, ]
    } else if (control == "non-target") {
      remainingModels <- condition_data[new_clust != targetClust & clust != targetClust, ]
    }
    paramSets$RemainingModels[i] <- nrow(remainingModels)
    
    rebelliousModels <- condition_data[new_clust != targetClust & clust == targetClust, ]
    paramSets$RebelliousModels[i] <- nrow(rebelliousModels)
    
    # Calculate log-fold change for each gene
    fc_values <- sapply(colnames(condition_data), function(gene) {
      avg_expr_converting <- mean(convertingModels[[gene]], na.rm = TRUE)
      avg_expr_remaining <- mean(remainingModels[[gene]], na.rm = TRUE)

      return(avg_expr_converting - avg_expr_remaining)
      
    })
    
    # Add the fold change values to the list
    fc_matrix[[i]] <- fc_values
  }
  
  # Convert the list to a data.frame and return
  fc_df <- do.call(rbind, fc_matrix)
  colnames(fc_df) <- colnames(condition_data)
  
  out_list <- list(fc_df=fc_df, resultSet=paramSets)
  
  return(out_list)
}

# Example usage:
# paramSets <- data.frame(ModelNo = 1:10, SetName = LETTERS[1:10]) # Define the dataframe with the required columns
# topoName <- "exampleTopo" # Replace with the actual topoName
# collectionName <- "exampleCollection" # Replace with the actual collectionName
# initialClust <- 1 # Replace with the actual initialClust
# targetClust <- 2 # Replace with the actual targetClust
# wt_data <- matrix(data = runif(100), nrow = 10) # Define the wild-type data as a dataframe or matrix
# clust <- sample(1:2, 10, replace = TRUE) # Define the vector of clusters

# fc_values <- genes_x_transitions(paramSets, topoName, collectionName, initialClust, targetClust, wt_data, clust)





filterModels <- function(racipe,
                         pca,
                         clust,
                         initialClust,
                         targetClust,
                         expName,
                         setName,
                         suffix=NA,
                         outDir,
                         desiredType = "Converting") {

  
  expDir <- file.path(outDir,expName)
  if(!dir.exists(expDir)) {
    dir.create(expDir, recursive = T)
  }
  
  plotDir <- file.path(expDir,"plots")
  if(!dir.exists(plotDir)) {
    dir.create(plotDir, recursive = T)
  }
  
  cbPalette <- palette.colors(palette = "Okabe-Ito")[2:9]
  names(cbPalette) <- c(1:8)
  
  # Pre-treat racipe results, save normalization moments
  simExp <- assay(racipe, 1)
  simExp <- log2(1+simExp)
  tmpMeans <- rowMeans(simExp)
  tmpSds <- apply(simExp,1,sd)
  
  racipeNorm <- sracipeNormalize(racipe)
  racipeData <- as.data.frame(t(assay(racipeNorm)))
  
  
  racipe_relax_fname <- file.path(expDir, paste0("racipe_",setName,"_relaxed.Rds"))
  racipe2 <- readRDS(racipe_relax_fname)
  racipe2Norm <- as.data.frame(t(assay(racipe2)))
  
  
  
  # Rescale data & assign clusters to perturbed data
  nCols <- ncol(racipe2Norm)
  racipe2Norm[,1:nCols] <- log2(1+racipe2Norm[,1:nCols]) # Log transform
  racipe2Norm[,1:nCols] <- sweep(racipe2Norm[,1:nCols], 2, tmpMeans, FUN = "-") # scale
  racipe2Norm[,1:nCols] <- sweep(racipe2Norm[,1:nCols], 2, tmpSds, FUN = "/") # scale
  newpca <- scale(racipe2Norm, pca$center, pca$scale) %*% pca$rotation 
  
  # newlabels <- assignClustersPC1(newpca) # PC1-based clustering replaced with KNN
  newlabels <- knn_classifier(racipe2Norm, racipeData, clust, k=25)
  
  
  # Identify converting and remaining models
  convertingModels <- which(newlabels == targetClust & clust != targetClust)
  
  remainingModels <- which(newlabels != targetClust & clust != targetClust)
  
  rebelliousModels <- which(newlabels != targetClust & clust == targetClust)
  
  
  if(desiredType == "Converting") {
    return(convertingModels)
  } else if(desiredType == "Remaining") {
    return(remainingModels)
  } else{
    return(rebelliousModels)
  }
  
  
}



plotCondition <- function(racipe,
                          pca,
                          clust,
                          initialClust,
                          targetClust,
                          expName,
                          setName,
                          plotDir=NA,
                          expDir=NA,
                          suffix=NA,
                          tmpMeans=NA,
                          tmpSds=NA,
                          outDir
                          ) {
  
  if(is.na(expDir)) {
    expDir <- file.path(outDir,expName)  
  }
  if(!dir.exists(expDir)) {
    dir.create(expDir, recursive = T)
  }
  
  if(is.na(plotDir)) {
    plotDir <- file.path(expDir,"plots")
  }
  if(!dir.exists(plotDir)) {
    dir.create(plotDir, recursive = T)
  }
  
  cbPalette <- palette.colors(palette = "Okabe-Ito")[2:9]
  names(cbPalette) <- c(1:8)
  
  # Pre-treat racipe results, save normalization moments
  simExp <- assay(racipe, 1)
  simExp <- log2(1+simExp)
  if(all(is.na(tmpMeans))) {
    tmpMeans <- rowMeans(simExp)  
  }
  if(all(is.na(tmpSds))) {
    tmpSds <- apply(simExp,1,sd)
  }
  
  
  
  racipeNorm <- sracipeNormalize(racipe)
  racipeData <- as.data.frame(t(assay(racipeNorm)))
  
  
  racipe_relax_fname <- file.path(expDir, paste0("racipe_",setName,"_relaxed.Rds"))
  racipe2 <- readRDS(racipe_relax_fname)
  racipe2Norm <- as.data.frame(t(assay(racipe2)))
  
  
  
  # Rescale data & assign clusters to perturbed data
  nCols <- ncol(racipe2Norm)
  racipe2Norm[,1:nCols] <- log2(1+racipe2Norm[,1:nCols]) # Log transform
  racipe2Norm[,1:nCols] <- sweep(racipe2Norm[,1:nCols], 2, tmpMeans, FUN = "-") # scale
  racipe2Norm[,1:nCols] <- sweep(racipe2Norm[,1:nCols], 2, tmpSds, FUN = "/") # scale
  newpca <- scale(racipe2Norm, pca$center, pca$scale) %*% pca$rotation 
  
  # newlabels <- assignClustersPC1(newpca) # PC1-based clustering replaced with KNN
  newlabels <- knn_classifier(racipe2Norm, racipeData, clust, k=25)
  
  
  
  # Plot
  label_df <- data.frame(x=c(-3,3),y=c(3,3),text=NA)
  label_df$text <- c(unname(table(newlabels)[1]), unname(table(newlabels)[2]))
  
  
  image <- ggplot() +
    geom_point(data=pca_df, aes(x=PC1,y=PC2),color="gray", alpha=0.6) +
    geom_point(data=as.data.frame(newpca), aes(x=PC1, y=PC2, color=as.factor(clust))) +
    ggtitle(paste0("Parameter set: ",setName)) +
    guides(color=guide_legend("Cluster")) +
    scale_color_manual(values=cbPalette) +
    geom_text(data=label_df, aes(x=x,y=y,label=text))
  
  if(is.na(suffix)) {
    plot_fname <- file.path(plotDir, paste0("pca_",setName,"_byOldClusters.pdf"))  
  } else {
    plot_fname <- file.path(plotDir, paste0("pca_",setName,"_byOldClusters_",suffix,".pdf"))
  }
  
  
  pdf(plot_fname, 10, 10)
  print(image)
  dev.off()
  
  
  
}




plotTransition <- function(racipe,
                           pca,
                           clust,
                           modelID,
                           initialClust,
                           targetClust,
                           expName,
                           setName,
                           suffix=NA,
                           outDir) {
  
  
  expDir <- file.path(outDir,expName)
  if(!dir.exists(expDir)) {
    dir.create(expDir, recursive = T)
  }
  
  plotDir <- file.path(expDir,"plots")
  if(!dir.exists(plotDir)) {
    dir.create(plotDir, recursive = T)
  }
  
  cbPalette <- palette.colors(palette = "Okabe-Ito")[2:9]
  names(cbPalette) <- c(1:8)
  
  # Pre-treat racipe results, save normalization moments
  simExp <- assay(racipe, 1)
  simExp <- log2(1+simExp)
  tmpMeans <- rowMeans(simExp)
  tmpSds <- apply(simExp,1,sd)
  
  racipeNorm <- sracipeNormalize(racipe)
  racipeData <- as.data.frame(t(assay(racipeNorm)))
  
  modelInit <- pca$x[modelID,c("PC1","PC2")]
  
  
  racipe_relax_fname <- file.path(expDir, paste0("racipe_",setName,"_relaxed.Rds"))
  racipe2 <- readRDS(racipe_relax_fname)
  racipe2Norm <- as.data.frame(t(assay(racipe2)))
  
  
  
  # Rescale data & assign clusters to perturbed data
  nCols <- ncol(racipe2Norm)
  racipe2Norm[,1:nCols] <- log2(1+racipe2Norm[,1:nCols]) # Log transform
  racipe2Norm[,1:nCols] <- sweep(racipe2Norm[,1:nCols], 2, tmpMeans, FUN = "-") # scale
  racipe2Norm[,1:nCols] <- sweep(racipe2Norm[,1:nCols], 2, tmpSds, FUN = "/") # scale
  newpca <- scale(racipe2Norm, pca$center, pca$scale) %*% pca$rotation 
  modelFinal <- newpca[modelID,c("PC1","PC2")]
  
  # newlabels <- assignClustersPC1(newpca) # PC1-based clustering replaced with KNN
  newlabels <- knn_classifier(racipe2Norm, racipeData, clust, k=25)
  
  # df for selected model
  modelDF <- cbind(modelInit, modelFinal)
  
  
  # Plot
  image <- ggplot() +
    geom_point(data=pca_df, aes(x=PC1,y=PC2),color="gray", alpha=0.6) +
    geom_point(data=as.data.frame(newpca), aes(x=PC1, y=PC2, color=as.factor(clust))) +
    geom_segment(aes(x=modelInit[1], y=modelInit[2], xend=modelFinal[1], yend=modelFinal[2]),
                 arrow = arrow(length=unit(0.5, 'cm')),size=2) +
    ggtitle(paste0("Model no. ",modelID," from parameter set: ",setName)) +
    guides(color=guide_legend("Cluster")) +
    scale_color_manual(values=cbPalette)
  
  if(is.na(suffix)) {
    plot_fname <- file.path(plotDir, paste0("pca_",setName,"_model=",modelID,".pdf"))  
  } else {
    plot_fname <- file.path(plotDir, paste0("pca_",setName,"_model=",modelID,"_",suffix,".pdf"))
  }
  
  
  pdf(plot_fname, 10, 10)
  print(image)
  dev.off()
  
  
}




replace_first_substring <- function(strings, replacement_map) {
  # Process each string in the list
  sapply(strings, function(s) {
    # Split the string into substrings
    substrings <- unlist(strsplit(s, "_"))
    
    # Replace the first substring if it exists in the replacement map
    if (substrings[1] %in% names(replacement_map)) {
      substrings[1] <- replacement_map[[substrings[1]]]
    }
    
    # Join the substrings back together
    paste(substrings, collapse = "_")
  })
}


# NB: Rows should be samples
normalizeExprMat <- function(exprMat, tmpMeans=NA, tmpSds=NA, giveMoments=F) {
  
  nCols <- ncol(exprMat)
  # Log-transform 
  exprMat <- log2(exprMat)
  
  # Standardize according to data features or provided moments
  if(all(c(is.na(tmpMeans), is.na(tmpSds)))) {
    tmpMeans <- rowMeans(log2(t(pyracipe_raw_states)))
    tmpSds <- apply(log2(t(pyracipe_raw_states)),1,sd)
  }
  
  exprMat[,1:nCols] <- sweep(exprMat[,1:nCols], 2, tmpMeans, FUN = "-") # scale
  exprMat[,1:nCols] <- sweep(exprMat[,1:nCols], 2, tmpSds, FUN = "/") # scale
  
  # Return normalized data
  if(giveMoments) {
    return(list(NormData=exprMt, SDs=tmpSds, Means=tmpMeans))
  } 
  return(exprMat)
  
}





# Classify transition between states
classify_transitions <- function(initClust, finalClust, diffs, targetClust, useDiffs=F, ...) {
  
  out <- rep("", length(diffs))
  if(useDiffs) {
    diff_verdicts <- diff_threshold(diffs, ...)  
  }
  
  
  #out[which(initClust == targetClust & finalClust == targetClust & diff_verdicts == 0)] <- "Target->Target (same)"
  #out[which(initClust == targetClust & finalClust == targetClust & diff_verdicts == 1)] <- "Target->Target (diff)"
  out[which(initClust == targetClust & finalClust == targetClust)] <- "Target->Target"
  out[which(initClust == targetClust & finalClust != targetClust)] <- "Rebellious"
  
  #out[which(initClust != targetClust & finalClust != targetClust & diff_verdicts == 0)] <- "Init->Init (same)"
  #out[which(initClust != targetClust & finalClust != targetClust & diff_verdicts == 1)] <- "Init->Init (diff)"
  out[which(initClust != targetClust & finalClust != targetClust)] <- "Init->Init"
  out[which(initClust != targetClust & finalClust == targetClust)] <- "Init->Target"
 
  return(out)
  
}


# Return 1 if the state is different above a threshold
diff_threshold <- function(diffs, threshold = 0.05) {
  out <- rep(1, length(diffs))
  out[which(diffs < threshold)] = 0
  return(out)
}



cells_x_signals <- function(paramSets, # dataframe w/ columns: ModelNo, SetName
                            topoName, # string
                            collectionName, # string
                            initialClust, # int
                            targetClust, # int
                            wt_data, # matrix/dataframe of size numSamples x numFeatures
                            clust, # vector of length numSamples containing integers
                            clust_all,
                            tmpMeans,
                            tmpSds,
                            InitRawStates,
                            ...) {
  # Initialize a matrix to store the cell fates
  out_matrix <- list()
  names <- list()
  
  # Iterate over each parameter set
  for (i in 1:nrow(paramSets)) {
    if(i %% 10 == 0) {
      print(paste0(i,"/",nrow(paramSets)))
    }
    
    # Extract the experiment name
    expName <- paramSets$SetName[i]
    
    # Extract noise level
    noise <- paramSets$Noise[i]
    
    # Read simulation data from file
    filepath <- file.path(getwd(), topoName, collectionName, paste0("racipe_", expName,"_relaxed.Rds"))
    if(!file.exists(filepath)) {
      next
    }
    condition_data <- readRDS(filepath)
    condition_data <- as.data.frame(t(assay(condition_data)))
    final_raw_states <- condition_data
    
    # Rescale data & assign clusters to perturbed data
    nCols <- ncol(condition_data)
    condition_data[,1:nCols] <- log2(1+condition_data[,1:nCols]) # Log transform
    condition_data[,1:nCols] <- sweep(condition_data[,1:nCols], 2, tmpMeans, FUN = "-") # scale
    condition_data[,1:nCols] <- sweep(condition_data[,1:nCols], 2, tmpSds, FUN = "/") # scale
    
    # Re-assign clusters
    new_clust <- knn_classifier(condition_data, wt_data, clust_all, k = 25)
    
    # Extract raw init and final states
    diffs <- rowSums(abs(as.matrix(final_raw_states) - as.matrix(InitRawStates)))
    verdicts <- classify_transitions(initClust = clust, finalClust = new_clust, 
                                    diffs = diffs, targetClust = targetClust, ...)
    
    
    # Add the fold change values to the list
    out_matrix[[i]] <- verdicts
    names[[i]] <- paramSets$SetName[i]
  }
  
  # Convert the list to a data.frame and return
  out_df <- as.data.frame(do.call(cbind, out_matrix))
  colnames(out_df) <- unlist(names)
  
  return(out_df)
}






##
calcTransitionRate <- function(paramSets,
                               setID,
                               racipe,
                               pca,
                               clust,
                               initialClust,
                               targetClust,
                               sigName,
                               outDir,
                               expName = NA,
                               plot=F,
                               noise = NA,
                               forceRerun = F,
                               forceRecompute = F,
                               anneal=F,
                               paramType = "G",
                               relax=T,
                               simTimes = c(10, 20, 50, 100, 200, 500, 1000, 2000),
                               simTimeRelax = 20,
                               save = F,
                               tcorr = NA,
                               scaledNoise = NA,
                               clamp=F) {
  
  # Setup directories
  if(!dir.exists(outDir)) {
    print("Error: outdir could not be found")
    return(NULL)
  } 
  
  if(is.na(expName)) {
    expName <- paste0("Param_Analysis_",Sys.Date())
    print(paste("Using experiment name",expName,"- please specify this to 
                consolidate future results for this analysis"))
  } 
  expDir <- file.path(outDir,expName)
  if(!dir.exists(expDir)) {
    dir.create(expDir, recursive = T)
  }
  
  plotDir <- file.path(expDir,"plots")
  if(!dir.exists(plotDir)) {
    dir.create(plotDir, recursive = T)
  }
  
  
  # Defaults
  cbPalette <- palette.colors(palette = "Okabe-Ito")[2:9]
  names(cbPalette) <- c(1:8)
  
  # Pre-treat racipe results, save normalization moments
  simExp <- assay(racipe, 1)
  simExp <- log2(1+simExp)
  tmpMeans <- rowMeans(simExp)
  tmpSds <- apply(simExp,1,sd)
  
  racipeNorm <- sracipeNormalize(racipe)
  racipeData <- as.data.frame(t(assay(racipeNorm)))
  
  
  # get gene info, out-degree
  geneinfo <- idGenes(clust, initialClust, targetClust, racipeData)
  geneinfo <- topoOutDegrees(geneinfo, 
                             topo = sracipeCircuit(racipe))
  
  # Get WT parameters
  params <- sracipeParams(racipe)  
  
  # Pre-allocate results data structure
  rs_list <- list()
  
  # Iterate over simulation times
  for(tID in seq_along(simTimes)) {
    simTime <- simTimes[tID]
    print(paste0("\nBeginning simulations with time ",simTime))
    
    # Modify parameters
    racipe2 <- racipe 
    sracipeParams(racipe2) <- params[,]
    num <- paramSets[setID,"NumGenes"]
    noisei <- paramSets[setID,"Noise"]
    if(!is.null(paramSets[setID,"Tau"])) {
      tcorr <- paramSets[setID,"Tau"]  
    }
    if(!is.null(paramSets[setID,"ScaledNoise"])) {
      scaledNoise <- paramSets[setID,"ScaledNoise"]  
    }
    
    # add signal
    if(clamp) {
      sig_clamp_genes <- getClampGenes(paramSets, setID)
      sig_clamp_gene_ids <- unlist(lapply(sig_clamp_genes, function(x) which(rownames(racipe) == x)))
      sig_clamp_df <- getClampDF(clamp_df, sig_clamp_genes, targetClust)
      colnames(sig_clamp_df) <- as.numeric(sig_clamp_gene_ids)
      sig_clamp_df <- as.matrix(sig_clamp_df)

      
    } else {
      for(sigID in 1:num) {
        sracipeParams(racipe2)[,paste0(paramType,"_",paramSets[setID,paste0("Species ",sigID)])] <- 
          params[,paste0(paramType,"_",paramSets[setID,paste0("Species ",sigID)])] * 
          paramSets[setID,paste0("Signal ",sigID)]
      }
    }

    # Set ICs to WT steady states
    sracipeIC(racipe2) <- assay(racipe)[,]
    
    
    # Simulate if not done already
    fname_racipe2 <- file.path(expDir, paste0("racipe_",paramSets[setID,"SetName"],"_t=",simTime,".Rds"))
    #}
    if(!file.exists(fname_racipe2) | forceRerun) {
      if(anneal) {
        racipe2 <- sracipeSimulate(racipe2, genIC = F, genParams = F, simulationTime = simTime, 
                                   initialNoise=noisei, nNoise=8, scaledNoise = T, anneal = T, 
                                   integrateStepSize = 0.2)  
      } else {
        # Check for the latest possible checkpoint we can use
        if(tID > 1) {
          # Step 1: Filter strings that start with "racipe" and end with a number followed by ".Rds"
          dirFiles <- dir(expDir)
          checkpointFiles <- grep("^racipe.*[0-9]\\.Rds$", dirFiles, value = TRUE)
          # Step 2: Extract the number that comes after "t=" and before an underscore
          timesSimulated <- as.numeric(regmatches(checkpointFiles, regexpr("(?<=t=)[0-9]+", checkpointFiles, perl = TRUE)))
          # Step 3: Continue from longest simulation not including present one
          useTime <- max(timesSimulated[which(timesSimulated < simTime)])
          
          # Read checkpoint data & insert state as initial conditions
          fname_checkpoint <- file.path(expDir, paste0("racipe_",paramSets[setID,"SetName"],"_t=",useTime,".Rds"))
          racipe_checkpoint <- readRDS(fname_checkpoint)
          
          print(paste0("Resuming from checkpoint at t=",useTime))
          
          sracipeIC(racipe2) <- assay(racipe_checkpoint, as.character(noisei))[,]
          
          # parameters should be equal, but check
          if(!all.equal(sracipeParams(racipe2), sracipeParams(racipe_checkpoint))) {
            print("Error: parameters don't match! Debug calcTransitionRate immediately")
          }
          
        }
        
        if(clamp) {
          racipe2 <- sracipeSimulate(racipe2, genIC = F, genParams = F, simulationTime = simTime,
                                     initialNoise=noisei, nNoise=1, scaledNoise = scaledNoise, nIC = 1,
                                     integrateStepSize = 0.2, simDet=F,
                                     stepper = "EM_Clamp", ouNoise_t = tcorr, clampGenes=sig_clamp_gene_ids,
                                     clampValues=sig_clamp_df)
          
          
        } else {
          racipe2 <- sracipeSimulate(racipe2, genIC = F, genParams = F, simulationTime = simTime,
                                     initialNoise=noisei, nNoise=1, scaledNoise = scaledNoise, 
                                     integrateStepSize = 0.2, simDet=F, stepper = "EM_OU", ouNoise_t=tcorr)
        }
        
      }
      
      saveRDS(racipe2, fname_racipe2)
    } else {
      racipe2 <- readRDS(fname_racipe2)
    }
    
    
    if(relax) {
      racipe_relax_fname <- file.path(expDir, paste0("racipe_",paramSets[setID,"SetName"],"_t=",simTime,"_relaxed.Rds"))
      if(!file.exists(racipe_relax_fname) | forceRerun) {
        
        # Remove signal
        racipe_relax <- racipe2 
        sracipeParams(racipe_relax) <- params[,]
        # Set ICs to perturbed steady states
        sracipeIC(racipe_relax) <- assay(racipe2, as.character(noisei))[,]
        
        
        racipe_relax <- sracipeSimulate(racipe_relax, genIC = F, genParams = F, simulationTime = simTimeRelax,
                                        integrateStepSize = 0.2, nNoise = 0, iNoise = 0, nCores = 1, simDet = T)
        saveRDS(racipe_relax, racipe_relax_fname)
      } else {
        racipe_relax <- readRDS(racipe_relax_fname)
      }
      
      racipe2 <- racipe_relax
      racipe2Norm <- as.data.frame(t(assay(racipe2)))
      
    } else {
      racipe2Norm <- as.data.frame(t(assay(racipe2,as.character(noisei))))
    }
    
    
    
    # Rescale data & assign clusters to perturbed data
    nCols <- ncol(racipe2Norm)
    racipe2Norm[,1:nCols] <- log2(1+racipe2Norm[,1:nCols]) # Log transform
    racipe2Norm[,1:nCols] <- sweep(racipe2Norm[,1:nCols], 2, tmpMeans, FUN = "-") # scale
    racipe2Norm[,1:nCols] <- sweep(racipe2Norm[,1:nCols], 2, tmpSds, FUN = "/") # scale
    newpca <- scale(racipe2Norm, pca$center, pca$scale) %*% pca$rotation 
    
    # newlabels <- assignClustersPC1(newpca) # PC1-based clustering replaced with KNN
    newlabels <- knn_classifier(racipe2Norm, racipeData, clust, k=25)
    

    # Compute & store # of transitions
    numTransitions <- length(which(newlabels == targetClust)) - length(which(clust == targetClust))
    rs_list[[tID]] <- list(NumTransitions = numTransitions,
                           NewClusters = newlabels,
                           Time = simTimes[tID])
    
    
  }
  
  
  if(save == F) {
    return(rs_list)
  } else {
    saveRDS(rs_list, file.path(expDir, paste0("transitionTimes_",paramSets[setID,"SetName"],".Rds")))
  }
  
  
  
  
  
}





theme_sticcc <- function() {
  font <- "Helvetica"   #assign font family up front
  
  theme_minimal() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      axis.ticks = element_blank(),          #strip axis ticks
      
      #since theme_minimal() already strips axis lines, 
      #we don't need to do that again
      
      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 28,                #set font size
        face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
        family = font,            #font family
        size = 16),               #font size
      
      plot.caption = element_text(           #caption
        family = font,            #font family
        size = 20,                 #font size
        hjust = 1),               #right align
      
      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 28),               #font size
      
      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 22),                #font size
      
      legend.title=element_text(size=18), 
      legend.text=element_text(size=16),
      
      
      
      #axis.text.x = element_text(            #margin for axis text
      #  margin=margin(5, b = 10))
      
      #since the legend often requires manual tweaking 
      #based on plot content, don't define it here
    )
}






params_x_transitions <- function(paramSets, # dataframe w/ columns: ModelNo, SetName
                                topoName, # string
                                collectionName, # string
                                initialClust, # int
                                targetClust, # int
                                wt_data, # matrix/dataframe of size numSamples x numFeatures
                                clust, # vector of length numSamples containing integers
                                tmpMeans,
                                tmpSds,
                                control = "non-target"
) {
  # Initialize a matrix to store the fold change values
  fc_matrix <- list()
  # Add results columns
  paramSets$InitialTargetStates <- NA
  paramSets$ConvertingModels <- NA
  paramSets$RemainingModels <- NA
  paramSets$RebelliousModels <- NA
  
  # Iterate over each parameter set
  for (i in 1:nrow(paramSets)) {
    # Extract the experiment name
    expName <- paramSets$SetName[i]
    
    # Extract noise level
    noise <- paramSets$Noise[i]
    
    # Read simulation data from file
    filepath <- file.path(getwd(), topoName, collectionName, paste0("racipe_", expName,".Rds"))
    filepath_final <- file.path(getwd(), topoName, collectionName, paste0("racipe_", expName,"_relaxed.Rds"))
    if(!file.exists(filepath) | !file.exists(filepath_final)) {
      next
    }
    condition_data <- readRDS(filepath)
    condition_data_final <- readRDS(filepath_final)
    condition_ss_final <- as.data.frame(t(assay(condition_data_final)))
    
    condition_data <- as.data.frame(sracipeParams(condition_data))
    condition_data_final <- as.data.frame(sracipeParams(condition_data_final))
    
    
    # Rescale data & assign clusters to perturbed data
    
    nCols <- ncol(condition_ss_final)
    condition_data_final[,1:nCols] <- log2(1+condition_data_final[,1:nCols]) # Log transform
    condition_data_final[,1:nCols] <- sweep(condition_data_final[,1:nCols], 2, tmpMeans, FUN = "-") # scale
    condition_data_final[,1:nCols] <- sweep(condition_data_final[,1:nCols], 2, tmpSds, FUN = "/") # scale
    

    # Re-assign clusters
    new_clust <- knn_classifier(condition_ss_final, wt_data, clust, k = 25)
    
    # Identify converting and remaining models
    paramSets$startingTargetPopulation[i] <- length(which(clust == targetClust))
    paramSets$startingInitPopulation[i] <- length(which(clust == initialClust))
    
    convertingModels <- condition_data[new_clust == targetClust & clust != targetClust, ]
    paramSets$ConvertingModels[i] <- nrow(convertingModels)
    
    if(control == "initial") {
      remainingModels <- condition_data[new_clust == initialClust & clust == initialClust, ]
    } else if (control == "non-target") {
      remainingModels <- condition_data[new_clust != targetClust & clust != targetClust, ]
    }
    paramSets$RemainingModels[i] <- nrow(remainingModels)
    
    rebelliousModels <- condition_data[new_clust != targetClust & clust == targetClust, ]
    paramSets$RebelliousModels[i] <- nrow(rebelliousModels)
    
    # Calculate log-fold change for each gene
    fc_values <- sapply(colnames(condition_data), function(gene) {
      avg_expr_converting <- mean(convertingModels[[gene]], na.rm = TRUE)
      avg_expr_remaining <- mean(remainingModels[[gene]], na.rm = TRUE)
      
      return(avg_expr_converting / avg_expr_remaining)
      
    })
    
    # Add the fold change values to the list
    fc_matrix[[i]] <- fc_values
  }
  
  # Convert the list to a data.frame and return
  fc_df <- do.call(rbind, fc_matrix)
  colnames(fc_df) <- colnames(condition_data)
  
  return(fc_df)
}












optimizeST_Clamp <- function(racipe, # NON NORMALIZED
                       pca, # same size as racipe
                       clust, # vector of length nrow(pca$x)
                       initialClust, # should correspond to clust
                       targetClust, # should correspond to clust
                       nSigGenes, # single value or vector of choices
                       clamp_df,
                       outDir,
                       expName = NA,
                       plot=F,
                       noise = NA, # single value or vector of choices
                       forceRerun = F, # whether to rerun signal simulations
                       forceRecompute = F, # whether to recompute final scores
                       checkpointSize=25, # how often to report progress to user
                       totalPerturbations = 500, # Only considered if randomParams is TRUE
                       relax=F, # whether to continue simulations for some time without a signal
                       simTime = 10,
                       simTimeRelax = 10,
                       onlyParams = F,
                       noise_tcorr = 10,
                       ...) {
  
  if(!dir.exists(outDir)) {
    print("Error: outdir could not be found")
    return(NULL)
  } 
  
  if(is.na(expName)) {
    expName <- paste0("Clamping_Analysis_",Sys.Date())
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
  if(all(is.na(nSigGenes))) {
    nSigGenes = c(1,2,3)
  }
  if(all(is.na(noise))) {
    noise = c(0.04, 0.1, 0.2, 0.5, 1)
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
  

  # Generate parameter sets
  paramSets <- optimizeSTParams_Clamp(targetClust,
                                nSigGenes,
                                geneinfo,
                                noise,
                                ...)   

  # Check for an existing paramSet file under this expName
  # If one is found, merge and remove duplicates
  old_pset_fname <- file.path(topoDir,expName,"result_summary.Rds")
  if(file.exists(old_pset_fname)) {
    paramSet_old <- readRDS(old_pset_fname)
    
    # By default, keep old results. If forcing parameter is true, discard these and recompute
    if(!forceRerun) {
      paramSets <- dplyr::bind_rows(paramSet_old, paramSets)
      paramSets <- paramSets[!duplicated(paramSets$SetName),]  
    } else {
      paramSets <- dplyr::bind_rows(paramSets, paramSet_old)
      paramSets <- paramSets[!duplicated(paramSets$SetName),] 
    }
    
  }
  
  paramSets <- paramSets[order(paramSets$TotalOutDegree, decreasing = T),]
  
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
    num <- paramSets[i,"NumGenes"]
    noisei <- paramSets[i,"Noise"]
    
    
    sig_clamp_genes <- getClampGenes(paramSets, i)
    sig_clamp_gene_ids <- unlist(lapply(sig_clamp_genes, function(x) which(rownames(racipe) == x)))
    sig_clamp_df <- getClampDF(clamp_df, sig_clamp_genes, targetClust)
    colnames(sig_clamp_df) <- as.numeric(sig_clamp_gene_ids)
    sig_clamp_df <- as.matrix(sig_clamp_df)
    

    # Set ICs to WT steady states
    sracipeIC(racipe2) <- assay(racipe)[,]
    
    # Simulate if not done already
    fname_racipe2 <- file.path(expDir, paste0("racipe_",paramSets[i,"SetName"],".Rds"))

    if(!file.exists(fname_racipe2) | forceRerun) {
        racipe2 <- sracipeSimulate(racipe2, genIC = F, genParams = F, simulationTime = simTime,
                                   initialNoise=noisei, nNoise=1, scaledNoise = T, nIC = 1,
                                   integrateStepSize = 0.2, simDet=F,
                                   stepper = "EM_Clamp", ouNoise_t = noise_tcorr, clampGenes=sig_clamp_gene_ids,
                                   clampValues=sig_clamp_df)
      
      saveRDS(racipe2, fname_racipe2)
    } else {
      racipe2 <- readRDS(fname_racipe2)
    }
    
    
    if(relax) {
      racipe_relax_fname <- file.path(expDir, paste0("racipe_",paramSets[i,"SetName"],"_relaxed.Rds"))
      if(!file.exists(racipe_relax_fname) | forceRerun) {
        
        
        racipe_relax <- racipe2 
        
        # Set ICs to perturbed steady states
        sracipeIC(racipe_relax) <- assay(racipe2, as.character(noisei))[,]
        
        # Simulate without signal
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
    
    # KNN classification of final states
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





optimizeSTParams_Clamp <- function(targetClust,
                             nSigGenes,
                             geneinfo,
                             noise,
                             ...) {
  paramSets = data.frame()
  for(n in nSigGenes) {
    paramSets_new <- as.data.frame(createParamSets_Clamp(targetClust=targetClust,
                                     nSigGenes = n,
                                     geneInfo = geneinfo,
                                     ...))
    
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
  relevant_columns <- grep("^Species", colnames(paramSets))
  sub_mat <- paramSets[, relevant_columns]
  
  # Apply the function to each row
  if(all(nSigGenes == 1)) {
    names_vector <- apply(as.data.frame(sub_mat), 1, create_name_Clamp, colnames=colnames(as.data.frame(sub_mat)))
  } else {
    names_vector <- apply(sub_mat, 1, create_name_Clamp, colnames=colnames(sub_mat))
  }
  
  paramSets$SetName <- paste0(names_vector,"_noise=",paramSets$Noise)
  
  
  # Leave empty columns to fill in
  #paramSets$OverlapLoss <- NA
  #paramSets$TargetLoss <- NA
  #paramSets$TotalLoss <- NA
  paramSets$Simulated <- FALSE
  paramSets$EndTargetStates <- NA
  
  # experimental
  #paramSets$SignalPower <- paramSets$TotalPageRank * paramSets$FoldChange * paramSets$Noise
  
  # reorder
  paramSets <- paramSets[order(paramSets$TotalPageRank, decreasing = T),]
  
  return(paramSets)
  
  
}


# Create the unified name function
create_name_Clamp <- function(row, colnames) {
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


# use hyperparameters and output from idGenes to create parameter sets
# return a nested list where each top-level element is a paramset, containing:
# sigGenes = list() of length nSigGenes
createParamSets_Clamp <- function(targetClust=2,
                            nSigGenes = 2,
                            geneInfo,
                            nPerturbations=NA,
                            order = T
) {
  
  
  k = nSigGenes
  nr = dim(geneInfo)[1]
  
  # Generate all combinations of k species from n total species
  species_combinations <- combn(nr, k)
  
  comb_df <- as.data.frame(t(species_combinations))
  comb_df$OutDegree <- 0
  comb_df$PageRank <- 0
  # sum the out-degree for each constituent gene in a signal
  for(c in 1:nSigGenes) {
    comb_df$OutDegree <- comb_df$OutDegree + geneInfo[comb_df[,paste0("V",c)],"outdegree"] 
    comb_df$PageRank <- comb_df$PageRank + geneInfo[comb_df[,paste0("V",c)],"PageRank"] 
  }
  
  if(!is.na(nPerturbations)) {
    ## Select top ranked combinations by total out-degree
    # subset ordered by out degree
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
    
  } else {
    out_deg_list <- comb_df[,"OutDegree"]
    pr_list <- comb_df[,"PageRank"]
  }
  
  # Initialize a list to hold the results
  all_combinations <- list()
  #clamp_dfs <- list()
  
  # Loop through each combination of species to perturb
  for (i in seq_len(ncol(species_combinations))) {
    current_combination <- species_combinations[, i]
    
    # Attach the species combination to the signaling combinations
    result <- matrix(current_combination, nrow=1, ncol=k, byrow=TRUE)
    colnames(result) <- c(paste("Species", seq_len(k)))
    
    
    # Store the result in the list
    all_combinations[[i]] <- result
  }
  
  # merge all combinations into one data frame
  final_result <- as.data.frame(do.call(rbind, all_combinations))
  
  
  # Rename signal genes
  for(r in 1:nSigGenes) {
    # Use match to find the index of each value in df1$int_column in df2$int_value
    index <- match(final_result[,paste0("Species ",r)], geneInfo$GeneID)
    
    # Replace df1$int_column with corresponding string values
    final_result[,paste0("Species ",r)] <- geneInfo$gene[index]
  }
  

  combinations_per_geneCombo = 1
  final_result$TotalOutDegree <- rep(out_deg_list, each=combinations_per_geneCombo)
  final_result$TotalPageRank <- rep(pr_list, each=combinations_per_geneCombo)
  
  return(final_result)
  
  
}



getClampGenes <- function(paramSets, 
                          idx) {
  out_genes <- c()
  for(i in 1:paramSets[idx,"NumGenes"]) {
    out_genes[i] <- paramSets[idx,paste0("Species ",i)]
  }
  return(out_genes)
  
}


getClampDF <- function(clamp_df, 
                       clamp_genes, 
                       targetClust) {
  numModels <- length(unique(clamp_df$Model))
  
  clamp_df_out <- as.data.frame(matrix(data=NA, nrow=numModels, ncol = length(clamp_genes)))
  for(gene in clamp_genes) {
    for(model in 1:numModels) {
      clamp_df_out[model, which(clamp_genes == gene)] <- clamp_df[which(clamp_df$ModelIndex == model &
                                                    clamp_df$Cluster == targetClust &
                                                    clamp_df$Gene == gene),"Expression"]
    }
  }
  
  return(clamp_df_out)
}



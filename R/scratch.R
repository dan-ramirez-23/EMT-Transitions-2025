## Splitting this function out into parts so I can process different signal sets on different cores. 

optimizeST_parallel <- function(racipe,
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
                       nCores = 4,
                       saveTrajectory = F,
                       printInterval = 10,
                       noise_tcorr = 5,
                       clamp=F,
                       clamp_df=NA,
                       tmpMeans, # means for genes in baseline experiment
                       tmpSds, # st devs for genes in baseline experiment
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
  #tmpMeans <- rowMeans(simExp)
  #tmpSds <- apply(simExp,1,sd)
  
  racipeNorm <- sracipeNormalize(racipe)
  racipeData <- as.data.frame(t(assay(racipeNorm)))
  
  
  # get gene info, out-degree
  geneinfo <- idGenes(clust, initialClust, targetClust, racipeData)
  geneinfo <- topoOutDegrees(geneinfo, 
                             topo = sracipeCircuit(racipe))
  
  # Get WT parameters
  params <- sracipeParams(racipe)  
  params_fname <- file.path(expDir, "params_WT.Rds")
  if(!file.exists(params_fname)) {
    saveRDS(params, params_fname)  
  }
  
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
  } else if(clamp){
    paramSets <- optimizeSTParams_Clamp(targetClust,
                                        nSigGenes,
                                        geneinfo,
                                        noise,
                                        ...)   
    
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
  
  #paramSets <- paramSets[order(paramSets$SignalPower, decreasing = T),]
  paramSets <- paramSets[order(paramSets$TotalOutDegree, decreasing = T),]

  
  if(onlyParams) {
    return(paramSets)
  }
  
  
  
  # simulate parameter sets
  currentMinLoss <- min(c(paramSets$TotalLoss, Inf), na.rm=T)
  newMinLoss <- F
  results = list()
  print(paste0("Beginning parameter set list of size ",nrow(paramSets)))
  
  
  ### The next ~150 lines need to be converted to a fully parameterized function. Then, process outputs & plot after
  ## Outline (only need to parallelize the starred items):
  
  ### PRE LOOP
  # Check for files, skip paramSet if existing
  # get signal traits, adjust parameters
  # save following in tmp files/lists:
  # racipe2 (with adjusted signal) --> racipeTmp, 
  # fname_racipe2,
  # fname_racipe_relax,
  # forceRerun = F,
  # simTime,
  # simTimeRelax,
  # noisei,
  # nCores = 4,
  # relax = T
  
  ### LOOP
  # simulate phase 1*
  # simulate phase 2*
  
  ### POST LOOP
  # compute metrics
  # plot
  skipList <- rep(FALSE, nrow(paramSets))
  racipeList <- character(length = nrow(paramSets))
  racipeFnameList <- character(length = nrow(paramSets))
  racipeRelaxFnameList <- character(length = nrow(paramSets))
  sig_clamp_gene_ids_list <- character(length = nrow(paramSets))
  sig_clamp_df_list <- character(length = nrow(paramSets))
  noiseList <- character(length = nrow(paramSets))
  
  for(i in 1:nrow(paramSets)) {
    
    
    fname_racipe2 <- file.path(expDir, paste0("racipe_",paramSets[i,"SetName"],".Rds"))
    racipe_relax_fname <- file.path(expDir, paste0("racipe_",paramSets[i,"SetName"],"_relaxed.Rds"))
    fname_racipe_tmp <- file.path(expDir, paste0("racipe_",paramSets[i,"SetName"],"_TMP.Rds"))
    fname_clamp_gene_ids <- file.path(expDir, paste0("clamp_gene_ids_",paramSets[i,"SetName"],".Rds"))
    fname_clamp_df <- file.path(expDir, paste0("clamp_df_",paramSets[i,"SetName"],".Rds"))
    
    racipeFnameList[[i]] <- fname_racipe2
    racipeRelaxFnameList[[i]] <- racipe_relax_fname
    racipeList[[i]] <- fname_racipe_tmp
    sig_clamp_gene_ids_list[[i]] <- fname_clamp_gene_ids
    sig_clamp_df_list[[i]] <- fname_clamp_df
    
    if(paramSets[i,"Simulated"] & !forceRecompute){
      print(paste0("Skipping set ", i))
      skipList[[i]] <- TRUE
      next
    }
    
    
    racipe2 <- racipe 
    sracipeParams(racipe2) <- params[,]
    num <- paramSets[i,"NumGenes"]
    noisei <- paramSets[i,"Noise"]
    noiseList[[i]] <- noisei
    
    
    # add signal
    if(clamp) {
      sig_clamp_genes <- getClampGenes(paramSets, i)
      sig_clamp_gene_ids <- unlist(lapply(sig_clamp_genes, function(x) which(rownames(racipe) == x)))
      sig_clamp_df <- getClampDF(clamp_df, sig_clamp_genes, targetClust)
      colnames(sig_clamp_df) <- as.numeric(sig_clamp_gene_ids)
      sig_clamp_df <- as.matrix(sig_clamp_df)
      
      saveRDS(sig_clamp_gene_ids, fname_clamp_gene_ids)
      saveRDS(sig_clamp_df, fname_clamp_df)
      
    } else {
      sig_clamp_gene_ids <- NA
      sig_clamp_df <- NA
      for(sigID in 1:num) {
        sracipeParams(racipe2)[,paste0(paramType,"_",paramSets[setID,paste0("Species ",sigID)])] <- 
          params[,paste0(paramType,"_",paramSets[setID,paste0("Species ",sigID)])] * 
          paramSets[setID,paste0("Signal ",sigID)]
      }
    }
    
    # Set ICs to WT steady states
    sracipeIC(racipe2) <- assay(racipe)[,]
    
    # Save tmp files for parallel processing
    
    if(!file.exists(fname_racipe_tmp)) {
      saveRDS(racipe2, fname_racipe_tmp)  
    }
    
    
  }
  print("Phase 1 done - preparing to compute transitions in parallel")
    
    
  ##### PREP PARALLEL LOOP #####
  requireNamespace("doFuture")
  registerDoFuture()
  plan(multisession, workers = nCores)
  maxIdx <- length(racipeList)
  
  

  ##### SIMULATE PHASES 1 & 2 #####
  x <- foreach(idx = as.list(c(1:length(racipeList))),
               racipeTmpFname = racipeList,
               racipe_fname = racipeFnameList,
               racipe_relax_fname = racipeRelaxFnameList,
               noisei = noiseList,
               fname_clamp_gene_ids = sig_clamp_gene_ids_list,
               fname_clamp_df = sig_clamp_df_list,
               .export = c("forceRerun","simTime","simTimeRelax","noisei","nCores","relax", 
                           "maxIdx", "params_fname", "skipList", "saveTrajectory", "printInterval",
                           "clamp","noise_tcorr")) %dorng% {
                             
                 library(sRACIPE, lib.loc = "/Users/danramirez/localR/4.2.2-arm")
                 
                 if(skipList[[idx]]) {
                   print(paste0("Skipping task no. ",idx," of ",maxIdx,": ", racipeTmpFname))
                 } else {
                   print(paste0("Beginning task no. ",idx," of ",maxIdx,": ", racipeTmpFname))
                   
                   sig_clamp_gene_ids <- readRDS(fname_clamp_gene_ids)
                   sig_clamp_df <- readRDS(fname_clamp_df)
                   
                   simulateST(racipeTmp = racipeTmpFname, 
                              paramsTmp = params_fname,
                              fname_racipe2 = racipe_fname,
                              racipe_relax_fname = racipe_relax_fname,
                              forceRerun = forceRerun,
                              simTime = simTime,
                              simTimeRelax = simTimeRelax,
                              noisei = as.numeric(noisei),
                              #nCores = nCores,
                              relax = relax,
                              noise_tcorr = noise_tcorr,
                              clamp=clamp,
                              sig_clamp_gene_ids = sig_clamp_gene_ids,
                              sig_clamp_df = sig_clamp_df
                              )
                 }
                 
                 

                 
                 #p()
                 
                 
               }
  #})
  
  print("Phase 2 done - preparing post-processing results")
  
    
    
  
    
  ##### PROCESS SIMULATION RESULTS #####
  for(i in 1:nrow(paramSets)) {
    racipe_tmp_fname <- racipeList[[i]]
    racipe_relax_fname <- racipeRelaxFnameList[[i]]
    racipe2 <- readRDS(racipe_relax_fname)
    racipe2Norm <- racipe2Norm <- as.data.frame(t(assay(racipe2)))
    
    # Remove tmp files
    file.remove(racipe_tmp_fname)
    
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
    numTransitions <- length(which(newlabels == targetClust)) - length(which(clust == targetClust))
    
    
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



simulateST <- function(racipeTmp, 
                       paramsTmp,
                       fname_racipe2,
                       racipe_relax_fname,
                       forceRerun = F,
                       simTime,
                       simTimeRelax,
                       noisei,
                       #nCores = 4,
                       relax = T,
                       anneal = F,
                       noise_tcorr = 5,
                       clamp = F,
                       sig_clamp_gene_ids = NA,
                       sig_clamp_df = NA) {
  
  
  racipe2 <- readRDS(racipeTmp)
  params <- readRDS(paramsTmp)
  
  if(!file.exists(fname_racipe2) | forceRerun) {
    if(anneal) {
      racipe2 <- sracipeSimulate(racipe2, genIC = F, genParams = F, simulationTime = simTime, 
                                 initialNoise=noisei, nNoise=8, scaledNoise = T, anneal = T, 
                                 integrateStepSize = 0.2)  
    } else if(clamp) {
      
      racipe2 <- sracipeSimulate(racipe2, genIC = F, genParams = F, simulationTime = simTime,
                                 initialNoise=noisei, nNoise=1, scaledNoise = T, nIC = 1,
                                 integrateStepSize = 0.2, simDet=F,
                                 stepper = "EM_Clamp", ouNoise_t = noise_tcorr, clampGenes=sig_clamp_gene_ids,
                                 clampValues=sig_clamp_df)
      
    } else {
      racipe2 <- sracipeSimulate(racipe2, genIC = F, genParams = F, simulationTime = simTime,
                                 initialNoise=noisei, nNoise=1, scaledNoise = T, 
                                 integrateStepSize = 0.2, simDet=F, stepper = "EM_OU", ouNoise_t=noise_tcorr)
    }
    
    
    tryCatch({
      saveRDS(racipe2, fname_racipe2)
      TRUE  # Return TRUE on successful save
    }, error = function(e) {
      message(sprintf("Error saving file %s: %s", filename, e$message))
      FALSE  # Return FALSE on error
    })
    
  } else {
    racipe2 <- readRDS(fname_racipe2)
    
  }
  
  
  if(relax) {
    if(!file.exists(racipe_relax_fname) | forceRerun) {
    
      # Remove signal
      racipe_relax <- racipe2 
      sracipeParams(racipe_relax) <- params[,]
      # Set ICs to perturbed steady states
      sracipeIC(racipe_relax) <- assay(racipe2, as.character(noisei))[,]
      
      racipe_relax <- sracipeSimulate(racipe_relax, genIC = F, genParams = F, simulationTime = simTimeRelax,
                                      integrateStepSize = 0.2, nNoise = 0, iNoise = 0, simDet = T)
                                      #, timeSeries = saveTrajectory, 
                                      #printStart = printInterval, printInterval = printInterval) #nCores = nCores, 
      
      
      tryCatch({
        saveRDS(racipe_relax, racipe_relax_fname)
        TRUE  # Return TRUE on successful save
      }, error = function(e) {
        message(sprintf("Error saving file %s: %s", filename, e$message))
        FALSE  # Return FALSE on error
      })
      
    } 
    
  } else {
    #racipe2Norm <- as.data.frame(t(assay(racipe2,as.character(noisei))))
  }
}






st_ensemble_analysis <- function(racipe,
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
                                nCores = 4,
                                clamp=F,
                                clamp_df=NA,
                                tmpMeans, # means for genes in baseline experiment
                                tmpSds, # st devs for genes in baseline experiment
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
  #tmpMeans <- rowMeans(simExp)
  #tmpSds <- apply(simExp,1,sd)
  
  racipeNorm <- sracipeNormalize(racipe)
  racipeData <- as.data.frame(t(assay(racipeNorm)))
  
  
  # get gene info, out-degree
  geneinfo <- idGenes(clust, initialClust, targetClust, racipeData)
  geneinfo <- topoOutDegrees(geneinfo, 
                             topo = sracipeCircuit(racipe))
  
  # Get WT parameters
  params <- sracipeParams(racipe)  
  params_fname <- file.path(expDir, "params_WT.Rds")
  if(!file.exists(params_fname)) {
    saveRDS(params, params_fname)  
  }
  
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
  } else if(clamp){
    paramSets <- optimizeSTParams_Clamp(targetClust,
                                        nSigGenes,
                                        geneinfo,
                                        noise,
                                        ...)   
    
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
  
  #paramSets <- paramSets[order(paramSets$SignalPower, decreasing = T),]
  #if(all(noise) == 0) {
  paramSets <- paramSets[order(paramSets$TotalOutDegree, decreasing = T),]
  #}
  
  if(onlyParams) {
    return(paramSets)
  }
  
  
  
  # simulate parameter sets
  currentMinLoss <- min(c(paramSets$TotalLoss, Inf), na.rm=T)
  newMinLoss <- F
  results = list()
  print(paste0("Beginning parameter set list of size ",nrow(paramSets)))
  
  
  ### The next ~150 lines need to be converted to a fully parameterized function. Then, process outputs & plot after
  ## Outline (only need to parallelize the starred items):
  
  ### PRE LOOP
  # Check for files, skip paramSet if existing
  # get signal traits, adjust parameters
  # save following in tmp files/lists:
  # racipe2 (with adjusted signal) --> racipeTmp, 
  # fname_racipe2,
  # fname_racipe_relax,
  # forceRerun = F,
  # simTime,
  # simTimeRelax,
  # noisei,
  # nCores = 4,
  # relax = T
  
  ### LOOP
  # simulate phase 1*
  # simulate phase 2*
  
  ### POST LOOP
  # compute metrics
  # plot
  skipList <- rep(FALSE, nrow(paramSets))
  racipeList <- character(length = nrow(paramSets))
  racipeFnameList <- character(length = nrow(paramSets))
  racipeRelaxFnameList <- character(length = nrow(paramSets))
  sig_clamp_gene_ids_list <- character(length = nrow(paramSets))
  sig_clamp_df_list <- character(length = nrow(paramSets))
  
  
  for(i in 1:nrow(paramSets)) {
    
    
    fname_racipe2 <- file.path(expDir, paste0("racipe_",paramSets[i,"SetName"],".Rds"))
    racipe_relax_fname <- file.path(expDir, paste0("racipe_",paramSets[i,"SetName"],"_relaxed.Rds"))
    fname_racipe_tmp <- file.path(expDir, paste0("racipe_",paramSets[i,"SetName"],"_TMP.Rds"))
    fname_clamp_gene_ids <- file.path(expDir, paste0("clamp_gene_ids_",paramSets[i,"SetName"],".Rds"))
    fname_clamp_df <- file.path(expDir, paste0("clamp_df_",paramSets[i,"SetName"],".Rds"))
    
    racipeFnameList[[i]] <- fname_racipe2
    racipeRelaxFnameList[[i]] <- racipe_relax_fname
    racipeList[[i]] <- fname_racipe_tmp
    sig_clamp_gene_ids_list[[i]] <- fname_clamp_gene_ids
    sig_clamp_df_list[[i]] <- fname_clamp_df

    
    if(paramSets[i,"Simulated"] & !forceRecompute){
      print(paste0("Skipping set ", i))
      skipList[[i]] <- TRUE
      next
    }
    
    
    racipe2 <- racipe 
    sracipeParams(racipe2) <- params[,]
    num <- paramSets[i,"NumGenes"]
    noisei <- paramSets[i,"Noise"]
    
    # add signal
    if(clamp) {
      sig_clamp_genes <- getClampGenes(paramSets, i)
      sig_clamp_gene_ids <- unlist(lapply(sig_clamp_genes, function(x) which(rownames(racipe) == x)))
      sig_clamp_df <- getClampDF(clamp_df, sig_clamp_genes, targetClust)
      colnames(sig_clamp_df) <- as.numeric(sig_clamp_gene_ids)
      sig_clamp_df <- as.matrix(sig_clamp_df)
      
      if(!file.exists(fname_clamp_gene_ids)) {
        saveRDS(sig_clamp_gene_ids, fname_clamp_gene_ids)
        saveRDS(sig_clamp_df, fname_clamp_df)
      }
      
    } else {
      sig_clamp_gene_ids <- NA
      sig_clamp_df <- NA
      for(sigID in 1:num) {
        sracipeParams(racipe2)[,paste0(paramType,"_",paramSets[setID,paste0("Species ",sigID)])] <- 
          params[,paste0(paramType,"_",paramSets[setID,paste0("Species ",sigID)])] * 
          paramSets[setID,paste0("Signal ",sigID)]
      }
    }
    
    # Set ICs to WT steady states
    sracipeIC(racipe2) <- assay(racipe)[,]
    
    # Save tmp files for parallel processing
    
    if(!file.exists(fname_racipe_tmp)) {
      saveRDS(racipe2, fname_racipe_tmp)  
    }
    
    
  }
  print("Phase 1 done - preparing to compute transitions in parallel")
  
  
  
  ##### PROCESS SIMULATION RESULTS #####
  for(i in 1:nrow(paramSets)) {
    
    racipe_relax_fname <- racipeRelaxFnameList[[i]]
    if(!file.exists(racipe_relax_fname)) {
      print(paste0("Skipping set ", i))
      next
    }
    racipe2 <- readRDS(racipe_relax_fname)
    racipe2Norm <- racipe2Norm <- as.data.frame(t(assay(racipe2)))
    
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
    numTransitions <- length(which(newlabels == targetClust)) - length(which(clust == targetClust))
    
    
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





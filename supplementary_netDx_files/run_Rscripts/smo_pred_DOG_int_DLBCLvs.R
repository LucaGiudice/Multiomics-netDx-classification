rm(list=ls())
rootData="/home/wr/%s"
#Load libraries ----
require(netDx)
require(netDx.examples)
library(matrixStats)
require(data.table)
#Function that takes a matrix: genes (measures) x patients
#Returns a simmetric similarity matrix of the weighted jaccard between all the patients
#It has been made to work with smo_pred2.R

source("/home/wr/sim_PEA.R")
source("/home/wr/sim_WJ2_v3.R")

#Build the path of the input and output files ----
setwd(sprintf(rootData,"DOG/env/data_multi/"))
rda_files = list.files(pattern="*.rda")
rda_files=rda_files[c(10,13,1)] #dog_gex_DLBCLvsFL.rda,dog_meth_DLBCLvsFL.rda,dog_cn_DLBCLvsFL.rda
pred_dir=sprintf(rootData,"DOG/predictions/CN_int_GEX_int_METH_%s")

datList <- list()
for(f_count in 1:length(rda_files)){
 load(rda_files[f_count])
 if(f_count==1){name_l="gex"}
 if(f_count==2){name_l="meth"}
 if(f_count==3){name_l="cn"}
 datList[[name_l]] <- geno;
}

  outDir=sprintf(pred_dir,"DLBCLvsFL")
  print("I am running with this outDir:")
  print(outDir)
  #Load data ----
  
  pheno[]=lapply(pheno,as.character)
  #PSN function passed to netDx ----
  makeNets <- function(dataList,groupList,netDir,numCores,...) {
    netList <- c(); netList2 <- c(); netList3 <- c()

    # create genetic nets
    if (!is.null(groupList[["gex"]])) {netList <- makePSN_NamedMatrix(dataList[["gex"]],
                                                                          rownames(dataList[["gex"]]),
                                                                          groupList[["gex"]],
                                                                          netDir,writeProfiles=FALSE,
                                                                          sparsify=TRUE,verbose=TRUE,
                                                                          append=TRUE,numC=numCores,sparsify_edgeMax=.Machine$integer.max,sparsify_maxInt=.Machine$integer.max,cutoff=.Machine$double.eps,simMetric="custom",customFunc=sim_WJ)
    }

    netList <- unlist(netList)
    cat(sprintf("\t%i gex-pathway nets\n", length(netList)))



    # create genetic nets
    if (!is.null(groupList[["cn"]])) {netList2 <- makePSN_NamedMatrix(dataList[["cn"]],
                                                                          rownames(dataList[["cn"]]),
                                                                          groupList[["cn"]],
                                                                          netDir,writeProfiles=FALSE,
                                                                          sparsify=TRUE,verbose=TRUE,
                                                                          append=TRUE,numC=numCores,sparsify_edgeMax=.Machine$integer.max,sparsify_maxInt=.Machine$integer.max,cutoff=.Machine$double.eps,
simMetric="custom",customFunc=sim_PEA)
    }
    netList2 <- unlist(netList2)
    cat(sprintf("\t%i cn-pathway nets\n", length(netList2)))



   # create genetic nets
    if (!is.null(groupList[["meth"]])) {netList3 <- makePSN_NamedMatrix(dataList[["meth"]],
                                                                          rownames(dataList[["meth"]]),
                                                                          groupList[["meth"]],
                                                                          netDir,writeProfiles=FALSE,
                                                                          sparsify=TRUE,verbose=TRUE,
                                                                          append=TRUE,numC=numCores,sparsify_edgeMax=.Machine$integer.max,sparsify_maxInt=.Machine$integer.max,cutoff=.Machine$double.eps,simMetric="custom",customFunc=sim_PEA)
    }

    netList3 <- unlist(netList3)
    cat(sprintf("\t%i meth-pathway nets\n", length(netList3)))

    netList <- c(netList,netList2,netList3) 
    cat(sprintf("Total of %i nets\n", length(netList)))

    return(netList)
  }

  #Setup to run predictor -----
  groupList	<- list()

  #pathFile	<- sprintf("%s/extdata/Human_160124_AllPathways.gmt",path.package("netDx.examples"))
  pathFile="/home/wr/Human_December_2018.gmt"
  pathwayList <- readPathways(pathFile,MIN_SIZE=5L,MAX_SIZE=300L)

  #For each pathway
  #Keep the index of the genes in geno that belong to the gene set
  gene_names=rownames(datList[["gex"]])
  netList <- foreach(k=1:length(pathwayList)) %do% {
    idx <- which(gene_names %in% pathwayList[[k]])
  }

  #Remove from the set of pathways and from netList those that do not have any gene in geno matrix
  n_genes_in_pathway=sapply(netList,FUN=length)
  pathwayList=pathwayList[n_genes_in_pathway>1]
  groupList[["gex"]] <- pathwayList

  pathwayList <- readPathways(pathFile,MIN_SIZE=5L,MAX_SIZE=300L)
  #For each pathway
  #Keep the index of the genes in geno that belong to the gene set
  gene_names=rownames(datList[["cn"]])
  netList <- foreach(k=1:length(pathwayList)) %do% {
    idx <- which(gene_names %in% pathwayList[[k]])
  }

  #Remove from the set of pathways and from netList those that do not have any gene in geno matrix
  n_genes_in_pathway=sapply(netList,FUN=length)
  pathwayList=pathwayList[n_genes_in_pathway>1]
  groupList[["cn"]] <- pathwayList

  pathwayList <- readPathways(pathFile,MIN_SIZE=5L,MAX_SIZE=300L)
  #For each pathway
  #Keep the index of the genes in geno that belong to the gene set
  gene_names=rownames(datList[["meth"]])
  netList <- foreach(k=1:length(pathwayList)) %do% {
    idx <- which(gene_names %in% pathwayList[[k]])
  }

  #Remove from the set of pathways and from netList those that do not have any gene in geno matrix
  n_genes_in_pathway=sapply(netList,FUN=length)
  pathwayList=pathwayList[n_genes_in_pathway>1]
  groupList[["meth"]] <- pathwayList

  #Run the predictor ----
  if (file.exists(outDir)){
    cat("the directory of output already exists")
    x <- readline("Do you want to overwrite the directory: y|n?")
    if(x=="y"){unlink("outDir", recursive=TRUE)}
    else {
      cat("Stop the execution")
      stop(call. = FALSE)
    }
  }

  runPredictor_nestedCV(pheno,
                        dataList=datList,groupList=groupList,
                        makeNetFunc=makeNets, ### custom network creation function
                        outDir=outDir, ## absolute path
                        numCores=18L,nFoldCV=10L, CVcutoff=8L,numSplits=30L, CVmemory=16L, trainProp=0.99
  )

  Sys.sleep(120)
  print("I finished a prediction")



rm(list=ls())
rootData="/home/wr/%s"
#Load libraries ----
require(netDx)
require(netDx.examples)
library(matrixStats)
require(data.table)
#Function that takes a matrix: genes (measures) x patients
#Returns a simmetric similarity matrix of the weighted jaccard between all the patients
#It has been made to work with smo_pred2.R

source("/home/wr/sim_PEA.R")
source("/home/wr/sim_WJ2_v3.R")

#Build the path of the input and output files ----
setwd(sprintf(rootData,"DOG/env/data_multi/"))
rda_files = list.files(pattern="*.rda")
rda_files=rda_files[c(11,14,2)] ##dog_gex_DLBCLvsMZL.rda,dog_meth_DLBCLvsMZL.rda,dog_cn_DLBCLvsMZL.rda
pred_dir=sprintf(rootData,"DOG/predictions/CN_int_GEX_int_METH_%s")

datList <- list()
for(f_count in 1:length(rda_files)){
 load(rda_files[f_count])
 if(f_count==1){name_l="gex"}
 if(f_count==2){name_l="meth"}
 if(f_count==3){name_l="cn"}
 datList[[name_l]] <- geno;
}

  outDir=sprintf(pred_dir,"DLBCLvsMZL")
  print("I am running with this outDir:")
  print(outDir)
  #Load data ----
  
  pheno[]=lapply(pheno,as.character)
  #PSN function passed to netDx ----
  makeNets <- function(dataList,groupList,netDir,numCores,...) {
    netList <- c(); netList2 <- c(); netList3 <- c()

    # create genetic nets
    if (!is.null(groupList[["gex"]])) {netList <- makePSN_NamedMatrix(dataList[["gex"]],
                                                                          rownames(dataList[["gex"]]),
                                                                          groupList[["gex"]],
                                                                          netDir,writeProfiles=FALSE,
                                                                          sparsify=TRUE,verbose=TRUE,
                                                                          append=TRUE,numC=numCores,sparsify_edgeMax=.Machine$integer.max,sparsify_maxInt=.Machine$integer.max,cutoff=.Machine$double.eps,simMetric="custom",customFunc=sim_WJ)
    }

    netList <- unlist(netList)
    cat(sprintf("\t%i gex-pathway nets\n", length(netList)))



    # create genetic nets
    if (!is.null(groupList[["cn"]])) {netList2 <- makePSN_NamedMatrix(dataList[["cn"]],
                                                                          rownames(dataList[["cn"]]),
                                                                          groupList[["cn"]],
                                                                          netDir,writeProfiles=FALSE,
                                                                          sparsify=TRUE,verbose=TRUE,
                                                                          append=TRUE,numC=numCores,sparsify_edgeMax=.Machine$integer.max,sparsify_maxInt=.Machine$integer.max,cutoff=.Machine$double.eps,
simMetric="custom",customFunc=sim_PEA)
    }
    netList2 <- unlist(netList2)
    cat(sprintf("\t%i cn-pathway nets\n", length(netList2)))



   # create genetic nets
    if (!is.null(groupList[["meth"]])) {netList3 <- makePSN_NamedMatrix(dataList[["meth"]],
                                                                          rownames(dataList[["meth"]]),
                                                                          groupList[["meth"]],
                                                                          netDir,writeProfiles=FALSE,
                                                                          sparsify=TRUE,verbose=TRUE,
                                                                          append=TRUE,numC=numCores,sparsify_edgeMax=.Machine$integer.max,sparsify_maxInt=.Machine$integer.max,cutoff=.Machine$double.eps,simMetric="custom",customFunc=sim_PEA)
    }

    netList3 <- unlist(netList3)
    cat(sprintf("\t%i meth-pathway nets\n", length(netList3)))

    netList <- c(netList,netList2,netList3) 
    cat(sprintf("Total of %i nets\n", length(netList)))

    return(netList)
  }

  #Setup to run predictor -----
  groupList	<- list()

  #pathFile	<- sprintf("%s/extdata/Human_160124_AllPathways.gmt",path.package("netDx.examples"))
  pathFile="/home/wr/Human_December_2018.gmt"
  pathwayList <- readPathways(pathFile,MIN_SIZE=5L,MAX_SIZE=300L)

  #For each pathway
  #Keep the index of the genes in geno that belong to the gene set
  gene_names=rownames(datList[["gex"]])
  netList <- foreach(k=1:length(pathwayList)) %do% {
    idx <- which(gene_names %in% pathwayList[[k]])
  }

  #Remove from the set of pathways and from netList those that do not have any gene in geno matrix
  n_genes_in_pathway=sapply(netList,FUN=length)
  pathwayList=pathwayList[n_genes_in_pathway>1]
  groupList[["gex"]] <- pathwayList

  pathwayList <- readPathways(pathFile,MIN_SIZE=5L,MAX_SIZE=300L)
  #For each pathway
  #Keep the index of the genes in geno that belong to the gene set
  gene_names=rownames(datList[["cn"]])
  netList <- foreach(k=1:length(pathwayList)) %do% {
    idx <- which(gene_names %in% pathwayList[[k]])
  }

  #Remove from the set of pathways and from netList those that do not have any gene in geno matrix
  n_genes_in_pathway=sapply(netList,FUN=length)
  pathwayList=pathwayList[n_genes_in_pathway>1]
  groupList[["cn"]] <- pathwayList

  pathwayList <- readPathways(pathFile,MIN_SIZE=5L,MAX_SIZE=300L)
  #For each pathway
  #Keep the index of the genes in geno that belong to the gene set
  gene_names=rownames(datList[["meth"]])
  netList <- foreach(k=1:length(pathwayList)) %do% {
    idx <- which(gene_names %in% pathwayList[[k]])
  }

  #Remove from the set of pathways and from netList those that do not have any gene in geno matrix
  n_genes_in_pathway=sapply(netList,FUN=length)
  pathwayList=pathwayList[n_genes_in_pathway>1]
  groupList[["meth"]] <- pathwayList

  #Run the predictor ----
  if (file.exists(outDir)){
    cat("the directory of output already exists")
    x <- readline("Do you want to overwrite the directory: y|n?")
    if(x=="y"){unlink("outDir", recursive=TRUE)}
    else {
      cat("Stop the execution")
      stop(call. = FALSE)
    }
  }

  runPredictor_nestedCV(pheno,
                        dataList=datList,groupList=groupList,
                        makeNetFunc=makeNets, ### custom network creation function
                        outDir=outDir, ## absolute path
                        numCores=18L,nFoldCV=10L, CVcutoff=8L,numSplits=30L, CVmemory=16L, trainProp=0.99
  )

  Sys.sleep(120)
  print("I finished a prediction")




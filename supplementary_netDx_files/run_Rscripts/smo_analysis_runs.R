rm(list=ls())

#Load libraries ----
suppressMessages(require(netDx))
suppressMessages(require(netDx.examples))
require("Biobase")

#Root path of the dataset ----
rootData="/home/wr/DOG/%s"
#Build the pathes of the input and output files ----
#set: path of the directory in which the result will be saved
outDir=sprintf(rootData,"analysis/")
#set: path of the directory containing all the results of the prediction
PredDir=sprintf(rootData,"predictions/")
n_runs=30
date=Sys.Date()
file_rda_2save=paste("30-int-analysis-",date,".rda",sep="")

#Check if outDir already exists ----
#false: create directory
if (!file.exists(outDir)) dir.create(outDir)
setwd(outDir)

#Load the name of directories (predictions) of the each case ----
all_preds <- list.dirs(PredDir, recursive = FALSE)
all_preds=all_preds[c(1,2,3)]
print(head(basename(all_preds)))
#WITHOUT BIN
#all_preds=all_preds[c(-1,-4)]


list_featScore=list()  
list_featSelNet=list()
#Prepare matrix for AUROC, AUPR, accuracy
res_auroc=matrix(0, n_runs, length(all_preds))
res_aupr=matrix(0, n_runs, length(all_preds))
res_acc=matrix(0, n_runs, length(all_preds))
colnames(res_auroc)=basename(all_preds)
colnames(res_aupr)=basename(all_preds)
colnames(res_acc)=basename(all_preds)


#For each directory prediction ----
for(k in 1:length(all_preds)) {
  #Set inDir
  inDir=all_preds[k]
  #set: path of the file "inputNets.txt" that belongs to the result
  netInfoFile=sprintf("%s/inputNets.txt",inDir)
  #load the name of directories (results) of the each prediction
  all_rngs <- list.dirs(inDir, recursive = FALSE)
  print(head(basename(all_rngs)))
  #Find the classes of the prediction
  dir(all_rngs[1])
  predClasses=setdiff(dir(all_rngs[1]),c("ids.txt","predictionResults.txt"))
  #Extract the performances
  predFiles <- unlist(lapply(all_rngs, function(x) paste(x, "predictionResults.txt", sep = "/")))
  #Check if there is a result for all the runs
  if(length(predFiles)==n_runs){
    print("I AM DOING:")
    print(inDir)
    
    png_file=substr(inDir,27,nchar(inDir))
    png_file=paste(png_file,".png",sep="")
    png(png_file, width=1000, height=1000)
    predPerf <- plotPerf(inDir, predClasses=predClasses)
    dev.off()
    #Create the matrix s.t. each column contains the measures of the runs
    res=matrix(0, n_runs, 3)
    colnames(res)=c("auroc","aupr","accuracy")
    rownames(res)=seq(1:n_runs)
    
    res[,1]=round(subListExtract(predPerf, name="auroc", simplify = TRUE, keep.names = TRUE),digits = 3)
    res[,2]=round(subListExtract(predPerf, name="aupr", simplify = TRUE, keep.names = TRUE),digits = 3)
    res[,3]=round(subListExtract(predPerf, name="accuracy", simplify = TRUE, keep.names = TRUE),digits = 3)
    
    res_auroc[,k]=res[,1]
    res_aupr[,k]=res[,2]
    res_acc[,k]=res[,3]
    #Feature-level scores are one of the outputs of feature-selection. ----
    #Note that in netDx, feature selection occurs once per patient class so each feature gets a predictive score for each class.
    featScores <- getFeatureScores(inDir,predClasses=predClasses)
    list_featScore[[inDir]]=featScores
    
    #Let us define feature-selected networks as those that score 7 out of 10 in at least 70% of the train/test splits
    featSelNet <- lapply(featScores, function(x) {
      callFeatSel(x, fsCutoff=9, fsPctPass=1)
    })
    
    n_1class_pathways=length(featSelNet[[1]])
    n_2class_pathways=length(featSelNet[[2]])
    pathways=unlist(featSelNet)
    names(pathways)=c(rep(predClasses[1],n_1class_pathways),rep(predClasses[2],n_2class_pathways))
    list_featSelNet[[inDir]]=pathways
    # index_cl1=!is.na(match(names(pathways),predClasses[1]))
    
  } else {
    #There isn't a result for each run so I don't put the performances in the table
    res_auroc[,k]=rep(0,n_runs)
    res_aupr[,k]=rep(0,n_runs)
    res_acc[,k]=rep(0,n_runs)
  }
  print("I FINISHED THE RUN:")
  print(k)
}
idx_wo_emptiness=colSums(res_acc)!=0
res_auroc=res_auroc[,idx_wo_emptiness]
res_aupr=res_aupr[,idx_wo_emptiness]
res_acc=res_acc[,idx_wo_emptiness]
save(res_auroc,res_aupr,res_acc,list_featSelNet,list_featScore,file=file_rda_2save)

# Supplementary repository to replicate the netDx classification in the paper: Integrated analysis of gene expression profiling, copy number variations and DNA methylation in canine B-cell indolent lymphomas

1. ***/log_netDx_files/*** contains the log files produced automatically by netDx itself to replicate the parameters used, the functions defines and so the classification themselves.

2. ***/input_data/*** contains the datasets and the gmt file of the pathways given in input to netDx:
   - which datasets are given in input for a specific classification are specified in the corresponding Rscript (e.g. the dataset dog_gex_DLBCLvsFL.rda is the first dataset given to netDx for the classification DLBCL vs FL in the Rscript smo_pred_DOG_int_DLBCLvs.R)
   - the order of the dataset given in input for a specific classification are specified in the corresponding Rscript (e.g. in the Rscript smo_pred_DOG_int_DLBCLvs.R for the classification DLBCL vs FL the order is dog_gex_DLBCLvsFL.rda,dog_meth_DLBCLvsFL.rda,dog_cn_DLBCLvsFL.rda)
   - the gmt file is taken in input by R with the command: pathFile="/home/wr/Human_December_2018.gmt"

3. ***/similarity_functions/*** contains the R functions which define the similarity measures
   - the functions are called in the Rscripts
   - the functions are used by netDx to compute the patient similarity networks during the classification
   
4. ***/run_Rscripts/*** contains the Rscripts used to run netDx for each classification and the Rscript to gain the performances from the output directories
   - smo_pred_DOG_int_DLBCLvs is the Rscript used to run netDx for the classifications DLBCLvsFL and DLBCLvsMZL
   - smo_pred_DOG_int_FLvsMZL is the Rscript used to run netDx for the classification FLvsMZL
   - In order to replicate the classifications, the user has to change all the paths defined in these scripts based on its computer and on where he has put the files in its operating system.
   
5. ***/output_performances/*** contains the plots of the AUROC performances

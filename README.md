# Supplementary repository to replicate the netDx classification in the paper: Integrated analysis of gene expression profiling, copy number variations and DNA methylation in canine B-cell indolent lymphomas

********************************
***CONTENT OF THIS REPOSITORY***
********************************
- Description of netDx
- Description of the parameters used, functions defined, input and output about our case in study
- Description of each file in the repository

********************************
***DESCRIPTION OF NETDX***
********************************
netDx is a supervised machine learning method for patient classification that exploits the recent patient similarity network (PSN) paradigm. In a PSN of netDx, each node is a patient and an edge between two of them is weighed with the corresponding pairwise similarity. The weight is computed based on the input patient’s values in a specific data matrix (gene expression, methylation or CNV) featuring only genes belonging to a pathway. As first consequence, netDx creates multiple PSNs from one dataset and each one specifically represents one biological process (feature). Second, converting all the different patient data types in similarity networks/pathways/features, netDx can integrate multiple omics. This enables the identification of strongly connected patient subgroups in the same pathway over all the input data. Briefly, in case FL and MZL are in contrast, the idea is the following one:  Given the 3 PSNs built from the 3 omics (gene expression, methylation, cnv) representing the gap junction pathway, if the FL dogs are very well connected between themselves, the FL dogs are not connected with MZL dogs and the MZL dogs are not very well connected between themselves, then gap junction pathway is a predictive biological process for the FL class.

Apart from the PSN paradigm, netDx starts processing input data in a cross-validation technique. Each class is randomly partitioned producing a training and testing set. The train samples are used to create the similarity networks, one for every pathway with existing genes in the input matrix related to each omic. The PSNs are sparsified for reducing the amount of memory usage. Next, netDx uses GeneMANIA. This algorithm scores each PSN based on how well it classifies one patient class (i.e. a PSN is scored optimal for the FLs if it connects them in a clique without outsiders). Once the PSNs are ranked in each group, netDx starts the validation using the testing set. It predicts each tester’s class based on how much is similar to the train patients in the pathways/PSNs with the highest score of each group (FL or MZL). Finally, netDx returns the AUROC and the most predictive pathways for class. The threshold value for considering AUROC good or bad is 0.5. This threshold is the same defined in the original Zhou et al. paper [1] about statistical methods in diagnostic medicine. More the classifier’s performances are close to 1 and more the classifier is considered able to distinguish between the two groups of interest. Otherwise, it is not considered better than a random classification of the patients. 

********************************
***DESCRIPTION OF OUR CASE***
********************************
In our study, the three data matrices have been pre-processed in the following way. Gene expression have been normalized in LogCPM, methylation levels were kept as counts and copy number data have been converted in LogRatio comparing the original values with respect the controls. Next, in each matrix we kept only the rows associated to a dog gene symbol (restriction given by how the pathways are described). And in case of multiple rows mapping the same gene, we summarised them computing the median by samples (one gene symbol is unique in a pathway). An exception regards the copy number data for classifying FL vs MZL. In fact, for this specific comparison and omic we did not sum up the rows.

After the data preprocessing, netDx v1.0.23 package has been downloaded from github [2]. For each pair of groups in contrast (DLBCL-FL, DLBCL-MZL, FL-MZL) , we gave in input to netDx: the corresponding gene expression dataset, methylation dataset, copy number variation dataset, a matrix containing the patient classes, a list of pathways with their members in HGCN symbol, two values specifying the range of the dimension of the usable pathways, a set of similarity functions for the computation of the patient similarity networks, three values to set the cutoff for the netDx edge removal operation and the parameters to perform 30 runs of cross-validation leave-one-out.  The made classifications based on the input data, user-defined parameters and functions can be replicated using the data provided with the ad-hoc github repository (https://github.com/LucaGiudice/suppl_dog) and following the log files (suppl. files *) produced automatically as output by the algorithm when it finished each classification.

The list of pathways is the one updated until December 2018 by the BaderLab research group and it has been downloaded from the laboratory corresponding site [3]. We used these gene sets because they have been collected, curated and annotated in the same way as the ones used in the netDx paper. The two next input values are 5 and 300. These are used by netDx to filter out all the pathways with a number of genes lower than 5 and greater than 300. As suggested in the netDx paper, this step decreases the chance to consider similar two patients for only one gene and the probability to discriminate the groups for very generic pathways as the apoptosis. About the similarity functions, the Weighted Jaccard has been used to compute similarity between patients based on gene expression data, while the Pearson Correlation has been used with the methylation and the copy number variation data. These two functions have been chosen empirically because they provided the best classification performances. Apart from the Pearson Correlation which is a well-known similarity function for genetic data, the Weighted Jaccard assesses how much two profiles are similar comparing the differences of their values in the same positions. In other words, lower is the difference between the expression of the same genes and greater is the similarity between two patients. The last three input values set thresholds which are used by netDx to remove edges. In fact, following the standard procedure, after the construction of the patient similarity networks, netDx removes all the edges with a similarity lower than 0.3 and, if this is not enough, keeps only 50 edges for each patient. These default operations are useful when the input datasets are large (300 patients or more) because otherwise the classifier is not able to finish the classification for memory issues. In our case we had a number of patients lower than 100 so we gave in input the largest values usable by a computer. In this way, the classifier has been able to use entirely the similarity networks without applying any kind of filtering.

1. ZHOU, X.‐H., OBUCHOWSKI, N. A. and MCCLISH, D. K. Statistical Methods in Diagnostic Medicine
2. https://github.com/BaderLab/netDx
3. http://download.baderlab.org/EM_Genesets/

********************************
***DESCRIPTION OF EACH FILE AND SUBDIRECTORY***
********************************

- ***/log_netDx_files/*** contains the log files produced automatically by netDx itself to replicate the parameters used, the functions defines and so the classification themselves.

- ***/input_data/*** contains the datasets and the gmt file of the pathways given in input to netDx:
   - which datasets are given in input for a specific classification are specified in the corresponding Rscript (e.g. the dataset dog_gex_DLBCLvsFL.rda is the first dataset given to netDx for the classification DLBCL vs FL in the Rscript smo_pred_DOG_int_DLBCLvs.R)
   - the order of the dataset given in input for a specific classification are specified in the corresponding Rscript (e.g. in the Rscript smo_pred_DOG_int_DLBCLvs.R for the classification DLBCL vs FL the order is dog_gex_DLBCLvsFL.rda,dog_meth_DLBCLvsFL.rda,dog_cn_DLBCLvsFL.rda)
   - the gmt file is taken in input by R with the command: pathFile="/home/wr/Human_December_2018.gmt"

- ***/similarity_functions/*** contains the R functions which define the similarity measures
   - the functions are called in the Rscripts
   - the functions are used by netDx to compute the patient similarity networks during the classification
   
- ***/run_Rscripts/*** contains the Rscripts used to run netDx for each classification and the Rscript to gain the performances from the output directories
   - smo_pred_DOG_int_DLBCLvs is the Rscript used to run netDx for the classifications DLBCLvsFL and DLBCLvsMZL
   - smo_pred_DOG_int_FLvsMZL is the Rscript used to run netDx for the classification FLvsMZL
   - In order to replicate the classifications, the user has to change all the paths defined in these scripts based on its computer and on where he has put the files in its operating system.
   
- ***/output_performances/*** contains the plots of the AUROC performances

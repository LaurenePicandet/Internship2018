#########################
##Worflow for performing feature selection + classification on tcga data
##Classification is made according to IDH_codel status 
##Author : Laurène Picandet, using Petr Nazarov scripts and FastICA algorithm implementation
##Many thanks to T.Kaoma 
##Date : 03.05.2018
##############################




##############################
#0. Libraries and scripts loading
##############################
source("functions_for_workflow.R")



##############################
#1. Importing data 
##############################


  ##############################
  #1a. Data type selection
  ##############################
data_type=menu(c("Gene expression rsem", "Gene expression count", "Junction expression"), title="Select data type for analysis")


if(data_type==0){  
  print("Please select a data type !")
  data_type=menu(c("Gene expression rsem", "Gene expression count", "Junction expression"), title="Select data type for analysis")
} else if (data_type==1){
  file_name="GBMLGG.rsem_180203.RData"
} else if (data_type==2){
  file_name="GBMLGG.count_180204.RData"
} else if (data_type==3){
  file_name="GBMLGG.junc_180306.RData"
}



#switching file path, according to the place where is run the script : local computer on Windows or HCP server on Linux
if(Sys.info()["sysname"]=="Windows"){
  prefix_path = "path on local computer" # insert path on local computer 
  nb_cores=1
} else if(Sys.info()["sysname"]=="Linux"){
  prefix_path = "path on remote HCP server" # insert pathon pcp server 
  nb_cores=7                                 # adjustnumber of cores 
}


  ##############################
  #1b. Data loading
  ##############################

rm(Data)                                          #removing pre-existing data 
print(paste0("loading ",prefix_path,file_name))
load(paste0(prefix_path,file_name))               #referred later in the script as "Data"
stopifnot(exists("Data"))                         #checking if data has been correctly loaded 
Y = Data$var$IDH_codel.subtype                    #creating vector with patients classes, will be used in classification



if(data_type==1){
  input_data=Data$rsem
  data_name="rsem"
} else if (data_type==2){
  input_data=Data$count 
  data_name="count"
} else if (data_type==3){
  input_data=Data$junc
  data_name="junction"
}

#cleaning memory
rm(Data,prefix_path,file_name)
gc()

##############################
#2. Pre-process before feature selection
##############################

x=as.matrix(input_data)
x=log2(1+x)

#checking data format
if (substr(colnames(x)[1], 1, 3) == "LGG" |substr(colnames(x)[1], 1, 3) == "GBM" ){ ##if samples in columns
  a=1
} else {
  a=2
}


x=filtering(x,3,a)



##############################
#3. Feature selection
##############################


feature_selection_method=menu(c("None", "PCA", "ICA", "top AUC", "top low FDR"), title="Which feature selection method would you like to use ? ")


#cleaning memory 
#rm(input_data)  # the original expression matrix is needed later, because we use it after DEG analysis to retrieve selected features from it 
#gc()


if (feature_selection_method == 1){ 
  method_name = "No feature selection"
  result_fs = noFeatureSelection_preparation(input_data) ## all pre-processing will be done in this function, this is why I take 'input_data' (which is unprocessed), and not x
  
  } else if (feature_selection_method == 2){ #PCA
  method_name = "PCA"
  start_time_fs = Sys.time()
  result_fs = PCA(x)
  end_time_fs = Sys.time()
  matrix_for_classification = result_fs$x[,1:100] 
  
} else if (feature_selection_method == 3){ #ICA
  method_name = "ICA"
  start_time_fs = Sys.time()
  result_fs = runICA(x, ncomp = 100, ntry = 20, ncores = nb_cores) ## ask for numbers of components and numbers of tries
  end_time_fs = Sys.time()
  matrix_for_classification = result_fs$M
  
} else if (feature_selection_method == 4){ #Differential expression analysis, selection of features according to their AUC
  method_name = "DEG_top_AUC"
  
  status_pair=menu(c("IDHmut-codel/IDHmut-non-codel", "IDHmut-codel/IDHwt", "IDHmut-non-codel/IDHwt"), title="Between which status pairs should the differential expression analysis be performed ? ")
  start_time_fs = Sys.time()
  if(status_pair==1){
    result_fs = DEA.limma(x, group = Y, key1="IDHmut-codel", key0="IDHmut-non-codel", counted=FALSE, return.auc=TRUE) #counted=TRUE if no log transformation
    matrix_for_classification = AUC_analysis(result_fs, input_data)
    
    
  } else if (status_pair==2){
    result_fs = DEA.limma(x, group = Y, key1="IDHmut-codel", key0="IDHwt", counted=FALSE, return.auc=TRUE) #counted=TRUE if no log transformation
    matrix_for_classification = AUC_analysis(result_fs, input_data)
    
    
  } else if (status_pair==3){
    result_fs = DEA.limma(x, group = Y, key1="IDHmut-non-codel", key0="IDHwt", counted=FALSE, return.auc=TRUE) #counted=TRUE if no log transformation
    matrix_for_classification = AUC_analysis(result_fs, input_data)
    
  }
  
  end_time_fs = Sys.time()
}


# maybe useless to store method name because we can guess it from the result format ? 

save(result_fs, file = paste0("result_",data_name,"_",method_name,"_",Sys.Date()))
cat(method_name, " calculation time : ")
print(end_time_fs - start_time_fs)




##############################
#4. Pre-process before classification
##############################
#insert here matrix and vector for classification if you don't do the whole workflow


# transpose if necessary 
if (substr(colnames(matrix_for_classification)[1], 1, 3) == "LGG" |substr(colnames(matrix_for_classification)[1], 1, 3) == "GBM" ){ ##if samples in columns
  m = t(matrix_for_classification)
  print("transposing matrix before classification ")
} else {
  m = matrix_for_classification
}

#remove Na values in labels vector 
id_na=is.na(Y)
if(any(id_na)){
  m = m[!id_na,]
  Y = Y[!id_na]
}
stopifnot(nrow(m) == length((Y)))
stopifnot(all(!is.na(Y)))


#shuffle the training and testing data 
irand = sample(nrow(m)) #vector with the rows numbers of expression matrix  randomly shuffled
m = m[irand,] # matrix with samples in rows randomly shuffled and junctions in columns
Y = Y[irand] #vector with class labels randomly shuffled

stopifnot(nrow(m) == length((Y)))
stopifnot(all(!is.na(Y)))


#cleaning memory 
rm(matrix_for_classification, result_fs)
gc()


##############################
#5 . Classification
##############################



classification_start_time = Sys.time()
class_result = cross_validation_classif_rf(m, Y, 10)  # change function if you also want svm results 
print(Sys.time() - classification_start_time)


#class_result = list with 10 factors of predicted subtypes for all patients
##############################
#6 . Analysing classification results 
##############################



accuracies = lapply(class_result, 
              function(x){

                # Y = vector with observed subtype
                # x = factor with predicted subtype
                str(x)
                class(x)
                CM.rf = getConfusionMatrix(Y,x) #creation of a confusion matrix
                print(CM.rf)
                #print(getMCError(CM.rf))
                accuracy_RF=1-getMCError(CM.rf) #obtaining the misclassification error and deducting accuracy 
              
              }
)


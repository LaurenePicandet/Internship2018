#comparing the classification accuracies obtained by RF between the different conditions 

#loading the data (not progvided on GitHub)
#R object will all accuracies : "ACCUR_RF"

str(ACCUR_RF) #observation of the structure of the object, result below 

# 'data.frame':	10 obs. of  8 variables:
# $ IC_junction  : num  0.969 0.97 0.969 0.966 0.969 ...
# $ PC_junction  : num  0.941 0.945 0.941 0.942 0.947 ...
# $ DE_junction  : num  0.945 0.945 0.948 0.948 0.948 ...
# $ rand_junction: num  0.884 0.892 0.895 0.88 0.892 ...
# $ IC_genes     : num  0.945 0.941 0.936 0.945 0.941 ...
# $ PC_genes     : num  0.913 0.917 0.911 0.91 0.91 ...
# $ DE_genes     : num  0.984 0.984 0.984 0.981 0.982 ...
# $ rand_genes   : num  0.953 0.956 0.951 0.944 0.954 ...

#the 10 observations are produced during the 10-fold cross validation 

#calculation of the mean accuracy of the 10 folds for each condition :
means = lapply(ACCUR_RF,mean)
str(means) 
# List of 8
# $ IC_junction  : num 0.969
# $ PC_junction  : num 0.942
# $ DE_junction  : num 0.947
# $ rand_junction: num 0.887
# $ IC_genes     : num 0.941
# $ PC_genes     : num 0.911
# $ DE_genes     : num 0.983
# $ rand_genes   : num 0.952



#########################################################
################Difference gene/junction data ###########
#########################################################



#calculation of the mean classification accuracy using gene expression

mean_gene = mean(unlist(means[5:8])) #unlist() to transform the list in a vector, required for mean calculation
# 0.9464815


#calculation of the mean classification accuracy using junction expression
mean_junction = mean(unlist(means[1:4]))
# 0.9362222

#difference 

diff = mean_gene - mean_junction 
 # 0.01025926



#########################################################
#########Accuracies by methods###########################
#########################################################

mean_ICA  = mean(unlist(means[c(1,5)])) #0.954963
mean_PCA  = mean(unlist(means[c(2,6)])) #0.9260741
mean_DEA  = mean(unlist(means[c(3,7)])) #0.9647407
mean_random = mean(unlist(means[c(4,8)])) #0.9196296

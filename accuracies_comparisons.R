#comparing the classification accuracies obtained by RF between the different conditions 

load("accuracies_rf_idh_2018-05-23")
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



#########################################################
################Difference gene/junction data ###########
#########################################################



#calculation of the mean classification accuracy using gene expression
gene = c(
  ACCUR_RF$IC_genes,     
  ACCUR_RF$PC_genes,   
  ACCUR_RF$DE_genes,     
  ACCUR_RF$rand_genes
)

mean_accuracies_gene = mean(gene)
confidence_interval_gene = qt(0.025, 19) * sd(gene) / sqrt(20) 

# 0.9464815


#calculation of the mean classification accuracy using junction expression
junction = c(
  ACCUR_RF$IC_junction,     
  ACCUR_RF$PC_junction,   
  ACCUR_RF$DE_junction,     
  ACCUR_RF$rand_junction
)

mean_accuracies_junction = mean(junction)
confidence_interval_junction = qt(0.025, 19) * sd(junction) / sqrt(20) 


#difference 

diff = mean_accuracies_gene - mean_accuracies_junction
# 0.01025926



#########################################################
#########Accuracies by methods + data type###############
#########################################################

means = colMeans(ACCUR_RF)
# IC_junction  PC_junction   DE_junction  rand_junction      
# 0.9691852    0.9416296     0.9466667    0.8874074     
# 
# 
# IC_genes     PC_genes      DE_genes     rand_genes 
# 0.9407407    0.9105185     0.9828148    0.9518519 


conf_int = function(x){
  n = length(x)
  CI = qt(0.025, n-1) * sd(x) / sqrt(n)
  return(CI)
}

confidence_intervals= apply(ACCUR_RF,2, conf_int)
# IC_junction      PC_junction   DE_junction    rand_junction      
# -0.001563961     -0.002605537  -0.001413051   -0.003351344 

#  IC_genes        PC_genes      DE_genes       rand_genes 
#  -0.002177657   -0.002745467   -0.001023852   -0.002559631 


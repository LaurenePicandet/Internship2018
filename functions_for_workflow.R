#################
###Functions and libraries used in worflow feature selection + classification
################

library(e1071) ## for SVM
library(randomForest) ## for randomForest
library(parallel)   ## for script parallelisation



PCA = function(x1,pc=NULL,scale=TRUE,ipc=1:2,imputed = NULL){
  #first make sure that the data is in a matrix form 	
  x1 = as.matrix(x1)
  #then make sure there are no missing values
  ## if there are NA - impute : here the missing values are replaced by the median value of the data set 
  if (sum(is.na(x1))>0){
    if (is.null(imputed)) imputed = median(x1,na.rm=TRUE)
    x1[is.na(x1)]=imputed
  }
  ## and make sure none of the variables has a null variance. 
  if (sum(apply(x1,1,var)==0)>0) #if one of the variables shows a null variance
    x1 = x1[apply(x1,1,var)>0,] #select only the lines of the matrix with a variance superior to 0 
  
  ##calculate PC
  if (is.null(pc)) pc=prcomp(t(x1),scale = scale)  #the matrix is transposed so that the genes (features) are in columns and the patients in rows 
  
  return(pc) ## here are the principal components 
}



###############################################################################
## Warp-ups for stable ICA: multi-run ICA providing consensus S and M matrices
## 
## GNU GPL (c) P.Nazarov, LIH, 2016-09-05 -> 2018-01-24
## Thanks: T.Kaoma and entire BIOMOD team
## Supported by FNR CORE grant 'DEMICS' 2018-2019
###############################################################################

##=============================================================================
## runICA - runs ICA several times using parallelization
##	X - data matix (samples - columns, features - rows)
##	ncomp - number of components
##	ntry - number of runs
##	show.every - how often do we show the progress messages 
##	ncores - number of cores to be set for calculation
##	Returns: list with X, M, S, etc
##-----------------------------------------------------------------------------
runICA = function(X,ncomp=3, ntry = 1, show.every=1, filter.thr = NULL,ncores=1){ #keep=NULL,
  
  ## install packages if absent
  if (!"fastICA" %in% rownames(installed.packages())){
    print("Cannot find `fastICA` - installing it from Bioconductor")
    source("https://bioconductor.org/biocLite.R")
    biocLite("fastICA")
  }
  if (ncores>1 &.Platform$OS.type=="unix" & !"doMC" %in% rownames(installed.packages())){
    print("Cannot find `doMC` - installing it from Bioconductor")
    source("https://bioconductor.org/biocLite.R")
    biocLite("doMC")
  }
  if (ncores>1 &.Platform$OS.type=="windows" & !"doSNOW" %in% rownames(installed.packages())){
    print("Cannot find `doSNOW` - installing it from Bioconductor")
    source("https://bioconductor.org/biocLite.R")
    biocLite("doSNOW")
  }
  
  require(fastICA)
  X = as.matrix(X)
  if (!is.null(filter.thr)) X = X[apply(X,1,max)>filter.thr,]
  Res = list()
  Res$X = X
  S = list()
  M = list()
  S[[1]] = matrix(nrow=nrow(X),ncol=ncomp)
  rownames(S[[1]]) = rownames(X)
  colnames(S[[1]]) = sprintf("ic.%d",1:ncomp)
  ## mixing matrix
  M[[1]] = matrix(nrow=ncomp,ncol = ncol(X))
  colnames(M[[1]]) = colnames(X)
  rownames(M[[1]]) = sprintf("ic.%d",1:ncomp)
  itry=1
  Res$S = S[[1]]
  Res$M = M[[1]]
  Res$S.best = S[[1]]
  Res$M.best = M[[1]]
  Res$mse = NA  ## mean square error bw X and S*M
  Res$mr2 = NA  ## mean correlation bw mixing profiles
  Res$n3 = NA  ## mean number of elements in |S| over 3
  ## do multiple tries
  idx.excludeSamples = sample(1:ncol(X),ntry, replace = (ntry > ncol(X)))
  if (ntry==1) idx.excludeSamples = integer(0)
  itry = 1
  
  ###############################
  ## Parallel section starts
  ##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if (ncores > 1) {
    require(foreach)
    if(.Platform$OS.type=="unix") {
      require(doMC)
      registerDoMC(ncores)
    }
    if(.Platform$OS.type=="windows") {
      require(doSNOW)
      cl = makeCluster(ncores) 
      registerDoSNOW(cl)
    }
  }
  
  cat("*** Starting",ifelse(ncores>1,"parallel",""),"calculation on",ncores,"core(s)...\n")
  cat("*** System:",.Platform$OS.type,"\n")
  cat("***",ncomp,"components,",ntry,"runs,",nrow(X),"features,",ncol(X),"samples.\n")
  cat("*** Start time:",as.character(Sys.time()),"\n")
  t0=Sys.time()
  flush.console()
  ## multi run ICA
  if (ncores > 1) {
    MRICA = foreach(itry=1:ntry) %dopar% {
      require(fastICA)
      SP = Res$S + NA
      MP = Res$M + NA
      if (length(idx.excludeSamples) == 0){
        ic = fastICA(X, n.comp = ncomp, alg.typ ="deflation")
        SP[,] = ic$S
        MP[,] = ic$A
      }else{
        x = X[,-idx.excludeSamples[itry]]
        ic = fastICA(x, n.comp = ncomp, alg.typ ="deflation")
        SP[,]= ic$S
        MP[,-idx.excludeSamples[itry]] = ic$A
        MP[is.na(MP)] = 0
      }
      return(list(S=SP,M=MP))
    }
  }else{
    require(fastICA)
    MRICA = list()
    for(itry in 1:ntry){
      MRICA[[itry]] = list()
      MRICA[[itry]]$S = Res$S + NA
      MRICA[[itry]]$M = Res$M + NA
      if (length(idx.excludeSamples) == 0){
        ic = fastICA(X, n.comp = ncomp, alg.typ ="deflation")
        MRICA[[itry]]$S[,] = ic$S
        MRICA[[itry]]$M[,] = ic$A
      }else{
        x = X[,-idx.excludeSamples[itry]]
        ic = fastICA(x, n.comp = ncomp, alg.typ ="deflation")
        MRICA[[itry]]$S[,]= ic$S
        MRICA[[itry]]$M[,-idx.excludeSamples[itry]] = ic$A
        MRICA[[itry]]$M[is.na(MRICA[[itry]]$M)] = 0
      }
      if (itry%%show.every == 0) {
        cat("try #",itry,"of",ntry,"\n")
        flush.console()
      }
    }
  }
  
  if(.Platform$OS.type=="windows" & ncores>1)  stopCluster(cl)
  
  cat("*** Done!", "\n")
  cat("*** End time:",as.character(Sys.time()),"\n")
  print(Sys.time()-t0)
  
  ##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ## ...end of parallel section
  ###############################
  
  S = lapply(MRICA,function(x)x$S)
  M = lapply(MRICA,function(x)x$M)
  rm(MRICA)
  
  for (itry in 1:ntry){
    Res$mse[itry] = sum(X - mean(X) - S[[itry]]%*%M[[itry]])^2 / ncol(X) / nrow(X)
    Res$mr2[itry] = mean((cor(t(M[[itry]]))^2)[upper.tri(matrix(0,nrow=ncomp,ncol=ncomp))])
    Res$n3[itry] = sum(apply(abs(S[[itry]])>3,2,sum))/ncomp
  }
  ## who is the best: min error and min correlated 
  Res$i.best = which.min(scale(Res$mse)+scale(Res$mr2))
  ## simply smallest error? #Res$i.best = which.min(Res$mse)
  if (length(Res$i.best) == 0 ) Res$i.best = 1
  Res$S.best = S[[Res$i.best]]
  Res$M.best = M[[Res$i.best]]
  ## if only one try - return
  if (ntry == 1) {
    Res$S = S[[1]]
    Res$M = M[[1]]
    return(Res)
  }
  
  ## correlate results
  ## s.cor - to which ic of the BEST decomposition we should address? 
  s.cor = matrix(nrow=ntry,ncol=ncomp)
  s.cor[Res$i.best,] = 1:ncomp
  itry = 1
  for (itry in (1:ntry)[-Res$i.best]) {
    r = cor(S[[itry]],S[[Res$i.best]])
    s.cor[itry,] = apply((r)^2,2,which.max)
    for (ic in 1:ncomp)
      s.cor[itry,ic] = s.cor[itry,ic] * sign(r[s.cor[itry,ic],ic])
  }
  ## build consensus S, M
  Res$S[,] = S[[1]]
  Res$M[,] = M[[1]]
  itry=2
  for (itry in 2:ntry){
    for (ic in 1:ncomp) {
      Res$S[,ic] = Res$S[,ic] + S[[itry]][,abs(s.cor[itry,ic])]* sign(s.cor[itry,ic])
      Res$M[ic,] = Res$M[ic,] + M[[itry]][abs(s.cor[itry,ic]),]* sign(s.cor[itry,ic])
    }
  }
  Res$S = Res$S / ntry
  Res$M = Res$M / ntry
  
  ## use consensus S, M to analyze stability
  Res$stab = s.cor + NA
  for (itry in (1:ntry)) {
    Res$stab[itry,] = diag(cor(Res$S,S[[itry]][,abs(s.cor[itry,])])^2)
  }
  colnames(Res$stab) = colnames(Res$S) 
  
  return(Res)
}





###############################################################################
## Collection of various warp-ups for DEA methods: RNAseq and microarrays
## (c) GNU GPL Petr Nazarov, Luxembourg Institute of Health, petr.nazarov[at]lih.lu
## last modification 2017-02-14
###############################################################################
library(compiler)
## I can use dir.create() !!!
requireFolder = function(folder=NULL,relative=TRUE){
  if (is.null(folder)) return(FALSE)
  if (folder=="") return(FALSE)
  if (!file.exists(folder)) {
    cat(sprintf("Folder '%s' is not found. Creating...\n",folder))
    if (Sys.info()[1] == "Windows"){
      if (relative){
        #try(shell(paste("mkdir",gsub("/","\\\\",file.path(getwd(),folder)))))
        try(shell(paste("mkdir \"",gsub("/","\\\\",file.path(getwd(),folder)),"\"",sep="" )))
        ## ToDo: try here parameter tranaslate=TRUE for slash->backslash
      }else{
        try(shell(paste("mkdir",gsub("/","\\\\",folder))))
      }
    }
    if (Sys.info()[1] == "Linux")
      try(system(  paste("mkdir",file.path(getwd(),folder)) ))    
    if (!file.exists(folder))
      stop("Cannot create folder to store dowloaded data!\n")
  }
  #return(TRUE)
}

getTopIdx=function(x,n){
  #return(sort(x, index.return=TRUE,decreasing=TRUE, na.last =TRUE)$ix[1:n])
  return(order(x,na.last=TRUE,decreasing=TRUE)[1:n])
}
getConcord = function(list1,list2){
  return( (length(intersect(list1,list2))/length(list1) + length(intersect(list1,list2))/length(list2))/2)
}
getJaccard = function(list1,list2){
  return( length(intersect(list1,list2))/(length(union(list1,list2))))
}
getMeanCI = function(x,alpha=0.05,do.print=FALSE){
  n=length(x[!is.na(x)])
  out = list()
  out$m = mean(x,na.rm=TRUE)
  out$me = -qt(alpha/2,n-1)*sd(x,na.rm=TRUE)/sqrt(n)
  
  if (do.print) print(sprintf("Mean(x[%d]) = %g +/- %g (CI %g%%)",n,out$m,out$me,100*(1-alpha)))
  return(out)
}

##=============================================================================
## calc.AUC - calculate Area Under roc-Curve 
##	x - data vector
##	group - classes vector
##	key0 - id for class 0 in group
##	key1 - id for class 1 in group
##=============================================================================
calc.AUC=function(x,group,key0,key1){
  require("caTools")
  idx0 = which(group == key0)
  idx1 = which(group == key1)
  group = group[c(idx0,idx1)]
  if (class(x) =="numeric" || class(x) =="integer"){
    auc = colAUC(t(x[c(idx0,idx1)]), factor(as.character(group)), plotROC=F)[1,1]
  }else{
    auc = colAUC(t(x[,c(idx0,idx1)]), factor(as.character(group)), plotROC=F)
  }
  return(auc)
}
##=============================================================================
## Normalization
## method = "DESeq","edgeR","voom"
##=============================================================================
if (FALSE){
  method="DESeq"
  doLog=TRUE
  group=NULL
}
norm.Counts=function(count, group=NULL, method="DESeq", doLog = FALSE){
  if (!method %in% c("DESeq","edgeR","voom")) stop("in norm.Counts {LibDEA}: unknown method",method,"\n")
  require(edgeR)
  require(DESeq2)
  if (is.null(group)) group = rep("one",ncol(count))
  group = factor(group)
  
  if (method =="edgeR"){
    dge = DGEList(counts=round(count),group=group)
    dge = calcNormFactors(dge,method = "TMM")
    x = cpm(dge, normalized.lib.sizes=TRUE)
    if (doLog) x = log2(0.5 + x)
  }
  if (method =="DESeq"){
    colData = data.frame(condition=group)
    if (nlevels(group)==1) {
      des = DESeqDataSetFromMatrix(round(count),colData, formula(~1))
    }else{
      des = DESeqDataSetFromMatrix(round(count),colData,formula(~ condition))
    }
    #des = estimateSizeFactors(des,type="iterate")
    des = estimateSizeFactors(des)
    x = counts(des,normalize=T)
    if (doLog) x = log2(1 + x)
  }
  return(x)
}

Count2FPKM = function(X, len){
  FPKM = X * 0
  for (i in 1:ncol(X)) FPKM[,i] = X[,i] / len / sum(X[,i]) * 1e9
  return(FPKM)
}

Count2TPM= function(X, len){
  TPM = X * 0
  for (i in 1:ncol(X)){
    TPM[,i] = X[,i] / len
    TPM[,i] = TPM[,i] /sum(TPM[,i]) * 1e6
  }
  return(TPM)
}
##=============================================================================
## DESeq
##=============================================================================
if (FALSE){
  count = NGS$counts
  group = paste(meta$Type,meta$State,sep=".")
  key0 = paste(type,"Normal",sep=".")
  key1 = paste(type,"Tumor",sep=".")
  name="zzz"
  folder="DEA"
  pair=NULL
}
DEA.DESeq = function(count,group,key0,key1,pair=NULL,name=NULL,folder="",return.auc=FALSE,adjust = "BH"){
  require(DESeq2)
  requireFolder(folder)
  if (ncol(count)!=length(group)) 
    stop("Length of 'group' should corresponds to the number of data columns")
  if (!is.null(pair))
    if (ncol(count)!=length(pair)) 
      stop("Length of 'pair' should corresponds to the number of data columns")
  idx0 = which(key0==group)
  idx1 = which(key1==group)
  group = factor(group[c(idx0,idx1)])
  group = relevel(group,key0)
  if (!is.null(pair)) {
    pat = factor(pair[c(idx0,idx1)])
    colData = data.frame(condition=group,pair=pat)
    des = formula(~ condition + pair)
  }else{
    colData = data.frame(condition=group)
    des = formula(~ condition)
  }
  des = DESeqDataSetFromMatrix(count[,c(idx0,idx1)],colData,design = des)
  des = DESeq(des)
  resultsNames(des)
  #v. <1.3 - prev.version!!!!!!! check versions!
  #res = results(des,sprintf("condition_%s_vs_%s",levels(group)[2],levels(group)[1]),pAdjustMethod = "BH")
  #v. >=1.3
  res = results(des,contrast=c("condition",levels(group)[2],levels(group)[1]),pAdjustMethod = adjust)
  #str(res)
  Tab = data.frame(id=res@rownames,res@listData,stringsAsFactors=F)
  names(Tab)=sub("log2FoldChange","logFC",sub("padj","FDR",names(Tab)))
  rownames(Tab)=Tab[,1]
  Tab$FDR[is.na(Tab$FDR)]=1
  res = Tab
  cat(sprintf("\nDESeq,%s%s: %d DEG (FDR<0.05), %d DEG (FDR<0.01)\n\n",
              ifelse(is.null(pair),"unpaired","paired"),
              ifelse(is.null(name),"",paste(" on",name)),
              sum(res$FDR<0.05),sum(res$FDR<0.01)))
  if (return.auc){
    if (require(caTools)){
      res$AUC = calc.AUC(x = count[,c(idx0,idx1)],group = group, key0,key1)[1,]
      ## count[,..] not all count as gtoup is reassigned
    }
  }
  if (!is.null(name)){
    if (!is.null(pair)){
      file.name = file.path(folder,sprintf("DEA(DESeq,paired)=%s_%s-%s.txt",name,key1,key0))
    }else{
      file.name = file.path(folder,sprintf("DEA(DESeq)=%s_%s-%s.txt",name,key1,key0))
    }
    write.table(res,file=file.name,sep="\t",quote=F,row.names=F)
  }
  
  res$logFC[is.na(res$logFC)] = 0
  return(res)
}


##=============================================================================
## edgeR
##=============================================================================
if (FALSE){
  
  #	ResP = DEA.edgeR(X,group=group,key0=key0,key1=key1,pair=xmeta$Batch,name=NULL,folder="DEA")
  count = X
  group = group
  key0 = key0
  key1 = key1
  pair=xmeta$Batch
  prior.df=10
  name=NULL
  folder=""
  return.auc=FALSE
  return.model=TRUE
  prefer.GLM=TRUE
}
##--------------------------------
DEA.edgeR=function(count,group,key0,key1,pair=NULL,prior.df=10,name=NULL,folder="",return.auc=FALSE,return.model=FALSE, prefer.GLM=FALSE,adjust="BH"){
  require(edgeR)
  requireFolder(folder)
  idx0 = which(group == key0)
  idx1 = which(group == key1)
  group = factor(group[c(idx0,idx1)])
  group = relevel(group,key0)
  dge = DGEList(counts=count[,c(idx0,idx1)],group=group)
  if (!is.null(pair) | prefer.GLM){ ## use GLM
    if (!is.null(pair)){
      pair = factor(pair[c(idx0,idx1)])
      design = model.matrix(~pair+group)  ## important factor order for glmLRT - (1)pat (2)group
    }else{
      design = model.matrix(~group)
    }
    dge = calcNormFactors(dge)
    dge = estimateGLMCommonDisp(dge,design)
    dge = estimateGLMTrendedDisp(dge,design)
    dge = estimateGLMTagwiseDisp(dge,design)
    fit = glmFit(dge, design)
    lrt = glmLRT(fit)
    res = topTags(lrt,n=nrow(count),adjust=adjust,sort.by="none")[[1]]
    res = data.frame(id=rownames(res),res,stringsAsFactors=F)
  }else{ ## use exactTest
    dge = calcNormFactors(dge)
    dge = estimateCommonDisp(dge)
    dge = estimateTagwiseDisp(dge,prior.df=prior.df)
    res= exactTest(dge)[[1]]
    res$FDR=p.adjust(res$PValue,"fdr")
    res = data.frame(id=rownames(res),res,stringsAsFactors=F)
  }
  cat(sprintf("\nedgeR,%s%s: %d DEG (FDR<0.05), %d DEG (FDR<0.01)\n\n",
              ifelse(is.null(pair),"unpaired","paired"),
              ifelse(is.null(name),"",paste(" on",name)),
              sum(res$FDR<0.05),sum(res$FDR<0.01)))	
  
  if (return.auc){
    if (require(caTools)){
      res$AUC = calc.AUC(x = count[,c(idx0,idx1)],group = group, key0,key1)[1,]
      ## count[,..] not all count as group is reassigned
    }
  }
  
  if (!is.null(name)){
    if (!is.null(pair)){
      file.name = file.path(folder,sprintf("DEA(edgeR,paired)=%s_%s-%s.txt",name,key1,key0))
    }else{
      file.name = file.path(folder,sprintf("DEA(edgeR)=%s_%s-%s.txt",name,key1,key0))
    }
    write.table(res,file=file.name,sep="\t",quote=F,row.names=F)
  }
  
  if (return.model){
    TMP = res
    res=list()
    res$name = ifelse(is.null(pair),"unpaired","paired")
    if (exists("fit")){
      res$dge = dge
      res$lmFit = fit
      res$lmLRT = lrt
    }else{
      res$model = "exactTest was used. No model."
    }
    res$topTable = TMP
  }
  
  return(res)
}

##=============================================================================
## LIMMA
##=============================================================================

## ToDo - correct & add q-q norm
DEA.limma=function(data=NULL,group,key0,key1,pair=NULL,name=NULL,folder="",counted=FALSE,norm="none",return.model=FALSE,return.auc=FALSE,adjust="BH",silent = FALSE){
  require(limma)
  requireFolder(folder)
  #idx0 = grep(key0,group)
  #idx1 = grep(key1,group)
  idx0 = which(group == key0)
  idx1 = which(group == key1)
  group = factor(group[c(idx0,idx1)])
  group = relevel(group,key0)
  
  if (is.null(rownames(data))) rownames(data) = sprintf("row%05d",1:nrow(data))
  
  if (!is.null(pair)){ ## paired
    pair = factor(pair[c(idx0,idx1)])
    design = model.matrix(~pair+group)  ## important factor order for glmLRT - (1)pat (2)group
    x = data[,c(idx0,idx1)]
    if (counted) {
      x = voom(x,design=design,plot=FALSE,normalize.method=norm)
    }else{
      if (norm=="scale") x = scale(x)
      if (norm=="quantile") x = normalizeQuantiles(x, ties=TRUE)
    }
    fit = lmFit(x,design=design)
    EB = eBayes(fit)
    res = data.frame(id=rownames(data),topTable(EB,coef=ncol(design),number=nrow(data),adjust=adjust,sort.by="none"),stringsAsFactors=F)
  }else{ ## unpaired
    design=model.matrix(~ group)
    #design=model.matrix(~ -1 + group)
    x = data[,c(idx0,idx1)]
    if (counted) {
      x = voom(x,design=design,plot=FALSE,normalize.method=norm)
    }else{
      if (norm=="scale") x = scale(x)
      if (norm=="quantile") x = normalizeQuantiles(x, ties=TRUE)
    }
    fit = lmFit(x,design=design)
    EB = eBayes(fit)
    res = data.frame(id=rownames(data),topTable(EB,coef=2,number=nrow(data),adjust=adjust,sort.by="none"),stringsAsFactors=F)
    #res = data.frame(id=rownames(data),topTable(EB,coef=1,number=nrow(data),adjust=adjust,sort.by="none"),stringsAsFactors=F)
  }
  names(res)=sub("adj.P.Val","FDR",names(res))	
  
  if (!silent)
    cat(sprintf("\nLimma,%s%s: %d,%d,%d,%d DEG (FDR<0.05, 1e-2, 1e-3, 1e-4)\n\n",
                ifelse(is.null(pair),"unpaired","paired"),
                ifelse(is.null(name),"",paste(" on",name)),
                sum(res$FDR<0.05),sum(res$FDR<0.01),sum(res$FDR<0.001),sum(res$FDR<0.0001)))
  
  if (return.auc){
    if (require(caTools)){
      res$AUC = calc.AUC(x = data[,c(idx0,idx1)],group = group, key0,key1)[1,]
      ## count[,..] not all count as gtoup is reassigned
    }
  }
  if (!is.null(name)){
    if (!is.null(pair)){
      file.name = file.path(folder,sprintf("DEA(limma,paired)=%s_%s-%s.txt",name,key1,key0))
    }else{
      file.name =file.path(folder,sprintf("DEA(limma)=%s_%s-%s.txt",name,key1,key0))
    }
    write.table(res,file=file.name,sep="\t",quote=F,row.names=F)
  }
  
  if (return.model){
    TMP = res
    res=list()
    res$name = ifelse(is.null(pair),"unpaired","paired")
    res$lmFit = fit
    res$eBayes = EB
    res$topTable = TMP
  }
  
  return(res)
}






######################################
####Pre-process before classification
####################################

## x = expression matrix, genes, junctions ( or matrix with level of methylations in beta values for example)
## thr = min threshold that the maximum expression value of a gene/junction 
## must reach in order to be kept 
## a = 1 if features in rows
## a = 2 if features in columns 

filtering=function(x,thr,a){ 
  ikeep = apply(x,a,max) > thr
  sum(ikeep)
  b = nrow(x)
  X = x[ikeep,]
  c = nrow(X)
  print(paste0(b-c," junctions/genes with low expression removed. Threshold = ",thr))
  return(X)
  
}

noFeatureSelection_preparation=function(x){
  x=as.matrix(x)
  if ("GBM" %in% colnames(input_data)[1]|"LGG" %in% colnames(input_data)[1]){
    x=t(x)  ## we want the samples in rows and the features in colums 
  }
  x=log2(1+x)

  
  
}

######################################
####Functions for classification
####################################

getConfusionMatrix = function(gr, gr.pred) {  ## gr = real classes  gr.pred= predicted classes
  if (class(gr) == "character") gr = factor(gr); gr.pred = factor(gr.pred)   
  if (class(gr) == "factor") gr = as.integer(gr); gr.pred = as.integer(gr.pred)  
  #--> transformation of the vectors containing real and predicted classes into integer factors with 3 levels
  # exemple : Y = 3 3 2 3 3 3 3 2 3 3 ... with : 
  #  1 = IDHmut-codel
  #  2 = IDHmut-non-codel
  #  3 = IDHwt  
  
  n = max(c(gr, gr.pred)) 
  
  print(paste0("confusion matrix of size ",n))
  Tab = matrix(nc = n, nr = n) #creation of a squared matrix that will become the confusion matrix
  rownames(Tab) = paste("pred", 1:n, sep = ".")
  colnames(Tab) = paste("group", 1:n, sep = ".")
  for (i in 1:n)
    for (j in 1:n){
      Tab[i,j] = sum((gr.pred == i) & (gr== j))
    }
  return(Tab)
}





### 3B Function that computes the misclassification error from a confusion matrix.  
## accuracy can be deduced from misclassification error : accuracy=1-getMCerror 
getMCError = function(CM) {
  x = 0
  for (i in 1:ncol(CM))
    x = x + CM[i,i]
  return(1 - x / sum(CM))
}



####
#function to apply classification with x-fold cross_validation, with RF and SVM algo
##
cross_validation_classif = function(X, Y, nb_break) {
  Pred.svm = Y #just to create container with the proper levels
  Pred.rf = Y
  Pred.svm[]=NA ## and now - clean it (we do not want to cheat ;)
  Pred.rf[]=NA
  
  
  #Create x equally size folds 
  folds = cut(seq(1,nrow(X)),breaks=nb_break,labels=FALSE)
  
  print(paste0("Classification start time : ",Sys.time(),"\n"))
  start_time = Sys.time()
  
  fold = 1
  for (fold in unique(folds)){ ## cross-validation loop (10 times)
    cat("Starting Cross-validation:",fold,"of ", nb_break," time = ",Sys.time(),"\n" )
    idx.test = which(fold == folds)    #vector with the indexes of one of the 10 parts of the data, that will be used as the training set
    idx.train = (1:nrow(X))[-idx.test] # the other 9 parts of the data will create the training set
    mod.svm = svm(x=X[idx.train,],y=Y[idx.train])   #creation of a classifier with the svm method
    mod.rf = randomForest(x=X[idx.train,],y=Y[idx.train], do.trace=TRUE) #creation of a classifier with the random forest method
    #do.trace = TRUE for more verbose output 
    Pred.svm[idx.test] = predict(mod.svm,X[idx.test,])  #applying the svm classifier on the data. output = vector with predicted classes of the data
    Pred.rf[idx.test] = predict(mod.rf,X[idx.test,])   #applying the RF classifier on the data. output = vector with predicted classes of the data
    
    flush.console()
  }
  end_time = Sys.time()
  print("Calculation time classification : ")
  print(end_time - start_time)
  
  result=list(
    idx.test = idx.test,
    pred.svm = Pred.svm,
    pred.rf  = Pred.rf,
    mod.svm  = mod.svm,
    mod.rf   = mod.rf
  )
  

  
  save(result, file=paste0("classification_result_",fold,".RData")) ##saving intermediate result
  return(result)
  #return(calc_time)
  
}


## addition of the saving of the accuracy at each iteration of 10 fold cross validation 
cross_validation_classif_rf = function(X, Y, nb_break) {
  Pred.rf = Y
  Pred.rf[]=NA
  
  
  #Create x equally size folds 
  folds = cut(seq(1,nrow(X)),breaks=nb_break,labels=FALSE)  #integer vector [1 1 1 ...2 2 2[..]10 10 10...]
  
  list_fold  = unique(folds) # integer vector : [1 2 3 4 5 6 7 8 9 10]
 
  
  ## cross-validation loop 
  ## predictions contains a list of nb_break factors
      # each factor is the result of one validation loop and contains subtype prediction for all patients
  
  predictions = lapply(list_fold, 
  
         function(x) { 

            cat("Starting Cross-validation:",x,"of ", nb_break," time = ",Sys.time(),"\n" )
            idx.test = which(x == folds)    #vector with the indexes of one of the 10 parts of the data, that will be used as the training set
            idx.train = (1:nrow(X))[-idx.test] # the other 9 parts of the data will create the training set
            mod.rf = randomForest(x=X[idx.train,],y=Y[idx.train], do.trace=TRUE) #creation of a classifier with the random forest method
            #do.trace = TRUE for more verbose output 
            
            Pred.rf = predict(mod.rf,X[idx.test,])   #applying the RF classifier on the data. output = vector with predicted classes of the data
            
            flush.console()
          
        
            save(Pred.rf, file=paste0("classification_result.RData")) ##saving intermediate result
            return(Pred.rf) # list with prediction for each fold of 10-fold cross-validation 
          
          
          }
  )
  
  return(predictions)

}


#http://sablab.net/scripts/sortDataFrame.r 
sortDataFrame <- function(x, key, ...) {
  if (missing(key)) {
    rn <- rownames(x)
    if (all(rn %in% 1:nrow(x))) rn <- as.numeric(rn)
    x[order(rn, ...), , drop=FALSE]
  } else {
    x[do.call("order", c(x[key], ...)), , drop=FALSE]
  }
}






####Analysing classification results 

# res1=output_classification$x1
# res2=output_classification$x2
# res3=output_classification$x3


AUC_analysis = function(output_classification, input_expression_matrix){
  
  ## selecting top 100 items with the highest AUC
  selection=tail(sortDataFrame(output_classification,'AUC'),100)$id   


  #selecting these top 100 genes in the original expression matrix
  top_100_AUC=input_expression_matrix[which(rownames(input_expression_matrix)%in%selection),]
 
  return(top_100_AUC)
  
  
}


FDR_analysis = function(output_classification, input_expression_matrix){
  
  ## selecting top 100 items with the lowest FDR
  selection=head(sortDataFrame(output_classification,'FDR'),100)$id   
  
  
  #selecting these top 100 genes in the original expression matrix
  low_100_FDR=input_expression_matrix[which(rownames(input_expression_matrix)%in%selection),]
  
  return(low_100_FDR)
  
  
}



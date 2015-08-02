#' Two Sample Test with Distribution-Valued Data
#' @param cases name of case group data (matrix)
#' 
#' @param control names of control group data (matrix)
#' 
#' @param rm.mean standarize each rows of cases and control to mean=0.
#' 
#' @param rm.var standarize each rows of cases and control to var=1.
#' 
#' @param paranum the number of quatile discretization + 1. Default is discretized by 1 \%.
#' 
#' @param q power of Wasserstein metric. Default is q = 2.
#' 
#' @param bsn the number of resampling. Default is bsn = 5000.
#' 
#' @param seed seed for random number generator.
#' 
#' @details this function is designed for two sample test based on Wasserstein metric. The function computes the the p-values based Wasserstein metric and resampling method. If rm.mean=F and rm.var=F, then statistical test is performed only based on more than 3rd order moments. 
#' 
#' @examples 
#' nrep <-12
#' cases <- Map(rbeta,rep(30,nrep),rep(1,nrep),rep(5,nrep)); cases <- do.call("rbind",cases)
#' control <- Map(rbeta,rep(30,nrep),rep(1,nrep),rep(5,nrep)); control <- do.call("rbind",control)
#' d3m(cases,control,paranum = 101, q = 2, bsn = 1000,seed = 100)
#' 
#' @return  pval p-value.
#' @return test.stat test statistic.
#' @return cases case group data used in the statistical test.
#' @return control control group data used in the statistical test.
#' @author Yusuke Matsui & Teppei Shimamura
#' 
#' @export
#' 

d3m <- function(cases, control, rm.mean = FALSE, rm.var = FALSE, paranum = 101, q = 2, bsn = 5000, seed = 100){
   
  #sourceCpp("./src/funcs.cpp")
  #library(Rcpp)
  
  if(rm.mean & rm.var){
    cases <- scale(x = cases,center = T,scale = T)
    control <- scale(x = control,center = T,scale = T)
  }else if(rm.mean & !rm.var){
    control <- scale(x = control,center = T,scale = F)
  }else if(!rm.mean & rm.var){
    cases <- scale(x = cases,center = F,scale = T)
  }

  d <- wasserCpp_mat(cases,control,paranum = paranum, q = q)
  #d <- .Call("wasserCpp_mat",cases,control,paranum,q,PACKAGE = "D3M")
  #d <- .Call("wasserCpp_mat",PACKAGE = "D3M")
  
  set.seed(seed)
  
  res <- permCpp(casesMat = cases,controlMat = control,bsn = bsn,qn = paranum, d = d)  
  #res <- .Call("permCpp",cases,control,bsn,paranum,d,PACKAGE = "D3M")  
  
  return(list(pval = res[[1]], test.stat = res[[2]], cases = cases, control = control))
  #pval <- res[[1]]
  
  #pval
}


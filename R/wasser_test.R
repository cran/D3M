#' Statistical Test with Wasserstein Metric
#' @param cases name of case group data (matrix sample * feature)
#' 
#' @param control names of control group data (matrix sample * feature)
#' 
#' @param test.stat test statistic
#' 
#' @param paranum the number of quatile discretization + 1. Default is discretized by 1 \%.
#' 
#' @param bsn the number of resampling. Default is bsn = 5000.
#' 
#' @param q power of Wasserstein metric. Default is q = 2.
#' 
#' @examples
#' 
#' nrep <- 12
#' cases <- Map(rbeta,rep(100,nrep),rep(1,nrep),rep(5,nrep))
#' cases <- do.call("rbind",cases)
#' control <- Map(rbeta,rep(100,nrep),rep(1,nrep),rep(5,nrep))
#' control <- do.call("rbind",control)
#' d <- wasserMetric(cases,control)
#' testRes <- wasser.test(cases = cases,control = control,test.stat = d)

#' @author Yusuke Matsui & Teppei Shimamura
#' 
#' @return  list of p-value and test statistics.
#' 
#' @export

wasser.test <- function(cases, control, test.stat, paranum = 101, bsn = 5000, q = 2){
  #library(D3M)
  
  #library(Rcpp)
  
  res <- permCpp(casesMat = cases, controlMat = control, bsn = bsn, qn = paranum, d = test.stat, q = q)

  o <- list(pval = res[[1]], tets.stat = test.stat)
  
  return(o)
}

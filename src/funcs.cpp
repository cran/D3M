//' @useDynLib D3M
//' @importFrom Rcpp sourceCpp


#include <Rcpp.h>
#include <stdlib.h>
#include <math.h>

// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
using namespace std;

//[[Rcpp::export]]
NumericVector quantileCpp(NumericVector x, NumericVector probs){
  
  NumericVector xs(x.size());
  
  xs = clone(x);
  
  std::sort(xs.begin(), xs.end());
  
  
  double bin = (double)(1.0 / xs.size());
  
  double prev = 0.0;
  
  int flag = probs.size();
  
  int cnt = 0;
  
  NumericVector result(probs.size());
  
  for(int i = 0; i < result.size(); i++){
    
    result[i] = 0;
    
  }
  
  
  double xs_min = 0.0;
  
  for(int i = 0; i < xs.size(); i++){
    
    if((prev + bin) >= probs[cnt]){
      
      if(i == 0){
        
        xs_min = min(xs);
        
        result[cnt] = xs_min  + ((probs[cnt] - prev) / bin) * (xs[i] - xs_min);
        
        flag--;
        
        cnt++;
        
      }else{
        
        result[cnt] = xs[ i - 1 ] + ((probs[cnt] - prev) / bin) * (xs[i] - xs[i - 1]);

        flag--;
        
        cnt++;
        
      }
      
      if(flag > 0){
        
        prev += bin;
        
        continue;
        
      }else{
        
        return result;
      }
      
    }else{
      
      prev +=  bin;
      
    }
  }  
  return result;
}


// [[Rcpp::export]]
double wasserCpp(NumericVector x, NumericVector y, int paranum = 101, int q = 2){
  //int len = x.size();   
  double d = 0;
  
  NumericVector pseq(paranum);
  NumericVector xq(paranum), yq(paranum);
  
  for(int i = 0; i < paranum; i++){
    pseq[i] = (double)i / (paranum - 1);
    //printf("%.3f\t",pseq[i]);
  }
  //printf("\n");
  
  //x = x[!is_na(x)];
  //y = y[!is_na(y)];
  
  xq = quantileCpp(x, pseq);
  yq = quantileCpp(y, pseq);
  
  for(int i = 0; i < paranum; i++){
    //printf("%.3f\t",xq[i]);
    //printf("%3.f\t",yq[i]);
    d += pow(xq[i] - yq[i],q);
  }
  return d;
}

// [[Rcpp::export]]
NumericVector wasserCpp_mat(NumericMatrix xMat, NumericMatrix yMat, int paranum = 101, int q = 2){
  
  int nr_xMat = xMat.nrow();
  int nc_xMat = xMat.ncol();
  int nr_yMat = yMat.nrow();
  int nc_yMat = yMat.ncol();
  //printf("x %d %d y %d %d\n",nr_xMat,nc_xMat,nr_yMat,nc_yMat);
  
  if(nr_xMat != nr_yMat){
    
    //printf("Error: size of matrix is different");
    
    double d = NA_REAL;
    
    return d;
  }
  
  int npos = nr_xMat;
  
  //printf("npos %d\n",npos);
  
  NumericVector x(nc_xMat);
  NumericVector y(nc_yMat);
  
  //printf("hoge1\n");
  
  NumericVector d_vec(npos);
  
  
  for(int pos = 0; pos < npos; pos++){
    
    //printf("hoge2\n");
    
    for(int j = 0; j < nc_xMat; j++){
      
      //printf("hoge3\n");
      
      x[j] = xMat(pos,j);
      
      //printf("x%d %.3f\t",j,x[j]);
      
    }
    
    for(int j = 0; j < nc_yMat; j++){

      y[j] = yMat(pos,j);
      
      //printf("y%d %.3f\t",j,y[j]);
    }
    
    //int len = x.size();
    double d = 0;
    
    NumericVector pseq(paranum);
    NumericVector xq(paranum), yq(paranum);
    
    for(int i = 0; i < paranum; i++){
      
      pseq[i] = (double)i / (paranum - 1);
      
    }
    
    //x = x[!is_na(x)];
    //y = y[!is_na(y)];
    
    xq = quantileCpp(x, pseq);
    yq = quantileCpp(y, pseq);
    
    for(int i = 0; i < paranum; i++){
      //printf("%.3f\t",xq[i]);
      d += pow(xq[i] - yq[i],q);
      
    }
    
    d_vec[pos] = d;
    
  }
  
  return d_vec;
}


inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// [[Rcpp::export]]
Rcpp::NumericVector randomShuffle(Rcpp::NumericVector a) {
  
  // clone a into b to leave a alone
  Rcpp::NumericVector b = Rcpp::clone(a);
  
  std::random_shuffle(b.begin(), b.end(), randWrapper);
  
  return b;
}

// [[Rcpp::export]]
List permCpp(NumericMatrix casesMat, NumericMatrix controlMat, NumericMatrix shuffleID, NumericVector d, int bsn = 10000, int qn = 101, int q = 2){ 
  
  int nr_casesMat = casesMat.nrow();
  int nc_casesMat = casesMat.ncol();
  //int nr_controlMat = controlMat.nrow();
  int nc_controlMat = controlMat.ncol();
  
  //printf("nr_caseMat %d\nnc_caseMat %d\nnr_controlMat %d\nnc_controlMat %d\n",nr_casesMat,nc_casesMat,nr_controlMat,nc_controlMat);
  
  int npos = nr_casesMat;
  NumericVector pval_vec(npos);
  
  NumericVector cases(nc_casesMat);
  NumericVector control(nc_controlMat);
  
  for(int pos = 0; pos < npos; pos++){
    
    for(int j = 0; j < nc_casesMat; j++){
      cases[j] = casesMat(pos,j);
      
    }
    
    for(int j = 0; j < nc_controlMat; j++){
      control[j] = controlMat(pos,j);
      // printf("i,j %d %d: %f\n",pos,j,control[j]);
    }
    
    //cases = cases[!is_na(cases)];
    //control = control[!is_na(control)];
    
    int ncases = cases.size();
    int ncontrol = control.size();
    int nsample = ncases + ncontrol;
    
    //printf("sample size:\ncase %d control %d total %d\n",ncases,ncontrol,nsample);
    
    NumericVector data(nsample);
    
    int cnt = 0;
    
    for(int i = 0; i < ncases; i++){
      
      data[i] = cases[i];
      
      cnt++;
    }
    
    for(int i = 0; i < ncontrol; i++){
      data[cnt] = control[i];
      
      cnt++;
    }
    
    
    NumericVector index(nsample);
    
    for(int i = 0; i < nsample; i++){
      
      index[i] = i;
    }
    
    NumericVector bootd(bsn);
    
    //NumericVector index2;
    
    NumericVector casedata(ncases);
    
    NumericVector controldata(ncontrol);
    
    double a;
    
    int nonaN = 0;
    
    for(int i = 0; i < bsn; i++){
      
      //index2 = randomShuffle(index);
      
      int k;
      
      cnt = 0;
      
      for(int j = 0; j < ncases; j++){
        
        //k = index2[cnt];
        
        k = shuffleID(cnt,i);
        
        casedata[j] = data[k];
        
        //printf("casedata[%d] = %.3f = data[%d]\n",j,casedata[j],k);
        
        cnt++;
        
      }
      
      for(int j = 0; j < ncontrol; j++){
        
        //k = index2[cnt];
        
        k = shuffleID(cnt,i);
        
        controldata[j] = data[k];
      
        //printf("controldata[%d] = %.3f = data[%d]\n",j,controldata[j],k);
         
        cnt++;
        
      }
      
      a = wasserCpp(casedata, controldata, qn, q);
      
      //printf("a[%d] = %.3f\n",i,a);
      
      if(!R_FINITE(a)) {
        
        bootd[i] = NAN;
        
      }else{
        
        bootd[i] = a;
        
        nonaN++;
      }
      
    }    
    
    NumericVector newbootd(nonaN);
    
    cnt = 0;
    
    for(int l = 0; l < bootd.size(); l++){
      
      if(R_FINITE(bootd[l])) {
        
        newbootd[cnt] = bootd[l];
        
        cnt++;
      }      
      
    }
    
    NumericVector temp = newbootd[newbootd >= d[pos]];
    
    double pval = (double)temp.size() / (double)newbootd.size();
    
    //printf("pval=temp.size()/newbootd.size()=%.3f/%.3f= %.3f\n",(double)temp.size(), (double)newbootd.size(),pval);
    
    if(pval < (double) (1 / bsn)){
      
      //int threshold_ind = (int)bsn * 0.995;
      
      //sort(newbootd.begin(),newbootd.end());
      
      //double threshold = newbootd[threshold_ind];
      
      NumericVector thresholdv = quantileCpp(newbootd, 0.995);
      
      double threshold = thresholdv[0];
      
      NumericVector ptemp = newbootd[newbootd > threshold];
      
      for(int j = 0; j < ptemp.size(); j++){
        
        ptemp[j] = ptemp[j] - threshold;
        
      }
      
      
      double lambda = 0;
      
      for(int j =0; j < ptemp.size(); j++){
        
        lambda += ptemp[j];
        
      }
      
      lambda = 1 / (lambda / ptemp.size());
      
      //printf("Estimated lambda: %f\n",lambda);
      
      
      //double param = -1 * lambda * (d[pos] - threshold);
      double param = -1 * lambda * d[pos];
      
      double param2 = -1 * lambda * threshold;
      
      double pval_threshold = exp(param2);
      
      double r = pval_threshold / 0.005;
      
      pval = exp(param) / r;
      
      //printf("semi-parametric p-value (lambda = %.3f) %.3e\n",lambda,pval);
      
    }
    
    
    pval_vec[pos] = pval;
        
  }
  
  return List::create(pval_vec,d,casesMat, controlMat);
    
}



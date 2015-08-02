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
  
  /*
  for(int i = 0; i < xs.size(); i++){
  printf("xs[%d] = %.3f\n",i,xs[i]);
  }
  */
  
  double xs_min = 0.0;
  
  for(int i = 0; i < xs.size(); i++){
    
    if((prev + bin) >= probs[cnt]){
      
      if(i == 0){
        
        xs_min = min(xs);
        
        result[cnt] = xs_min  + ((probs[cnt] - prev) / bin) * (xs[i] - xs_min);
        
        /*
        printf("xs_min = %.3f\n",xs_min);
        
        printf("(probs[%d] - prev) / bin = %.2f - %.2f / %.2f = %.2f\n",cnt, probs[cnt], prev, bin,(probs[cnt] - prev) / bin);
        
        printf("result[%d] = %.3f\n",cnt,result[cnt]);
        */
        
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
  int len = x.size();   
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
  
  //xq = percentile_rcpp(x, pseq);
  //yq = percentile_rcpp(y, pseq);
  
  xq = quantileCpp(x, pseq);
  yq = quantileCpp(y, pseq);
  
  for(int i = 0; i < len; i++){
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
      
      //printf("hoge4\n");
      
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
    
    //xq = percentile_rcpp(x, pseq);    
    //yq = percentile_rcpp(y, pseq);
    
    xq = quantileCpp(x, pseq);
    yq = quantileCpp(y, pseq);
    
    for(int i = 0; i < paranum; i++){
      
      d += pow(xq[i] - yq[i],q);
      
    }
    
    d_vec[pos] = d;
    
  }
  
  return d_vec;
}


// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// [[Rcpp::export]]
Rcpp::NumericVector randomShuffle(Rcpp::NumericVector a) {
  
  // clone a into b to leave a alone
  Rcpp::NumericVector b = Rcpp::clone(a);
  
  std::random_shuffle(b.begin(), b.end(), randWrapper);
  
  return b;
}

// [[Rcpp::export]]
List permCpp(NumericMatrix casesMat, NumericMatrix controlMat, NumericVector d, int bsn = 10000, int qn = 101, int q=2){ 
  
  
  
  //List permCpp(NumericVector cases, NumericVector control, double d, NumericMatrix shuffleID, int bsn = 10000, int qn = 101){  
  
  int nr_casesMat = casesMat.nrow();
  int nc_casesMat = casesMat.ncol();
  //int nr_controlMat = controlMat.nrow();
  int nc_controlMat = controlMat.ncol();
  
  //printf("nr_caseMat %d\nnc_caseMat %d\nnr_controlMat %d\nnc_controlMat %d\n",nr_casesMat,nc_casesMat,nr_controlMat,nc_controlMat);
  
  /*
  if(nr_casesMat != nr_controlMat){
  //printf("Error: size of matrix is different");
  int pval = NA_INTEGER;
  int bootd = NA_INTEGER;
  return List::create(pval,bootd);
  }
  */
  int npos = nr_casesMat;
  NumericVector pval_vec(npos);
  //List bootd_list(npos);
  //List res(npos);
  
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
    
    NumericVector index2;
    
    NumericVector casedata(ncases);
    
    NumericVector controldata(ncontrol);
    
    double a;
    
    int nonaN = 0;
    
    for(int i = 0; i < bsn; i++){
      
      
      //srand(i);
      
      index2 = randomShuffle(index);
      
      int k;
      
      cnt = 0;
      
      for(int j = 0; j < ncases; j++){
        
        k = index2[cnt];
        
        //k = shuffleID(i,cnt);
        
        casedata[j] = data[k];
        
        //printf("casedata:%.3f\t",casedata[j]);
        
        cnt++;
        
      }
      
      for(int j = 0; j < ncontrol; j++){
        
        k = index2[cnt];
        
        //k = shuffleID(i,cnt);
        
        controldata[j] = data[k];
        
        //printf("controldata:%.3f",controldata[j]);
        
        cnt++;
        
      }
      
      a = wasserCpp(casedata, controldata, qn, q);
      
      
      if(isnan(a) || isinf(a)){
        
        bootd[i] = NAN;
        
      }else{
        
        bootd[i] = a;
        
        nonaN++;
      }
      //printf("bootd[%d] = %.3f\n",i,bootd[i]);
      //bootd[i] = wasserCpp(casedata, controldata, qn, q);
      
    }
    //printf("nonaN = %d\n",nonaN);
    
    
    /*
    for(int l = 0; l < bootd.size(); l++){
    
    printf("l:%d\n",l);
    if((!isnan(bootd[l])) || (!isinf(bootd[l]))){
    
    nonaN++;  
    
    }      
    }
    */
    
    
    NumericVector newbootd(nonaN);
    
    cnt = 0;
    
    for(int l = 0; l < bootd.size(); l++){
      
      if(!isnan(bootd[l]) && !isinf(bootd[l])){
        
        newbootd[cnt] = bootd[l];
        
        cnt++;
      }      
      
    }
    
    
    /*
    bool temp[newbootd.size()];
    int num = 0;
    for(int l = 0; l < newbootd.size(); l++){
    //temp[l] = (bool)((bootd[l] >= d[pos]) && (!isnan(bootd[l]) && (!isinf(bootd[l]))));
    temp[l] = (bool)((newbootd[l] >= d[pos]));
    //temp[l] = (bool)(isnan(bootd[l]));
    
    num += temp[l];
    
    if(temp[l]){
    printf("l: %d d[%d]: %.3f bootd: %.3f temp[%d]: %s\n",l, pos,d[pos],newbootd[l],l,temp[l]?"true":"false"); 
    }
    
    }
    
    
    NumericVector temp2(num);
    
    cnt = 0;
    
    for(int l = 0; l < newbootd.size(); l++){
    
    if(temp[l]){
    
    temp2[cnt] = newbootd[l];
    
    //printf("temp2[%d]=%.3f\n",cnt,temp2[cnt]);
    
    cnt++;
    }
    }
    */
    
    
    //printf("number of True: %d",num);
    
    
    
    
    
    NumericVector temp = newbootd[newbootd >= d[pos]];
    
    
    
    double pval = (double)temp.size() / (double)newbootd.size();
    
    //printf("bootstrap p-value (%d times) %.3f\n",bsn, pval);
    
    
    //pval = 0.000001;
    
    //printf("pval: %.5f 1 / bsn: %.5f\n",pval, 1/bsn);
    
    
    if(pval < (double) (1 / bsn)){
      
      int threshold_ind = (int)bsn * 0.995;
      
      sort(newbootd.begin(),newbootd.end());
      
      double threshold = newbootd[threshold_ind];    
      
      //NumericVector ptemp(estn); 
      
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
      
      
      double param = -1 * lambda * (d[pos] - threshold);
      
      pval = exp(param);
      
      //printf("semi-parametric p-value (lambda = %.3f) %.3e\n",lambda,pval);
      
      //return ptemp;
      //pval.names() = CharacterVector::create("pval");
      //bootd.names() = CharacterVector::create("bootd");
      //return List::create(pval,bootd);
      
    }
    
    
    pval_vec[pos] = pval;
    
    //bootd_list[pos] = newbootd;
    
  }
  
  
  //return List::create(pval_vec,bootd_list,casesMat, controlMat);
  return List::create(pval_vec,d,casesMat, controlMat);
  
  //return 0;
  
}



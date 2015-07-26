//' @useDynLib D3M
//' @importFrom Rcpp evalCpp


#include <Rcpp.h>
#include <stdlib.h>
#include <math.h>

// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector sort_rcpp(Rcpp::NumericVector& x)
{
  std::vector<double> tmp = Rcpp::as<std::vector<double> > (x);    // or NumericVector tmp = clone(x);

  std::sort(tmp.begin(), tmp.end());
  
  return wrap(tmp);
}

// [[Rcpp::export]]
NumericVector percentile_rcpp(NumericVector& x, NumericVector& percentile)
{
  NumericVector tmp_sort = sort_rcpp(x);
  
  int size_per = percentile.size();
  
  NumericVector percentile_vec = no_init(size_per);
  
  for (int ii = 0; ii < size_per; ii++)
  
  {
  
	  double size_per = tmp_sort.size() * percentile[ii];
      double size_per_round;
    
	  if (size_per < 1.0)
    
	  {
      size_per_round = 1.0;
    }
    else
    {
      size_per_round = round(size_per);
    }
    percentile_vec[ii] = tmp_sort[size_per_round-1];  // For extreme case such as size_per_round == tmp_sort.size() to avoid overflow
  }
  return percentile_vec;
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
  
  xq = percentile_rcpp(x, pseq);
  yq = percentile_rcpp(y, pseq);
    
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
    
    xq = percentile_rcpp(x, pseq);    
    yq = percentile_rcpp(y, pseq);
    
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
  int nr_controlMat = controlMat.nrow();
  int nc_controlMat = controlMat.ncol();
  
  //printf("nr_caseMat %d\nnc_caseMat %d\nnr_controlMat %d\nnc_controlMat %d\n",nr_casesMat,nc_casesMat,nr_controlMat,nc_controlMat);
  
  if(nr_casesMat != nr_controlMat){
    //printf("Error: size of matrix is different");
    int pval = NA_INTEGER;
    int bootd = NA_INTEGER;
    return List::create(pval,bootd);
  }
  
  int npos = nr_casesMat;
  NumericVector pval_vec(npos);
  List bootd_list(npos);
  List res(npos);
  
  NumericVector cases(nc_casesMat);
  NumericVector control(nc_controlMat);
  
  for(int pos = 0; pos < npos; pos++){
    
    for(int j = 0; j < nc_casesMat; j++){
      cases[j] = casesMat(pos,j);
      //printf("%f\n",cases[j]);
    }
    
    for(int j = 0; j < nc_controlMat; j++){
      control[j] = controlMat(pos,j);
    }
    
    cases = cases[!is_na(cases)];

    control = control[!is_na(control)];
    
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
    
    for(int i = 0; i < bsn; i++){
      
      //srand(i);
      
      index2 = randomShuffle(index);
      
      int k;
      cnt = 0;
      for(int j = 0; j < ncases; j++){
        
        k = index2[cnt];
        
        //k = shuffleID(i,cnt);
        
        casedata[j] = data[k];
        
        cnt++;
        
      }
      
      for(int j = 0; j < ncontrol; j++){
        
        k = index2[cnt];
        
        //k = shuffleID(i,cnt);
        
        controldata[j] = data[k];
        
        cnt++;
        
      }
      
      bootd[i] = wasserCpp(casedata, controldata, qn, q);
    }
    
    NumericVector temp = bootd[bootd >= d[pos]];
    
    
    double pval = (double)temp.size() / (double)bootd.size();
    
    //printf("bootstrap p-value (%d times) %.3f\n",bsn, pval);
    
    if(pval == 0){
      
      int threshold_ind = (int)bsn * 0.995;
      
      sort(bootd.begin(),bootd.end());
      
      double threshold = bootd[threshold_ind];
      
      NumericVector ptemp = bootd[bootd > threshold];
      
      for(int j = 0; j < ptemp.size();j++){
        
        ptemp[j] = ptemp[j] - threshold;
        
      }
      
      double lambda = 0;
      
      for(int j =0; j <= ptemp.size(); j++){
        
        lambda += ptemp[j];
        
      }
      
      lambda = 1 / (lambda / ptemp.size());
      
      //printf("Estimated lambda: %f",lambda);
      
      double param = -1 * lambda * (d[pos] - threshold);
      
      pval = exp(param);
      
      //printf("semi-parametric p-value (lambda = %.3f) %.3e\n",lambda,pval);
      
      //return ptemp;
      //pval.names() = CharacterVector::create("pval");
      //bootd.names() = CharacterVector::create("bootd");
      //return List::create(pval,bootd);
    }
    
    pval_vec[pos] = pval;
    
    bootd_list[pos] = bootd;
    
  }
  
  
  //return List::create(pval_vec,bootd_list,casesMat, controlMat);
  return List::create(pval_vec,d,casesMat, controlMat);
  
}



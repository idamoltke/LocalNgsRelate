
#include <cstdio>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <cassert>
#include <cstring>
#include <cfloat>
#include <vector>
#include "analysisfunctions.h"
#include "bfgs.h"


//####################################
// Global/internal constants and definitions
//####################################

extern double alim[];
extern int runoldrelateV; // Remove when done?
double klim[2]= {DBL_EPSILON,1-DBL_EPSILON};//{0,1};//
bool verbose = false;
bool testrandomstart = true;

// For indexing into vector with transitions
/*
  0  , 1 , 2 , 3 , 4 , 5 ,6  , 7 , 8
  T00,T01,T02,T10,T11,T12,T20,T21,T22
*/
int T00 =0;
int T01 =1;
int T02 =2;
int T10 =3;
int T11 =4;
int T12 =5;
int T20 =6;
int T21 =7;
int T22 =8;


//####################################
// Internal help structs
//####################################

// Struct for 
typedef struct toOptim2_t{
  const genome &g;
  const std::vector<perChr> &pc;
  const para &p;
  double *tsk;
  toOptim2_t(const genome &gg,const std::vector<perChr> &pcc,const para &pp) : g(gg),pc(pcc), p(pp) {};
}toOptim2;


//####################################
// Internal help functions
//####################################


// Function getting max number in a vector of length 3
// ------------------------------------------------------------
double max(double a[3]){
  double m=a[0];
  for(int i=1;i<2;i++)
    if(a[i]>m)
      m=a[i];
  return m;
}

// Function getting max number out of three doubles
// ------------------------------------------------------------
double max(double a,double b,double c){
  return std::max(a,std::max(b,c));
}

// Function getting the index the max of 3 doubles
// ------------------------------------------------------------

int whichmax(double a,double b,double c){
  double aa[] = {a,b,c};
  int r =0;
  for(int i=1;i<3;i++)
    if(aa[i]>aa[r])
      r=i;
  return r;
}


// Function sampling a random real mumber between to numbers
// ------------------------------------------------------------
double myRand(double low,double up){
  assert(up>low);
  double r = (drand48());//r=drand48();
  double range = up-low;
  return r*range+low;
}

double myRand2(double low,double up){
  assert(up>low);
  double r = (double)rand()/(double)(RAND_MAX);
  double range = up-low;
  return r*range+low;
}


// Functions for adding while protecting for underflow
// ------------------------------------------------------------
// NB: currently also in NGSrelate.cpp

double addProtect2(double a,double b){
  //function does: log(exp(a)+exp(b)) while protecting for underflow
  if(std::isinf(a)&&std::isinf(b))
    return log(0);
  if(std::isinf(a))
    return b;
  if(std::isinf(b))
    return a;
  double maxVal;// = std::max(a,b));
  if(a>b)
    maxVal=a;

  else
    maxVal=b;
  double sumVal = exp(a-maxVal)+exp(b-maxVal);
  return log(sumVal) + maxVal;
}

double addProtect3(double a,double b, double c){
  //function does: log(exp(a)+exp(b)+exp(c)) while protecting for underflow
  double maxVal;// = std::max(a,std::max(b,c));
  if(a>b&&a>c)
    maxVal=a;
  else if(b>c)
    maxVal=b;
  else
    maxVal=c;
  double sumVal = exp(a-maxVal)+exp(b-maxVal)+exp(c-maxVal);
  return log(sumVal) + maxVal;
}

// Function for fixing underflow... (setting all numbers in an array a below "tole" to 0) 
// ---------------------------------------------------------------------------------------

void fixunderflow(double *a,int l, double tole){
  for(int i=0;i<l;i++)
    if(a[i]<tole)
      a[i]=0.0;

}

// Function for getting distance between positions of neighbouring loci
// -----------------------------------------------------------------------

double *diffPos(const int *pos,int len){
  double *ret = new double [len+1];
  ret[0] = .0;
  for(int i=0;i<len-1;i++)
    ret[i+1] = (pos[i+1]-pos[i])/1e6;
  ret[len] = .0;
  return ret;

}

// Function for transposing a in_x x in_y matrix called in
// ------------------------------------------------------------

template <typename T>
void transpose2(T***in,size_t in_x,size_t in_y){
  int nInd = in_y;
  int nSites = in_x;

  T **tmp = new T*[nInd];
  for(size_t i=0;i<in_y;i++)
    tmp[i] = new T [3*nSites];
  
  for(int i=0;i<nInd;i++)
    for(int j=0;j<nSites;j++)
      for(int o=0;o<3;o++)
	tmp[i][j*3+o] = (*in)[j][i*3+o];
  
  (*in) = tmp;
}

// Function for writing a set of parameter value to a file
// ------------------------------------------------------------

void printPars(FILE *fp,para p){
  fprintf(fp,"\tpair=(%d,%d) alpha=%f k0=%f k1=%f k2=%f\n",p.pair[0],p.pair[1],p.a,p.k0,p.k1,p.k2);
}

//void printStuff(const std::vector<perChr>& pc ){
//  for(uint i=0;i<pc.size();i++)
//    fprintf(stderr,"i=%d chr=%s nSites=%d (%d,%d) (%f,%f)\n",i,pc[i].name,pc[i].nSites,pc[i].pos[0],pc[i].pos[pc[i].nSites-1],min<double>(pc[i].freq,pc[i].nSites),max<double>(pc[i].freq,pc[i].nSites));

//}



//####################################
// Functions needed for inference
//####################################

// Function for calculating a (from Albrechtsen et al. 2009)                                                                      
// ------------------------------------------------------------

double calculateA(double k0,double k1, double k2,double phi){
  double ma,mb,xa,xb,m,a,sq,pw;
  if(k2==0){
    ma = 1-log(k1)/log(2);
    mb = 0;
  }
  else{
    pw = pow(k1+2*k2,2);
    sq = sqrt(pw-4*k2);
    xa = (k1+2*k2+sq)/2.0;
    xb = k2/xa;
    ma = 1-log(xa)/log(2);
    mb = 1-log(xb)/log(2);
  }
  m = ma+mb;
  a = -m*log(1-phi);
  
  if(std::isnan(a)){//there is a bug in the bfgs algortihm, that makes the optim algo try outside of parameter space
    return DBL_MAX;
  }

  else if(std::isinf(a)){
    if(verbose){
      fprintf(stderr,"std::isinf in calculate.a -> k0=%f,k1=%f,k2=%f,phi=%f\n",k0,k1,k2,phi);
      if(k2==0)
	fprintf(stderr,"k2=0,\t");
      fprintf(stderr,"ma=%f,mb=%f,m=%f\n",ma,mb,m);
    }
  }
  return a;
}

// Function for calculating emissions
// --------------------------------------     
/* Input:   p1 and p2 are the indeces of the pair of individuals being analysed  
            pc is struct with all relavant info  (nSites, log freqs and log gls for all sites and all inds)
            res is pointer to a matrix 3xnSites where the emissions will be stored (should be allocated beforehand)
   Output:  ret is filled with log emissions so ret[j][i] contains 
            log(P(Data obs in site i|X=j) = sum over G_real P(G_real|X=j)P(Data obs|G_real) 
	    where the latter are the gls of the two individuals multiplied 
*/

void calclogemis(int p1,int p2,const perChr &pc,double **ret){

  // Extracting an array of 3*nSites of log GLs in log for each of the two individuals
  double *loggl1 = pc.loggl[p1];
  double *loggl2 = pc.loggl[p2];
  
  // Extracting arrays of the frequencies in log
  double *logf1,*logf2;
  logf1 = pc.logqerf;//freqA
  logf2 = pc.logfreq;//freqa

  // Calculating log emission for each site for each of the 3 possible IBD states by summing over all possible genotypes 
  if(!runoldrelateV){ // Remove when done?
  //AA,AA
    for(int s=0;s<pc.nSites;s++){
    double c1 = loggl1[3*s]+loggl2[3*s];
    ret[0][s] = 4*logf1[s]+c1;
    ret[1][s] = 3*logf1[s]+c1;
    ret[2][s] = 2*logf1[s]+c1;
  }
  //AA,Aa
  for(int s=0;s<pc.nSites;s++){
    double c1 = loggl1[3*s]+loggl2[3*s+1];
    ret[0][s] = addProtect2(ret[0][s],M_LN2+3*logf1[s]+logf2[s]+c1);
    ret[1][s] = addProtect2(ret[1][s],2*logf1[s]+logf2[s]+c1); 
    //no third
  }
  //AA,aa
  for(int s=0;s<pc.nSites;s++) {
    double c1 = loggl1[3*s]+loggl2[3*s+2];
    ret[0][s] = addProtect2(ret[0][s],2*(logf1[s]+logf2[s])+c1);
    //no sec
    //no third
  }
  //Aa,AA
  for(int s=0;s<pc.nSites;s++) {
    double c1 = loggl1[3*s+1]+loggl2[3*s];
    ret[0][s] = addProtect2(ret[0][s],M_LN2+3*logf1[s]+logf2[s]+c1);
    ret[1][s] = addProtect2(ret[1][s],2*logf1[s]+logf2[s]+c1);
    //no third
  }
  //Aa,Aa
  for(int s=0;s<pc.nSites;s++) {
    double c1 = loggl1[3*s+1]+loggl2[3*s+1];
    ret[0][s] = addProtect2(ret[0][s],2*M_LN2+2*(logf1[s]+logf2[s])+c1);
    ret[1][s] = addProtect2(ret[1][s],logf1[s]+logf2[s]+c1);
    ret[2][s] = addProtect2(ret[2][s],M_LN2+logf1[s]+logf2[s]+c1);
  }
  //Aa,aa
  for(int s=0;s<pc.nSites;s++) {
     double c1 = loggl1[3*s+1]+loggl2[3*s+2];
    ret[0][s] = addProtect2(ret[0][s],M_LN2+logf1[s]+3*logf2[s]+c1);
    ret[1][s] = addProtect2(ret[1][s],logf1[s]+2*logf2[s]+c1);
    //no third
  }
  //aa,AA
  for(int s=0;s<pc.nSites;s++) {
    double c1 = loggl1[3*s+2]+loggl2[3*s];
    ret[0][s] = addProtect2(ret[0][s],2*(logf1[s]+logf2[s])+c1);
    //no sec
    //no third
  }
  //aa,Aa
  for(int s=0;s<pc.nSites;s++) {
    double c1 = loggl1[3*s+2]+loggl2[3*s+1];
    ret[0][s] = addProtect2(ret[0][s],M_LN2+logf1[s]+3*logf2[s]+c1);
    ret[1][s] = addProtect2(ret[1][s],logf1[s]+2*logf2[s]+c1);
    //no third
  }
  for(int s=0;s<pc.nSites;s++) {
    double c1 = loggl1[3*s+2]+loggl2[3*s+2];
    ret[0][s] = addProtect2(ret[0][s],4*logf2[s]+c1);
    ret[1][s] = addProtect2(ret[1][s],3*logf2[s]+c1);
    ret[2][s] = addProtect2(ret[2][s],2*logf2[s]+c1);
  }

  
  }else{ // Remove this and everything below when done?

  //AA,AA (same)    
  for(int s=0;s<pc.nSites;s++){
    double c1 = loggl1[3*s]+loggl2[3*s];
    ret[0][s] = 4*logf1[s]+c1;
    ret[1][s] = 3*logf1[s]+c1;
    ret[2][s] = 2*logf1[s]+c1;
  }
  //AA,Aa
  for(int s=0;s<pc.nSites;s++){
    double c1 = loggl1[3*s]+loggl2[3*s+1];
    ret[0][s] = addProtect2(ret[0][s],M_LN2+M_LN2+3*logf1[s]+logf2[s]+c1); // added an extra LN2
    ret[1][s] = addProtect2(ret[1][s],M_LN2+2*logf1[s]+logf2[s]+c1);       // added an extra LN2
    //no third
  }
  //AA,aa
  for(int s=0;s<pc.nSites;s++) {
    double c1 = loggl1[3*s]+loggl2[3*s+2];
    ret[0][s] = addProtect2(ret[0][s],M_LN2+2*(logf1[s]+logf2[s])+c1);     // added an extra LN2
    //no sec
    //no third
  }
  //Aa,AA
  for(int s=0;s<pc.nSites;s++) {
    double c1 = loggl1[3*s+1]+loggl2[3*s];
    ret[0][s] = addProtect2(ret[0][s],M_LN2+M_LN2+3*logf1[s]+logf2[s]+c1); // added an extra LN2
    ret[1][s] = addProtect2(ret[1][s],M_LN2+2*logf1[s]+logf2[s]+c1);       // added an extra LN2
    //no third
  }
  //Aa,Aa
  for(int s=0;s<pc.nSites;s++) {
    double c1 = loggl1[3*s+1]+loggl2[3*s+1];
    ret[0][s] = addProtect2(ret[0][s],2*M_LN2+2*(logf1[s]+logf2[s])+c1);
    ret[1][s] = addProtect2(ret[1][s],logf1[s]+logf2[s]+c1);
    ret[2][s] = addProtect2(ret[2][s],M_LN2+logf1[s]+logf2[s]+c1);
  }
  //Aa,aa
  for(int s=0;s<pc.nSites;s++) {
     double c1 = loggl1[3*s+1]+loggl2[3*s+2];
    ret[0][s] = addProtect2(ret[0][s],M_LN2+M_LN2+logf1[s]+3*logf2[s]+c1);  // added an extra LN2  
    ret[1][s] = addProtect2(ret[1][s],M_LN2+logf1[s]+2*logf2[s]+c1);        // added an extra LN2  
    //no third
  }
  //aa,AA
  for(int s=0;s<pc.nSites;s++) {
    double c1 = loggl1[3*s+2]+loggl2[3*s];
    ret[0][s] = addProtect2(ret[0][s],M_LN2+2*(logf1[s]+logf2[s])+c1);      // added an extra LN2
    //no sec
    //no third
  }
  //aa,Aa
  for(int s=0;s<pc.nSites;s++) {
    double c1 = loggl1[3*s+2]+loggl2[3*s+1];
    ret[0][s] = addProtect2(ret[0][s],M_LN2+M_LN2+logf1[s]+3*logf2[s]+c1);  // added an extra LN2 
    ret[1][s] = addProtect2(ret[1][s],M_LN2+logf1[s]+2*logf2[s]+c1);        // added an extra LN2 
    //no third
  }
  for(int s=0;s<pc.nSites;s++) {
    double c1 = loggl1[3*s+2]+loggl2[3*s+2];
    ret[0][s] = addProtect2(ret[0][s],4*logf2[s]+c1);
    ret[1][s] = addProtect2(ret[1][s],3*logf2[s]+c1);
    ret[2][s] = addProtect2(ret[2][s],2*logf2[s]+c1);
  }
  }
}


// Function for calculating transitions
// --------------------------------------     
/* Input:   pars is an array with a, k0, k1 and k2
            dpos is an array with distances between neightboring sites in Mb
	    l is nSites+1 so number of transitions
            res is a matrix 9xnSites where the transitions will be stored (should be allocated beforehand)
   Output:  ret is filled with transitions so ret[j][i] contains 
            P(X_i|X_i-1) for each j of nine possible combinations of X_i, X_i-1 in {0,1,2} 
*/

void calctrans(const double *pars,double *dpos,int l,double **ret) {

  // Set parameter values
  double a = pars[0];
  double k0= pars[1];
  double k1= pars[2];
  double k2= pars[3];

  // Precalc at=a*dpos
  double *at = new double[l];
  for(int i=0;i<l;i++)
    at[i] = a*dpos[i];

  // Calc transition probs for  X=0 -> 
  // T0->1 
  for(int i=0;i<l;i++)
    ret[T01][i] = (1.0-exp(-at[i]))*k1;
  // T0->2 
  if(k2==0){
    for(int i=0;i<l;i++)
      ret[T02][i] = 0;
  }else{
    for(int i=0;i<l;i++)
      ret[T02][i] = exp(-at[i]*k1)*k2/(k1-1)+exp(-at[i])*k1+exp(-at[i])*k0*k1/(k1-1)+k2;
    fixunderflow(ret[T02],l,1e-15); // NB. necessary?
  } // NB changed slightly 
  
  //T0->0
  for(int i=0;i<l;i++)
    ret[T00][i] = 1-ret[T01][i]-ret[T02][i];

  // Calc transition probs for  X=1 -> 
  //T1->0 and T1->2
  if(k1==0){
    memcpy(ret[T10],ret[T01],sizeof(double)*l); // not relevant since x!=1 so set to 0
    memcpy(ret[T12],ret[T01],sizeof(double)*l); // not relevant since x!=1 so set to 0
  }else{
    for(int i=0;i<l;i++){
      ret[T10][i] = ret[T01][i]*k0/k1;
      ret[T12][i] = ret[T01][i]*k2/k1;
    }
  }
  
  //T1->1
  for(int i=0;i<l;i++)
    ret[T11][i] = 1-ret[T10][i]-ret[T12][i];

  // Calc transition probs for  X=2 -> 
  //T2->1
  memcpy(ret[T21],ret[T01],sizeof(double)*l);

  //T2->0;
  if(k1==1){
    for(int i=0;i<l;i++)
      ret[T20][i] = 0;
  }else{
    for(int i=0;i<l;i++)
      ret[T20][i] = exp(-at[i]*k1)*k0/(k1-1.0)+exp(-at[i])*k1+exp(-at[i])*k2*k1/(k1-1.0)-(1-k1)*k0/(k1-1.0);
    fixunderflow(ret[T20],l,1e-15); // NB. necessary?
  }// NB changed slightly 
  
 //T2->2
  for(int i=0;i<l;i++)
    ret[T22][i] = 1-ret[T21][i]-ret[T20][i];
  fixunderflow(ret[T22],l,1e-15); // NB. necessary?
  
}



// Function for calculating log likelihood fast
// ------------------------------------------  
/* Input:   pars is an array with k0, k1 and k2 // NB function is called with allpars+1 so c(a,k0,k1,k2)[2:4]=c(k0,k1,k2) 
            logeprob is a 3xnSites matrix with log emission probabilities
            tprob is a 9x(nSites+1) matrix with transition probabilities 
            nSites is numbers of sites available for analysis
   Output:  returns log likelihood for one chr as double 
*/

double calcLoglikeChr(double *pars,double **logeprob,double **tprob,int nSites) {
  
  // First check if data is inconsistent with pars, if so return log(0) // NB added
  bool skip = false;
  double res = 0;
  for(int i=1;i<nSites;i++){
    if((pars[1]>(1-1e-15))&&(std::isinf(-logeprob[1][i]))){
      skip = true;
    }
    if((pars[2]>(1-1e-15))&&(std::isinf(-logeprob[2][i]))){
      skip = true;
    }
  }
  if(skip){
    res = log(0);
  }else{
  // If there are no inconsistencies log like is calculated 
  double logres[3]={log(pars[0]),log(pars[1]),log(pars[2])}; // c(log(k0),log(k1),log(k2))

  double lp[3];

  double m=0;
  for(int i=1;i<nSites;i++){
    lp[0] = exp(logres[0]+logeprob[0][i-1]-m);
    lp[1] = exp(logres[1]+logeprob[1][i-1]-m);
    lp[2] = exp(logres[2]+logeprob[2][i-1]-m);

    logres[0] = m+log(tprob[T00][i]*lp[0]+tprob[T10][i]*lp[1]+tprob[T20][i]*lp[2]);
    logres[1] = m+log(tprob[T01][i]*lp[0]+tprob[T11][i]*lp[1]+tprob[T21][i]*lp[2]);
    logres[2] = m+log(tprob[T02][i]*lp[0]+tprob[T12][i]*lp[1]+tprob[T22][i]*lp[2]);

    for(int ii=0;ii<3;ii++)
      lp[ii] = log(lp[ii])+m;
    m = max(lp);

    if(std::isnan(m)||std::isnan(logres[0])){
      fprintf(stderr,"Problems in likelihood at site:%d %f %f\n",i,m,logres[0]);
      exit(0);
    }
  }
  
  lp[0] = logeprob[0][nSites-1]+logres[0];
  lp[1] = logeprob[1][nSites-1]+logres[1];
  lp[2] = logeprob[2][nSites-1]+logres[2];

  m=max(lp);
  res = m+log(exp(lp[0]-m) + exp(lp[1]-m)+ exp(lp[2]-m));

  }
  return res;
}


// Function for calculating likelihoods for a full genome
// --------------------------------------     
/* Input:   pars is an array with a, k0, k1 and k2
            g is a struct that contains info (dpos, nSites, log emissions, transitions) for each of g.results.size() chromosomes 
   Output:  the total likelihood across chromosomes is returned as a double
*/

double calcLoglikeGenome(double *pars,const genome &g){

  assert(pars[0]>0);
  if(fabs(pars[1]+pars[2]+pars[3]-1)>1e-6)
    fprintf(stderr,"alpha=%f k0=%f k1=%f k2=%f sum:%e\n",pars[0],pars[1],pars[2],pars[3],pars[1]+pars[2]+pars[3]-1);

  double perchrloglike[g.results.size()];
  double totloglike =0;
  for(size_t i=0;i<g.results.size();i++){
    hmmRes results = g.results[i];
    calctrans(pars,results.dpos,results.nSites+1,results.trans);
    perchrloglike[i] = calcLoglikeChr(pars+1,results.logemis,results.trans,results.nSites); 
    totloglike += perchrloglike[i];
  }

  return totloglike;
}


// Function for doing forward calculations 
// -----------------------------------------------------  
/* Input:   pars is an array with k0, k1 and k2 // NB function is called with allpars+1 so c(a,k0,k1,k2)[2:4]=c(k0,k1,k2) 
            logeprob is 3xnSites matrix with log emission probabilities
            tprob is 9x(nSites+1) matrix with transition probabilities 
            nSites is numbers of sites available for analysis
	    loglike is where the final log likelihood is stored
   Output:  a 3 x nSites matrix res with forward values is returned
*/

double **forward(double *pars,double **logeprob,double **tprob,int nSites,double &loglike){
  
  // Make 3 x nSites matrix to store forward values and fill in these values 
  double **logf = new double*[3];

  // -  init
  double logks[3]={log(pars[0]),log(pars[1]),log(pars[2])}; // c(log(k0),log(k1),log(k2))  

  for(int k=0;k<3;k++){
    logf[k] = new double[nSites];
    logf[k][0] = logks[k]+logeprob[k][0]; 
  }

  double fm[3];
  double m=0;

  // - recursion
  for(int i=1;i<nSites;i++){
    fm[0] = exp(logf[0][i-1]-m);
    fm[1] = exp(logf[1][i-1]-m);
    fm[2] = exp(logf[2][i-1]-m);
    logf[0][i] = m+log(tprob[T00][i]*fm[0]+tprob[T10][i]*fm[1]+tprob[T20][i]*fm[2])+logeprob[0][i];
    logf[1][i] = m+log(tprob[T01][i]*fm[0]+tprob[T11][i]*fm[1]+tprob[T21][i]*fm[2])+logeprob[1][i];
    logf[2][i] = m+log(tprob[T02][i]*fm[0]+tprob[T12][i]*fm[1]+tprob[T22][i]*fm[2])+logeprob[2][i];
    m = max(logf[0][i],logf[1][i],logf[2][i]);
    
    if(std::isnan(m)||std::isnan(logks[0])){
      fprintf(stderr,"\t-> Problems in forward at site:%d\n",i);
      exit(0);
    }
  }
  // - termination
  m = max(logf[0][nSites-1],logf[1][nSites-1],logf[2][nSites-1]);
  loglike = m+log(exp(logf[0][nSites-1]-m) + exp(logf[1][nSites-1]-m)+ exp(logf[2][nSites-1]-m));
  // NB: the -m+m is a trick (add protect) to avoid underflow

  return logf;
}


// Function for doing backward calculations 
// -----------------------------------------------------  
/* Input:   pars is an array with k0, k1 and k2 // NB function is called with allpars+1 so c(a,k0,k1,k2)[2:4]=c(k0,k1,k2) 
            logeprob is 3xnSites matrix with log emission probabilities
            tprob is 9x(nSites+1) matrix with transition probabilities 
            nSites is numbers of sites available for analysis
	    loglike is where the final log likelihood is stored
   Output:  a 3 x nSites matrix res with backward values is returned
*/

double **backward(double *pars,double **logeprob,double **tprob,int nSites,double &loglike){

  // Make 3 x nSites matrix to store forward backward and fill in these values 
  double **logb = new double*[3];

  // - init  
  double logks[3]={log(pars[0]),log(pars[1]),log(pars[2])}; // c(log(k0),log(k1),log(k2))  
  for(int i=0;i<3;i++){
    logb[i] = new double[nSites];
    logb[i][nSites-1] = 0.0 ; 
  }
  double bm[3];
  double m=0;

  // - recursion  
  for(int i=nSites-2;i>=0;i--){
    bm[0] = exp(logeprob[0][i+1]+logb[0][i+1]-m);
    bm[1] = exp(logeprob[1][i+1]+logb[1][i+1]-m);
    bm[2] = exp(logeprob[2][i+1]+logb[2][i+1]-m);

    logb[0][i] = m+log(tprob[T00][i+1]*bm[0]+tprob[T01][i+1]*bm[1]+tprob[T02][i+1]*bm[2]);
    logb[1][i] = m+log(tprob[T10][i+1]*bm[0]+tprob[T11][i+1]*bm[1]+tprob[T12][i+1]*bm[2]);
    logb[2][i] = m+log(tprob[T20][i+1]*bm[0]+tprob[T21][i+1]*bm[1]+tprob[T22][i+1]*bm[2]);

    m = max(logb[0][i],logb[1][i],logb[2][i]);

    if(std::isnan(m)||std::isnan(logks[0])){
      fprintf(stderr,"Probs in backward at site:%d\n",i);
      exit(0);
    }
  }

  // - termination 
  double logp[3]={logb[0][0]+logeprob[0][0]+logks[0],
		  logb[1][0]+logeprob[1][0]+logks[1],
		  logb[2][0]+logeprob[2][0]+logks[2]};
  m=max(logp);
  loglike = m+log(exp(logp[0]-m) + exp(logp[1]-m)+ exp(logp[2]-m));
  // NB: the -m+m is a trick (add protect) to avoid underflow  

  return logb;
}

//####################################
// Init functions
//####################################


// Function that takes data from  beagle files and frequencies and makes a vector with all relevant data split up by chr
// ------------------------------------------------------ 
/* Input:   bgl is data read from beagle file 
            freq is vectors with frequencies
   Output:  returns a vector with nchr perChr structs with both gls and freqs (after maf filtering!) 
*/

std::vector<perChr> makeDat(const bgl& in,std::vector<double>& freq,double minfreq,int switchmaf){

  int start=0;
  int stop=1;

  // Make vector that tracks switches between chrs and counts sites in each chr that passes maf filter
  std::vector<int> splits;
  std::vector<int> keepcnts;
  int keepcnt = 0;
  int totalkeep = 0;

  if((freq[0]>minfreq)&&((1-freq[0])>minfreq)){
    keepcnt++;
  }
  for(int i=1;i<in.nSites;i++){
    if(strcmp(in.chr[i],in.chr[i-1])!=0){
      splits.push_back(i);
      totalkeep +=keepcnt;
      keepcnts.push_back(keepcnt);
      keepcnt=0;
    }
    if((freq[i]>minfreq)&&((1-freq[i])>minfreq)){
      keepcnt++;
    }

  }
  totalkeep +=keepcnt;
  keepcnts.push_back(keepcnt);
  splits.push_back(in.nSites);

  fprintf(stderr,"\t-> Data divided into %lu chromosomes.\n",splits.size());
  std::vector<perChr> ret;

  // Then fill in info for each chr for those SNPs that pass the maf filter
  int freqcnt=0;
  for(size_t i=0;i<splits.size();i++){
    perChr pc;
    pc.name = strdup(in.chr[splits[i]-1]);
    pc.nInd=in.nInd;
    pc.nSites = keepcnts[i];
    if(pc.nSites==0){
      fprintf(stderr,"\t-> Warning: The number of sites left on chromosome %zu after frequency filtering is 0, so this chromosome is excluded\n", i+1);
    }else{

    pc.pos = new int[pc.nSites];
    pc.loggl = new double*[pc.nSites];
    pc.keepind = new int [pc.nSites];
    pc.logfreq = new double [pc.nSites];
    pc.logqerf = new double [pc.nSites];
    pc.dpos = new double [pc.nSites-1];
    
    start = i!=0?splits[i-1]:0;
    stop = splits[i];
    int cnt =0;
    for(int j=start;j<stop;j++){
      if((freq[freqcnt]>minfreq)&&(1-freq[freqcnt])>minfreq){ 
	pc.pos[cnt] = in.pos[j];
	pc.loggl[cnt] = in.genos[j];
	if(switchmaf>0){
	  pc.logfreq[cnt] = log(1-freq[freqcnt]);
	  pc.logqerf[cnt] = log(freq[freqcnt]);
	}else{
	  pc.logfreq[cnt] = log(freq[freqcnt]);
	  pc.logqerf[cnt] = log(1-freq[freqcnt]);
	}
	cnt++;
      }
      freqcnt++;
    }
    transpose2<double>(&pc.loggl,pc.nSites,pc.nInd);
    ret.push_back(pc);
  }
  }
  fprintf(stderr,"\t-> Frequency filter removed %lu sites out of %lu\n",freq.size()-totalkeep,freq.size());
  if(totalkeep==0){
    fprintf(stderr,"\t-> Error: The number of sites left after frequency filtering is 0, so running an analysis makes no sense\n");
    exit(0);
  }

  return ret;
}



// Function for allocating and precomputing the data structures for given pair that donw depend on the a,k0,k1,k2doing viterbi decoding
// ------------------------------------------------------ 
/* Input:   pc is a vector of perChr structs with info about positions and nSites in each chr 
            p is a structure that contains information about which pairs to analyse
   Output:  returns a genome structure with precalculated log emissions and distances between loci for each chr 
            plus allocated space for transitions and other info (e.g. decoding results) for later
*/

genome mkGenome(const std::vector<perChr> &pc,const para &p){
  genome g;
  // For each chromosome make relevant datastructures and precalculate log emission and transition probs
  for(size_t i=0;i<pc.size();i++){
    hmmRes res;
    res.pos = pc[i].pos;
    res.nSites = pc[i].nSites;
    res.trans = new double*[9];
    for(int j=0;j<9;j++)
      res.trans[j] = new double[pc[i].nSites+1];
    res.logemis = new double*[3];
    for(int j=0;j<3;j++){
      res.logemis[j] = new double[pc[i].nSites];
    }

    calclogemis(p.pair[0],p.pair[1],pc[i],res.logemis);
    res.dpos = diffPos(pc[i].pos,pc[i].nSites);
    g.results.push_back(res);
  }
  return g;
}


//####################################
// Decode functions
//####################################

// Function for doing viterbi decoding
// -----------------------------------------------------  
/* Input:   pars is a vector with k0, k1 and k2 // NB function is called with allpars+1 so c(a,k0,k1,k2)[2:4]=c(k0,k1,k2) 
            logeprob is 3xnSites matrix with log emission probabilities
            tprob is 9x(nSites+1) matrix with transition probabilities 
            nSites is numbers of sites available for analysis
   Output:  a nSites long vector with the inferred viterbi path is returned
*/

char *viterbi(double *pars,double **logeprob,double **tprob,int nSites){

  // - make container for tmp result, res and ptr
  double logres[3]={log(pars[0]),log(pars[1]),log(pars[2])}; // c(log(k0),log(k1),log(k2))  
  double **res = new double*[3];
  char **ptr = new char*[3];

  // - init res (as P(X0=x)P(g0|X0=x) for the 3 diff possible x vals) and ptr  
  for(int i=0;i<3;i++){
    ptr[i] = new char[nSites];
    res[i] = new double[nSites];
    res[i][0] = logres[i]+logeprob[i][0]; 
  }

  // - recursion    
  for(int i=1;i<nSites;i++){
    res[0][i] = logeprob[0][i]+max(log(tprob[T00][i])+res[0][i-1],log(tprob[T10][i])+res[1][i-1],log(tprob[T20][i])+res[2][i-1]);
    res[1][i] = logeprob[1][i]+max(log(tprob[T01][i])+res[0][i-1],log(tprob[T11][i])+res[1][i-1],log(tprob[T21][i])+res[2][i-1]);
    res[2][i] = logeprob[2][i]+max(log(tprob[T02][i])+res[0][i-1],log(tprob[T12][i])+res[1][i-1],log(tprob[T22][i])+res[2][i-1]);
    ptr[0][i] =           whichmax(log(tprob[T00][i])+res[0][i-1],log(tprob[T10][i])+res[1][i-1],log(tprob[T20][i])+res[2][i-1]);
    ptr[1][i] =           whichmax(log(tprob[T01][i])+res[0][i-1],log(tprob[T11][i])+res[1][i-1],log(tprob[T21][i])+res[2][i-1]);
    ptr[2][i] =           whichmax(log(tprob[T02][i])+res[0][i-1],log(tprob[T12][i])+res[1][i-1],log(tprob[T22][i])+res[2][i-1]);
  }

  // - termination 
  double loglike = max(res[0][nSites-1],res[1][nSites-1],res[2][nSites-1]);
  char *vit = new char[nSites];
  vit[nSites-1] = whichmax(res[0][nSites-1],res[1][nSites-1],res[2][nSites-1]); 

  // - traceback
  for(int i=(nSites-1);i>0;i--)
    if(vit[i]==2)
      vit[i-1] = ptr[2][i];
    else if(vit[i]==1)
      vit[i-1] = ptr[1][i];
    else 
      vit[i-1] = ptr[0][i];

  // - return most probable path 
  return vit;
}


// Function for doing posterior decoding 
// -----------------------------------------------------  
/* Input:   fw result of running the forward algorithm
            bw result of running the backward algorithm
            nSites is numbers of sites available for analysis
	    loglike the log of the likelihood obtained e.g. using the backward algorithm
   Output:  a nSites long vector with the inferred viterbi path is returned
*/

double **post(double **fw,double **bw,int nSites,double loglike){
  double **res = new double*[3];
  
  for(int i=0;i<3;i++){
    res[i] = new double[nSites];
  }
  for(int i=0;i<3;i++)
    for(int s=0;s<nSites;s++){
      res[i][s] = exp(fw[i][s]+bw[i][s]-loglike);
    }

  return res;
}

// Function for doing all decoding across all chromosomes
// --------------------------------------------------------  
/* Input:   pars is a vector with a, k0, k1 and k2 
            g is a struct that contains info (dpos, nSites, log emissions, transitions) for each of g.results.size() chromosomes 
   Output:  a nSites long vector with the inferred viterbi path is returned
*/

void forward_backward_decode_viterbi(double *pars,genome &g){
  assert(pars[0]>0);
  if(fabs(pars[1]+pars[2]+pars[3]-1)>1e-6)
    fprintf(stderr,"Alpha=%f k0=%f k1=%f k2=%f sum of ks:%f\n",pars[0],pars[1],pars[2],pars[3],pars[1]+pars[2]+pars[3]);
  double loglikefw,loglikebw;
  loglikefw=loglikebw=0;
  for(size_t i=0;i<g.results.size();i++){
    hmmRes &res = g.results[i];
    double tmp =0;
    res.forward = forward(pars+1,res.logemis,res.trans,res.nSites,tmp);
    loglikefw += tmp;
    tmp=0;
    res.backward = backward(pars+1,res.logemis,res.trans,res.nSites,tmp);
    loglikebw += tmp;
    res.post = post(res.forward,res.backward,res.nSites,tmp);
    res.viterbi = viterbi(pars+1,res.logemis,res.trans,res.nSites);
  }

}



//####################################
// Optim functions
//####################################


// Functions for doing optim when CalcA==TRUE and k2==0. So only one free parameter
// --------------------------------------------------------  

// 1) Function that returns -loglike for parameter values provided CalcA==TRUE and k2==0 (input pars={k0})
double bfgs_call_k2zero_calcA2(const double* pars,const void *dats){

  toOptim2 *to = (toOptim2*)dats;

  // Make array inV with all parameters including a calculated from the given ks
  double *inV = to->tsk;
  inV[1]=pars[0]; // k0
  inV[2]=1-inV[1];// k1
  inV[3]=0;       // k2 
  inV[0]= calculateA(inV[1],inV[2],inV[3],PHI); // a calculated

  // If calculated a is outside alim a very high -loglike is returned
  // (to make sure bfgs does not use this)
  if(inV[0]<alim[0]||inV[0]>alim[1]){ // changed
    if(verbose){
      fprintf(stderr,"CalcA is outside bounds: a=%f k0=%f\n",inV[0],inV[1]);
    }
    return DBL_MAX; 
  }

  // If not then -loglike is calculated and returned 
  double minusloglike=-calcLoglikeGenome(inV,to->g);
  if(verbose){
    fprintf(stderr,"\t   [%s] a=%f k0=%f ",__FUNCTION__,inV[0],pars[0]);
    fprintf(stderr," loglike=%f\n",-minusloglike);
  }
  return minusloglike;
}

// 2) function that uses the above function+random start parameters values as input to bfgs
//    (which then minimizes the above function, and thus -loglike and returns -min(-loglike)=max(loglike))
double run_optim_k2zero_calcA2(const genome &g,const std::vector<perChr>&pc,para&p){

  // Init optim 
  toOptim2 *to = new toOptim2(g,pc,p);
  double tsk[4];//tmparray used of saving stackpointers
  to->tsk = tsk;

  // Find random starting value for k0 and specify limits to the parameter
  double pars[1] = {myRand2(0,1)};//{(double)rand()/(double)(RAND_MAX)};//{myRand(0,1)};
  double k1 = 1-pars[0];
  double k2 = 0;
  p.a = calculateA(pars[0],k1,k2,PHI);
  if(testrandomstart){
    while((p.a<alim[0])||(p.a>alim[1])||(pow((k1+2*k2),2.0)<4*k2)){
      pars[0]=myRand(0,1);
      k1=1-pars[0];
      k2=0;
      p.a=calculateA(pars[0],k1,k2,PHI);
    }
  }
  
  double lbd[1]={klim[0]};
  double ubd[1]={klim[1]};
  int nbd[1]={2};

  // Run bfgs with the random start values to find pars that maximizes loglike
  // (by minimizing bfgs_call_k2zero_calcA2 which returns -loglike) 
  double opt = findmax_bfgs(1,pars,(void *)to,bfgs_call_k2zero_calcA2,NULL,lbd, ubd,nbd, -1);

  // Set pars according to bfgs and user defined values
  p.k0=pars[0];
  p.k1=1-p.k0;
  p.k2=0;
  p.a=calculateA(p.k0,p.k1,p.k2,PHI);

  // Return loglike
  return opt;
}


// Functions for doing optim when CalcA==TRUE. So two free parameters (k0 and k1)
// --------------------------------------------------------  

// 1) Function that returns -loglike for parameter values provided and CalcA==TRUE (input pars={k0,k1})
double bfgs_call_full_calcA2(const double* pars,const void *dats){

  toOptim2 *to = (toOptim2*)dats;

  // Make array inV with all parameters including a calculated from the given ks
  double *inV = to->tsk;
  inV[1]=pars[0];
  inV[2]=pars[1];
  inV[3]=1-inV[1]-inV[2];
  inV[0]= calculateA(inV[1],inV[2],inV[3],PHI);

  // If the parameter values are outside their limits a very large -loglike values is returned (to make sure bfgs does not use these)
  if((pars[0]+pars[1])>1)
    return DBL_MAX;
  if(std::isnan(inV[0]))
    return DBL_MAX;
  if((inV[0]<alim[0])||(inV[0]>alim[1])){ // NB added
    return DBL_MAX;
  }
  // If not then -loglike is calculated and returned 
  double minusloglike=-calcLoglikeGenome(inV,to->g);
  if(verbose){
    fprintf(stderr,"[%s] k0=%f k1=%f\n",__FUNCTION__,pars[0],pars[1]);
    fprintf(stderr,"loglike=%f\n",-minusloglike);
  }
  return minusloglike;
}

// 2) function that uses the above function+random start parameters values as input to bfgs
//    (which then minimizes the above function, and thus -loglike and returns -min(-loglike)=max(loglike))
double run_optim_full_calcA2(const genome &g,const std::vector<perChr>&pc,para&p){

  // Init optim 
  toOptim2 *to = new toOptim2(g,pc,p);
  double tsk[4];//tmparray used of saving stackpointers
  to->tsk = tsk;

  // Find random starting value for k0 and k1 and specify limits to the parameters
  double pars[2];
  double k0 = myRand2(0,1);
  double k1 = myRand2(0,1-k0);
  double k2 = 1-k0-k1;
  if(testrandomstart){
    while((p.a<alim[0])||(p.a>alim[1])||(pow((k1+2*k2),2.0)<4*k2)){
      k0 = myRand2(0,1);
      k1 = myRand2(0,1-k0);
      k2 = 1-k0-k1;
      p.a = calculateA(k0,k1,k2,PHI);
    }
    pars[0]=k0;
    pars[1]=k1;
  }
    
  double lbd[2]={klim[0],klim[0]};
  double ubd[2]={klim[1],klim[1]};
  int nbd[2]={2,2};

  // Run bfgs with the random start values to find pars that maximizes loglike
  // (by minimizing bfgs_call_full_calcA2 which returns -loglike) 
  double opt= findmax_bfgs(2,pars,(void *)to, bfgs_call_full_calcA2,NULL,lbd, ubd,nbd, -1);
  
  // Set pars according to bfgs and user defined values
  p.k0=pars[0];
  p.k1=pars[1];
  p.k2=1-pars[0]-pars[1];
  p.a=calculateA(p.k0,p.k1,p.k2,PHI);

  // Return loglike
  return opt;
}


// Functions for doing optim when CalcA==FALSE but fixk2to0==TRUE. So two free parameters (a and k0)
// --------------------------------------------------------  

// 1) Function that returns -loglike for parameter values provided and CalcA==TRUE (input pars={k0,k1})
double bfgs_call_full_k2zero2(const double* pars,const void *dats){

  toOptim2 *to = (toOptim2*)dats;

  // Make array inV with all parameters including a calculated from the given ks
  double *inV = to->tsk;
  inV[1]=pars[1];
  inV[2]=1-pars[1];
  inV[3]=0;
  inV[0]=pars[0];

  // If the parameter values are outside their limits a very large -loglike values is returned
  // (to make sure bfgs does not use these)
  if((inV[1]>1)||(inV[1]<0))
    return DBL_MAX;
  if((inV[0]<alim[0])||(inV[0]>alim[1])){
    // NB added since there is a bug in the bfgs algorithm that makes the optim algo try outside of parameter space!
    return DBL_MAX;
  }

  // If not then -loglike is calculated and returned 
  double minusloglike=-calcLoglikeGenome(inV,to->g);
  if(verbose){
    fprintf(stderr,"\t   [%s] a=%f k0=%f ",__FUNCTION__,pars[0],pars[1]);
    fprintf(stderr,"loglike=%f\n",-minusloglike);
  }
  return minusloglike;
}

// 2) function that uses the above function+random start parameters values as input to bfgs
//    (which then minimizes the above function, and thus -loglike and returns -min(-loglike)=max(loglike))
double run_optim_full_k2zero2(const genome &g,const std::vector<perChr>&pc,para&p){

  // Init optim 
  toOptim2 *to = new toOptim2(g,pc,p);
  double tsk[4];//tmparray used of saving stackpointers
  to->tsk = tsk;

  // Find random starting value for a and k0 and specify limits to the parameters
  double pars[2];
  pars[0]=myRand(alim[0],alim[1]);
  pars[1]=myRand(0,1);
  double lbd[2]={alim[0],klim[0]};
  double ubd[2]={alim[1],klim[1]};
  int nbd[2]={2,2};

  // Run bfgs with the random start values to find pars that maximizes loglike
  // (by minimizing bfgs_call_k2zero2 which returns -loglike) 
  double opt= findmax_bfgs(2,pars,(void *)to, bfgs_call_full_k2zero2,NULL,lbd, ubd,nbd, -1);

  // Set pars according to bfgs and user defined values
  p.k0=pars[1];
  p.k1=1-pars[1];
  p.k2=0;
  p.a=pars[0];

  // Return loglike
  return opt;
}


// Functions for doing optim when CalcA==FALSE. So three free parameters (a, k0 and k1)
// --------------------------------------------------------  

// 4 parameter version

// 1) Function that returns -loglike for parameter values provided (input pars={a,k0,k1,k2})
double bfgs_call_full3(const double* pars,const void *dats){

  toOptim2 *to = (toOptim2*)dats;

  // Make array inV with all parameters including a calculated from the given ks
  double *inV = to->tsk;
  inV[0]=pars[0];
  inV[1]=pars[1];
  inV[2]=pars[2];
  inV[3]=pars[3];
  double tsum=inV[1]+inV[2]+inV[3];
  inV[1] /= tsum;
  inV[2] /= tsum;
  inV[3] /= tsum;

  // If the parameter values are outside their limits a very large -loglike values is returned
  // (to make sure bfgs does not use these)
  if((inV[0]<alim[0])||(inV[0]>alim[1])) 
    return DBL_MAX;

  // If not then -loglike is calculated and returned   
  double minusloglike=-calcLoglikeGenome(inV,to->g);
  if(verbose){
    fprintf(stderr,"[%s] %f %f %f %f\n",__FUNCTION__,pars[0],pars[1],pars[2],pars[3]);
    fprintf(stderr,"loglike=%f\n",-minusloglike);
  }
  return minusloglike;
}


// 2) function that uses the above function+random start parameters values as input to bfgs
//    (which then minimizes the above function, and thus -loglike and returns -min(-loglike)=max(loglike))
double run_optim_full3(const genome &g,const std::vector<perChr>&pc,para&p){

  // Init optim
  toOptim2 *to = new toOptim2(g,pc,p);
  double tsk[4];//tmparray used of saving stackpointers
  to->tsk = tsk;

  // Find random starting value for a, k0, k1 and k2 and specify limits to the parameters  
  double pars[4];
  pars[0]=myRand(alim[0],alim[1]);
  pars[1]=myRand(0,1);
  pars[2]=myRand(0,1-pars[1]); //Changed: was originally myRand(0,pars[1]); 
  pars[3]=1-pars[1]-pars[2];
  //  double pars[3]={0.045854, 0.823746,0.176254};
  double lbd[4]={alim[0],klim[0],klim[0],klim[0]};
  double ubd[4]={alim[1],klim[1],klim[1],klim[1]};
  int nbd[4]={2,2,2,2};

  // Run bfgs with the random start values to find pars that maximizes loglike
  // (by minimizing bfgs_call_full3 which returns -loglike) 
  double opt= findmax_bfgs(4,pars,(void *)to, bfgs_call_full3,NULL,lbd, ubd,nbd, -1);

  // Set pars according to bfgs and user defined values
  p.a=pars[0];
  double tsum = pars[1]+pars[2]+pars[3];
  p.k0=pars[1]/tsum;
  p.k1=pars[2]/tsum;
  p.k2=pars[3]/tsum;

  // Return loglike
  return opt;
}

// 3 parameter version
// 1) Function that returns -loglike for parameter values provided (input pars={a,k0,k1})
double bfgs_call_full2(const double* pars,const void *dats){
  toOptim2 *to = (toOptim2*)dats;
  double *inV = to->tsk;

  // Make array inV with all parameters including a calculated from the given ks
  inV[0]=pars[0];
  inV[1]=pars[1];
  inV[2]=pars[2];
  inV[3]=1-pars[1]-pars[2];

  // If the parameter values are outside their limits a very large -loglike values is returned
  // (to make sure bfgs does not use these)
  if((inV[1]+inV[2])>1){
    //fprintf(stderr,"[%s] a=%f k0=%f k1=%f ",__FUNCTION__,pars[0],pars[1],pars[2]);
    //fprintf(stderr,"loglike=%f\n",-DBL_MAX);
    return DBL_MAX;
  }
  // If not then -loglike is calculated and returned   
  double minusloglike=-calcLoglikeGenome(inV,to->g);
  if(verbose){
    fprintf(stderr,"[%s] a=%f k0=%f k1=%f ",__FUNCTION__,pars[0],pars[1],pars[2]);
    fprintf(stderr,"loglike=%f\n",-minusloglike);
  }
  return minusloglike;
}

// 2) function that uses the above function+random start parameters values as input to bfgs
//    (which then minimizes the above function, and thus -loglike and returns -min(-loglike)=max(loglike))
double run_optim_full2(const genome &g,const std::vector<perChr>&pc,para&p){

  // Init optim
  toOptim2 *to = new toOptim2(g,pc,p);
  double tsk[4];//tmparray used of saving stackpointers
  to->tsk = tsk;

  // Find random starting value for a, k0 and k1 and specify limits to the parameters  
  double pars[3];
  pars[0]=myRand(alim[0],alim[1]);
  pars[1]=myRand(0,1);
  pars[2]=myRand(0,1-pars[1]); 
  double lbd[3]={alim[0],klim[0],klim[0]};
  double ubd[3]={alim[1],klim[1],klim[1]};
  int nbd[3]={2,2,2};
  
  // Run bfgs with the random start values to find pars that maximizes loglike
  // (by minimizing bfgs_call_full2 which returns -loglike) 
  double opt= findmax_bfgs(3,pars,(void *)to, bfgs_call_full2,NULL,lbd, ubd,nbd, -1);
  p.a=pars[0];
  p.k0=pars[1];
  p.k1=pars[2];
  p.k2=1-pars[1]-pars[2];

  // Return loglike
  return opt;
}



// Function for running the parameter optim according to input parameters
// ------------------------------------------------------ 
/* Input:   p is a structure that contains information about input values of a, k0, k1, k2 and which pairs to analyse
            g is a genome structure with precalcalted logemission and allocated space for inference results
            pc is a vector of perChr structs with info about nInds, positions, frequencies and nSites in each chr 
            calcA a boolean that states whether to calculate a (as opposed to estimate it)
   Output:  returns a genome structure with precalculated log emissions and distances between loci for each chr 
            plus allocated space for transitions and other info (e.g. decoding results) for later
*/
double doOptim(para &p,const genome &g,const std::vector<perChr>&pc,int calcA,int seed=0,int nopti=1){

  srand48((unsigned) seed);
  double loglike;
  double curmaxloglike=-DBL_MAX;
  double curmaxpars[4];

  if(calcA==1&&p.k0==-1&&p.k1==-1&&p.k2==0&&p.a==-1){
    fprintf(stderr,"\t-> [%s] 1) Optimizing k0 (with k2=zero, k1=1-k0 and calcA=TRUE)\n",__FUNCTION__);
    for(int i=0;i<nopti;i++){
      loglike = run_optim_k2zero_calcA2(g,pc,p);
      fprintf(stderr,"\t   Obtained optimum %d: a=%f k0=%f k1=%f k2=%f loglike=%.16f",i, p.a,p.k0,p.k1,p.k2,loglike);
      if(loglike>curmaxloglike){
	curmaxpars[0]=p.a;
	curmaxpars[1]=p.k0;
	curmaxpars[2]=p.k1;
	curmaxpars[3]=p.k2;
	curmaxloglike = loglike;
	fprintf(stderr," (new optimum)\n");
      }else{
	fprintf(stderr,"\n");
      }
    }
  }else if(calcA==1&&p.k0==-1&&p.k1==-1&&p.k2==-1&&p.a==-1){
    fprintf(stderr,"\t-> [%s] 2) Optimizing k0 and k1 (with k2=1-k0-k1 and calcA=TRUE)\n",__FUNCTION__);
    for(int i=0;i<nopti;i++){
      loglike = run_optim_full_calcA2(g,pc,p);
      fprintf(stderr,"\t   Obtained optimum %d: a=%f k0=%f k1=%f k2=%f loglike=%.16f",i, p.a,p.k0,p.k1,p.k2,loglike);
      if(loglike>curmaxloglike){
	curmaxpars[0]=p.a;
	curmaxpars[1]=p.k0;
	curmaxpars[2]=p.k1;
	curmaxpars[3]=p.k2;
	curmaxloglike = loglike;
	fprintf(stderr," (new optimum)\n");
      }else{
	fprintf(stderr,"\n");
      }
    }
  }else if(calcA!=1&p.k0==-1&&p.k1==-1&&p.k2==0&&p.a==-1){
    fprintf(stderr,"\t-> [%s] 3) Optimizing a and k0 (with k1=1-k0 and k2=0)\n",__FUNCTION__);
    for(int i=0;i<nopti;i++){
      loglike = run_optim_full_k2zero2(g,pc,p);
      fprintf(stderr,"\t   Obtained optimum %d: a=%f k0=%f k1=%f k2=%f loglike=%.16f",i, p.a,p.k0,p.k1,p.k2,loglike);
      if(loglike>curmaxloglike){
	curmaxpars[0]=p.a;
	curmaxpars[1]=p.k0;
	curmaxpars[2]=p.k1;
	curmaxpars[3]=p.k2;
	curmaxloglike = loglike;
	fprintf(stderr," (new optimum)\n");
      }else{
	fprintf(stderr,"\n");
      }
    }
  }else if(calcA!=1&p.k0==-1&&p.k1==-1&&p.k2==-1&&p.a==-1){
    fprintf(stderr,"\t-> [%s] 4) Optimizing a, k0 and k1 (with k2=1-k0-k1)\n",__FUNCTION__);
    for(int i=0;i<nopti;i++){
      loglike = run_optim_full2(g,pc,p);
      fprintf(stderr,"\t   Obtained optimum %d: a=%f k0=%f k1=%f k2=%f loglike=%.16f",i, p.a,p.k0,p.k1,p.k2,loglike);
      if(loglike>curmaxloglike){
	curmaxpars[0]=p.a;
	curmaxpars[1]=p.k0;
	curmaxpars[2]=p.k1;
	curmaxpars[3]=p.k2;
	curmaxloglike = loglike;
	fprintf(stderr," (new optimum)\n");
      }else{
	fprintf(stderr,"\n");
      }
    }
  }else{
    fprintf(stderr,"\t-> [%s] Optimization not implemented for this combination of options\n",__FUNCTION__);
    exit(0);
  }
  p.a=curmaxpars[0];
  p.k0=curmaxpars[1];
  p.k1=curmaxpars[2];
  p.k2=curmaxpars[3];
  loglike = curmaxloglike;
  return loglike;
} 








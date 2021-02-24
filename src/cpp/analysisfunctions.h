#define PHI 0.013

double calculateA(double k0,double k1, double k2,double phi);

typedef struct{
  char *name;
  int nInd;
  int nSites;
  int *pos;
  double *logfreq;//in log
  double *logqerf;//in log 1-freq
  double **loggl;//in log
  double *dpos;// distacne between positions
}perChr;

typedef struct{
  double a;
  double k0;
  double k1;
  double k2;
  int pair[2];
  int calcA;
}para;


//for each chr
typedef struct{
  int nSites;
  int *pos;
  double **post;
  double **forward;
  double **backward;
  double **logemis;
  double **trans;
  char *viterbi;
  double *dpos;//<-only used internally, this is not something that should be used for results // maybe move to perChr
}hmmRes;

typedef struct{
  double pars[4];
  double like;
  double pologlike;
  double uloglike;
  std::vector<hmmRes> results;
  std::vector<int> keepsite;
}genome;


// Struct for data (all the data from the beagle file)
typedef struct{
  double **genos;
  char **chr;
  int *pos;
  int nSites;
  int nInd;
}bgl;


//hmm analysis(const perChr &pc,double *freq,para p,int calcA);
//double addProtect2(double a,double b);
//double addProtect3(double a,double b, double c);
void printPars(FILE *fp,para p);
//std::vector<perChr> makeDat(const bgl& in);
std::vector<perChr> makeDat(const bgl& in,std::vector<double>& freq,double minfreq,int switchmaf);
genome mkGenome(const std::vector<perChr> &pd,const para &p);
double doOptim(para &p,const genome &g,const std::vector<perChr>&pc,int calcA,int seed,int nopti);
void forward_backward_decode_viterbi(double *pars,genome &g);
double calcLoglikeGenome(double *pars,const genome &g);




#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h> 
#include <libgen.h> 
#include <time.h>
#include <vector>
#include <zlib.h>  
#include <cmath>
#include <cassert>
#include "analysisfunctions.h"
#include "filehandlingfunctions.h"


typedef struct{
  const char *freqfile;
  const char *glbinfile;
  const char *glbeaglefile;
  const char *outname;
  int fixk2to0;
  int nInd;
  para p;
  float minMaf;
  int switchMaf;
  int seed;
  int nOpti;
}cArg;


void print_info(FILE *fp){
  fprintf(fp, "\n\n");
  fprintf(fp, "You are using LocalNgsRelate version 0.999 (build time: %s:%s)\n",__DATE__,__TIME__);
  fprintf(fp, "\nUsage: ./localngsrelate  [options] \n");
  fprintf(fp, "\nRequired options:\n");
  fprintf(fp, "   -f         <filename>    Name of file with frequencies\n");
  //  fprintf(fp, "   -g         <fileprefix>  Prefix of files with genotype likelihoods and positions (suffix .glf.gz and gls.pos.gz)\n"); # Not implemented
  fprintf(fp, "   -gbeagle   <fileprefix>  Prefix of files with genotype likelihoods (assumed suffix is .beagle.gz)\n");
  fprintf(fp, "   -n         <INT>         Number of samples in genotype likelihood file\n");
  fprintf(fp, "   -a         <INT>         First individual used for analysis? (zero offset). NB if not specified it is set to 0!\n");
  fprintf(fp, "   -b         <INT>         Second individual used for analysis? (zero offset). NB. if not specified it is set to 1!\n\n");

  fprintf(fp, "Optional options:\n");
  fprintf(fp, "   -calcA     <INT>         Calculate the parameter a instead of estimating it (if set to 1, default is 0)\n");
  fprintf(fp, "   -fixA      <FLOAT>       Set parameter a (rate of change between IBD states) to this value. NB only works when k0, k1 and k2 are also provided\n"); 
  fprintf(fp, "   -minA      <FLOAT>       Set lower limit to parameter a (default 0.001). NB. it has to be > 0\n");
  fprintf(fp, "   -maxA      <FLOAT>       Set upper limit to parameter a (default 0.15)\n");
  fprintf(fp, "   -k0        <FLOAT>       Set k0 to this value (must be used with -k1 and -k2)\n");
  fprintf(fp, "   -k1        <FLOAT>       Set k1 to this value (must be used with -k0 and -k2)\n");
  fprintf(fp, "   -k2        <FLOAT>       Set k2 to this value (must be used with -k0 and -k1)\n");
  fprintf(fp, "   -fixk2to0  <INT>         Set k2 to 0 and estimate only k0 and k1\n");
  fprintf(fp, "   -l         <INT>         minor allele frequency filter (default 0.05)\n");
  fprintf(fp, "   -s         <INT>         Should you swich the freq with 1-freq? (default 0)\n");
  fprintf(fp, "   -r         <FLOAT>       Seed for rand (default 0)\n");
  fprintf(fp, "   -N         <INT>         Number of times to start parameter estimation for each pair with random seed (default 1)\n");
  fprintf(fp, "   -O         <fileprefix>  Prefix for name of outputfiles\n");
  fprintf(fp, "\n");
  exit(0);
}

double alim[2];
int runoldrelateV; 

int main(int argc, char **argv){

  // Start run time clock
  clock_t t=clock();
  time_t t2=time(NULL);
  
  // Set default parameters
  cArg ca;
  ca.freqfile=NULL;
  ca.glbinfile =NULL;
  ca.glbeaglefile=NULL;
  ca.outname=NULL;
  ca.nInd=-1;
  ca.p.calcA=-1;
  ca.p.pair[0]=0;
  ca.p.pair[1]=1;
  ca.p.a=-1;
  ca.p.k0=-1;
  ca.p.k1=-1;
  ca.p.k2=-1;
  ca.fixk2to0 = 0;
  ca.minMaf = 0.05;
  ca.switchMaf = 0;
  ca.seed = 20;
  ca.nOpti = 20;
  alim[0]=0.00000001;
  alim[1]=0.10;
  runoldrelateV = 0; // Remove?
  
  // Specify the expected options for parsing purposes
  static struct option long_options[] = {
      {"f",       required_argument, 0,  'f' },
      {"g"     ,  required_argument, 0, 'g' },
      {"gbeagle", required_argument, 0,  'G' },
      {"n",       required_argument, 0,  'n' },
      {"a",       required_argument, 0,  'a' },
      {"b",       required_argument, 0,  'b' },
      {"calcA",   required_argument, 0,  'c' },
      {"fixA",    required_argument, 0,  'A' },
      {"minA",    required_argument, 0,  'h' },
      {"maxA",    required_argument, 0,  'j' },
      {"k0",      required_argument, 0,  'x' },
      {"k1",      required_argument, 0,  'y' },
      {"k2",      required_argument, 0,  'w' },
      {"fixk2to0",required_argument, 0,  'v' },
      {"l",       required_argument, 0,  'l' },
      {"s",       required_argument, 0,  's' },
      {"r",       required_argument, 0,  'r' },
      {"N",       required_argument, 0,  'N' },
      {"O",       required_argument, 0,  'O' },
      {"old",     required_argument, 0,  'o' }, // Remove when done?
      {0,        0,                0,  0   }
    };
  
  // Check if any options are specified if not then print calling instruction
  if(argc==1){// if no arguments, print info on program
    print_info(stderr);     
    return 0;
  }
  
  // If there are any options specified then read them
  int opt= 0;
  int long_index = 0;
  while ((opt = getopt_long_only(argc, argv,"", 
				 long_options, &long_index )) != -1) {
      switch (opt) {
	// Reading in arguments and setting parameter values accordingly
      case 'f': ca.freqfile = strdup(optarg); break;
      case 'g': ca.glbinfile = strdup(optarg); break;
      case 'G': ca.glbeaglefile = strdup(optarg); break;
      case 'n': ca.nInd = atoi(optarg); break;
      case 'a': ca.p.pair[0] = atoi(optarg); break;
      case 'b': ca.p.pair[1] = atoi(optarg); break;
      case 'c': ca.p.calcA = atoi(optarg); break;
      case 'A': ca.p.a = atof(optarg); break;
      case 'h': alim[0] = atof(optarg); break;
      case 'j': alim[1] = atof(optarg); break;
      case 'x': ca.p.k0 = atof(optarg); break;
      case 'y': ca.p.k1 = atof(optarg); break;
      case 'w': ca.p.k2 = atof(optarg); break;
      case 'v': ca.fixk2to0 = atoi(optarg); ca.p.k2 = 0; break;
      case 'l': ca.minMaf = atof(optarg); break;
      case 's': ca.switchMaf = atoi(optarg); break;
      case 'r': ca.seed = atoi(optarg); break;
      case 'N': ca.nOpti = atoi(optarg); break;
      case 'O': ca.outname = strdup(optarg); break;
      case 'o': runoldrelateV = atoi(optarg); break; // Remove when done?
      default: return 0;//{fprintf(stderr,"unknown arg:\n");return 0;}
	print_info(stderr);
      }
  }
      
  // Start log file and write command to screen
  int nooutname = 0;
  if(ca.outname==NULL){
    nooutname = 1;
    if(ca.glbeaglefile!=NULL)
      ca.outname = ca.glbeaglefile;
    else
      ca.outname = ca.glbinfile;
  }

  std::vector<char *> dumpedFiles;
  FILE *flog=openFile(ca.outname,".log",dumpedFiles);
  fprintf(stderr,"\n\t-> You are using LocalNgsRelate version 0.999 (build time: %s:%s)\n",__DATE__,__TIME__);
  fprintf(stderr,"\t-> Command running is:\n\t   ");
  fprintf(flog,"\t-> You are using LocalNgsRelate version 0.999 (build time: %s:%s)\n",__DATE__,__TIME__);
  fprintf(flog,"\t-> Command running is:\n\t   ");
  for(int i=0;i<argc;i++){
    fprintf(flog,"%s ",argv[i]);
    fprintf(stderr,"%s ",argv[i]);
    if(i==6||i==12){
      fprintf(stderr,"\n\t   ");
    }
  }
  fprintf(stderr,"\n");
  fprintf(flog,"\n");
  if(nooutname){
     fprintf(stderr,"\t->NB. Since no output prefix was provided %s will be used as prefix for all output\n\n",ca.outname);
  }
  
  // Check if all necessary arguments have been provided
  bool stop = false;

  if(ca.p.a!=-1 && ca.p.calcA==1){
    fprintf(stderr,"\n\t## Error: you can't supply both -fixA and -calcA \n");
    stop=true;
  }
  if(ca.glbeaglefile==NULL&&ca.glbinfile==NULL){
    fprintf(stderr,"\n\t## Error: please supply a genotype likelihood file using the option -gbeagle\n");
    stop=true;
  }else if(ca.glbeaglefile!=NULL&&ca.glbinfile!=NULL){
    fprintf(stderr,"\n\t## Error: please supply only 1 input data file using -gbeagle\n");
    stop=true;
  }
  if(ca.freqfile==NULL){
    fprintf(stderr,"\n\t## Error: please supply both a genotype likelihood file and a frequency file\n");
    stop=true;
  }
  if(ca.nInd<1){
    fprintf(stderr,"\n\t## Error: number of individuals in the beagle file needs to be specified and it has to be positive\n");
    exit(0);;
  }
  if(ca.p.pair[0]<0){
    fprintf(stderr,"\n\t## Error: the specified -a value cannot be negative\n");
    stop=true;
  }
  if((ca.p.pair[0]>(ca.nInd-1))){
    fprintf(stderr,"\n\t## Error: the specified -a value is incompatible with the number of individuals specified\n");
    stop=true;
  }
  if((ca.p.pair[1]>(ca.nInd-1))){
    fprintf(stderr,"\n\t## Error: the specified -b value is incompatible with the number of individuals specified\n");
    stop=true;
  }
  if(ca.p.pair[1]<0){
    fprintf(stderr,"\n\t## Error: the specified -b value cannot be negative\n");
    stop=true;
  }
  if((ca.p.a > 0)&&((ca.p.k0<0)||(ca.p.k1<0)||(ca.p.k2<0))){
    fprintf(stderr,"\n\t## Error: you cannot use the option -fixA without also using -k0, -k1 and -k2\n");
    stop=true;
  }
  if((ca.p.a < 0)&&(ca.p.calcA<0)&&((ca.p.k0>-1)&&(ca.p.k1>-1)&&(ca.p.k2>-1))){
    fprintf(stderr,"\n\t## Error: you cannot use the options -k0, -k1 and -k2 without either using -fixA or -calcA\n");
    stop=true;
  }
  if((ca.p.calcA>0)&&((ca.p.k0==1)&&(ca.p.k1==0)&&(ca.p.k2==0))){
    fprintf(stderr,"\n\t## Error: if you set -k0 1 aCalc cannot be used (only works for related pairs)\n");
    stop=true;
  }

  // Check if user provided values are valid
  if((ca.p.k0>-1)||(ca.p.k1>-1)||(ca.p.k2>-1)){
    if(((ca.p.k0<0)||(ca.p.k1<0)||(ca.p.k2<0))&ca.fixk2to0==0){
      fprintf(stderr,"\n\t## Error: -k0, -k1 and -k2 must all be specified together, values have to be non-negative and sum to 1 (perhaps consider --fixk2to0)\n");
      stop=true;
    }else{
    double sumoffixed = 0;
    if(ca.p.k0>0)
      if(ca.p.k0>1){
	fprintf(stderr,"\n\t## Error: provided k0 value was above 1 (%f)\n",ca.p.k0);
	stop=true;
      }else{
	sumoffixed+=ca.p.k0;
      }
    if(ca.p.k1>0)
      if(ca.p.k1>1){
	fprintf(stderr,"\n\t## Error: provided k1 value was above 1 (%f)\n",ca.p.k1);
	stop=true;
      }else{
	sumoffixed+=ca.p.k1;
      }
    if(ca.p.k2>0)
      if(ca.p.k2>1){
	fprintf(stderr,"\n\t## Error: provided k2 value was above 1 (%f)\n",ca.p.k2);
	stop=true;
      }else{
	sumoffixed+=ca.p.k2;
      }
    if((sumoffixed!=1)&ca.fixk2to0==0){
      fprintf(stderr,"\n\t## Error: provided k values do not sum to 1\n");
      stop=true;
    }
    }
    if((ca.p.a > -1)&&!(ca.p.a>0)){
      fprintf(stderr,"\n\t## Error: provided a is not above 0 as required\n");
      stop=true;
    }
  }
  if(!(alim[0]>0)){
    fprintf(stderr,"\n\t## Error: provided minA value is not above 0 as it should be\n");
    stop=true;
  }
  if(!(alim[1]>0)){
    fprintf(stderr,"\n\t## Error: provided maxA value is not above 0 as it should be\n");
    stop=true;
  }
  if(alim[0]>alim[1]){
    fprintf(stderr,"\n\t## Error: provided maxA value is smaller than the provided minA\n");
    stop=true;
  }
  if(stop)
    exit(0);
  

  // Read in freqs
  std::vector<double> freq;
  int nSitesInFreqfile = readFrequencyFile(ca.freqfile,freq);
  if((nColInFile(ca.freqfile)!=1)){
    fprintf(stderr, "\n\t## Error: there is more than one column (%d) in frequency file: %s\n", nColInFile(ca.freqfile), ca.freqfile);
    exit(0);
  }
  
  // Read in gls
  std::vector<perChr> pd;
  if(ca.glbeaglefile!=NULL){
    bgl d=readBeagle(ca.glbeaglefile);
    if(d.nSites!=nSitesInFreqfile){
      fprintf(stderr,"\n\t## Error: The number of sites in the beagle file does not fit with number of sites in freqfile: %d vs. %d\n",d.nSites,nSitesInFreqfile);
      exit(0);
    }
    if(d.nInd!=ca.nInd){
      fprintf(stderr,"\n\t## Error: The number of individuals in the beagle file does not fit with number of individuals specified with -n: %d vs. %d\n",d.nInd,ca.nInd);
      exit(0);
    }
    pd = makeDat(d,freq,ca.minMaf,ca.switchMaf);
  }else{
    printf("\n\t## Error the -g (using binary glf files) is not implemented yet\n");
    exit(0);
  }
  

  // Print args to a log file 
  genome g = mkGenome(pd,ca.p);
  double loglike;

  if(ca.p.k0!=-1&&ca.p.k1!=-1&&ca.p.k2!=-1){
    double ts = ca.p.k0+ca.p.k1+ca.p.k2;
    ca.p.k0 /= ts;
    ca.p.k1 /= ts;
    ca.p.k2 /= ts;

    // If ks are proved by the user and calcA is set: calc a and the corresponding log like
    if(ca.p.calcA==1 && ca.p.a==-1){
      ca.p.a=calculateA(ca.p.k0,ca.p.k1,ca.p.k2,PHI);
      fprintf(stderr,"\t-> Calculating a to ");
      fprintf(stderr,"%f\n",ca.p.a);
      double pars[] = {ca.p.a,ca.p.k0,ca.p.k1,ca.p.k2};
      loglike =  calcLoglikeGenome(pars,g);
    }
    // If ks AND a are user provided: then calc the corresponding log like
    if(ca.p.calcA==-1 && ca.p.a!=-1 && ca.p.a>0.000000000000000000000001){
      fprintf(stderr,"\t-> Setting a to %f\n",ca.p.a);
      double pars[] = {ca.p.a,ca.p.k0,ca.p.k1,ca.p.k2};
      loglike =  calcLoglikeGenome(pars,g);
    }
    // If one or more parameters need to be estimated this is done
  }
  else{
    fprintf(stderr,"\t-> Performing parameter estimation now...\n");
    loglike = doOptim(ca.p,g,pd,ca.p.calcA,ca.seed,ca.nOpti);
  }
  fprintf(stderr,"\t-> Parameters estimated/set to: a=%.16f k0=%.16f k1=%.16f k2=%.16f loglike=%.16f\n",ca.p.a,ca.p.k0,ca.p.k1,ca.p.k2,loglike);
  fprintf(flog,"\t-> Parameters estimated/set to: a=%.16f k0=%.16f k1=%.16f k2=%.16f loglike=%.16f\n",ca.p.a,ca.p.k0,ca.p.k1,ca.p.k2,loglike);
  
  FILE *fres=openFile(ca.outname,".parameters",dumpedFiles);
  fprintf(fres,"%.16f %.16f %.16f %.16f %.16f\n",ca.p.a,ca.p.k0,ca.p.k1,ca.p.k2,loglike);
  fclose(fres);
  
  // Now do the forward, backward, posterior decoding and viterbi
  if(std::isinf(-loglike)){
    fprintf(stderr,"\t-> Likelihood is 0 with the current parameter values, so no inference was performed.\n");
  }else{
  fprintf(stderr,"\t-> Performing IBD inference along the genome now...\n");
  double pars[] = {ca.p.a,ca.p.k0,ca.p.k1,ca.p.k2};
  forward_backward_decode_viterbi(pars,g);
  
  // Write results to a file
  gzFile add = openFileGz(ca.outname,".IBDtractinference.gz","wb",dumpedFiles);
  gzprintf(add,"Chr\tPos\tFreq\tViterbi\tPost0\tPost1\tPost2\n");
  for(int i=0;i<pd.size();i++){
    perChr en = pd[i];
    hmmRes to = g.results[i];
    assert(en.nSites==to.nSites);
    for(int j=0;j<en.nSites;j++){
      gzprintf(add,"%s\t%d\t%f\t",en.name,to.pos[j],exp(en.logfreq[j]));
      gzprintf(add,"%d\t",to.viterbi[j]);
      gzprintf(add,"%.16e\t%.16e\t%.16e\n",to.post[0][j],to.post[1][j],to.post[2][j]);
    }
  }
  gzclose(add);

  }

 // Calc log like for unrelated and parent offspring
  double tmppars[] = {ca.p.a, 1,0,0};
  g.uloglike = calcLoglikeGenome(tmppars,g);
  tmppars[1] = 0; tmppars[2] = 1;
  g.pologlike = calcLoglikeGenome(tmppars,g);
  fprintf(stderr,"\t-> With the same parameters loglike for unrelated is: %.16f and loglike for parent offspring is: %.16f\n",g.uloglike,g.pologlike);
  fprintf(flog,"\t-> With the same parameters loglike for unrelated is: %.16f and loglike for parent offspring is: %.16f\n",g.uloglike,g.pologlike);

  fprintf(flog,"\t-> All analyses are done.\n");
  fprintf(stderr,"\t-> All analyses are done.\n");
  fprintf(flog,"\t   [ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(flog,"\t   [ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  fprintf(stderr,"\t   [ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr,"\t   [ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  fclose(flog); 
  fprintf(stderr,"\n\t   Files created are:\n");
  for(int i=0;1&&i<dumpedFiles.size();i++){
    fprintf(stderr,"\t   %s\n",dumpedFiles[i]);
    free(dumpedFiles[i]);
  }


  

  /*
  for(int j=0;j<pd.size();j++){
    fprintf(stderr,"Info about chr %s\n",pd[j].name);
    for(int i=0;i<pd[j].nSites;i++){
      fprintf(stderr,"freq: %f 1-freq:%f pos: %d\n",
	      exp(pd[j].logfreq[i]),exp(pd[j].logqerf[i]),pd[j].pos[i]);
    }
  }
  */
  

  // Get rid of sites with keep==0?
  
  //delete [] freq ;
  
  
  
  return 0;
      
}

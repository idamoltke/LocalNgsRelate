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
#include "version.h"
#include "analysisfunctions.h"
#include "filehandlingfunctions.h"
#include "vcf.h"
#include "parse_args.h"

std::vector<char *> dumpedFiles;

double alim[2];

void print_info(FILE *fp){
  fprintf(fp, "\n\n");
  fprintf(fp, "You are using LocalNgsRelate version %s  (build time: %s:%s)\n",LocalNgsRelate_version,__DATE__,__TIME__);
  fprintf(fp, "\nUsage: ./localngsrelate  [options] \n");
  fprintf(fp, "\nRequired options:\n");
  fprintf(fp, "   -f         <filename>    Name of file with frequencies\n");
  fprintf(fp, "   -B         <filename>    Name of bcffile\n");
  //  fprintf(fp, "   -g         <fileprefix>  Prefix of files with genotype likelihoods and positions (suffix .glf.gz and gls.pos.gz)\n"); # Not implemented
  fprintf(fp, "   -gbeagle   <fileprefix>  Prefix of files with genotype likelihoods (assumed suffix is .beagle.gz)\n");
  fprintf(fp, "   -n         <INT>         Number of samples in genotype likelihood file\n");
  fprintf(fp, "   -a         <INT>         First individual used for analysis? (zero offset). NB if not specified it is set to 0!\n");
  fprintf(fp, "   -b         <INT>         Second individual used for analysis? (zero offset). NB. if not specified it is set to 1!\n\n");

  fprintf(fp, "Optional options:\n");
  fprintf(fp, "   -calcA     <INT>         Calculate the parameter a instead of estimating it (if set to 1, default is 0)\n");
  fprintf(fp, "   -fixA      <FLOAT>       Set parameter a (rate of change between IBD states) to this value. NB only works when k0, k1 and k2 are also provided\n"); 
  fprintf(fp, "   -minA      <FLOAT>       Set lower limit to parameter a (default 0.001). NB. it has to be > 0\n");
  fprintf(fp, "   -maxA      <FLOAT>       Set upper limit to parameter a (default 0.1)\n");
  fprintf(fp, "   -k0        <FLOAT>       Set k0 to this value (must be used with -k1 and -k2)\n");
  fprintf(fp, "   -k1        <FLOAT>       Set k1 to this value (must be used with -k0 and -k2)\n");
  fprintf(fp, "   -k2        <FLOAT>       Set k2 to this value (must be used with -k0 and -k1)\n");
  fprintf(fp, "   -fixk2to0  <INT>         Set k2 to 0 and estimate only k0 and k1\n");
  fprintf(fp, "   -l         <INT>         minor allele frequency filter (default 0.05)\n");
  fprintf(fp, "   -s         <INT>         Should you swich the freq with 1-freq? (default 0)\n");
  fprintf(fp, "   -r         <FLOAT>       Seed for rand (default 20)\n");
  fprintf(fp, "   -N         <INT>         Number of times to start parameter estimation for each pair with random seed (default 1)\n");
  fprintf(fp, "   -O         <fileprefix>  Prefix for name of outputfiles\n");
  fprintf(fp, "   -T         <STRING>      For -B bcf use PL (default) or GT tag\n");
  fprintf(fp, "   -F         <STRING>      For -h vcf use allele frequency TAG e.g. AFngsrelate (default)\n");  
  fprintf(fp, "\n");
  exit(0);
}

int runoldrelateV;


int main(int argc, char **argv){

  // Check if any options are specified if not then print calling instruction
  if(argc==1){// if no arguments, print info on program
    print_info(stderr);     
    return 0;
  }
  
  // Start run time clock
  clock_t t=clock();
  time_t t2=time(NULL);

  cArg ca = get_pars(argc,argv);
  

 
  // Read in gls
  std::vector<perChr> pd;
  if(ca.glbeaglefile!=NULL){
    // Read in freqs
    std::vector<double> freq;
    int nSitesInFreqfile = readFrequencyFile(ca.freqfile,freq);
    if((nColInFile(ca.freqfile)!=1)){
      fprintf(stderr, "\n\t## Error: there is more than one column (%d) in frequency file: %s\n", nColInFile(ca.freqfile), ca.freqfile);
      exit(0);
    }
  
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
  }else if(ca.vcffile!=NULL){
    pd = readbcfvcf_bgl(ca.vcffile,0,ca.minMaf,ca.vcf_format_field,ca.vcf_allele_field,NULL,ca.switchMaf);

  }


  else{
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
  fprintf(ca.flog,"\t-> Parameters estimated/set to: a=%.16f k0=%.16f k1=%.16f k2=%.16f loglike=%.16f\n",ca.p.a,ca.p.k0,ca.p.k1,ca.p.k2,loglike);
  
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
  fprintf(ca.flog,"\t-> With the same parameters loglike for unrelated is: %.16f and loglike for parent offspring is: %.16f\n",g.uloglike,g.pologlike);

  fprintf(ca.flog,"\t-> All analyses are done.\n");
  fprintf(stderr,"\t-> All analyses are done.\n");
  fprintf(ca.flog,"\t   [ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(ca.flog,"\t   [ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  fprintf(stderr,"\t   [ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr,"\t   [ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  fclose(ca.flog); 
  fprintf(stderr,"\n\t   Files created are:\n");
  for(int i=0;1&&i<dumpedFiles.size();i++){
    fprintf(stderr,"\t   %s\n",dumpedFiles[i]);
    free(dumpedFiles[i]);
  }
  
  return 0;
}

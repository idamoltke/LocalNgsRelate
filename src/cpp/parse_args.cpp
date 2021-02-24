#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <zlib.h>
#include "parse_args.h"
#include "filehandlingfunctions.h"
#include "version.h"

extern std::vector<char *> dumpedFiles;

extern double alim[2];

cArg get_pars(int argc,char **argv){
  // Set default parameters
  cArg ca;
  ca.freqfile=NULL;
  ca.glbeaglefile=NULL;
  ca.vcffile=NULL;
  ca.outname=NULL;
  ca.vcf_format_field = strdup("PL"); // can take PL or GT
  ca.vcf_allele_field = strdup("AFngsrelate"); // can take any tag value e.g. AF AF1 etc
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
  ca.runoldrelateV = 0; // Remove?
  
  // Specify the expected options for parsing purposes
  static struct option long_options[] = {
      {"f",       required_argument, 0,  'f' },
      {"beagle"     ,  required_argument, 0, 'g' },
      {"bcffile", required_argument, 0,  'B' },
      {"T",       required_argument, 0,  'T' },
      {"F",       required_argument, 0,  'F' },
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
  
  // If there are any options specified then read them
  int opt= 0;
  int long_index = 0;
  while ((opt = getopt_long(argc, argv,"f:g:B:T:F:n:a:b:c:A:h:j:x:y:w:v:l:s:r:N:O:o:", 
				 long_options, &long_index )) != -1) {
    switch (opt) {
	// Reading in arguments and setting parameter values accordingly
      case 'f': ca.freqfile = strdup(optarg); break;
      case 'g': ca.glbeaglefile = strdup(optarg); break;
      case 'B': ca.vcffile = strdup(optarg); break;
      case 'T': ca.vcf_format_field = strdup(optarg); break;
      case 'F': ca.vcf_allele_field = strdup(optarg); break;	
      case 'n': ca.nInd = atoi(optarg); break;
      case 'a': ca.p.pair[0] = atoi(optarg); break;
      case 'b': ca.p.pair[1] = atoi(optarg); break;
      case 'c': ca.p.calcA = atoi(optarg); break;
      case 'A': ca.p.a = atof(optarg); break;
      case 'h': ca.alim[0] = atof(optarg); break;
      case 'j': ca.alim[1] = atof(optarg); break;
      case 'x': ca.p.k0 = atof(optarg); break;
      case 'y': ca.p.k1 = atof(optarg); break;
      case 'w': ca.p.k2 = atof(optarg); break;
      case 'v': ca.fixk2to0 = atoi(optarg); ca.p.k2 = 0; break;
      case 'l': ca.minMaf = atof(optarg); break;
      case 's': ca.switchMaf = atoi(optarg); break;
      case 'r': ca.seed = atoi(optarg); break;
      case 'N': ca.nOpti = atoi(optarg); break;
      case 'O': ca.outname = strdup(optarg); break;
      case 'o': ca.runoldrelateV = atoi(optarg); break; // Remove when done?
      default: {
	fprintf(stderr,"\t-> Unknown arg: %d %s\n",opt,optarg);
	exit(0);
      }
    }
  }
  fprintf(stderr,"-B = %s\n",ca.vcffile);
  // Start log file and write command to screen
  int nooutname = 0;
  if(ca.outname==NULL){
    nooutname = 1;
    if(ca.glbeaglefile!=NULL)
      ca.outname = ca.glbeaglefile;
    if(ca.vcffile!=NULL)
      ca.outname = ca.vcffile;
  }

  ca.flog = openFile(ca.outname,".log",dumpedFiles);
  fprintf(stderr,"\n\t-> You are using LocalNgsRelate version %s (build time: %s:%s)\n",LocalNgsRelate_version,__DATE__,__TIME__);
  fprintf(stderr,"\t-> Command running is:\n\t   ");
  fprintf(ca.flog,"\t-> You are using LocalNgsRelate version Ts (build time: %s:%s)\n",LocalNgsRelate_version,__DATE__,__TIME__);
  fprintf(ca.flog,"\t-> Command running is:\n\t   ");
  for(int i=0;i<argc;i++){
    fprintf(ca.flog,"%s ",argv[i]);
    fprintf(stderr,"%s ",argv[i]);
    if(i==6||i==12){
      fprintf(stderr,"\n\t   ");
    }
  }
  fprintf(stderr,"\n");
  fprintf(ca.flog,"\n");
  if(nooutname){
     fprintf(stderr,"\t->NB. Since no output prefix was provided %s will be used as prefix for all output\n\n",ca.outname);
  }
  
  // Check if all necessary arguments have been provided
  bool stop = false;

  if(ca.p.a!=-1 && ca.p.calcA==1){
    fprintf(stderr,"\n\t## Error: you can't supply both -fixA and -calcA \n");
    stop=true;
  }
  if(ca.glbeaglefile==NULL&&ca.vcffile==NULL){ //
    fprintf(stderr,"\n\t## Error: please supply a genotype likelihood file using the option -gbeagle or vcffile with -B\n");
    stop=true;
  }else if(ca.glbeaglefile!=NULL&&ca.vcffile!=NULL){
    fprintf(stderr,"\n\t## Error: please supply only 1 input data file using -gbeagle -B\n");
    stop=true;
  }
  if(ca.freqfile==NULL&&ca.vcffile==NULL){
    fprintf(stderr,"\n\t## Error: please supply both a genotype likelihood file and a frequency file\n");
    stop=true;
  }
  if(ca.nInd<1&&ca.vcffile==NULL){
    fprintf(stderr,"\n\t## Error: number of individuals in the beagle file needs to be specified and it has to be positive\n");
    exit(0);;
  }
  if(ca.p.pair[0]<0){
    fprintf(stderr,"\n\t## Error: the specified -a value cannot be negative\n");
    stop=true;
  }
  if(ca.vcffile==NULL){
    if((ca.p.pair[0]>(ca.nInd-1))){
      fprintf(stderr,"\n\t## Error: the specified -a value is incompatible with the number of individuals specified\n");
      stop=true;
    }
    if((ca.p.pair[1]>(ca.nInd-1))){
      fprintf(stderr,"\n\t## Error: the specified -b value is incompatible with the number of individuals specified\n");
      stop=true;
    }
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
    if(((ca.p.k0<0)||(ca.p.k1<0)||(ca.p.k2<0))&&ca.fixk2to0==0){
      fprintf(stderr,"\n\t## Error: -k0, -k1 and -k2 must all be specified together, values have to be non-negative and sum to 1 (perhaps consider --fixk2to0)\n");
      stop=true;
    }else{
    double sumoffixed = 0;
    if(ca.p.k0>0){
      if(ca.p.k0>1){
	fprintf(stderr,"\n\t## Error: provided k0 value was above 1 (%f)\n",ca.p.k0);
	stop=true;
      }else{
	sumoffixed+=ca.p.k0;
      }
    }
    if(ca.p.k1>0){
      if(ca.p.k1>1){
	fprintf(stderr,"\n\t## Error: provided k1 value was above 1 (%f)\n",ca.p.k1);
	stop=true;
      }else{
	sumoffixed+=ca.p.k1;
      }
    }
    if(ca.p.k2>0){
      if(ca.p.k2>1){
	fprintf(stderr,"\n\t## Error: provided k2 value was above 1 (%f)\n",ca.p.k2);
	stop=true;
      }else{
	sumoffixed+=ca.p.k2;
      }
    }
    if((sumoffixed!=1)&&ca.fixk2to0==0){
      fprintf(stderr,"\n\t## Error: provided k values do not sum to 1\n");
      stop=true;
    }
    }
    if((ca.p.a > -1)&&!(ca.p.a>0)){
      fprintf(stderr,"\n\t## Error: provided a is not above 0 as required\n");
      stop=true;
    }
  }
  if(!(ca.alim[0]>0)){
    fprintf(stderr,"\n\t## Error: provided minA value is not above 0 as it should be\n");
    stop=true;
  }
  if(!(ca.alim[1]>0)){
    fprintf(stderr,"\n\t## Error: provided maxA value is not above 0 as it should be\n");
    stop=true;
  }
  if(ca.alim[0]>ca.alim[1]){
    fprintf(stderr,"\n\t## Error: provided maxA value is smaller than the provided minA\n");
    stop=true;
  }
  if(stop)
    exit(0);
  return ca;
}

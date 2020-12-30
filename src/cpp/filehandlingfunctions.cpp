#include <cstring>
#include <map>
#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <zlib.h>
#include <vector>
#include <cmath>
#include <libgen.h>
#include <sys/stat.h>  
#include "analysisfunctions.h"
#include "filehandlingfunctions.h"

#define LENS  4096


int SIG_COND =1;//if we catch signal then quit program nicely


/*
  Check number of columns in a file
*/

int nColInFile(const char *fname){
  const char *delims = "\t \n";
  FILE* fp = fopen(fname, "r");
  if (fp == NULL)
    exit(EXIT_FAILURE);

  char* line = NULL;
  size_t len = 0;
  int l = getline(&line, &len, fp);
  if (l == -1){
    fprintf(stderr, "\n Something is wrong with %s\n\n", fname);
    exit(0);
  }
  int ncol=1;
  strtok(line,delims);
  while(strtok(NULL,delims)){
    ncol++;
  }

  fclose(fp);

  if (line)
    free(line);

  return ncol;
}

 
/*
  Reads frequencies from the file "fname" into the vector ret
*/

size_t readFrequencyFile(const char *fname,std::vector<double> &ret) {
  assert(ret.size()==0);
  gzFile gz = Z_NULL;
  if (((gz = gzopen(fname, "r"))) == Z_NULL) {
    fprintf(stderr, "\n\t## Error: problem opening the file %s\n",fname);
    exit(0);
  }
  int nbytes=10000;
  char *buf = new char[nbytes];
  while (gzgets(gz, buf, nbytes)) {
    // reads line by line
    ret.push_back(atof(buf));
    if(atof(buf)<0){
      fprintf(stderr, "\n\t## Error: problem with content of frequency file %s. At least 1 frequency is below 0!\n",fname);
      exit(0);
    }
    if(atof(buf)>1){
      //fprintf(stderr, "[%s]\t-> Problem with content of frequency file %s: at least 1 frequency is above 1!\n",__FUNCTION__, fname);
      fprintf(stderr, "\n\t## Error: problem with content of frequency file %s. At least 1 frequency is above 1!\n",fname);
      exit(0);
    }
  }
  fprintf(stderr, "\t-> Frequency file \'%s\' successfully parsed.\n\t   It contained frequencies from %lu sites\n",
          fname, ret.size());
  gzclose(gz);
  delete[] buf;
  return ret.size();
}


/*
  Returns the bgl struct containing all data from a beagle file.

  It finds the nsamples from counting the header
  It finds the number of sites by quering every line in a std::vector
  After the file has been read in total it reloops over the lines in the vector and parses data
 */

bgl readBeagle(const char* fname) {
  const char *delims = "\t \n";
  char *suffix = strdup(".beagle.gz");
  char *fullfname = new char[strlen(fname)+strlen(suffix)+1];
  strcpy(fullfname,fname);
  strncat(fullfname,suffix,strlen(suffix));
  gzFile fp = NULL;
  if(Z_NULL==(fp=gzopen(fullfname,"r"))){
    fprintf(stderr,"\n\t## Error: problem opening the file %s\n",fullfname);
    exit(0);
  }
  
  bgl ret;
  char buf[LENS];

  //find number of columns
  gzgets(fp,buf,LENS);
  strtok(buf,delims);
  int ncols=1;
  while(strtok(NULL,delims))
    ncols++;
  if(0!=(ncols-3) %3 ){
    fprintf(stderr,"\n\t##Problem parsing beagle file: ncols=%d will exit\n",ncols);
    exit(0);
  }
  ret.nInd = (ncols-3)/3;//this is the number of samples
  
  //read every line into a vector
  std::vector<char*> tmp;
  while(gzgets(fp,buf,LENS))
    tmp.push_back(strdup(buf));
  
  //now we know the number of sites
  ret.nSites=tmp.size();
  ret.chr = new char*[ret.nSites];
  ret.pos = new int[ret.nSites];
  ret.genos= new double*[ret.nSites];

  //then loop over the vector and parsing every line
  for(int s=0;SIG_COND&& (s<ret.nSites);s++){
    ret.chr[s] = strdup(strtok(tmp[s],"_"));
    ret.pos[s] = atoi(strtok(NULL,delims));
    strtok(NULL,delims);//major
    strtok(NULL,delims);//minor
    ret.genos[s]= new double[3*ret.nInd];
    for(int i=0;i<ret.nInd*3;i++){
      ret.genos[s][i] = atof(strtok(NULL,delims));
      if(ret.genos[s][i]<0){
	fprintf(stderr,"Likelihoods must be positive\n");
	fprintf(stderr,"site %d ind %d geno %d has value %f\n",s,int(i*1.0/3),i%3,ret.genos[s][i]);
	exit(0);
      }
    }
    for(int i=0;i<ret.nInd;i++){
      double tmpS = 0.0;
      for(int g=0;g<3;g++)
	tmpS += ret.genos[s][i*3+g];
      if(!(tmpS>0)){
	fprintf(stderr,"The sum of likelihoods for a genotypes must be positive\n");
	fprintf(stderr,"individual %d site %d has sum %f\n",i,s,tmpS);
	exit(0);
      } 
      for(int g=0;1&&g<3;g++)
	ret.genos[s][i*3+g] = log(ret.genos[s][i*3+g]);
    }
    free(tmp[s]);
  }
  
  fprintf(stderr, "\t-> Beagle file \'%s\' successfully parsed.\n\t   It contained data from  %d individuals at %d sites.\n",
	  fullfname,ret.nInd,ret.nSites);

  gzclose(fp); //clean up filepointer
  return ret;
}

int fexists(const char* str){///@param str Filename given as a string.
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}


gzFile getGz(const char*fname,const char* mode){

  //  fprintf(stderr,"\t-> opening: %s\n",fname);
  gzFile fp=Z_NULL;
  if(NULL==(fp=gzopen(fname,mode))){
    fprintf(stderr,"\n\t## Error opening gzFile handle for the file %s\n",fname);
    exit(0);
  }
  return fp;
}

gzFile openFileGz(const char* a,const char* b,const char *mode,std::vector<char *>& dumpedFiles){
  if(0)
    fprintf(stderr,"[%s] %s%s\n",__FUNCTION__,a,b);
  char *c = new char[strlen(a)+strlen(b)+1];
  strcpy(c,a);
  strncat(c,b,strlen(b));
  //  fprintf(stderr,"\t-> Dumping file: %s\n",c);
  dumpedFiles.push_back(strdup(c));
  gzFile fp = getGz(c,mode);
  delete [] c;
  return fp;
}

FILE *openFile(const char* a,const char* b,std::vector<char *>& dumpedFiles){
  if(0)
    fprintf(stderr,"[%s] %s %s",__FUNCTION__,a,b);
  char *c = new char[strlen(a)+strlen(b)+1];
  strcpy(c,a);
  strncat(c,b,strlen(b));
  //  fprintf(stderr,"\t-> Dumping file: %s\n",c);
  if(0&&fexists(c)){//ANDERS DAEMON DRAGON HATES THIS
    fprintf(stderr,"File: %s exists will exist\n",c);
    fflush(stderr);
    exit(0);
  }
  dumpedFiles.push_back(strdup(c));
  FILE *fp = fopen(c,"w");
  delete [] c;
  return fp;
}



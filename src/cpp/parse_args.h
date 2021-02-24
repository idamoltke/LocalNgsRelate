#pragma once
#include "analysisfunctions.h"

typedef struct{
  char *freqfile;
  char *glbeaglefile;
  char *vcffile;
  char *outname;
  char *vcf_format_field;
  char *vcf_allele_field;
  int fixk2to0;
  int nInd;
  para p;
  float minMaf;
  int switchMaf;
  int seed;
  int nOpti;
  FILE *flog;
  int runoldrelateV;
}cArg;

cArg get_pars(int argc,char **argv);

void print_carg(FILE *fp,cArg ca);

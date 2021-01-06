# LocalNgsRelate

## Brief description
This page contains information about the program called LocalNgsRelate, which can be used to infer IBD sharing along the genomes of two individuals from low-depth Next Generation Sequencing (NGS) data by using genotype likelihoods instead of called genotypes. To be able to infer the IBD sharing you will need to know the population allele frequencies and have genotype likelihoods. This can be obtained e.g. using the program ANGSD as shown in example 1 below. For more information about ANGSD see here: http://popgen.dk/angsd/index.php/Quick_Start. 

## How to download and install
On a linux or mac system with git and g++ installed LocalNgsRelate can be downloaded and installed as follows:

```
git clone https://github.com/idamoltke/LocalNgsRelate
cd LocalNgsRelate/src/cpp/;make
```

## Help and run options
To see how the program can be run and to see all options, you can type the following command in the folder where the src code is:

```
./localngsrelate 
```
Or from a different folder with full path "path/" to the src code folder

```
path/localngsrelate 
```

## Input file format
### Formal description
LocalNgsRelate takes two files as input: a file with genotype likelihoods (-gbeagle) and a file with population allele frequencies (-f) for the sites there are genotype likelihoods for. The genotype likelihood file needs to contain a line for each site with 3 values for each individual (one log transformed genotype likelihood for each of the 3 possible genotypes encoded as 'double's) and it needs to be in beagle format and gz compressed (see e.g. http://www.popgen.dk/angsd/index.php/Beagle_input). 
The frequency file needs to contain a line per site with the allele frequency of the site in it. 
For examples of the two types of input files see the files in the folder exampledata which are described below.


### Example input files
The example files included here (in the folder exampledata) are made based on LWK 1000G low-depth samples. There are two files:

1) LWK.beagle.gz  a file with genotype likelihoods for 101 individuals at 130709 sites
2) LWK.freq       a file with allele frequency estimates for these sites (based on data from 101 individuals)

Among the LWK individuals there at least 4 pairs of relatives:



| Index ind1 | Index ind2   | 1000G IDs ind1 | 1000G IDs ind2   |Relationships |
| ------------ | ------------ | ------------ | ------------------------ | :------------- | 
|   0   |       1    | "NA19027" | "NA19042"       | Half-sibs |
|   2   |  3    | "NA19313" | "NA19331"       | Parent-offspring |
|   3   |  4    | "NA19331" | "NA19334"       | Full siblings |
|   5   |  6    | "NA19451" | "NA19452"       | First cousins |

For a description of how this dataset was made see "Making input data" below.

## Output format
Successfully running the program should lead to 3 output files and if run the program with "–o exampleoutput" these will be called

1) exampleoutput.log
2) exampleoutput.parameters
3) exampleoutput.IBDtractinference.gz

The .log file just logs the command called.
The .parameter file contains the obtained ML extimates of a, k0, k1 and k2 followed by the log likelihood for those parameter values.
The IBDtractinference.gs file contains the IBD tract inference results of the following format

```
Chr     Pos     Freq            Viterbi Post0                   Post1                    Post2
1       30860   0.184148        1       4.7034837882467595e-01  5.2965162117352105e-01  0.0000000000000000e+00
1       52238   0.362546        1       4.6995311651972155e-01  5.3004688347834739e-01  0.0000000000000000e+00
1       73822   0.058160        1       4.6941249114032702e-01  5.3058750885676254e-01  0.0000000000000000e+00
1       94986   0.073169        1       4.6888250067835235e-01  5.3111749931621532e-01  0.0000000000000000e+00
1       118617  0.318185        1       4.6829086584325724e-01  5.3170913415127108e-01  0.0000000000000000e+00
```

If you as shown below run the program with “2>&1 | tee exampleoutput.fulllog” added to the end of your command your will get an additional output file (exampleoutput.fulllog), which contains everything that would otherwise be written to the screen. This can be helpful if you want to assess convergence (for details see “Run examples” below). To get this to work you need the program tee installed.


## Run examples

### Making input data
NB We describe how the data in the folder exampledata was made as an example of how input data to LocalNgsRelate can be made. Running the commands your self will take som time so feel free to skip that and simply jump to the description of how to run LocalNgsRelate below.

To make the example input files in exampldata we first downloaded bamfiles for 101 LWK samples sequenced to about 6x as follows:

xxx

Then we ran the program ANGSD to get a genotype likelihood files as well as frequency estimates:

```
./ANGSD -b LWK.bamlist -GL 1 -doMajorMinor 3 -doMaf 1 -sites LWK.sites -rf LWK.chrs -out LWK -minMapQ 30 -minQ 20 -P 4
```

In that command: 

1) LWK.bamlist is a file that contains a line for each bamfile and this line contains the full path of the given bamfile on your system

2) LWK.sites contains a line for each site to include in the dataset and this line contains 4 pieces of information tab separated, namely chromosome number, position, major allele and minor allele. Note that this needs to be indexed by ANGSD for details see http://www.popgen.dk/angsd/index.php/Sites     

3) LWK.chr contains a line for each chromosome the sites file have sites on. So if the sites file contains sites on all 22 autosomes it will contain 22 lines, one for each chromosome. 

For details on what the different options in ANGSD mean see http://www.popgen.dk/angsd/index.php/Genotype_Likelihoods and http://www.popgen.dk/angsd/index.php/Allele_Frequencies. 

The above command will create a file LWK.beagle.gz and a file LWK.mafs.gz. The first one contains genotype likelihoods and is in exactly the format that is needed for LocalNgsRelate. The other one, which contains allele frequencies, needs to be altered a bit to work as an input file for LocalNgsRelate. But all we did to make the final input file fromit was to run the following command:

```
zcat LWK.mafs.gz | cut -f5 |sed 1d > LWK.freq
```


###	Analysing selected pairs
Let’s here try to infer IBD tracks for a few pairs of individuals in the example dataset. We can do this from the command line by first creating a folder for the results and moving into this folder (from the folder called src):

```
mkdir ../../exampleruns
cd ../../exampleruns
```

And then we can run the following commands to analyse the 4 listed related pairs (see “Example input files”):

```
../src/cpp/localngsrelate -a 0 -b 1 -gbeagle ../exampledata/LWK -f ../exampledata/LWK.freq -n 101 -fixk2to0 1 -O exampleoutput_ind0_ind1 2>&1 | tee exampleoutput_ind0_ind1.fulllog
../src/cpp/localngsrelate -a 2 -b 3 -gbeagle ../exampledata/LWK -f ../exampledata/LWK.freq -n 101 -fixk2to0 1 -O exampleoutput_ind2_ind3 2>&1 | tee exampleoutput_ind2_ind3.fulllog
../src/cpp/localngsrelate -a 5 -b 6 -gbeagle ../exampledata/LWK -f ../exampledata/LWK.freq -n 101 -fixk2to0 1 -O exampleoutput_ind5_ind6 2>&1 | tee exampleoutput_ind5_ind6.fulllog
../src/cpp/localngsrelate -a 3 -b 4 -gbeagle ../exampledata/LWK -f ../exampledata/LWK.freq -n 101 -O exampleoutput_ind3_ind4 2>&1 | tee exampleoutput_ind3_ind4.fulllog
```

Note that we here specify the input files with –gbeagle and –f, the number of individuals in the data with –n and the output name prefix with the option –O. Also, note that we are here running the program with default values except for in the first 3 runs where we set the option -fixk2to0 1 since we only expect k2 to be above 0 for full siblings or twins. Finally note that we ended each command line with “2>&1 | tee exampleoutput_indx_indy.fulllog” where x and y are the indices of the analysed individuals. This produces a long log file which can be used to assess convergence (see below for how).


### Assessing convergence
If you ran the program as above with “2>&1 | tee exampleoutput_indx_indy.fulllog” you can use the script getlikes.sh in the folder src/scripts like this. E.g. for for individuals 0 and 1 we can run the command:
 
```
../src/scripts/getlikes.sh exampleoutput_ind0_ind1.fulllog 
```

This should produce a sorted list of log likelihoods for the 20 rounds of optimization that was run (20 is the default number) like this:

```
17	-178820.6488210122624878
13	-178820.6488210122915916
7	-178820.6488210122915916
16	-178820.6488210123206954
12	-178820.6488210123497993
19	-178820.6488210123497993
8	-178820.6488210123497993
15	-178820.6488210123789031
0	-178820.6488210124080069
3	-178820.6488210124080069
6	-178820.6488210124080069
10	-178820.6488210124371108
14	-178820.6488210124371108
18	-178820.6488210124371108
2	-178820.6488210124371108
11	-178820.6488210124662146
1	-178820.6488210124662146
5	-178820.6488210124662146
9	-178820.6488210124953184
4	-178820.6488210125535261
```

As can be seen the top 5 log likelihoods differ very little suggesting that convergence was likely reached (and thus the maximum likelihood solution was found by the program). If this had not beent he case you might consider rerunning the program with a higher number of optimization rounds, say 100 (by rerunning the same command as above but with “–N 100” added).

```
### Plotting the results
We can plot the final results opening R and using a script in the folder src/scripts called plotLocalNgsRelateOutput.R. E.g. we can plot the results for the full sibling pair on chromosome 2 as follows:


# Read in plotting script
source("../src/scripts/plotLocalNgsRelateOutput.R")

## Read in results        
nam = "exampleoutput_ind3_ind4"                                                 
res = read.table(paste(nam,".IBDtractinference.gz",sep=""),header=T)
      
## Plot posterior 
pdf(paste(nam,"_posteriorandviterbi_chr2.pdf",sep=""),h=4,w=16)                                               
plot.localngsrelate.posterior(res,chr=2)
plot.localngsrelate.viterbi(res,chr=2)
graphics.off()               
```

Or the results for the parent offspring pair on all chromosomes 

```
## Read in plotting script
source("../src/scripts/plotLocalNgsRelateOutput.R")

## Read in results        
nam = "exampleoutput_ind2_ind3"                                                 
res = read.table(paste(nam,".IBDtractinference.gz",sep=""),header=T)
      
## Plot posterior 
pdf(paste(nam,"_posteriorandviterbi_allchrs.pdf",sep=""),h=4,w=16)                                               
plot.localngsrelate.posterior(res)
plot.localngsrelate.viterbi(res)
graphics.off()               

```
The results can be seen in here:

Plots for chromosome 2 (ind 3 and 4):  
http://popgen.dk/ida/LocalNgsRelate/exampleoutput_ind3_ind4_posteriorandviterbi_chr2.pdf

Plots for all chromosomes (ind 2 and 3):  
http://popgen.dk/ida/LocalNgsRelate/exampleoutput_ind2_ind3_posteriorandviterbi_allchrs.pdf


## Citing and references
For questions contact: ida@binf.ku.dk

If you use the program in a publication please cite: xxx


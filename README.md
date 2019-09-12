# GEM

Building DockerFile 
## Execution
### Steps to be followed - 
#### Building and running docker images
Execute in your terminal - 
 - Clone this git repository to your local machine <br> 
 ```
   git clone https://github.com/large-scale-gxe-methods/GEM.git
 ```
 - Change your working directory to GEM <br>
 ```
   cd GEM/
 ```
 - Build the base docker image <br> 
 ```
   docker build -t large-scale-gxe-methods/gem:v0.7 . 
 ```
 
#### Testing installation
 - Run GEM within a contrainter with the example files <br> 
 ```
$ cat example/GEM_Input.param
SAMPLE_ID_HEADER  // Column header name of sample ID in pheno data file
sampleid
PHENOTYPE  // 0 for continous data using linear regression, 1 for binary data using logistic regression
1
PHENO_HEADER  // column header name of phenotype data in pheno data file
pheno2
COVARIATES_HEADERS  // column header names of the slected covariates in the pheno data file
cov1 cov2 cov3
MISSING // missing value key of phenotype file
NaN
ROBUST  // 0 for non-robust analysis, 1 for robust analysis
1
STREAM_SNPS  // SNP numbers for each GWAS analysis
20
NUM_OF_INTER_COVARIATE  // Number of interactive covariate data columns with Geno (1 means the first column of covariate data, 2 means the first two selected columns, ...)
1
LOGISTIC_CONVERG_TOL  // Convergence tolerance for logistic regression
0.000001
DELIMINATOR
,
GENO_FILE_PATH  // Path of geno file
/GEMexample/example.bgen
PHENO_FILE_PATH  // Path of pheno data file
/GEMexample/example.pheno
SAMPLE_FILE_PATH // Path of sample file. Executed only when bgen file does not provide sample order. Format follows UK BioBank.

OUTPUT_PATH // Path of output file
/GEMexample/my_example.out

 $ docker run --rm -v /PATH/TO/RESPOSITORY/large-scale-gxe-methods/GEM/example/:/GEMexample large-scale-gxe-methods/gem:v0.7 /GEM/GEM /GEMexample/GEM_Input.param
*********************************************************
Number of command-line arguments: 2
Parameter input file is: /GEMexample/GEM_Input.param
*********************************************************
The Selected Phenotype Data Header Name is: pheno2
Linear or Binary? Binary
Robust or Non-Robust Analysis? Robust
The Total Number of Selected Covariates is: 3
The Selected Covariate Column Header Names are:  cov1   cov2   cov3
*********************************************************
Before ID Matching and checking missing values:
Size of the pheno vector is: 500 X 1
Size of the selected covariate matrix (including first column for interception values) is: 500 X 4
End of reading pheno and covariate data.
*********************************************************
General information of bgen file.
offset: 6028
L_H:
BGEN snpBlocks (Mbgen): 1000
BGEN samples (Nbgen): 500
CompressedSNPBlocks: 1
Layout: 2
SampleIdentifiers: 1
****************************************************************************
After processes of sample IDMatching and checking missing values, the sample size changes from 500 to 250.
****************************************************************************
Sample IDMatching and checking missing values processes have been completed.
New pheno and covariate data vectors with the same order of sample ID sequence of geno data are updated.
****************************************************************************
Starting GWAS.
Logistic regression reaches convergence after 5 steps.
*********************************************************
Streaming SNPs for speeding up GWAS analysis in parallel.
SNP numbers for each GWAS analysis are : 20
*********************************************************
Total Wall Time = 0.221223  Seconds
Total CPU Time = 0.278795  Seconds
Execution Wall Time = 0.0424184  Seconds
*********************************************************
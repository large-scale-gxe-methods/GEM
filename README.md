# GEM  

GEM (Gene-Environment interaction analysis for Millions of samples) is a software program for large-scale gene-environment interaction testing in samples from unrelated individuals. It enables genome-wide association studies in up to millions of samples while allowing for multiple exposures, control for genotype-covariate interactions, and robust inference.  


<br />
Current version: 1.4   

<br />
Additional documentation:  
https://large-scale-gxe-methods.github.io/GEM-website/index.html  

<br />

## Contents
- [Installation](#installation)
- [Usage](#usage)
- [Contact](#contact)
- [License](#license)

<br />

## Installation  
Library Dependencies:  
* C++11 compiler.  
* BLAS/LAPACK. For Intel processors, we recommend that GEM is compiled with an optimized math routine library such as the Math Kernal Library for top performance.  
* Boost C++ libraries. GEM links the following Boost libraries:  ```boost_program_options, boost_thread, boost_system, and boost_filesystem```.  

<br />

To install GEM, run the following lines of code.
 ```
 git clone https://github.com/large-scale-gxe-methods/GEM
 cd GEM/
 cd src/  
 make  
 ```

<br />
<br />
<br />

## Usage

### Running GEM

1. [Command Line Options](#command-line-options)  
1. [Input Files](#input-files)
1. [Output File Format](#output-file-format)
1. [Examples](#examples)

<br /> 
<br />

### Command Line Options

Once GEM is installed, the executable ```./GEM``` can be used to run the program.  
For a list of options, use ```./GEM --help```.  

<details>
     <summary> <b>List of Options</b> </summary>

```
General Options:

--help
    Prints the available options of GEM and exits.  
   
--version
    Prints the version of GEM and exits.
  
  
  
Input/Output File Options:  

--pheno-file  
     Path to the phenotype file.  

--bgen  
     Path to the BGEN file.  

--sample  
     Path to the sample file.  
     Required when the BGEN file does not contain sample identifiers.  
  
--pfile  
     Path and prefix to the .pgen, .pvar, and .psam files.  
     If this flag is used, then --pgen/--pvar/--psam don't need to be specified.
  
--pgen  
     Path to the pgen file.  
  
--pvar  
     Path to the pvar file.  
  
--psam  
     Path to the psam file.  
  
--bfile  
     Path and prefix to the .bed, .bim and .fam files.  
     If this flag is used, then --bed/--bim/--fam don't need to be specified.
  
--bed  
     Path to the bed file.  
  
--bim  
     Path to the bim file.  
  
--fam  
     Path to the fam file.  
  
--out  
     Full path and extension to where GEM output results.  
     Default: gem.out
  
--output-style  
     Modifies the output of GEM. Must be one of the following:
	minimum: Output the summary statistics for only the GxE and marginal G terms.
        meta: 'minimum' output plus additional fields for the main G and any GxCovariate terms.
               For a robust analysis, additional columns for the model-based summary statistics will be included.
        full: 'meta' output plus additional fields needed for re-analyses of a subset of interactions.
	Default: minimum   



Phenotype File Options:

--sampleid-name  
     Column name in the phenotype file that contains sample identifiers.  

--pheno-name  
     Column name in the phenotype file that contains the phenotype of interest. 
     If the number of levels (unique observations) is 2, the phenotype is treated as binary; otherwise it is assumed to be continuous.

--exposure-names  
     One or more column names in the phenotype file naming the exposure(s) to be included in interaction tests.  
     If no exposures are included, GEM will only perform the marginal test.  

--int-covar-names  
     Any column names in the phenotype file naming the covariate(s) for which interactions should be included 
     for adjustment (do not include with --exposure-names).  

--covar-names  
     Any column names in the phenotype file naming the covariates for which only main effects should be included
     for adjustment (do not include with --exposure-names or --int-covar-names).  

--robust
     0 for model-based standard errors and 1 for robust standard errors.
        Default: 0  

--tol 
     Convergence tolerance for logistic regression.  
        Default: 0.0000001  

--delim  
     Delimiter separating values in the phenotype file. Tab delimiter should be represented as \t and space delimiter as \0.
        Default: , (comma-separated)  

--missing-value  
     Indicates how missing values in the phenotype file are stored.  
        Default: NA
  
--center 
     0 for no centering to be done and 1 to center ALL exposures and covariates.
        Default: 1

--scale
     0 for no scaling to be done and 1 to scale ALL exposures and covariates by the standard deviation.
        Default: 0

--categorical-names
     Names of the exposure or interaction covariate that should be treated as categorical.
        Default: None

--cat-threshold
     A cut-off to determine which exposure or interaction covariate not specified using --categorical-names 
     should be automatically treated as categorical based on the number of levels (unique observations).
        Default: 20
   


Filtering Options:  

--maf
     Minimum threshold value [0, 0.5] to exclude variants based on the minor allele frequency.
        Default: 0.001
  
--miss-geno-cutoff
     Maximum threshold value [0, 1.0] to filter variants based on the missing genotype rate.  
        Default: 0.05  
  
--include-snp-file  
     Path to file containing a subset of variants in the specified genotype file to be used for analysis. 
     The first line in this file is the header that specifies which variant identifier in the genotype file  
     is used for ID matching. This must be 'snpid' (PLINK or BGEN) or 'rsid' (BGEN only). There should be one variant 
     identifier per line after the header.
  
  

Performance Options:  

--threads
     Set number of compute threads.
    	  Default: ceiling(detected threads / 2)  

--stream-snps 
     Number of SNPs to analyze in a batch. Memory consumption will increase for larger values of stream-snps.  
        Default: 1

```
</details>

<br /> 

### Input Files

* #### Phenotype File
    A file which should contain a sample identifier column and columns for the phenotypes, exposures, and covariates. The ordering of the columns does not matter.
    All inputs should be coded numerically (e.g., males/females as 0/1)

* #### Genotype Files
[*BGEN genotype file*](https://www.well.ox.ac.uk/~gav/bgen_format/spec/latest.html).  
Variants that are non-biallelic should be filtered from the BGEN file. Note that since there are no indication of a REF/ALT allele in the BGEN file, the second allele is the effect allele counted in association testing.   

A [.sample file](https://www.well.ox.ac.uk/~gav/qctool_v2/documentation/sample_file_formats.html) is required as input when the .bgen file does not contain a sample identifier block.
     

* [*Plink BED files*](https://www.cog-genomics.org/plink/1.9/)  
[**.fam**](https://www.cog-genomics.org/plink/2.0/formats#fam) - The .fam file can be space or tab-delimited and must contain at least 2 columns where the first column is the family ID (FID) and the second column is the individual ID (IID). GEM will use the IID column for sample identifier matching with the phenotype file.  

[**.bim**](https://www.cog-genomics.org/plink/2.0/formats#bim) - The .bim file can also be space or tab-delimited and should be in the following order: the chromosome, variant id, cM (optional), base-pair coordinate, ALT allele, and REF allele.  

[**.bed**](https://www.cog-genomics.org/plink/2.0/formats#bed) - A bed file must be stored in variant-major form. The ALT allele specified in the .bim file is the effect allele counted in association testing.   


* [*Plink2.0 PGEN files*](https://www.cog-genomics.org/plink/2.0/)  
[***.psam***](https://www.cog-genomics.org/plink/2.0/formats#psam) - The .psam file is a tab-delimited text file containing the sample information. If header lines are present, the last header line should contain a column with the name #IID (if the first column is not #FID) or IID (if the first column is #FID) that holds the individual ID for sample identifier matching with the phenotype file. All previous header lines will be ignored. If no header line beginning with #IID or #FID is present, then the columns are assumed to be in .fam file order.  

[**.pvar**](https://www.cog-genomics.org/plink/2.0/formats#pvar) - The .pvar file is a tab-delimited text file containing the variant information. If header lines are present, the last header line should start with #CHROM. If #CHROM is present, then the columns POS, ID, REF, and ALT must also be present. All previous header lines will be ignored. If the .pvar file contain no header lines beginning with #CHROM, it is assumed that the columns are in .bim file order.  

[**.pgen**](https://www.cog-genomics.org/plink/2.0/formats#pgen) - The .pgen file should be filtered for non-biallelic variants. The ALT allele specified in the .pvar file is the effect allele counted in association testing.
     

    
<br /> 

### Output File Format  

GEM will write results to the output file specified with the --out parameter (or 'gem.out' if no output file is specified).  
Below are details of the possible column headers in the output file.  

```diff 
SNPID              - The SNP identifier as retrieved from the genotype file.
RSID               - The reference SNP ID number. (BGEN only)
CHR                - The chromosome of the SNP.
POS                - The physical position of the SNP.
Non_Effect_Allele  - The allele not counted in association testing.  
Effect_Allele      - The allele that is counted in association testing.  
N_samples          - The number of samples without missing genotypes.
AF                 - The allele frequency of the effect allele.
N_catE_*           - The number of non-missing samples in each combination of strata for all of the categorical exposures and interaction covariates.
AF_catE_*          - The allele frequency of the effect allele for each combination of strata for all of the catgorical exposure or interaction covariate.

Beta_Marginal           - The coefficient estimate for the marginal genetic effect (i.e., from a model with no interaction terms).
SE_Beta_Marginal        - The model-based SE associated with the marginal genetic effect estimate.
robust_SE_Beta_Marginal - The robust SE associated with the marginal genetic effect estimate.

Beta_G             - The coefficient estimate for the genetic main effect (G).
Beta_G-*           - The coefficient estimate for the interaction or interaction covariate terms.
SE_Beta_G          - Model-based SE associated with the the genetic main effect (G).  
SE_Beta_G-*        - Model-based SE associated with any GxE or interaction covariate terms.
robust_SE_Beta_G   - Robust SE associated with the the genetic main effect (G).  
robust_SE_Beta_G-* - Robust SE associated with any GxE or interaction covariate terms.
Cov_Beta_G_G-*          - Model-based covariance between the genetic main effect (G) and any GxE or interaction covariate terms.
Cov_Beta_G-*_G-*        - Model-based covariance between any GxE or interaction covariate terms.
robust_Cov_Beta_G_G-*   - Robust covariance between the genetic main effect (G) and any GxE or interaction covariate terms.
robust_Cov_Beta_G-*_G-* - Robust covariance between any GxE or interaction covariate terms.

P_Value_Marginal           - Marginal genetic effect p-value from model-based SE.
P_Value_Interaction        - Interaction effect p-value from model-based SE.
P_Value_Joint              - Joint test p-value (K+1 degrees of freedom test of genetic and interaction effect) from model-based SE.
robust_P_Value_Marginal    - Marginal genetic effect p-value from robust SE.
robust_P_Value_Interaction - Interaction effect p-value from robust SE.
robust_P_Value_Joint       - Joint test p-value (K+1 degrees of freedom test of genetic and interaction effect) from robust SE.
```

<br />
The --output-style flag can be used to specify which columns should be included in the output file.  

#### minimum:  
The minimum option includes the variant information, Beta_Marginal, SE_Beta_Marginal, coefficient estimates for only the GxE terms, and depending on the --robust option, SE and covariance for only the GxE terms.

#### meta:  
Includes each of the possible outputs listed above when applicable. For a model-based analysis (--robust 0), the columns containing the "robust" prefix (robust_*) are excluded in the output file.

#### full:  
The "full" option provides, in addition to "meta", an initial header line with the residual variance estimate necessary for re-analysis of a subset of interactions using only summary statistics (for example, switching an exposure and interaction covariate).”

<br />
<br />

### Examples  
<br />

To run GEM using the example data, execute GEM with the following code.
```unix
./GEM --bgen example.bgen --sample example.sample --pheno-file example.pheno --sampleid-name sampleid --pheno-name pheno2 --covar-names cov2 cov3 --exposure-names cov1 --robust 1 --missing-value NaN --out my_example.out
```
The results should look like the following output file [my_example.out](https://github.com/large-scale-gxe-methods/GEM/blob/master/example/my_example.out).  

<br />
<br />


## Contact 
For comments, suggestions, bug reports and questions, please contact Han Chen (Han.Chen.2@uth.tmc.edu), Alisa Manning (AKMANNING@mgh.harvard.edu), Kenny Westerman (KEWESTERMAN@mgh.harvard.edu), or Duy Pham (Duy.T.Pham@uth.tmc.edu). For bug reports, please include an example to reproduce the problem without having to access your confidential data.

<br />
<br />

## License 

 ```
 GEM : Gene-Environment interaction analysis for Millions of samples
 Copyright (C) 2018-2021  Liang Hong, Han Chen, Duy Pham
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ```

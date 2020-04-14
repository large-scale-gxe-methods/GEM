# GEM  

GEM (Gene-Environment interaction analysis for Millions of samples) is a software program for large-scale gene-environment interaction testing in samples from unrelated individuals. It enables genome-wide association studies in up to millions of samples while allowing for multiple exposures, control for genotype-covariate interactions, and robust inference.

<br />

## Contents
- [Installation](#installation)
- [Usage](#usage)
- [License](#license)

<br />

## Installation  
Library Dependencies:  
* BLAS/LAPACK. For top performance, we recommend that the BLAS/LAPACK libraries are linked with optimized math routine libraries such as the Math Kernal Library (MKL) on Intel processors. For AMD processors, ATLAS or OPENBLAS may be better alternatives.  
* Boost C++ libraries. GEM links the following Boost libraries  ```boost_program_options boost_thread boost_system boost_filesystem boost_iostreams``` that will need to be installed prior to executing the makefile.  
* Zlib for BGEN files with zlib compression type.  
* ZSTD for BGEN v1.3 formats.

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

Once GEM is installed, the executable ```./GEM``` can be used to run the program.  
For a list of options, use ```./GEM --help```.  

#### Command Line Options

```

General Options:

--help
    Prints the available options of GEM and exits.  
   
   
--verison
    Prints the version of GEM and exits.
  
  
  
Input File Options:  

--pheno-file  
     Path to the phenotype file.  

--bgen  
     Path to the BGEN file.  

--sample  
     Path to the sample file.  
     Required when the BGEN file does not contain sample identifiers.

--out  
     Full path and extension to where GEM output results.  
     Default: gem.out
  
  
  
Phenotype File Options:

--sampleid-name  
     Column name in the phenotype file that contains sample identifiers.  

--pheno-name  
     Column name in the phenotype file that contains the phenotype of interest.  

--exposure-names  
     One or more column names in the phenotype file naming the exposure(s) to be included in interaction tests.  

--int-covar-names  
     Any column names in the phenotype file naming the covariate(s) for which interactions should be included 
     for adjustment (mutually exclusive with --exposure-names).  

--covar-names  
     Any column names in the phenotype file naming the covariates for which only main effects should be included
     for adjustment (mutually exclusive with both --exposure-names and --int-covar-names).  

--pheno-type
     0 indicates a continuous phenotype and 1 indicates a binary phenotype.  

--robust  
     0 for model-based standard errors and 1 for robust standard errors.
        Default: 0  

--tol 
     Convergence tolerance for logistic regression.  
        Default: 0.0000001  

--delim  
     Delimiter separating values in the phenotype file.  
        Default: , (comma-separated)  

--missing-value  
     Indicates how missing values in the phenotype file are stored.  
        Default: NA
  
  
   
Filtering Options:  

--maf <value>
    Threshold value [0, 1.0] to exclude variants based on the minor allele frequency.
       Default: 0.001
  
  
  
Performance Options:  

--threads <value>
    Set number of compute threads.
    	  Default: ceiling(detected threads / 2)  

--stream-snps  
    Number of SNPs to analyze in a batch. Memory consumption will increase for larger values of stream-snps.  
          Default: 1

```

<br /> 

#### Input Files

* ##### Phenotype File
    A file which should contain a sample identifier column and columns for phenotypes and covariates.  

* ##### BGEN File
    GEM can support v1.1, v1.2 or v1.3 .bgen genotype files described here [BGEN Format](https://www.well.ox.ac.uk/~gav/bgen_format/spec/latest.html).

* ##### Sample File
    A .sample file is required when the .bgen file does not contain sample identifiers.  
    Formats for .sample files should follow QCTOOL v2 ([.sample example](https://www.well.ox.ac.uk/~gav/qctool_v2/documentation/sample_file_formats.html)) format.
    
<br /> 

#### Output File Format  

GEM will write results to the output file specified to the --out paramater (or 'gem.out' if no output file is specified).  
Below are details of the column header in the output file.  

```diff
# SNP Info
SNPID   - The SNP identifier as retrieved from the BGEN file.
rsID    - The reference SNP ID number.
CHR     - The chromosome of the SNP.
POS     - The physical position of the SNP.
Allele1 - The first allele in the BGEN file.
Allele2 - The second allele in the BGEN file that is counted in association testing.
AF      - The allele frequency of the second allele in the BGEN file.


# Betas and Variances
Beta_Marginal      - The coefficient estimate for the marginal genetic effect
Var_Beta_Marginal  - The variance associated with the marginal genetic effect estimate
Beta_Interaction_k - The coefficient estimate for the kth interction term, 
                     where k = {1..length(--exposure-names)}
Var_Beta_Interaction_k_j - The variance associated with the kth and jth interaction term, 
                           where j = {1..length(--exposure-names)}  

# P-values
P_Value_Marginal    - Marginal genetic effect p-value
P_Value_Interaction - Interaction effect p-value
P_Value_Joint       - Joint test p-value (K+1 degrees of freedom test of genetic effect)
```

<br />
<br />

### Examples  
<br />

To run GEM using the example data, execute GEM with the following code.
```unix
./GEM --bgen example.bgen --pheno-file example.pheno --sampleid-name sampleid --pheno-name pheno2 --covar-names cov2 cov3 --exposure-names cov1 --pheno-type 1 --robust 1 --missing-value NaN 
```
The results should look like the following output file [example.out](https://github.com/large-scale-gxe-methods/GEM/blob/master/example/example.out).  

<br />

GEM allows interaction covariate adjustments by specifying the ```--int-covar-names``` parameter.     
```unix
./GEM --bgen example.bgen --pheno-file example.pheno --sampleid-name sampleid --pheno-name pheno2 --covar-names cov2 --exposure-names cov1 --int-covar-names cov3 --pheno-type 1 --robust 1 --missing-value NaN 
```

<br />
<br />



## License 

 ```
 GEM : Gene-Environment interaction analysis for Millions of samples
 Copyright (C) 2018,2019  Liang Hong, Han Chen
 
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

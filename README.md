# GEM

## Contents
- [Installation](#installation)
- [Usage](#usage)
- [License](#license)

<br />

## Installation  
Dependencies:  
* C++11 compiler
* Intel Math Kernel Library (MKL)
* Boost C++ libraries

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

--help
    Prints the available options of GEM and exits.  
   
   
--verison
    Prints the version of GEM and exits.
   

--param <param_file_path>
    Path to the GEM paramater file. (REQUIRED)
   
   
--maf <value>
    Threshold value [0, 1.0] to exclude variants based on the minor allele frequency.
       Default: 0.001


--threads <value>
    Set number of compute threads.
    	  Default: ceiling(detected threads / 2)
```

<br /> 

#### Input Files

* ##### Parameter File Format (Required)
    At minimum, GEM requires a parameter file (.param) as input.  
    For an example of the parameter file format, see the [.param](https://github.com/large-scale-gxe-methods/GEM/blob/master/example/GEM_Input.param) file example.  

* ##### Genotype Files
    Currently, GEM can only process v1.1, v1.2 or v1.3 .bgen genotype files described here [BGEN Format](https://www.well.ox.ac.uk/~gav/bgen_format/spec/latest.html).  
    Future updates will allow different genotype files as input.  

* ##### .sample Files
    A .sample file is required when the .bgen file does not contain sample identifiers.  
    Formats for .sample files should follow QCTOOL v2 ([.sample example](https://www.well.ox.ac.uk/~gav/qctool_v2/documentation/sample_file_formats.html)) format.
    
<br /> 

#### Output File Format  

GEM will write results to the output file specified in the parameter file ([.param](https://github.com/large-scale-gxe-methods/GEM/blob/master/example/GEM_Input.param)).  
Below are details of the column header in the output file.  

```diff
### SNP Info
SNPID   - The SNP identifier.
rsID    - The reference identifier of the SNP.
CHR     - The chromosome of the SNP.
POS     - The position of the SNP of its CHR.
Allele1 - The reference allele (REF) assuming that the REF allele is first.
Allele2 - The alternative allele (ALT) assuming that the ALT allele is second.
AF      - The allele frequency of the ALT allele.


### Beta and Variance
Beta_Main     - The beta value for the main effect
Var_Beta_Main - The variance of the main effect

To complete..


### P-values
P_Value_Main        - Main effect p-value
P_Value_Interaction - Interaction effect p-value
P_Value_Joint       - Joint p-value
```

<br />
<br />

### Examples  
<br />

To run GEM using the example data, execute GEM with the following code.
```unix
./GEM --param ../example/GEM_Input.param
```
The results should look like the following output file [example.out](https://github.com/large-scale-gxe-methods/GEM/blob/master/example/example.out).

<br />
<br />

To exclude variants with minor allele frequency (MAF) < 0.2, use the --maf parameter.
```
./GEM --param ../example/GEM_Input.param --maf 0.2
```

<br />

To spawn 3 threads for multithreading, use the --threads parameter
```
./GEM --param ../example/GEM_Input.param --maf 0.2 --threads 3
```

<br />
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

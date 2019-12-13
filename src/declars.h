#include <vector>
#include <cstring>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <utility>
#include <cctype>
#include <sstream>
#include <cstdlib>
#include <cassert>
#include <ctime>
#include <chrono>
#include <omp.h>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <exception>


#include <boost/algorithm/string/replace.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/program_options.hpp>
#include "zlib.h"

using std::string;
using std::vector;
using std::cin;
using std::cout;
using std::endl;
using std::cerr;
using std::getline;
using std::find;
using std::unordered_map;
using std::unordered_set;
namespace po = boost::program_options;


typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long uint64;
typedef unsigned short ushort;


typedef struct {
  /*** Pheno type data information ***/
  int samSize; /* Sample Size */
  int phenoTyp; /* 0 for continous data using linear regression, 1 for binary data using logistic regression */
  char phenoHeader[300]; /* column name of phenotype data in the pheno file */
  char samIDHeader[300]; /* column name of sample ID in the pheno file */
  int numSelCol; /* total columns of covariate data */
  vector <string> covSelHeaders; /* column header names of the slected covariate data in the pheno data file */
  char MissingKey[300];

  /*** GWAS setting parameter ***/
  int robust; /* 0 for non-robust analysis, 1 for robust analysis */
  int IDMatching; /* 0 for not checking IDMatching, 1 for checking */
  char outputfile[300]; // output file name

  /*** Streaming SNPs ***/
  int stream_snps; /* SNP numbers for each GWAS analysis */
  int Sq; /* Number of interactive covariate data columns with Geno */

  /*** Logistic regression ***/
  double epsilon; /* Convergence tolerance */

  /*** deliminator ***/
  char delim_pheno[300]; // deliminator in phenotype file
  char delim_sample[300]; // deliminator in sample file

} PARAMETERS;

#include "MatrixUtils.h"
#include "ReadParameters.h"

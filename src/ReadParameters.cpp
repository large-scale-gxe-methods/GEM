/****************************************************************************
  ReadParameters.cpp reads parameters from a file
****************************************************************************/

#include "declars.h"

void ReadParameters(char *paramfile, char *genofile, char *phenofile, char *samplefile, PARAMETERS *prm) {

  char lineinfo[500], keyword[200];
  FILE *fparams;
  char tmp[300];

  fparams = fopen(paramfile, "rb");
  if (fparams == NULL) {
    cout << "Cannot open parameter file " << paramfile << endl;
    return;
  }

  while (!feof(fparams)) {
    fgets(lineinfo, 500, fparams);
    /*** read in the keyword about parameters in the next line ***/
    sscanf(lineinfo, "%s", keyword);
    /*** read in the next line of text ***/
    fgets(lineinfo, 500, fparams);

    /*** parse the line text based on the keyword ***/
    if (keyword[0] >= 'A' && keyword[0] <= 'Z') {
    /* the first letter of the keyword must be an upper case letter; otherwise skip the next data line */
      if (!strcmp(keyword, "SAMPLE_SIZE"))  // Sample size
        sscanf(lineinfo, "%d", &prm->samSize); 
      else if (!strcmp(keyword, "SAMPLE_ID_HEADER")) // Header name of sample ID
        sscanf(lineinfo, "%s", prm->samIDHeader);
      else if (!strcmp(keyword, "PHENOTYPE"))  // Type of pheno data (continuous or binary)
        sscanf(lineinfo, "%d", &prm->phenoTyp);
      else if (!strcmp(keyword, "PHENO_HEADER")) // Header name of pheno data
        sscanf(lineinfo, "%s", prm->phenoHeader);
      else if (!strcmp(keyword, "COVARIATE_TOT_COLMUN_NUMS"))  // Total columns of covariate data
        sscanf(lineinfo, "%d", &prm->numSelCol);
      else if (!strcmp(keyword, "ROBUST"))  // Robust analysis or not
        sscanf(lineinfo, "%d", &prm->robust);
      else if (!strcmp(keyword, "SAMPLE_ID_MATCHING"))  // Check sample ID matching order in geno and pheno files or not
        sscanf(lineinfo, "%d", &prm->IDMatching);
      else if (!strcmp(keyword, "COVARIATES_HEADERS")) { // Column numbers of selected covariate data
//        sscanf(lineinfo, "%d %d %d", &prm->colSelVec[0], &prm->colSelVec[1], &prm->colSelVec[2]);
	std::istringstream iss(lineinfo);
	string value;
	vector <string> values;
	while (getline(iss, value, ' ')) values.push_back(value);
	for (int i = 0; i < values.size(); i++) {
	  sscanf(values[i].c_str(), "%s", tmp);
	  prm->covSelHeaders.push_back(string(tmp));
	}
      }
      else if (!strcmp(keyword, "STREAM_SNPS"))  // SNP numbers for each GWAS analysis
        sscanf(lineinfo, "%d", &prm->stream_snps);
      else if (!strcmp(keyword, "NUM_OF_INTER_COVARIATE"))  // Number of interactive covariate data columns with Geno
        sscanf(lineinfo, "%d", &prm->Sq);
      else if (!strcmp(keyword, "LOGISTIC_CONVERG_TOL"))  // Convergence tolerance for logistic regression
        sscanf(lineinfo, "%lf", &prm->epsilon);
      else if (!strcmp(keyword, "GENO_FILE_PATH"))  // Path of geno file
        sscanf(lineinfo, "%s", genofile);
      else if (!strcmp(keyword, "PHENO_FILE_PATH"))  // Path of pheno data file
        sscanf(lineinfo, "%s", phenofile);
      else if (!strcmp(keyword, "SAMPLE_FILE_PATH"))  // Path of sample file
        sscanf(lineinfo, "%s", samplefile);
      else if (!strcmp(keyword, "DELIMINATOR")) // define deliminators
	sscanf(lineinfo, "%s", prm->delim_pheno);
      else if (!strcmp(keyword, "MISSING"))
	sscanf(lineinfo, "%s", prm->MissingKey);
      else if (!strcmp(keyword, "OUTPUT_PATH"))
	sscanf(lineinfo, "%s", prm->outputfile);
      else
        cout << "Does not recognize keyword: " << keyword << endl;
    }
  }
  fclose(fparams);
}

/****************************************************************************
  ReadParameters.cpp reads parameters from a file
****************************************************************************/


#include "declars.h"


#define VERSION "8.1"
void print_help();

void ReadParameters(char *paramfile, char *genofile, char *phenofile, char *samplefile, PARAMETERS *prm) {

  char lineinfo[500], keyword[200];
  FILE *fparams;
  char tmp[300];

  fparams = fopen(paramfile, "rb");
  if (fparams == 0) {
    cout << "\nError: Cannot open parameter file " << paramfile << "\n\n";
    exit(1);
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





// Function to process command line arguments
void CommandLine::processCommandLine(int argc, char* argv[]) {


	// GEM parameters. Details are printed from the print_help() function below.

	// General commands
	po::options_description general("General options");
	general.add_options()
		("help", "")
		("version", "");


	// Input file options
	po::options_description files("Input file options");
	files.add_options()
		("param", po::value<string>(), "");


	// Filtering options
	po::options_description filter("Filtering options");
	filter.add_options()
		("maf", po::value<double>()->default_value(0.001), "");


	//Performance options
	po::options_description performance("Performance options");
	performance.add_options()
		("threads", po::value<int>()->default_value(ceil((boost::thread::hardware_concurrency() / 2))), "");




	// Combine all options together
	po::options_description all("Options");
	all.add(general).add(files).add(filter).add(performance);

	po::variables_map out;


	// QC
	try {
		po::store(po::command_line_parser(argc, argv)
			.options(all)
			.style(po::command_line_style::unix_style
				| po::command_line_style::allow_long_disguise)
			.run(),
			out);

	}
	catch (po::error const& e) {
		std::cerr << e.what() << endl;
		exit(EXIT_FAILURE);
	}
	po::notify(out);




	// Prints the help option and exits
	if (out.count("help")) {
		print_help();
		exit(1);
	}

	// Prints the version and exits
	if (out.count("version")) {
		cout << "\nGEM version: " << VERSION << "\n" << endl;
		exit(1);
	}


	// Read parameter file
	if (out.count("param")) {
		pFile = out["param"].as<string>();
	}
	else {
		cerr << "\nERROR: parameter file (-param) is needed." << endl;
		exit(1);
	}



	// Minor allele frequency threshold
	if (out.count("maf")) {
		MAF = out["maf"].as<double>();
	}


	if (out.count("threads")) {
		threads = out["threads"].as<int>();
	}
}






void print_help() {


	// Welcome and version output
	cout << "\n\nWelcome to GEM" << endl;
	cout << "Version: " << VERSION << endl << endl << endl;
	cout << "Usage: GEM <options>" << endl << endl;



	cout << "General Options: " << endl
		 << "   -help \t\t Prints available options and exits." << endl
		 << "   -version \t\t Prints the version of GEM and exits." << endl;
	cout << endl << endl;



	cout << "Input Options: " << endl
		 << "   -param \t\t Path to GEM parameter file." << endl;
	cout << endl << endl;



	cout << "Filtering Options: " << endl
		 << "   -maf \t\t Threshold to filter variants based on the minor allele frequency.\n \t\t\t    Default: 0.001" << endl;
	cout << endl << endl;



	cout << "Performance Options:" << endl
	     << "   -threads \t\t Set number of compute threads \n \t\t\t    Default: ceiling(detected threads / 2)" << endl;
    cout << endl << endl;

	     
}
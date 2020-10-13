/****************************************************************************
  ReadParameters.cpp reads parameters from a file
****************************************************************************/


#include "declars.h"
#define VERSION "1.2"

void print_help();


// Function to process command line arguments
void CommandLine::processCommandLine(int argc, char* argv[]) {


	// GEM parameters. Details are printed from the print_help() function below.

	// General commands
	po::options_description general("General options");
	general.add_options()
		("help, h", "")
		("version", "");

	// Input file options
	po::options_description files("Input file options");
	files.add_options()
		("bgen", po::value<std::string>(), "")
		("sample", po::value<std::string>(), "")
		("pfile", po::value<std::string>(), "")
		("pgen", po::value<std::string>(), "")
		("pvar", po::value<std::string>(), "")
		("psam", po::value<std::string>(), "")
		("bfile", po::value<std::string>(), "")
		("bed", po::value<std::string>(), "")
		("bim", po::value<std::string>(), "")
		("fam", po::value<std::string>(), "")
		("pheno-file", po::value<std::string>(), "")
		("out", po::value<std::string>()->default_value("gem.out"), "");

	// Phenotype file
	po::options_description phenofile("Phenotype file options");
	phenofile.add_options()
		("sampleid-name", po::value<std::string>(), "")
		("pheno-type", po::value<int>(), "")
		("pheno-name", po::value<std::string>(), "")
		("covar-names", po::value<std::vector<std::string>>()->multitoken(), "")
		("int-covar-names", po::value<std::vector<std::string>>()->multitoken(), "")
		("exposure-names", po::value<std::vector<std::string>>()->multitoken(), "")
		("delim", po::value<std::string>()->default_value(","), "")
		("missing-value", po::value<std::string>()->default_value("NA"), "")
		("robust", po::value<int>()->default_value(0), "")
		("tol", po::value<double>()->default_value(.0000001));

	// Filtering options
	po::options_description filter("Filtering options");
	filter.add_options()
		("maf", po::value<double>()->default_value(0.001), "")
		("miss-geno-cutoff", po::value<double>()->default_value(0.05), "")
		("include-snp-file", po::value<std::string>(), "");

	//Performance options
	po::options_description performance("Performance options");
	performance.add_options()
		("threads", po::value<int>()->default_value(ceil((boost::thread::hardware_concurrency() / 2))), "")
		("stream-snps", po::value<int>()->default_value(1), "");

	// Combine all options together
	po::options_description all("Options");
	all.add(general).add(files).add(phenofile).add(filter).add(performance);

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



	// General
	if (out.count("help")) {
		print_help();
		exit(1);
	}
	if (out.count("version")) {
		cout << "\nGEM version: " << VERSION << "\n" << endl;
		exit(1);
	}


	// Input Files
	if (out.count("pheno-file")) {
		phenoFile = out["pheno-file"].as<string>();
	}
	else {
		cerr << "\nERROR: Phenotype file (--pheno-file) is needed. \n\n";
		exit(1);
	}


	if (out.count("bgen")) {
		bgenFile = out["bgen"].as<string>();
		strcpy(genofile, bgenFile.c_str());
		useBgenFile = true;

		if (out.count("pgen") || out.count("pfile") || out.count("bed") || out.count("bfile")) {
			cerr << "\nERROR: Only one genotype file format can be used.\n\n";
			exit(1);
		}
	}
	else if (out.count("pfile")) {

		if (out.count("pgen")) {
			cerr << "\nERROR: --pfile and --pgen cannot be used simultaneously.\n\n";
			exit(1);
		}

		if (out.count("bed") || out.count("bfile")) {
			cerr << "\nERROR: --pfile and --bed/--bfile cannot be used simultaneously.\n\n";
			exit(1);
		}

		pgenFile = out["pfile"].as<string>() + ".pgen";
		psamFile = out["pfile"].as<string>() + ".psam";
		pvarFile = out["pfile"].as<string>() + ".pvar";
		strcpy(pgenfile, pgenFile.c_str());
		usePgenFile = true;

	}
	else if (out.count("pgen")) {
		if (out.count("pfile")) {
			cerr << "\nERROR: --pgen and --pfile cannot be used simultaneously.\n\n";
			exit(1);
		}
		if (out.count("bed") || out.count("bfile")) {
			cerr << "\nERROR: --pgen and --bed/--bfile cannot be used simultaneously.\n\n";
			exit(1);
		}
		if (!out.count("psam")) {
			cerr << "\nERROR: .psam file (--psam) is needed.\n\n";
			exit(1);
		}
		if (!out.count("pvar")) {
			cerr << "\nERROR: .pvar file (--pvar) is needed.\n\n";
			exit(1);
		}
		pgenFile = out["pgen"].as<string>() + ".pgen";
		psamFile = out["psam"].as<string>() + ".psam";
		pvarFile = out["pvar"].as<string>() + ".pvar";
		strcpy(pgenfile, pgenFile.c_str());
		usePgenFile = true;
	}
	else if (out.count("bfile")) {
		if (out.count("bed")) {
			cerr << "\nERROR: --bed and --bfile cannot be used simultaneously.\n\n";
			exit(1);
		}
		if (out.count("pgen") || out.count("pfile")) {
			cerr << "\nERROR: --bfile and --pfile/--pgen cannot be used simultaneously.\n\n";
			exit(1);
		}
		bedFile = out["bfile"].as<string>() + ".bed";
		famFile = out["bfile"].as<string>() + ".fam";
		bimFile = out["bfile"].as<string>() + ".bim";
		strcpy(pgenfile, pgenFile.c_str());
		useBedFile = true;
	}
	else if (out.count("bed")) {
		if (out.count("bfile")) {
			cerr << "\nERROR: --pgen and --pfile cannot be used simultaneously.\n\n";
			exit(1);
		}
		if (out.count("pgen") || out.count("pfile")) {
			cerr << "\nERROR: --bed and --pfile/--pgen cannot be used simultaneously.\n\n";
			exit(1);
		}
		if (!out.count("fam")) {
			cerr << "\nERROR: .fam file (--fam) is needed.\n\n";
			exit(1);
		}
		if (!out.count("bim")) {
			cerr << "\nERROR: .bim file (--bim) is needed.\n\n";
			exit(1);
		}
		bedFile = out["bed"].as<string>() + ".bed";
		famFile = out["fam"].as<string>() + ".fam";
		bimFile = out["bim"].as<string>() + ".bim";
		strcpy(pgenfile, pgenFile.c_str());
		useBedFile = true;
	}
	else {
		cerr << "\nERROR: Genotype file (--bgen) / (--pfile/--pgen) / (--bfile/--bed) is needed.\n\n";
		exit(1);
	}


	if (out.count("sample")) {
		if (out.count("pgen") || out.count("pfile") || out.count("bed") || out.count("bfile")) {
			cerr << "\nERROR: --sample flag should only be used with --bgen.\n\n";
			exit(1);
		}
		sampleFile = out["sample"].as<string>();
		strcpy(samplefile, sampleFile.c_str());
		useSampleFile = true;
	}
	if (out.count("out")) {
		outFile = out["out"].as<string>();

		std::ofstream results(outFile, std::ofstream::binary);
		if (!results) {
			cerr << "\nERROR: Output file could not be opened.\n\n";
			exit(1);
		}

		if (results.fail()) {
			cerr << "\nERROR: Output file could not be opened.\n\n";
			exit(1);
		}

		results << "test" << endl;
		if (results.fail()) {
			cout << "\nERROR: Cannot write to output file.\n\n";
		}
		results.close();
		boost::filesystem::remove(outFile.c_str());
	}



	// Phenotype file
	if (out.count("sampleid-name")) {
		sampleID = out["sampleid-name"].as<string>();
	}
	else {
		cerr << "\nERROR: Sample ID column name (--sample-col) is not specified. \n\n";
		exit(1);
	}
	if (out.count("pheno-name")) {
		phenoName = out["pheno-name"].as<string>();
	}
	else {
		cerr << "\nERROR: Phenotype column name (--pheno-name) is not specified.\n\n";
		exit(1);
	}
	if (out.count("pheno-type")) {
		phenoType = out["pheno-type"].as<int>();
		switch (phenoType) {
		case 0:
			break;
		case 1:
			break;
		default:
			cerr << "\nERROR: --pheno-type must be 0 (continuous) or 1 (binary). \n\n";
			exit(1);
		}
	}
	else {
		cerr << "\nERROR: --pheno-type is not specified. Must be 0 (continuous) or 1 (binary). \n\n";
		exit(1);
	}

	if (out.count("covar-names")) {
		cov = out["covar-names"].as<std::vector<std::string>>();
	}
	if (out.count("int-covar-names")) {
		icov = out["int-covar-names"].as<std::vector<std::string>>();
	}
	if (out.count("exposure-names")) {
		exp = out["exposure-names"].as<std::vector<std::string>>();
	}
	else {
		cerr << "\n ERROR: No exposures (--exposure-names) specified. \n\n";
		exit(1);
	}
	if (out.count("delim")) {
		delim = out["delim"].as<string>();
		strcpy(pheno_delim, delim.c_str());
	}
	if (out.count("missing-value")) {
		missing = out["missing-value"].as<std::string>();
	}
	if (out.count("robust")) {
		robust = out["robust"].as<int>();
		switch (robust) {
		case 0:
			break;
		case 1:
			break;
		default:
			cerr << "\nERROR: --robust must be 0 (non-robust) or 1 (robust). \n\n";
			exit(1);
		}
	}
	if (out.count("tol")) {
		tol = out["tol"].as<double>();
	}



	// Filtering options
	if (out.count("maf")) {
		MAF = out["maf"].as<double>();

		if (MAF < 0.0 || MAF > 1) {
			cerr << "\nERROR: Please specify --maf with a value greater than or equal to 0 and less than or equal to 1.\n\n";
			exit(1);
		}
	}
	if (out.count("miss-geno-cutoff")) {
		missGenoRate = out["miss-geno-cutoff"].as<double>();
		if (missGenoRate < 0.0 || missGenoRate > 1.0) {
			cerr << "\nERROR: Please specify --miss-geno-cutoff with a value greater than or equal to 0 and less than or equal to 1.\n\n";
			exit(1);
		}

	}
	if (out.count("include-snp-file")) {
		includeVariantFile = out["include-snp-file"].as<std::string>();
		if (out.count("pgen") || out.count("pfile") || out.count("bed") || out.count("bfile")) {
			cerr << "\nERROR: --include-snp-file currently unsupported for plink genotype files.\n\n";
			exit(1);
		}
	}

	if (includeVariantFile.empty()) {
		doFilters = false;
	}
	else {
		doFilters = true;
	}


	// Performance options
	if (out.count("threads")) {
		threads = out["threads"].as<int>();

		if (threads <= 0) {
			cerr << "\nERROR: Please specify --threads with a value greater than 0.\n\n";
			exit(1);
		}
	}
	if (out.count("stream-snps")) {
		stream_snps = out["stream-snps"].as<int>();

		if (stream_snps <= 0) {
			cerr << "\nERROR: Please specify --stream_snps with a value greater than 0.\n\n";
			exit(1);
		}
	}
}






void print_help() {

	cout << "\n\nWelcome to GEM" << endl;
	cout << "Version: " << VERSION << endl;
	cout << "(C) 2018-2020 Liang Hong, Han Chen, Duy Pham \n";
	cout << "GNU General Public License v3\n\n\n";


	cout << "General Options: " << endl
		<< "   --help \t\t Prints available options and exits." << endl
		<< "   --version \t\t Prints the version of GEM and exits." << endl;
	cout << endl << endl;



	cout << "Input File Options: " << endl
		<< "   --pheno-file \t Path to the phenotype file." << endl
		<< "   --bgen \t\t Path to the BGEN file." << endl
		<< "   --sample \t\t Path to the sample file. Required when the BGEN file does not contain sample identifiers." << endl
		<< "   --pfile \t\t Path and prefix to the .pgen, .pvar, and .psam files." << endl
		<< "   --pgen \t\t Path to the pgen file." << endl
		<< "   --pvar \t\t Path to the pvar file." << endl
		<< "   --psam \t\t Path to the psam file." << endl
		<< "   --bfile \t\t Path and prefix to the .bed, .bim and .fam files." << endl
		<< "   --bed \t\t Path to the bed file." << endl
		<< "   --bim \t\t Path to the bim file." << endl
		<< "   --fam \t\t Path to the fam file." << endl
		<< "   --out \t\t Full path and extension to where GEM output results. \n \t\t\t    Default: gem.out" << endl;
	cout << endl << endl;


	cout << "Phenotype File Options: " << endl
		<< "   --sampleid-name \t Column name in the phenotype file that contains sample identifiers." << endl
		<< "   --pheno-name \t Column name in the phenotype file that contains the phenotype of interest." << endl
		<< "   --exposure-names \t One or more column names in the phenotype file naming the exposure(s) to be included in interaction tests." << endl
		<< "   --int-covar-names \t Any column names in the phenotype file naming the covariate(s) for which interactions should\n \t\t\t   be included for adjustment (mutually exclusive with --exposure-names)." << endl
		<< "   --covar-names \t Any column names in the phenotype file naming the covariates for which only main effects should\n \t\t\t   be includedfor adjustment (mutually exclusive with both --exposure-names and --int-covar-names)." << endl
		<< "   --pheno-type \t 0 indicates a continuous phenotype and 1 indicates a binary phenotype." << endl
		<< "   --robust \t\t 0 for model-based standard errors and 1 for robust standard errors. \n \t\t\t    Default: 0" << endl
		<< "   --tol \t\t Convergence tolerance for logistic regression. \n \t\t\t    Default: 0.0000001" << endl
		<< "   --delim \t\t Delimiter separating values in the phenotype file. \n \t\t\t    Default: , (comma-separated)" << endl
		<< "   --missing-value \t Indicates how missing values in the phenotype file are stored. \n \t\t\t    Default: NA" << endl;
	cout << endl << endl;


	cout << "Filtering Options: " << endl
		<< "   --maf \t\t Threshold to filter variants based on the minor allele frequency.\n \t\t\t    Default: 0.001" << endl
		<< "   --miss-geno-cutoff \t Threshold to filter variants based on the missing genotype rate.\n \t\t\t    Default: 0.05" << endl
		<< "   --include-snp-file \t Path to file containing a subset of variants in the specified BGEN file to be used for analysis. The first\n \t\t\t   line in this file is the header that specifies which variant identifier in the BGEN file is used for ID\n \t\t\t   matching. This must be either 'snpid' or 'rsid'. There should be one variant identifier per line after the header.\n \t\t\t   Variants not listed in this file will be excluded from analysis." << endl;
	cout << endl << endl;



	cout << "Performance Options:" << endl
		<< "   --threads \t\t Set number of compute threads \n \t\t\t    Default: ceiling(detected threads / 2)" << endl
		<< "   --stream-snps \t Number of SNPs to analyze in a batch. Memory consumption will increase for larger values of stream-snps.\n \t\t\t    Default: 1" << endl;
	cout << endl << endl;

}
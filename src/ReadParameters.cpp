/****************************************************************************
  ReadParameters.cpp reads parameters from a file
****************************************************************************/


#include "declars.h"
#define VERSION "1.4"

void print_help();


// Function to process command line arguments
void CommandLine::processCommandLine(int argc, char* argv[]) {


    cout << "\n*********************************************************\n";
    cout << "Welcome to GEM v" << VERSION << "\n";
    cout << "(C) 2018-2021 Liang Hong, Han Chen, Duy Pham \n";
    cout << "GNU General Public License v3\n";
    cout << "*********************************************************\n";


    // GEM parameters. Details are printed from the print_help() function below.

    // General commands
    po::options_description general("General options");
    general.add_options()
        ("help, h", "")
        ("version", "");

    // Input/Output file options
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
        ("out", po::value<std::string>()->default_value("gem.out"), "")
        ("output-style", po::value<std::string>()->default_value("minimum"), "");

    // Phenotype file
    po::options_description phenofile("Phenotype file options");
    phenofile.add_options()
        ("sampleid-name", po::value<std::string>(), "")
        ("pheno-name", po::value<std::string>(), "")
        ("covar-names", po::value<std::vector<std::string>>()->multitoken(), "")
        ("int-covar-names", po::value<std::vector<std::string>>()->multitoken(), "")
        ("exposure-names", po::value<std::vector<std::string>>()->multitoken(), "")
        ("delim", po::value<std::string>()->default_value(","), "")
        ("missing-value", po::value<std::string>()->default_value("NA"), "")
        ("robust", po::value<int>()->default_value(0), "")
        ("tol", po::value<double>()->default_value(.000001))
        ("center", po::value<int>()->default_value(1))
        ("scale", po::value<int>()->default_value(0))
        ("categorical-names", po::value<std::vector<std::string>>()->multitoken(), "")
        ("cat-threshold", po::value<int>()->default_value(20));

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


    // Input/Output Files
    if (out.count("pheno-file")) {
        phenoFile = out["pheno-file"].as<string>();

    }
    else {
        cerr << "\nERROR: Phenotype file (--pheno-file) is needed. \n\n";
        exit(1);

    }

    if (out.count("bgen")) {
        if (out.count("pgen") || out.count("pfile") || out.count("bed") || out.count("bfile")) {
            cerr << "\nERROR: Only one genotype file format can be used.\n\n";
            exit(1);
        }

        bgenFile = out["bgen"].as<string>();
        useBgenFile = true;

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

        pgenFile = out["pgen"].as<string>();
        psamFile = out["psam"].as<string>();
        pvarFile = out["pvar"].as<string>();
        usePgenFile = true;

    }
    else if (out.count("bfile")) {
        if (out.count("bed")) {
            cerr << "\nERROR: --bed and --bfile cannot be used simultaneously.\n\n";
            exit(1);
        }

        bedFile = out["bfile"].as<string>() + ".bed";
        famFile = out["bfile"].as<string>() + ".fam";
        bimFile = out["bfile"].as<string>() + ".bim";
        useBedFile = true;

    }
    else if (out.count("bed")) {
        if (!out.count("fam")) {
            cerr << "\nERROR: .fam file (--fam) is needed.\n\n";
            exit(1);
        }
        if (!out.count("bim")) {
            cerr << "\nERROR: .bim file (--bim) is needed.\n\n";
            exit(1);
        }

        bedFile = out["bed"].as<string>();
        famFile = out["fam"].as<string>();
        bimFile = out["bim"].as<string>();
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
            cerr << "\nERROR: Cannot write to output file.\n\n";
            results.close();
            exit(1);
        }
        results.close();

        boost::filesystem::remove(outFile.c_str());
    }
    if (out.count("output-style")) {
        outStyle = out["output-style"].as<string>();
        if (outStyle.compare("minimum") != 0 && outStyle.compare("meta") != 0  && outStyle.compare("full") != 0 ) {
            cerr << "\nERROR: --output-style should be minimum, meta or full.\n\n";
            exit(1);
        }
    }

    // Phenotype file
    if (out.count("sampleid-name")) {
        sampleID = out["sampleid-name"].as<string>();

    }
    else {
        cerr << "\nERROR: Sample ID column name (--sampleid-name) is not specified. \n\n";
        exit(1);

    }
    if (out.count("pheno-name")) {
        phenoName = out["pheno-name"].as<string>();

    }
    else {
        cerr << "\nERROR: Phenotype column name (--pheno-name) is not specified.\n\n";
        exit(1);

    }


    if (out.count("exposure-names")) {
        exp = out["exposure-names"].as<std::vector<std::string>>();

        for (size_t i = 0; i < exp.size(); i++) {
            if (phenoName.compare(exp[i]) == 0) {
                cerr << "\nERROR: Exposure " << exp[i] << " is also specified as the phenotype (--pheno-name)." << "\n\n";
                exit(1);
            }

            expHM[exp[i]] += 1;
            if (expHM[exp[i]] > 1) {
                cerr << "\nERROR: Exposure " + exp[i] + " is specified more than once.\n\n";
                exit(1);
            }
        }
        numExpSelCol = exp.size();

    }
    if (out.count("int-covar-names")) {
        icov = out["int-covar-names"].as<std::vector<std::string>>();
        if (exp.size() == 0) {
            cerr << "\nERROR: --int-covar-names should not be included when there are no exposures. \n\n";
            exit(1);
        }
        for (size_t i = 0; i < icov.size(); i++) {
            if (expHM.find(icov[i]) != expHM.end()) {
                cerr << "\nERROR: Interactive covariate " << icov[i] << " is specified as an interaction covariate (--int-covar-names) and exposure (--exposure-names)." << "\n\n";
                exit(1);
            }
            if (phenoName.compare(icov[i]) == 0) {
                cerr << "\nERROR: Interactive covariate " << icov[i] << " is also specified as the phenotype (--pheno-name)." << "\n\n";
                exit(1);
            }

            intHM[icov[i]] += 1;
            if (intHM[icov[i]] > 1) {
                cerr << "\nERROR: Interactive covariate " + icov[i] + "is specified more than once.\n\n";
                exit(1);
            }
        }
        numIntSelCol = icov.size();

    }
    if (out.count("covar-names")) {
        cov = out["covar-names"].as<std::vector<std::string>>();

        for (size_t i = 0; i < cov.size(); i++) {
            if (expHM.find(cov[i]) != expHM.end()) {
                cerr << "\nERROR: Covariate " << cov[i] << " is specified as a covariate (--covar-names) and exposure (--exposure-names)." << "\n\n";
                exit(1);
            }
            if (intHM.find(cov[i]) != intHM.end()) {
                cerr << "\nERROR: Covariate " << cov[i] << " is specified as a covariate (--covar-names) and interaction covariate (--int-covar-names)." << "\n\n";
                exit(1);
            }
            if (phenoName.compare(cov[i]) == 0) {
                cerr << "\nERROR: Covariate " << cov[i] << " is also specified as the phenotype (--pheno-name)." << "\n\n";
                exit(1);
            }

            covHM[cov[i]] += 1;
            if (covHM[cov[i]] > 1) {
                cerr << "\nERROR: Covariate " + cov[i] + " is specified more than once.\n\n";
                exit(1);
            }
        }
        numSelCol = cov.size();

    }
    if (out.count("delim")) {
        string s_delim = out["delim"].as<string>();

        char delim[300];
        strcpy(delim, s_delim.c_str());
        if ((delim[0] == '\\' && delim[1] == 't') || delim[0] == 't') {
            pheno_delim = '\t';
        }
        else if ((delim[0] == '\\' && delim[1] == '0') || delim[0] == '0') {
            pheno_delim = ' ';
        }
        else {
            pheno_delim = delim[0];
        }

    }
    if (out.count("missing-value")) {
        missing = out["missing-value"].as<std::string>();

    }
    if (out.count("robust")) {
        robust = out["robust"].as<int>();

        if (robust != 0 && robust != 1) {
            cerr << "\nERROR: Please specify --robust with a value equal to 0 (false) or 1 (true). \n\n";
            exit(1);
        }

    }
    if (out.count("tol")) {
        tol = out["tol"].as<double>();

    }
    if (out.count("center")) {
        center = out["center"].as<int>();

        if (center != 0 && center != 1) {
            cerr << "\nERROR: Please specify --center with a value equal to 0 (false) or 1 (true).\n\n";
            exit(1);
        }

    }
    if (out.count("scale")) {
        scale = out["scale"].as<int>();

        if (scale != 0 && scale != 1) {
            cerr << "\nERROR: Please specify --scale with a value equal to 0 (false) or 1 (true).\n\n";
            exit(1);
        }
    }
    if (out.count("categorical-names")) {
        cat_names = out["categorical-names"].as<std::vector<std::string>>();
        if (cat_names.size() > 0) {
            for (size_t i = 0; i < cat_names.size(); i++) {
                if (expHM.find(cat_names[i]) != expHM.end()) {
                    continue;
                }
                else if (intHM.find(cat_names[i]) != intHM.end()) {
                    continue;
                } 
                else {
                    cerr << "\nERROR: " << cat_names[i] << " is not an exposure or interaction exposure.\n\n";
                    exit(1);
                }
            }
        }
    }
    if (out.count("cat-threshold")) {
        cat_threshold = out["cat-threshold"].as<int>();

        if (cat_threshold <= 1) {
            cerr << "\nERROR: Please specify --cat-threshold with a value greater than 1.\n\n";
            exit(1);
        }
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
            cerr << "\nERROR: Please specify --stream-snps with a value greater than 0.\n\n";
            exit(1);
        }
    }


    // Print parameter info
    cout << "The Phenotype File is: " << phenoFile << "\n";
    cout << "The Genotype File is: ";
    if (useBgenFile) {
        cout << bgenFile << "\n";
    } else if (usePgenFile) {
        cout << pgenFile << "\n";
    } else {
        cout << bedFile << "\n";
    }

    cout << "Model-based or Robust: "; robust == 0 ? cout << "Model-based \n\n" : cout << "Robust \n\n";

    if (numSelCol == 0) {
        cout << "No Covariates Selected" << "\n";
    }
    else {
        cout << "The Total Number of Selected Covariates is: " << numSelCol << '\n';
        cout << "The Selected Covariates are:  ";
        for (int i = 0; i < numSelCol; i++) {
            cout << cov[i] << "   ";
        }
        cout << "\n";
    }

    if (numIntSelCol == 0) {
        cout << "No Interaction Covariates Selected" << "\n";
    }
    else {
        cout << "The Total Number of Selected Interaction Covariates is: " << numIntSelCol << "\n";
        cout << "The Selected Interaction Covariates are:  ";
        for (int i = 0; i < numIntSelCol; i++) {
            cout << icov[i] << "   ";
        }
        cout << "\n";
    }

    if (numExpSelCol == 0) {
        cout << "No Exposures Selected" << "\n";
    }
    else {
        cout << "The Total Number of Exposures is: " << numExpSelCol << '\n';
        cout << "The Selected Exposures are:  ";
        for (int i = 0; i < numExpSelCol; i++) {
            cout << exp[i] << "   ";
        }
        cout << "\n\n";
    }

    cout << "Categorical Threshold: " << cat_threshold << "\n";
    cout << "Minor Allele Frequency Threshold: " << MAF << "\n";
    cout << "Number of SNPS in batch: " << stream_snps << "\n";
    cout << "Number of Threads: " << threads << "\n";
    cout << "Output File: " << outFile << "\n";
    cout << "*********************************************************\n";
}






void print_help() {

    cout << "General Options: " << endl
        << "   --help \t\t Prints available options and exits." << endl
        << "   --version \t\t Prints the version of GEM and exits." << endl;
    cout << endl << endl;



    cout << "Input/Output File Options: " << endl
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
        << "   --out \t\t Full path and extension to where GEM output results. \n \t\t\t    Default: gem.out" << endl
        << "   --output-style \t Modifies the output of GEM. Must be one of the following: \n\t\t\t    minimum: Output the summary statistics for only the GxE and marginal G terms. \n \t\t\t    meta: 'minimum' output plus additional fields for the main G and any GxCovariate terms \n \t\t\t    full: 'meta' output plus additional fields needed for re-analyses of a subset of interactions \n \t\t\t    Default: minimum" << endl;       
    cout << endl << endl;


    cout << "Phenotype File Options: " << endl
        << "   --sampleid-name \t Column name in the phenotype file that contains sample identifiers." << endl
        << "   --pheno-name \t Column name in the phenotype file that contains the phenotype of interest." << endl
        << "   --exposure-names \t One or more column names in the phenotype file naming the exposure(s) to be included in interaction tests." << endl
        << "   --int-covar-names \t Any column names in the phenotype file naming the covariate(s) for which interactions should\n \t\t\t   be included for adjustment (mutually exclusive with --exposure-names)." << endl
        << "   --covar-names \t Any column names in the phenotype file naming the covariates for which only main effects should\n \t\t\t   be included for adjustment (mutually exclusive with both --exposure-names and --int-covar-names)." << endl
        << "   --robust \t\t 0 for model-based standard errors and 1 for robust standard errors. \n \t\t\t    Default: 0" << endl
        << "   --tol \t\t Convergence tolerance for logistic regression. \n \t\t\t    Default: 0.0000001" << endl
        << "   --delim \t\t Delimiter separating values in the phenotype file. Tab delimiter should be represented as \\t and space delimiter as \\0. \n \t\t\t    Default: , (comma-separated)" << endl
        << "   --missing-value \t Indicates how missing values in the phenotype file are stored. \n \t\t\t    Default: NA" << endl
        << "   --center \t\t 0 for no centering to be done and 1 to center ALL exposures and covariates. \n \t\t\t    Default: 1" << endl
        << "   --scale \t\t 0 for no scaling to be done and 1 to scale ALL exposures and covariates by the standard deviation. \n \t\t\t    Default: 0" << endl
        << "   --categorical-names \t Names of the exposure or interaction covariate that should be treated as categorical. \n \t\t\t    Default: None" << endl
        << "   --cat-threshold \t A cut-off to determine which exposure or interaction covariate is categorical based on the number of unique elements. \n \t\t\t    Default: 20" << endl;
    cout << endl << endl;


    cout << "Filtering Options: " << endl
        << "   --maf \t\t Threshold to filter variants based on the minor allele frequency.\n \t\t\t    Default: 0.001" << endl
        << "   --miss-geno-cutoff \t Threshold to filter variants based on the missing genotype rate.\n \t\t\t    Default: 0.05" << endl
        << "   --include-snp-file \t Path to file containing a subset of variants in the specified genotype file to be used for analysis. The first\n \t\t\t   line in this file is the header that specifies which variant identifier in the genotype file is used for ID\n \t\t\t   matching. This must be 'snpid' or 'rsid' (BGEN only). There should be one variant identifier per line after the header.\n \t\t\t   Variants not listed in this file will be excluded from analysis." << endl;
    cout << endl << endl;



    cout << "Performance Options:" << endl
        << "   --threads \t\t Set number of compute threads \n \t\t\t    Default: ceiling(detected threads / 2)" << endl
        << "   --stream-snps \t Number of SNPs to analyze in a batch. Memory consumption will increase for larger values of stream-snps.\n \t\t\t    Default: 1" << endl;
    cout << endl << endl;

}

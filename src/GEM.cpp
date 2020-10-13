/*  GEM : Gene-Environment interaction analysis for Millions of samples
 *  Copyright (C) 2018-2020  Liang Hong, Han Chen, Duy Pham
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

 /*
   software description and change log

   version 1: GWAS of snp one by one
   version 2: streaming multi-snps and GWAS analysis of multi-snps simultaneously
   version bgen: multi-snps and bgen
   version Log:  logistic regression for binary phenodata
   7/27/18:  input column header name in pheno file not index number
   8/01/18:  sample ID matching
   8/30/18:  initial sample Size is not necessary equal
   8/30/18:  hash table of genoUnMatchID for significant unmatching numbers;
   12/10/18: aritrary stream_snps implemented, omp parallel with > 20 covariates (HC)
   2/7/19:   speed up reading genotype data without looking up in genoUnMatchID

   To-Do List:
   1. OOP
 */

#include "declars.h"
#define VERSION 1.2




int main(int argc, char* argv[]) {

	cout << "\n*********************************************************\n";
	cout << "Welcome to GEM v" << VERSION << "\n";
	cout << "(C) 2018-2020 Liang Hong, Han Chen, Duy Pham \n";
	cout << "GNU General Public License v3\n";
	cout << "*********************************************************\n";

	// Process command line
	CommandLine cmd;
	cmd.processCommandLine(argc, argv);



	/****************************************************
	 Parameters
	****************************************************/
	double epsilon = cmd.tol;

	int samSize;
	int phenoCol;
	int samIDCol;
	int robust = cmd.robust;
	int phenoTyp = cmd.phenoType;
	int stream_snps = cmd.stream_snps;
	char delim = cmd.pheno_delim[0];
	if (delim == '\0') delim = ' ';

	string phenoHeaderName = cmd.phenoName;
	string samIDHeaderName = cmd.sampleID;
	string phenoMissingKey = cmd.missing;
	string output = cmd.outFile;
	std::unordered_map<string, int> covHM;
	std::unordered_map<string, int> intHM;
	std::unordered_map<string, int> expHM;
	vector<string> covSelHeadersName;
	vector<string> intCovSelHeadersName;
	vector<string> expCovSelHeadersName;


	string expErrorString;
	int numExpSelCol = cmd.exp.size();
	for (int i = 0; i < numExpSelCol; i++) {

		if (phenoHeaderName.compare(cmd.exp[i]) == 0) {
			cerr << "\nERROR: Exposure " << cmd.exp[i] << " is also specified as the phenotype (--pheno-name)." << "\n\n";
			exit(1);
		}

		expHM[cmd.exp[i]] += 1;
		if (expHM[cmd.exp[i]] > 1) {
			expErrorString += cmd.exp[i] + " ";
			continue;
		}
		expCovSelHeadersName.push_back(cmd.exp[i]);
	}
	numExpSelCol = expCovSelHeadersName.size();

	if (!expErrorString.empty()) {
		cout << "WARNING: Exposure " + expErrorString + "is specified more than once.\n \t Removing the duplicate exposure names... \n\n";
	}




	string intCovErrorString;
	int numIntSelCol = cmd.icov.size();
	for (int i = 0; i < numIntSelCol; i++) {
		if (expHM.find(cmd.icov[i]) != expHM.end()) {
			cerr << "\nERROR: Interactive covariate " << cmd.icov[i] << " is specified as an interaction covariate (--int-covar-names) and exposure (--exposure-names)." << "\n\n";
			exit(1);
		}
		if (phenoHeaderName.compare(cmd.icov[i]) == 0) {
			cerr << "\nERROR: Interactive covariate " << cmd.icov[i] << " is also specified as the phenotype (--pheno-name)." << "\n\n";
			exit(1);
		}
		intHM[cmd.icov[i]] += 1;
		if (intHM[cmd.icov[i]] > 1) {
			intCovErrorString += cmd.icov[i] + " ";
			continue;
		}
		intCovSelHeadersName.push_back(cmd.icov[i]);
	}

	numIntSelCol = intCovSelHeadersName.size();
	int Sq = numIntSelCol + numExpSelCol;

	if (!intCovErrorString.empty()) {
		cout << "\nWARNING: Interactive covariates " + intCovErrorString + "is specified more than once.\n \t Removing the duplicate interaction covariate names... \n\n";
	}




	string covErrorString;
	int numSelCol = cmd.cov.size();
	for (int i = 0; i < numSelCol; i++) {

		if (expHM.find(cmd.cov[i]) != expHM.end()) {
			cerr << "\nERROR: Covariate " << cmd.cov[i] << " is specified as a covariate (--covar-names) and exposure (--exposure-names)." << "\n\n";
			exit(1);
		}
		if (intHM.find(cmd.cov[i]) != intHM.end()) {
			cerr << "\nERROR: Covariate " << cmd.cov[i] << " is specified as a covariate (--covar-names) and interaction covariate (--int-covar-names)." << "\n\n";
			exit(1);
		}
		if (phenoHeaderName.compare(cmd.cov[i]) == 0) {
			cerr << "\nERROR: Covariate " << cmd.cov[i] << " is also specified as the phenotype (--pheno-name)." << "\n\n";
			exit(1);
		}
		covHM[cmd.cov[i]] += 1;
		if (covHM[cmd.cov[i]] > 1) {
			covErrorString += cmd.cov[i] + " ";
			continue;
		}
		covSelHeadersName.push_back(cmd.cov[i]);
	}
	numSelCol = covSelHeadersName.size();

	if (!covErrorString.empty()) {
		cout << "\nWARNING: Covariates " + covErrorString + "is specified more than once.\n \t Removing the duplicate covariate names... \n\n";
	}








	/****************************************************
	 Print parameter info
	****************************************************/
	cout << "The Phenotype File is: " << cmd.phenoFile << "\n";
	cout << "The Selected Phenotype is: " << phenoHeaderName << '\n';
	cout << "Continuous or Binary? "; phenoTyp == 0 ? cout << "Continuous \n" : cout << "Binary \n";
	cout << "Robust or Non-Robust Analysis? "; robust == 0 ? cout << "Non-Robust \n\n" : cout << "Robust \n\n";

	if (numSelCol == 0) {
		cout << "No Covariates Selected" << "\n";
	}
	else {
		cout << "The Total Number of Selected Covariates is: " << numSelCol << '\n';
		cout << "The Selected Covariates are:  ";
		for (int i = 0; i < numSelCol; i++) {
			cout << covSelHeadersName[i] << "   ";
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
			cout << intCovSelHeadersName[i] << "   ";
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
			cout << expCovSelHeadersName[i] << "   ";
		}
		cout << "\n\n";
	}

	if (cmd.phenoType == 1) {
		cout << "Logistic Convergence Threshold: " << cmd.tol << "\n";
	}
	cout << "Minor Allele Frequency Threshold: " << cmd.MAF << "\n";
	cout << "Number of Threads: " << cmd.threads << "\n";
	cout << "Output File: " << cmd.outFile << "\n";
	cout << "*********************************************************\n";





	// Rearranging covSelHeaders
	numSelCol = numSelCol + numIntSelCol + numExpSelCol;
	vector<int> colSelVec(numSelCol);
	for (int i = numExpSelCol - 1; i >= 0; i--) { covSelHeadersName.insert(covSelHeadersName.begin(), expCovSelHeadersName[i]); }
	if (numIntSelCol != 0) {
		for (int i = numIntSelCol - 1; i >= 0; i--) { covSelHeadersName.insert(covSelHeadersName.begin(), intCovSelHeadersName[i]); }
	}




	// Start clock
	auto wall0 = std::chrono::system_clock::now();
	std::clock_t cpu0 = std::clock();

	/***************************************************
	 Reading phenotype file headers
	***************************************************/
	std::unordered_map<string, int> colNames;

	// Process header line: store column names and assign corresponding column numbers.
	string phenopath(cmd.phenoFile);
	long unsigned int  phenoncols;

	std::ifstream finph;
	finph.open(phenopath);
	if (!finph.is_open()) {
		cerr << "\nERROR: Cannot open phenotype file. \n\n" << endl;
		exit(1);
	}

	string line;
	getline(finph, line);
	std::istringstream issHead(line);
	string headerName;
	int header_i = 0;

	while (getline(issHead, headerName, delim)) {
		if (colNames.find(headerName) != colNames.end()) {
			cerr << "\nERROR: There are duplicate header names (" << headerName << ") in the phenotype file.\n\n";
			exit(1);
		}

		headerName.erase(std::remove(headerName.begin(), headerName.end(), '\r'), headerName.end());
		colNames[headerName] = header_i;
		++header_i;
	}


	phenoncols = colNames.size();
	if (colNames.find(phenoHeaderName) == colNames.end()) {
		cerr << "\nERROR: Cannot find phenotype column " << phenoHeaderName << " in phenotype file. \n\n";
		exit(1);
	}
	else {
		phenoCol = colNames[phenoHeaderName];
	}

	if (colNames.find(samIDHeaderName) == colNames.end()) {
		cerr << "\nERROR: Cannot find sample ID column " << samIDHeaderName << " in phenotype file. \n\n";
		exit(1);
	}
	else {
		samIDCol = colNames[samIDHeaderName];
	}

	for (int i = 0; i < numSelCol; i++) {
		if (colNames.find(covSelHeadersName[i]) == colNames.end()) {
			cerr << "\nERROR: Cannot find covariate column " << covSelHeadersName[i] << " in phenotype file. \n\n";
			exit(1);
		}
		else {
			colSelVec[i] = colNames[covSelHeadersName[i]];
		}
	}

	if (numIntSelCol != 0) {
		for (int i = 0; i < numIntSelCol; i++) {
			if (colNames.find(intCovSelHeadersName[i]) == colNames.end()) {
				cerr << "\nERROR: Cannot find interaction covariate column " << intCovSelHeadersName[i] << " in phenotype file. \n\n";
				exit(1);
			}
		}
	}

	for (int i = 0; i < numExpSelCol; i++) {
		if (colNames.find(expCovSelHeadersName[i]) == colNames.end()) {
			cerr << "\nERROR: Cannot find exposure column " << expCovSelHeadersName[i] << " in phenotype file. \n\n";
			exit(1);
		}
	}

	// Erase all elements and leaving it with a size of 0.
	colNames.clear();







	// count sample size
	int nrows = 0;
	while (getline(finph, line)) nrows++;
	samSize = nrows;
	finph.clear();
	finph.seekg(0, finph.beg);
	getline(finph, line);
	cout << "Before ID Matching and checking missing values... \n";
	cout << "Size of the phenotype vector is: " << samSize << " X 1\n";
	cout << "Size of the selected covariate matrix (including first column for interception values) is: " << samSize << " X " << numSelCol + 1 << '\n';


	// initialize data matrix
	vector <double> phenodata(samSize);
	vector <string> sampleIds(samSize);
	vector <double> covdata(samSize * (numSelCol + 1));
	//  vector <int> genoUnMatchID;// geno index with unmathcing ID
	unordered_set<int> genoUnMatchID;
	// A Hashmap phenodata for IDMatching process.
	// key is smapleID in phenotype file,
	// values are pheno-data, and the selected variables values.
	unordered_map<string, vector<string>> phenomap;
	//  if (IDMatching != 0 && IDMatching != 1) {
	//    cerr << "ERROR: ID Matching parameter is either 0 or 1.\n";
	//    exit(1);
	//  }
	  // read and store data
	for (int r = 0; r < samSize; r++) {
		getline(finph, line);
		std::istringstream iss(line);
		string value;
		string temvalue;
		vector <string> values;
		while (getline(iss, value, delim)) values.push_back(value);
		if (values.size() != phenoncols) {
			cerr << "ERROR: Wrong number of entries in data row:\n";
			cerr << line << '\n';
			cerr << "Expected " << phenoncols << " fields; parsed " << values.size() << '\n';
			exit(1);
		}
		sscanf(values[phenoCol].c_str(), "%lf", &phenodata[r]);
		sampleIds[r] = values[samIDCol];
		covdata[r * (numSelCol + 1)] = 1.0;
		for (int c = 0; c < numSelCol; c++)
			sscanf(values[colSelVec[c]].c_str(), "%lf", &covdata[r * (numSelCol + 1) + c + 1]);

		// IDMatching
	//    if (IDMatching == 1) {
		phenomap[values[samIDCol]] = { values[phenoCol] };
		for (int c = 0; c < numSelCol; c++)
			phenomap[values[samIDCol]].push_back(values[colSelVec[c]]);
		//    }
	}
	finph.close();
	if (nrows != samSize) {
		cerr << "ERROR: Wrong number of total row numbers:\n";
		cerr << "Expected " << samSize << " sample numbers; while reading row number is " << nrows << '\n';
		cerr << "Please also check the header line is at the top of the pheno data file! \n";
		exit(1);
	}
	cout << "End of reading phenotype and covariate data. \n";
	cout << "*********************************************************\n";





	if (cmd.usePgenFile) {

		Pgen pgen;

		pgen.processPgenHeader(cmd.pgenFile);
		pgen.processPvar(pgen, cmd.pvarFile);
		pgen.processPsam(pgen, cmd.psamFile, phenomap, phenoMissingKey, phenodata, covdata, numSelCol, samSize);

		sampleIds.clear();
		phenomap.clear(); 
		phenodata.clear();
		covdata.clear();

		samSize = pgen.new_samSize;
		double* phenoY = &pgen.new_phenodata[0];
		double* covX = &pgen.new_covdata[0];

		vector <double> residvec(samSize);
		// for logistic regression
		vector <double> miu(samSize);
		int Check = 1; // convergence condition of beta^(i+1) - beta^(i)
		int iter = 1;

		cout << "Precalculations and fitting null model..." << endl;
		auto start_time = std::chrono::high_resolution_clock::now();
		// transpose(X) * X
		double* XTransX = new double[(numSelCol + 1) * (numSelCol + 1)];
		matTmatprod(covX, covX, XTransX, samSize, numSelCol + 1, numSelCol + 1);
		// invert (XTransX)
		matInv(XTransX, numSelCol + 1);
		// transpose(X) * Y
		double* XTransY = new double[(numSelCol + 1)];
		matTvecprod(covX, phenoY, XTransY, samSize, numSelCol + 1);
		// beta = invert(XTransX) * XTransY
		double* beta = new double[(numSelCol + 1)];
		matvecprod(XTransX, XTransY, beta, numSelCol + 1, numSelCol + 1);

		while ((phenoTyp == 1) && (Check != (numSelCol + 1))) { // logistic regression
			iter++;
			// X * beta
			double* XbetaFL = new double[samSize];
			matvecprod(covX, beta, XbetaFL, samSize, numSelCol + 1);
			double* Yip1 = new double[samSize];
			// W * X and W * Y
			double* WX = new double[samSize * (numSelCol + 1)];
			double* WYip1 = new double[samSize];
			for (int i = 0; i < samSize; i++) {
				miu[i] = exp(XbetaFL[i]) / (1.0 + exp(XbetaFL[i]));
				Yip1[i] = XbetaFL[i] + (phenoY[i] - miu[i]) / (miu[i] * (1 - miu[i]));
				WYip1[i] = miu[i] * (1 - miu[i]) * Yip1[i];
				for (int j = 0; j < numSelCol + 1; j++) {
					WX[i * (numSelCol + 1) + j] = miu[i] * (1 - miu[i]) * covX[i * (numSelCol + 1) + j];
				}
			}
			// transpose(X) * WX
			matTmatprod(covX, WX, XTransX, samSize, numSelCol + 1, numSelCol + 1);
			// invert (XTransX)
			matInv(XTransX, numSelCol + 1);
			// transpose(X) * WYip1
			matTvecprod(covX, WYip1, XTransY, samSize, numSelCol + 1);
			// beta = invert(XTransX) * XTransY
			double* betaT = new double[(numSelCol + 1)];
			matvecprod(XTransX, XTransY, betaT, numSelCol + 1, numSelCol + 1);
			Check = 0;
			for (int i = 0; i < numSelCol + 1; i++) {
				if (std::abs(betaT[i] - beta[i]) <= epsilon) Check++;
				beta[i] = betaT[i];
			}

			delete[] Yip1;
			delete[] WYip1;
			delete[] WX;
			delete[] XbetaFL;
			delete[] betaT;
		}

		// X * beta
		double* Xbeta = new double[samSize];
		matvecprod(covX, beta, Xbeta, samSize, numSelCol + 1);

		// X*[invert (XTransX)]
		double* XinvXTX = new double[samSize * (numSelCol + 1)];

		if (phenoTyp == 1) {
			double* WX = new double[samSize * (numSelCol + 1)];
			for (int i = 0; i < samSize; i++) {
				miu[i] = exp(Xbeta[i]) / (1.0 + exp(Xbeta[i]));
				Xbeta[i] = miu[i];
				for (int j = 0; j < numSelCol + 1; j++) {
					WX[i * (numSelCol + 1) + j] = miu[i] * (1.0 - miu[i]) * covX[i * (numSelCol + 1) + j];
				}
			}
			// transpose(X) * WX
			matTmatprod(covX, WX, XTransX, samSize, numSelCol + 1, numSelCol + 1);
			// invert (XTransX)
			matInv(XTransX, numSelCol + 1);
			matmatprod(WX, XTransX, XinvXTX, samSize, numSelCol + 1, numSelCol + 1);
			delete[] WX;

			cout << "Logistic regression reaches convergence after " << iter << " steps...\n";
		}
		else {
			matmatprod(covX, XTransX, XinvXTX, samSize, numSelCol + 1, numSelCol + 1);
		}


		// residual = Y - X * beta
		double sigma2 = 0;
		for (int i = 0; i < samSize; i++) {
			residvec[i] = phenoY[i] - Xbeta[i];
			sigma2 += residvec[i] * residvec[i];
		}
		// sqr(sigma) = transpose(resid)*resid/[samSize-(numSelCol+1)]
		sigma2 = sigma2 / (samSize - (numSelCol + 1));
		if (phenoTyp == 1) sigma2 = 1.0;

		cout << "Execution time... ";
		auto end_time = std::chrono::high_resolution_clock::now();
		printExecutionTime(start_time, end_time);
		cout << "Done.\n";
		cout << "*********************************************************\n";

		double* resid = &residvec[0];

		delete[] XTransY;
		delete[] beta;
		delete[] Xbeta;

		cout << "Streaming SNPs for speeding up GWAS analysis in parallel. \n";
		cout << "Number of SNPs in each batch is: " << stream_snps << "\n\n";

		pgen.numSelCol = numSelCol;
		pgen.numIntSelCol = numIntSelCol;
		pgen.numExpSelCol = numExpSelCol;
		pgen.Sq = Sq;
		pgen.robust = robust;
		pgen.stream_snps = stream_snps;
		pgen.maf = cmd.MAF;
		pgen.missGeno = cmd.missGenoRate;
		pgen.miu = miu;
		pgen.phenoTyp = phenoTyp;
		pgen.covX = covX;
		pgen.XinvXTX = XinvXTX;
		pgen.resid = resid;
		pgen.sigma2 = sigma2;
		pgen.outFile = cmd.outFile;


		if (!cmd.doFilters) {
			//Preparing for parallelizing of BGEN file
			uint32_t threads = cmd.threads;
			pgen.threads = threads;
			if (pgen.raw_variant_ct < threads) {
				threads = pgen.raw_variant_ct;
				cout << "Number of variants (" << pgen.raw_variant_ct << ") is less than the number of specified threads (" << threads << ")...\n";
				cout << "Using " << threads << " thread(s) instead... \n\n";
			}

			pgen.begin.resize(threads);
			pgen.end.resize(threads);
			if (threads > 1) {
				cout << "The second allele in the PGEN file will be used for association testing.\n";
				cout << "Running multithreading...\n";
			}
			else {
				cout << "The second allele in the PGEN file will be used for association testing.\n";
				cout << "Running with single thread...\n";
			}

			for (uint32_t t = 0; t < threads; t++) {
				pgen.begin[t] = floor((pgen.raw_variant_ct / threads) * t);

				if ((t + 1) == (threads)) {
					pgen.end[t] = pgen.raw_variant_ct - 1;
				}
				else {
					pgen.end[t] = floor(((pgen.raw_variant_ct / threads) * (t + 1)) - 1);
				}
			}
		}
		else {
			
		}


		boost::thread_group thread_grp;
		start_time = std::chrono::high_resolution_clock::now();
		for (uint32_t i = 0; i < pgen.threads; ++i) {
			thread_grp.create_thread(boost::bind(&gemPGEN, pgen.begin[i], pgen.end[i], cmd.pgenFile, cmd.pvarFile, i, boost::ref(pgen)));
		}
		thread_grp.join_all();
		cout << "Joining threads... \n";
		end_time = std::chrono::high_resolution_clock::now();
		cout << "Execution time... ";
		printExecutionTime(start_time, end_time);
		cout << "Done. \n";
		cout << "*********************************************************\n";



		// Write all results from each thread to 1 file
		cout << "Combining results... \n";
		start_time = std::chrono::high_resolution_clock::now();
		std::ofstream results(output, std::ofstream::binary);
		results << "ID" << "\t" << "CHROM" << "\t" << "POS" << "\t" << "Non_Effect_Allele" << "\t" << "Effect_Allele" << "\t" << "N_Samples" << "\t" << "AF" << "\t" << "Beta_Marginal" << "\t" << "Var_Beta_Marginal" << "\t";
		for (int i = 1; i <= numExpSelCol; i++) {
			results << "Beta_Interaction" << "_" << i << "\t";
		}
		for (int i = 1; i <= numExpSelCol; i++) {
			for (int j = 1; j <= numExpSelCol; j++) {
				results << "Var_Beta_Interaction" << "_" << i << "_" << j << "\t";
			}
		}
		results << "P_Value_Marginal" << "\t" << "P_Value_Interaction" << "\t" << "P_Value_Joint\n";


		for (uint32_t i = 0; i < pgen.threads; i++) {
			std::string threadOutputFile = cmd.outFile + "_bin_" + std::to_string(i) + ".tmp";
			std::ifstream thread_output(threadOutputFile);
			results << thread_output.rdbuf();

			thread_output.close();
			boost::filesystem::remove(threadOutputFile.c_str());

		}
		results.close();
		end_time = std::chrono::high_resolution_clock::now();
		cout << "Execution time... ";
		printExecutionTime(start_time, end_time);
		cout << "Done. \n";



		// Finished
		cout << "*********************************************************\n";
		std::chrono::duration<double> wallduration = (std::chrono::system_clock::now() - wall0);
		double cpuduration = (std::clock() - cpu0) / (double)CLOCKS_PER_SEC;
		cout << "Total Wall Time = " << wallduration.count() << "  Seconds\n";
		cout << "Total CPU Time = " << cpuduration << "  Seconds\n";
		//cout << "Execution Wall Time = " << exetime << "  Seconds\n";
		cout << "*********************************************************\n";


		delete[] XTransX;
		delete[] XinvXTX;

	}


	if (cmd.useBedFile) {

		Bed bed;
		bed.processBed(cmd.bedFile, cmd.bimFile, cmd.famFile);
		bed.processFam(bed, cmd.famFile, phenomap, phenoMissingKey, phenodata, covdata, numSelCol, samSize);

		sampleIds.clear();
		phenomap.clear(); 
		phenodata.clear();
		covdata.clear();

		samSize = bed.new_samSize;
		double* phenoY = &bed.new_phenodata[0];
		double* covX = &bed.new_covdata[0];
		vector <double> residvec(samSize);

		// for logistic regression
		vector <double> miu(samSize);
		int Check = 1; // convergence condition of beta^(i+1) - beta^(i)
		int iter = 1;

		cout << "Precalculations and fitting null model..." << endl;
		auto start_time = std::chrono::high_resolution_clock::now();
		// transpose(X) * X
		double* XTransX = new double[(numSelCol + 1) * (numSelCol + 1)];
		matTmatprod(covX, covX, XTransX, samSize, numSelCol + 1, numSelCol + 1);
		// invert (XTransX)
		matInv(XTransX, numSelCol + 1);
		// transpose(X) * Y
		double* XTransY = new double[(numSelCol + 1)];
		matTvecprod(covX, phenoY, XTransY, samSize, numSelCol + 1);
		// beta = invert(XTransX) * XTransY
		double* beta = new double[(numSelCol + 1)];
		matvecprod(XTransX, XTransY, beta, numSelCol + 1, numSelCol + 1);

		while ((phenoTyp == 1) && (Check != (numSelCol + 1))) { // logistic regression
			iter++;
			// X * beta
			double* XbetaFL = new double[samSize];
			matvecprod(covX, beta, XbetaFL, samSize, numSelCol + 1);
			double* Yip1 = new double[samSize];
			// W * X and W * Y
			double* WX = new double[samSize * (numSelCol + 1)];
			double* WYip1 = new double[samSize];
			for (int i = 0; i < samSize; i++) {
				miu[i] = exp(XbetaFL[i]) / (1.0 + exp(XbetaFL[i]));
				Yip1[i] = XbetaFL[i] + (phenoY[i] - miu[i]) / (miu[i] * (1 - miu[i]));
				WYip1[i] = miu[i] * (1 - miu[i]) * Yip1[i];
				for (int j = 0; j < numSelCol + 1; j++) {
					WX[i * (numSelCol + 1) + j] = miu[i] * (1 - miu[i]) * covX[i * (numSelCol + 1) + j];
				}
			}
			// transpose(X) * WX
			matTmatprod(covX, WX, XTransX, samSize, numSelCol + 1, numSelCol + 1);
			// invert (XTransX)
			matInv(XTransX, numSelCol + 1);
			// transpose(X) * WYip1
			matTvecprod(covX, WYip1, XTransY, samSize, numSelCol + 1);
			// beta = invert(XTransX) * XTransY
			double* betaT = new double[(numSelCol + 1)];
			matvecprod(XTransX, XTransY, betaT, numSelCol + 1, numSelCol + 1);
			Check = 0;
			for (int i = 0; i < numSelCol + 1; i++) {
				if (std::abs(betaT[i] - beta[i]) <= epsilon) Check++;
				beta[i] = betaT[i];
			}

			delete[] Yip1;
			delete[] WYip1;
			delete[] WX;
			delete[] XbetaFL;
			delete[] betaT;
		}

		// X * beta
		double* Xbeta = new double[samSize];
		matvecprod(covX, beta, Xbeta, samSize, numSelCol + 1);

		// X*[invert (XTransX)]
		double* XinvXTX = new double[samSize * (numSelCol + 1)];

		if (phenoTyp == 1) {
			double* WX = new double[samSize * (numSelCol + 1)];
			for (int i = 0; i < samSize; i++) {
				miu[i] = exp(Xbeta[i]) / (1.0 + exp(Xbeta[i]));
				Xbeta[i] = miu[i];
				for (int j = 0; j < numSelCol + 1; j++) {
					WX[i * (numSelCol + 1) + j] = miu[i] * (1.0 - miu[i]) * covX[i * (numSelCol + 1) + j];
				}
			}
			// transpose(X) * WX
			matTmatprod(covX, WX, XTransX, samSize, numSelCol + 1, numSelCol + 1);
			// invert (XTransX)
			matInv(XTransX, numSelCol + 1);
			matmatprod(WX, XTransX, XinvXTX, samSize, numSelCol + 1, numSelCol + 1);
			delete[] WX;

			cout << "Logistic regression reaches convergence after " << iter << " steps...\n";
		}
		else {
			matmatprod(covX, XTransX, XinvXTX, samSize, numSelCol + 1, numSelCol + 1);
		}


		// residual = Y - X * beta
		double sigma2 = 0;
		for (int i = 0; i < samSize; i++) {
			residvec[i] = phenoY[i] - Xbeta[i];
			sigma2 += residvec[i] * residvec[i];
		}
		// sqr(sigma) = transpose(resid)*resid/[samSize-(numSelCol+1)]
		sigma2 = sigma2 / (samSize - (numSelCol + 1));
		if (phenoTyp == 1) sigma2 = 1.0;

		cout << "Execution time... ";
		auto end_time = std::chrono::high_resolution_clock::now();
		printExecutionTime(start_time, end_time);
		cout << "Done.\n";
		cout << "*********************************************************\n";

		double* resid = &residvec[0];

		delete[] XTransY;
		delete[] beta;
		delete[] Xbeta;

		cout << "Streaming SNPs for speeding up GWAS analysis in parallel. \n";
		cout << "Number of SNPs in each batch is: " << stream_snps << "\n\n";

		bed.numSelCol = numSelCol;
		bed.numIntSelCol = numIntSelCol;
		bed.numExpSelCol = numExpSelCol;
		bed.Sq = Sq;
		bed.robust = robust;
		bed.stream_snps = stream_snps;
		bed.maf = cmd.MAF;
		bed.missGeno = cmd.missGenoRate;
		bed.miu = miu;
		bed.phenoTyp = phenoTyp;
		bed.covX = covX;
		bed.XinvXTX = XinvXTX;
		bed.resid = resid;
		bed.sigma2 = sigma2;
		bed.outFile = cmd.outFile;

		if (!cmd.doFilters) {
			//Preparing for parallelizing of BGEN file
			uint32_t threads = cmd.threads;
			bed.threads = threads;
			if (bed.n_variants < threads) {
				threads = bed.n_variants;
				cout << "Number of variants (" << bed.n_variants << ") is less than the number of specified threads (" << threads << ")...\n";
				cout << "Using " << threads << " thread(s) instead... \n\n";
			}

			bed.begin.resize(threads);
			bed.end.resize(threads);
			if (threads > 1) {
				cout << "The SNP major allele in the BED file will be used for association testing.\n";
				cout << "Running multithreading...\n";
			}
			else {
				cout << "The SNP major allele in the BED file will be used for association testing.\n";
				cout << "Running with single thread...\n";
			}

			for (uint32_t t = 0; t < threads; t++) {
				bed.begin[t] = floor((bed.n_variants / threads) * t);

				if ((t + 1) == (threads)) {
					bed.end[t] = bed.n_variants - 1;
				}
				else {
					bed.end[t] = floor(((bed.n_variants / threads) * (t + 1)) - 1);
				}
			}
		}
		else {

		}


		boost::thread_group thread_grp;
		start_time = std::chrono::high_resolution_clock::now();
		for (uint32_t i = 0; i < bed.threads; ++i) {
			 thread_grp.create_thread(boost::bind(&gemBED, bed.begin[i], bed.end[i], cmd.bedFile, cmd.bimFile, i, boost::ref(bed)));
		}
		thread_grp.join_all();
		cout << "Joining threads... \n";
		end_time = std::chrono::high_resolution_clock::now();
		cout << "Execution time... ";
		printExecutionTime(start_time, end_time);
		cout << "Done. \n";
		cout << "*********************************************************\n";


		// Write all results from each thread to 1 file
		cout << "Combining results... \n";
		start_time = std::chrono::high_resolution_clock::now();
		std::ofstream results(output, std::ofstream::binary);
		results << "ID" << "\t" << "CHROM" << "\t" << "POS" << "\t" << "Non_Effect_Allele" << "\t" << "Effect_Allele" << "\t" << "N_Samples" << "\t" << "AF" << "\t" << "Beta_Marginal" << "\t" << "Var_Beta_Marginal" << "\t";
		for (int i = 1; i <= numExpSelCol; i++) {
			results << "Beta_Interaction" << "_" << i << "\t";
		}
		for (int i = 1; i <= numExpSelCol; i++) {
			for (int j = 1; j <= numExpSelCol; j++) {
				results << "Var_Beta_Interaction" << "_" << i << "_" << j << "\t";
			}
		}
		results << "P_Value_Marginal" << "\t" << "P_Value_Interaction" << "\t" << "P_Value_Joint\n";


		for (uint32_t i = 0; i < bed.threads; i++) {
			std::string threadOutputFile = cmd.outFile + "_bin_" + std::to_string(i) + ".tmp";
			std::ifstream thread_output(threadOutputFile);
			results << thread_output.rdbuf();

			thread_output.close();
			boost::filesystem::remove(threadOutputFile.c_str());

		}
		results.close();
		end_time = std::chrono::high_resolution_clock::now();
		cout << "Execution time... ";
		printExecutionTime(start_time, end_time);
		cout << "Done. \n";



		// Finished
		cout << "*********************************************************\n";
		std::chrono::duration<double> wallduration = (std::chrono::system_clock::now() - wall0);
		double cpuduration = (std::clock() - cpu0) / (double)CLOCKS_PER_SEC;
		cout << "Total Wall Time = " << wallduration.count() << "  Seconds\n";
		cout << "Total CPU Time = " << cpuduration << "  Seconds\n";
		//cout << "Execution Wall Time = " << exetime << "  Seconds\n";
		cout << "*********************************************************\n";


		delete[] XTransX;
		delete[] XinvXTX;

	}



	if (cmd.useBgenFile) {
		/***************************************************************
		  Read General information of bgen data.
		  Conduct Sample IDMatching process if necessary.
		***************************************************************/
		Bgen bgen;
		bgen.processBgenHeaderBlock(cmd.genofile);
		bgen.processBgenSampleBlock(bgen, cmd.samplefile, cmd.useSampleFile, phenomap, phenoMissingKey, phenodata, covdata, numSelCol, samSize);


		/******************************************************************
		  Genome-Wide Associate Study using Linear or Logistic regression
		  and Processing Geno Files (.bgen)
		******************************************************************/
		cout << "Starting GWAS... \n\n";
		samSize = bgen.new_samSize;
		double* phenoY = &bgen.new_phenodata[0];
		double* covX = &bgen.new_covdata[0];
		vector <double> residvec(samSize);

		// for logistic regression
		vector <double> miu(samSize);
		int Check = 1; // convergence condition of beta^(i+1) - beta^(i)
		int iter = 1;


		cout << "Precalculations and fitting null model..." << endl;
		auto start_time = std::chrono::high_resolution_clock::now();
		// transpose(X) * X
		double* XTransX = new double[(numSelCol + 1) * (numSelCol + 1)];
		matTmatprod(covX, covX, XTransX, samSize, numSelCol + 1, numSelCol + 1);
		// invert (XTransX)
		matInv(XTransX, numSelCol + 1);
		// transpose(X) * Y
		double* XTransY = new double[(numSelCol + 1)];
		matTvecprod(covX, phenoY, XTransY, samSize, numSelCol + 1);
		// beta = invert(XTransX) * XTransY
		double* beta = new double[(numSelCol + 1)];
		matvecprod(XTransX, XTransY, beta, numSelCol + 1, numSelCol + 1);

		
		while ((phenoTyp == 1) && (Check != (numSelCol + 1))) { // logistic regression
			iter++;
			// X * beta
			double* XbetaFL = new double[samSize];
			matvecprod(covX, beta, XbetaFL, samSize, numSelCol + 1);
			double* Yip1 = new double[samSize];
			// W * X and W * Y
			double* WX = new double[samSize * (numSelCol + 1)];
			double* WYip1 = new double[samSize];
			for (int i = 0; i < samSize; i++) {
				miu[i] = exp(XbetaFL[i]) / (1.0 + exp(XbetaFL[i]));
				Yip1[i] = XbetaFL[i] + (phenoY[i] - miu[i]) / (miu[i] * (1 - miu[i]));
				WYip1[i] = miu[i] * (1 - miu[i]) * Yip1[i];
				for (int j = 0; j < numSelCol + 1; j++) {
					WX[i * (numSelCol + 1) + j] = miu[i] * (1 - miu[i]) * covX[i * (numSelCol + 1) + j];
				}
			}
			// transpose(X) * WX
			matTmatprod(covX, WX, XTransX, samSize, numSelCol + 1, numSelCol + 1);
			// invert (XTransX)
			matInv(XTransX, numSelCol + 1);
			// transpose(X) * WYip1
			matTvecprod(covX, WYip1, XTransY, samSize, numSelCol + 1);
			// beta = invert(XTransX) * XTransY
			double* betaT = new double[(numSelCol + 1)];
			matvecprod(XTransX, XTransY, betaT, numSelCol + 1, numSelCol + 1);
			Check = 0;
			for (int i = 0; i < numSelCol + 1; i++) {
				if (std::abs(betaT[i] - beta[i]) <= epsilon) Check++;
				beta[i] = betaT[i];
			}

			delete[] Yip1;
			delete[] WYip1;
			delete[] WX;
			delete[] XbetaFL;
			delete[] betaT;
		}

		// X * beta
		double* Xbeta = new double[samSize];
		matvecprod(covX, beta, Xbeta, samSize, numSelCol + 1);

		// X*[invert (XTransX)]
		double* XinvXTX = new double[samSize * (numSelCol + 1)];

		if (phenoTyp == 1) {
			double* WX = new double[samSize * (numSelCol + 1)];
			for (int i = 0; i < samSize; i++) {
				miu[i] = exp(Xbeta[i]) / (1.0 + exp(Xbeta[i]));
				Xbeta[i] = miu[i];
				for (int j = 0; j < numSelCol + 1; j++) {
					WX[i * (numSelCol + 1) + j] = miu[i] * (1.0 - miu[i]) * covX[i * (numSelCol + 1) + j];
				}
			}
			// transpose(X) * WX
			matTmatprod(covX, WX, XTransX, samSize, numSelCol + 1, numSelCol + 1);
			// invert (XTransX)
			matInv(XTransX, numSelCol + 1);
			matmatprod(WX, XTransX, XinvXTX, samSize, numSelCol + 1, numSelCol + 1);
			delete[] WX;

			cout << "Logistic regression reaches convergence after " << iter << " steps...\n";
		}
		else {
			matmatprod(covX, XTransX, XinvXTX, samSize, numSelCol + 1, numSelCol + 1);
		}


		// residual = Y - X * beta
		double sigma2 = 0;
		for (int i = 0; i < samSize; i++) {
			residvec[i] = phenoY[i] - Xbeta[i];
			sigma2 += residvec[i] * residvec[i];
		}
		// sqr(sigma) = transpose(resid)*resid/[samSize-(numSelCol+1)]
		sigma2 = sigma2 / (samSize - (numSelCol + 1));
		if (phenoTyp == 1) sigma2 = 1.0;

		cout << "Execution time... ";
		auto end_time = std::chrono::high_resolution_clock::now();
		printExecutionTime(start_time, end_time);
		cout << "Done.\n";
		cout << "*********************************************************\n";

		double* resid = &residvec[0];

		delete[] XTransY;
		delete[] beta;
		delete[] Xbeta;

		cout << "Streaming SNPs for speeding up GWAS analysis in parallel. \n";
		cout << "Number of SNPs in each batch is: " << stream_snps << "\n\n";


		bgen.numSelCol = numSelCol;
		bgen.numIntSelCol = numIntSelCol;
		bgen.numExpSelCol = numExpSelCol;
		bgen.Sq = Sq;
		bgen.robust = robust;
		bgen.stream_snps = stream_snps;
		bgen.maf = cmd.MAF;
		bgen.missGeno = cmd.missGenoRate;
		bgen.miu = miu;
		bgen.phenoTyp = phenoTyp;
		bgen.covX = covX;
		bgen.XinvXTX = XinvXTX;
		bgen.resid = resid;
		bgen.sigma2 = sigma2;
		bgen.outFile = cmd.outFile;


		start_time = std::chrono::high_resolution_clock::now();
		bgen.getPositionOfBgenVariant(bgen, cmd);
		end_time = std::chrono::high_resolution_clock::now();
		cout << "Execution time... ";
		printExecutionTime(start_time, end_time);
		cout << "Done.\n";
		cout << "*********************************************************\n";


		//Preparing for parallelizing of BGEN file
		if (cmd.threads > 1) {
			cout << "The second allele in the BGEN file will be used for association testing.\n";
			cout << "Running multithreading...\n";
		}
		else {
			cout << "The second allele in the BGEN file will be used for association testing.\n";
			cout << "Running with single thread...\n";
		}


		boost::thread_group thread_grp;
		start_time = std::chrono::high_resolution_clock::now();
		for (uint i = 0; i < bgen.threads; ++i) {
			thread_grp.create_thread(boost::bind(&BgenParallelGWAS, bgen.Mbgen_begin[i], bgen.Mbgen_end[i], bgen.bgenVariantPos[i], bgen.keepVariants[i], cmd.genofile, bgen.filterVariants, i, boost::ref(bgen)));
		}
		thread_grp.join_all();
		cout << "Joining threads... \n";
		end_time = std::chrono::high_resolution_clock::now();
		cout << "Execution time... ";
		printExecutionTime(start_time, end_time);
		cout << "Done. \n";
		cout << "*********************************************************\n";





		// Write all results from each thread to 1 file
		cout << "Combining results... \n";
		start_time = std::chrono::high_resolution_clock::now();
		std::ofstream results(output, std::ofstream::binary);
		results << "SNPID" << "\t" << "rsID" << "\t" << "CHR" << "\t" << "POS" << "\t" << "Non_Effect_Allele" << "\t" << "Effect_Allele" << "\t" << "N_Samples" << "\t" << "AF" << "\t" << "Beta_Marginal" << "\t" << "Var_Beta_Marginal" << "\t";
		for (int i = 1; i <= numExpSelCol; i++) {
			results << "Beta_Interaction" << "_" << i << "\t";
		}
		for (int i = 1; i <= numExpSelCol; i++) {
			for (int j = 1; j <= numExpSelCol; j++) {
				results << "Var_Beta_Interaction" << "_" << i << "_" << j << "\t";
			}
		}
		results << "P_Value_Marginal" << "\t" << "P_Value_Interaction" << "\t" << "P_Value_Joint\n";


		for (uint i = 0; i < bgen.threads; i++) {
			std::string threadOutputFile = cmd.outFile + "_bin_" + std::to_string(i) + ".tmp";
			std::ifstream thread_output(threadOutputFile);
			results << thread_output.rdbuf();

			thread_output.close();
			boost::filesystem::remove(threadOutputFile.c_str());

		}
		results.close();
		end_time = std::chrono::high_resolution_clock::now();
		cout << "Execution time... ";
		printExecutionTime(start_time, end_time);
		cout << "Done. \n";






		// Finished
		cout << "*********************************************************\n";
		std::chrono::duration<double> wallduration = (std::chrono::system_clock::now() - wall0);
		double cpuduration = (std::clock() - cpu0) / (double)CLOCKS_PER_SEC;
		cout << "Total Wall Time = " << wallduration.count() << "  Seconds\n";
		cout << "Total CPU Time = " << cpuduration << "  Seconds\n";
		//cout << "Execution Wall Time = " << exetime << "  Seconds\n";
		cout << "*********************************************************\n";



		delete[] XTransX;
		delete[] XinvXTX;
	}

	return 0;
}


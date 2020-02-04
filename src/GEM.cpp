/*  GEM : Gene-Environment interaction analysis for Millions of samples
 *  Copyright (C) 2018,2019  Liang Hong, Han Chen
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


int main(int argc, char* argv[]) {

	char paramfile[300];
	char phenofile[300];
	char genofile[300];
	char samplefile[300];



	// Process command line
	CommandLine cmd;
	cmd.processCommandLine(argc, argv);




	/****************************************************
	  Call subroutine to read parameters from a file
	****************************************************/
	PARAMETERS prm;
	strcpy(paramfile, cmd.pFile.c_str());
	ReadParameters(paramfile, genofile, phenofile, samplefile, &prm);





	/****************************************************
	Parameters
	****************************************************/
	// int samSize = prm.samSize;
	// int IDMatching = prm.IDMatching;

	double exetime = 0;
	double epsilon = prm.epsilon;

	int samSize;
	int phenoCol;
	int samIDCol;
	int Sq = prm.Sq;
	int robust      = prm.robust;
	int phenoTyp    = prm.phenoTyp;
	int stream_snps = prm.stream_snps;
	int numSelCol   = prm.covSelHeaders.size();

	string phenoHeaderName(prm.phenoHeader);
	string samIDHeaderName(prm.samIDHeader);
	string output(prm.outputfile);

	vector<int> colSelVec(numSelCol);
	vector<string> covSelHeadersName;
	for (int i = 0; i < numSelCol; i++) covSelHeadersName.push_back(prm.covSelHeaders[i]);











	/***************************************************
	Processing pheno file
	***************************************************/
	auto wall0 = std::chrono::system_clock::now();
	std::clock_t cpu0 = std::clock();

	string phenoMissingKey(prm.MissingKey);
	std::unordered_map<string, int> colNames;
	int phenoncols;
	//  char delim = ' '; // deliminator of pheno data file
	char delim = prm.delim_pheno[0];
	if (delim == '\0') delim = ' ';

	// process header line: store column names
	// and assign corresponding column numbers.
	std::ifstream finph;
	string phenopath(phenofile);
	finph.open(phenopath);
	string line;
	getline(finph, line);
	std::istringstream issHead(line);
	string headerName;
	int header_i = 0;

	while (getline(issHead, headerName, delim)) {
		if (colNames.find(headerName) != colNames.end()) {
			cerr << "There has duplicate header names.\n";
			exit(1);
		}

		headerName.erase(std::remove(headerName.begin(), headerName.end(), '\r'), headerName.end());
		colNames[headerName] = header_i;
		++header_i;
	}


	phenoncols = colNames.size();
	if (colNames.find(phenoHeaderName) == colNames.end()) {
		cerr << "Pheno header name is wrong.\n";
		exit(1);
	}
	else
		phenoCol = colNames[phenoHeaderName];
	if (colNames.find(samIDHeaderName) == colNames.end()) {
		cerr << "Sample ID header name is wrong.\n";
		exit(1);
	}
	else
		samIDCol = colNames[samIDHeaderName];

	for (int i = 0; i < numSelCol; i++) {
		if (colNames.find(covSelHeadersName[i]) == colNames.end()) {
			cerr << "Covariate header name is wrong.\n";
			exit(i);
		}
		else
			colSelVec[i] = colNames[covSelHeadersName[i]];
	}

	// Erase all elements and leaving it with a size of 0.
	colNames.clear();

	// print out header names and select pheno columns
	cout << "*********************************************************\n";
	cout << "The Selected Phenotype Data Header Name is: " << phenoHeaderName << '\n';
	cout << "Linear or Binary? ";
	phenoTyp == 0 ? cout << "Linear \n" : cout << "Binary \n";
	cout << "Robust or Non-Robust Analysis? ";
	robust == 0 ? cout << "Non-Robust \n" : cout << "Robust \n";
	cout << "The Total Number of Selected Covariates is: " << numSelCol << '\n';
	cout << "The Selected Covariate Column Header Names are:  ";
	for (int i = 0; i < numSelCol; i++) {
		cout << covSelHeadersName[i] << "   ";
	}
	cout << '\n';
	cout << "MAF filtering threshold: " << cmd.MAF << endl;
	//  cout << "Check Missing Values in Pheno Data File And Match of Order Sequence of Sample IDs?\n";
	//  IDMatching == 0 ? cout << "No Chekcing! \n" : cout << "Yes, Checking Please! \n"; 

	  // count sample size
	int nrows = 0;
	while (getline(finph, line)) nrows++;
	samSize = nrows;
	finph.clear();
	finph.seekg(0, finph.beg);
	getline(finph, line);
	cout << "*********************************************************\n";
	cout << "Before ID Matching and checking missing values: \n";
	cout << "Size of the pheno vector is: " << samSize << " X 1\n";
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
	cout << "End of reading pheno and covariate data. \n";
	cout << "*********************************************************\n";











	/***************************************************************
	  Read General information of bgen data.
	  Conduct Sample IDMatching process if necessary.
	***************************************************************/
	Bgen bgen;
	bgen.processBgenHeaderBlock(genofile);



	uint maxLA  = 65536;
	char* samID = new char[maxLA + 1];


	/**** sample identifier block ********/
	/* exists when SampleIdentifiers = 1 */
	if (bgen.SampleIdentifiers == 1) {
		uint LS1; fread(&LS1, 4, 1, bgen.fin); // std::cout << "LS1: " << LS1 << std::endl; // LS1 + L_H <= offset
		uint Nrow; fread(&Nrow, 4, 1, bgen.fin); // cout << "Nrow: " << Nrow << " " << std::flush;
		if (Nrow != bgen.Nbgen) {
			cerr << "ERROR: Nrow = " << Nrow << " does not match Nbgen = " << bgen.Nbgen << '\n';
			exit(1);
		}


		int k = 0;
		for (uint m = 0; m < bgen.Nbgen; m++) {
			ushort LSID; fread(&LSID, 2, 1, bgen.fin); // std::cout << "LSID: " << LSID << " ";
			fread(samID, 1, LSID, bgen.fin); // std::cout << "samID: " << samID << " " << std::endl;

			// IDMatching
			string strtmp(samID);
			int itmp = k;
			if (phenomap.find(strtmp) != phenomap.end()) {
				auto tmp_valvec = phenomap[strtmp];
				if (find(tmp_valvec.begin(), tmp_valvec.end(), phenoMissingKey) == tmp_valvec.end()) {
					sscanf(tmp_valvec[0].c_str(), "%lf", &phenodata[k]);
					// covdata[k*(numSelCol+1) + 0] = 1.0;
					for (int c = 0; c < numSelCol; c++)
						sscanf(tmp_valvec[c + 1].c_str(), "%lf", &covdata[k * (numSelCol + 1) + c + 1]);
					k++;
				}
				// erase the used element in phenomap
				phenomap.erase(strtmp);
			}
			// save the index with unmatched ID into genoUnMatchID.
	//        if (itmp == k) genoUnMatchID.push_back(m);
			if (itmp == k) genoUnMatchID.insert(m);
		}
		// After IDMatching, resizing phenodata and covdata, and updating samSize;
		phenodata.resize(k);
		covdata.resize(k * (numSelCol + 1));
		samSize = k;
		cout << "****************************************************************************\n";
		if (genoUnMatchID.empty())
			cout << "After processes of sample IDMatching and checking missing values, the sample size does not change.\n";
		else
			cout << "After processes of sample IDMatching and checking missing values, the sample size changes from "
			<< samSize + genoUnMatchID.size() << " to " << samSize << ".\n";

		cout << "****************************************************************************\n";
		cout << "Sample IDMatching and checking missing values processes have been completed.\n";
		cout << "New pheno and covariate data vectors with the same order of sample ID sequence of geno data are updated.\n";
		cout << "****************************************************************************\n";
		//      } // end of if IDMatching == 1
	} // end SampleIdentifiers == 1
	else {
		std::ifstream fIDMat;
		fIDMat.open(samplefile);
		string IDline;
		getline(fIDMat, IDline);
		getline(fIDMat, IDline);

		int k = 0;
		for (uint m = 0; m < bgen.Nbgen; m++) {
			// IDMatching
			getline(fIDMat, IDline);
			std::istringstream iss(IDline);
			string strtmp;
			iss >> strtmp;
			//        string strtmp(samID);
			int itmp = k;
			if (phenomap.find(strtmp) != phenomap.end()) {
				auto tmp_valvec = phenomap[strtmp];
				if (find(tmp_valvec.begin(), tmp_valvec.end(), phenoMissingKey) == tmp_valvec.end()) {
					sscanf(tmp_valvec[0].c_str(), "%lf", &phenodata[k]);
					// covdata[k*(numSelCol+1) + 0] = 1.0;
					for (int c = 0; c < numSelCol; c++)
						sscanf(tmp_valvec[c + 1].c_str(), "%lf", &covdata[k * (numSelCol + 1) + c + 1]);
					k++;
				}
				// erase the used element in phenomap
				phenomap.erase(strtmp);
			}
			// save the index with unmatched ID into genoUnMatchID.
	  //        if (itmp == k) genoUnMatchID.push_back(m);
			if (itmp == k) genoUnMatchID.insert(m);
		}
		// After IDMatching, resizing phenodata and covdata, and updating samSize;
		fIDMat.close();

		phenodata.resize(k);
		covdata.resize(k * (numSelCol + 1));
		samSize = k;
		cout << "****************************************************************************\n";
		if (genoUnMatchID.empty())
			cout << "After processes of sample IDMatching and checking missing values, the sample size does not change.\n";
		else
			cout << "After processes of sample IDMatching and checking missing values, the sample size changes from "
			<< samSize + genoUnMatchID.size() << " to " << samSize << ".\n";

		cout << "****************************************************************************\n";
		cout << "Sample IDMatching and checking missing values processes have been completed.\n";
		cout << "New pheno and covariate data vectors with the same order of sample ID sequence of geno data are updated.\n";
		cout << "****************************************************************************\n";
	}
	sampleIds.clear(); // clear memory
	phenomap.clear(); // clear phenomap
	// end of sample identifier block










	/******************************************************************
	  Genome-Wide Associate Study using Linear or Logistic regression
	  and Processing Geno Files (.bgen)
	******************************************************************/
	cout << "Starting GWAS. \n";
	double* phenoY = &phenodata[0];
	double* covX   = &covdata[0];
	vector <double> residvec(samSize);

	// for logistic regression
	vector <double> miu(samSize);
	int Check = 1; // convergence condition of beta^(i+1) - beta^(i)
	int iter = 1;


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
		delete[] WX;

		cout << "Logistic regression reaches convergence after " << iter << " steps.\n";
		cout << "*********************************************************\n";
	}

	// X*[invert (XTransX)]
	double* XinvXTX = new double[samSize * (numSelCol + 1)];
	matmatprod(covX, XTransX, XinvXTX, samSize, numSelCol + 1, numSelCol + 1);
	// residual = Y - X * beta
	double sigma2 = 0;
	for (int i = 0; i < samSize; i++) {
		residvec[i] = phenoY[i] - Xbeta[i];
		sigma2 += residvec[i] * residvec[i];
	}
	// sqr(sigma) = transpose(resid)*resid/[samSize-(numSelCol+1)]
	sigma2 = sigma2 / (samSize - (numSelCol + 1));
	if (phenoTyp == 1) sigma2 = 1.0;

	double* resid = &residvec[0];

	delete[] XTransY;
	delete[] beta;
	delete[] Xbeta;

	cout << "Streaming SNPs for speeding up GWAS analysis in parallel. \n";
	cout << "Number of SNPs in each batch is: " << stream_snps << '\n';
	cout << "*********************************************************\n";
	vector <double> ZGSvec(samSize * (1 + Sq) * stream_snps);
	vector <double> ZGSR2vec(samSize * (1 + Sq) * stream_snps);

	// for logistic regression
	vector <double> WZGSvec(samSize * (1 + Sq) * stream_snps);
	double* WZGS = &WZGSvec[0];


	int stream_snps1 = stream_snps;

	vector<long int> include_idx(samSize);
	int ii = 0;
	for (uint i = 0; i < bgen.Nbgen; i++) {
		if (genoUnMatchID.find(i) == genoUnMatchID.end()) {
			include_idx[ii] = i;
			ii++;
		}
	}
	genoUnMatchID.clear();




	bgen.numSelCol = numSelCol;
	bgen.Sq = Sq;
	bgen.samSize = samSize;
	bgen.robust = robust;
	bgen.stream_snps = stream_snps;
    bgen.include_idx = include_idx;
	bgen.maf = cmd.MAF;
	bgen.miu = miu;
	bgen.phenoTyp = phenoTyp;
	bgen.covX = covX;
	bgen.XinvXTX = XinvXTX;
	bgen.resid  = resid;
	bgen.sigma2 = sigma2;




	// Identifying the start position of each BGEN variant block for parallizing.
	cout << "Detected " << boost::thread::hardware_concurrency() << " available core(s). Using " << cmd.threads << " for parallel processing..." << endl;
	cout << "Dividing BGEN file into " << cmd.threads << " blocks..." << endl;

	vector<int> Mbgen_begin(cmd.threads);
	vector<int> Mbgen_end(cmd.threads);
	vector<long int> bgenVariantPos(cmd.threads);

	for (int t = 0; t < cmd.threads; t++) {
		Mbgen_begin[t] = floor((bgen.Mbgen / cmd.threads) * t);

		if ((t + 1) == (cmd.threads)) {
			Mbgen_end[t] = bgen.Mbgen - 1;
		}else {
			Mbgen_end[t] = floor(((bgen.Mbgen / cmd.threads) * (t + 1)) - 1);
		}
	}

	cout << "Identifying start position of each block... Execution time: ";
	auto start_time = std::chrono::high_resolution_clock::now();
	bgenVariantPos = getPositionOfBgenVariant(bgen.fin, bgen.offset, bgen.Mbgen, bgen.Nbgen, bgen.CompressedSNPBlocks, bgen.Layout, Mbgen_begin);
	auto end_time = std::chrono::high_resolution_clock::now();

	auto execution_time_ms  = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
	auto execution_time_sec = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
	auto execution_time_min = std::chrono::duration_cast<std::chrono::minutes>(end_time - start_time).count();

	if (execution_time_min > 0)
		cout << "" << execution_time_min % 60 << "m, ";
	if (execution_time_sec > 0)
		cout << "" << execution_time_sec % 60 << "s, ";
	if (execution_time_ms > 0)
		cout << "" << execution_time_ms % long(1E+3) << " ms";
	cout << "\n";
	fclose(bgen.fin);











	// Preparing for parallelizing of BGEN file
	boost::asio::io_service io_service;
	boost::thread_group thread_grp;
	boost::asio::io_service::work work(io_service);
	std::vector<boost::shared_future<void> > pending_data; // vector of futures

	// Create threads
	for (int i = 0; i < cmd.threads; ++i) {
		thread_grp.create_thread(boost::bind(&boost::asio::io_service::run, &io_service));
	}

	// Assigning each thread a chunk of the BGEN file.
	for (int i = 0; i < cmd.threads; i++) {
		typedef boost::packaged_task<void> task_t;
		boost::shared_ptr<task_t> task = boost::make_shared<task_t>(
			boost::bind(&BgenParallelGWAS, Mbgen_begin[i], Mbgen_end[i], bgenVariantPos[i], genofile, i, boost::ref(bgen)));

		boost::shared_future<void> fut = task->get_future();
		pending_data.push_back(fut);

		io_service.post(boost::bind(&task_t::operator(), task));
	}

	boost::wait_for_all(pending_data.begin(), pending_data.end());
	//pending_data.clear();











	// Write all results from each thread to 1 file
	std::ofstream results(output, std::ofstream::binary);
	results << "SNPID" << "\t" << "rsID" << "\t" << "CHR" << "\t" << "POS" << "\t" << "Allele1" << "\t" << "Allele2" << "\t" << "AF" << "\t" << "Beta_Main" << "\t" << "Var_Beta_Main" << "\t";
	for (int i = 1; i <= Sq; i++) {
		results << "Beta_Interaction" << "_" << i << "\t";

		for (int i = 1; i <= Sq; i++) {
			for (int j = 1; j <= Sq; j++) {
				results << "Var_Beta_Interaction" << "_" << i << "_" << j << "\t";
			}
		}
		results << "P_Value_Main" << "\t" << "P_Value_Interaction" << "\t" << "P_Value_Joint\n";
	}

	for (int i = 0; i < cmd.threads; i++) {
		std::string threadOutputFile = "gem_bin_" + std::to_string(i) + ".tmp";
		std::ifstream thread_output(threadOutputFile);
		results << thread_output.rdbuf();

		thread_output.close();
		boost::filesystem::remove(threadOutputFile.c_str());

	}

	results.close();



	phenodata.clear();
	covdata.clear();
	delete[] XTransX;
	delete[] XinvXTX;
	delete[] samID;
	include_idx.clear();





	std::chrono::duration<double> wallduration = (std::chrono::system_clock::now() - wall0);
	double cpuduration = (std::clock() - cpu0) / (double)CLOCKS_PER_SEC;
	cout << "Total Wall Time = " << wallduration.count() << "  Seconds\n";
	cout << "Total CPU Time = "  << cpuduration << "  Seconds\n";
	//cout << "Execution Wall Time = " << exetime << "  Seconds\n";
	cout << "*********************************************************\n";

	return 0;
}


#include "declars.h"
#include "ReadBed.h"




void Bed::processBed(string bedFile, string bimFile, string famFile) {

	cout << "General information of BED file. \n";
	FILE* fin = fopen(bedFile.c_str(), "rb");
	if (fin == 0) {
		cerr << "\nERROR: Bed file could not be opened.\n\n";
		exit(1);
	}

	unsigned char magic[3];
	if (!fread(magic, 1, 3, fin)) {
		cerr << "\nERROR: Cannot read bed file magic byte numbers.\n\n"; 
		exit(1);
	}
	if (magic[0] != 0x6C || magic[1] != 0x1B) { 
		cerr << "ERROR: " << bedFile << " is not a bed file (incorrect 'magic numbers').\n\n";
		exit(1);
	}
	if (magic[2] != 0x01) {
		cerr << "ERROR: " << bedFile << " is not in SNP major format.\n\n";
		exit(1);
	}

	// Fam file
	std::ifstream readFam;
	string IDline;
	readFam.open(famFile);
	if (!readFam) { 
		cerr << "\nERROR: Cannot open .fam file " << famFile << ".\n\n";
		exit(1);
	}

	getline(readFam, IDline);
	std::istringstream iss(IDline);
	string value;
	vector <string> values;
	char tmpDelim = ' ';
	while (getline(iss, value, tmpDelim)) {
		value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
		values.push_back(value);
	}

	if (values.size() <= 1) {
		values.clear();
		readFam.clear();
		readFam.seekg(0, readFam.beg);
		getline(readFam, IDline);
		std::istringstream iss(IDline);
		while (getline(iss, value, '\t')) {
			value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
			values.push_back(value);
		}
		if (values.size() > 1) {
			tmpDelim = '\t';
		}
		else {
			cerr << "\nERROR: .fam file should be space or tab separated.\n\n";
		}
	}
	famDelim = tmpDelim;
	uint nsamples = 1;
	while (getline(readFam, IDline)) {
		   nsamples++;
	}
	readFam.close();


	// Bim File
	std::ifstream readBim;
	string var;
	readBim.open(bimFile);
	if (!readBim) {
		cerr << "\nERROR: Cannot open .bim file " << bimFile << ".\n\n";
		exit(1);
	}

	getline(readBim, IDline);
	std::istringstream iss2(IDline);
	values.clear();
	tmpDelim = ' ';
	while (getline(iss2, value, tmpDelim)) {
		value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
		values.push_back(value);
	}

	if ((values.size() != 5) && (values.size() != 6)) {
		values.clear();
		readBim.clear();
		readBim.seekg(0, readBim.beg);
		getline(readBim, IDline);
		std::istringstream iss2(IDline);
		while (getline(iss2, value, '\t')) {
			value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
			values.push_back(value);
		}
		if ((values.size() == 5) || (values.size() == 6)) {
			tmpDelim = '\t';
		}
		else {
			cerr << "\nERROR: .fam file should be space or tab separated.\n\n";
		}
	}
	bimDelim = tmpDelim;
	bimLast = values.size() - 1;
	uint nvars = 1;
	while (getline(readBim, var)) {
		nvars++;
	}
	readBim.close();

	n_samples = nsamples;
	n_variants = nvars;

	cout << "Number of variants: " << nvars << endl;
	cout << "Number of samples: " << nsamples << endl;

	fclose(fin);
}



void Bed::processFam(Bed bed, string famFile, unordered_map<string, vector<string>> phenomap, string phenoMissingKey, vector<double> phenodata, vector<double> covdata, int numSelCol, int samSize, double center, double scale) {


	unordered_set<int> genoUnMatchID;

	std::ifstream fIDMat;
	string IDline, FID, strtmp;
	fIDMat.open(famFile);
	if (!fIDMat.is_open()) {
		cerr << "\nERROR: .fam file could not be opened.\n\n";
		exit(1);
	}

	int k = 0;
	for (uint m = 0; m < bed.n_samples; m++) {
		getline(fIDMat, IDline);
		std::istringstream iss(IDline);

		string value;
		vector <string> values;
		while (getline(iss, value, bed.famDelim)) {
			value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
			values.push_back(value);
		}
		string strtmp = values[1];
		int itmp = k;
		if (phenomap.find(strtmp) != phenomap.end()) {
			auto tmp_valvec = phenomap[strtmp];
			if (find(tmp_valvec.begin(), tmp_valvec.end(), phenoMissingKey) == tmp_valvec.end()) {
				sscanf(tmp_valvec[0].c_str(), "%lf", &phenodata[k]);
				for (int c = 0; c < numSelCol; c++) {
					sscanf(tmp_valvec[c + 1].c_str(), "%lf", &covdata[k * (numSelCol + 1) + c + 1]);
				}
				k++;
			}
			phenomap.erase(strtmp);
		}

		if (itmp == k) genoUnMatchID.insert(m);
	}
	fIDMat.close();

	//After IDMatching, resizing phenodata and covdata, and updating samSize;
	phenodata.resize(k);
	covdata.resize(k * (numSelCol + 1));
	samSize = k;

	if (samSize == 0) {
		cout << "\nERROR: sample size changed from " << samSize + genoUnMatchID.size() << " to " << samSize << ".\n";
		cout << "       check if sample ids are consistent between the phenotype file and .fam file, or check\n";
		cout << "       if (--sampleid-name) is specified correctly. \n\n";
		exit(1);
	}


	int ii = 0;
	include_idx.resize(samSize);
	for (uint i = 0; i < bed.n_samples; i++) {
		if (genoUnMatchID.find(i) == genoUnMatchID.end()) {
			include_idx[ii] = i;
			ii++;

		}
	}

	cout << "****************************************************************************\n";
	if (genoUnMatchID.empty()) {
		cout << "After processes of sample IDMatching and checking missing values, the sample size does not change.\n\n";
	}
	else {
		cout << "After processes of sample IDMatching and checking missing values, the sample size changes from "
			<< samSize + genoUnMatchID.size() << " to " << samSize << ".\n\n";
	}
	cout << "Sample IDMatching and checking missing values processes have been completed.\n";
	cout << "New pheno and covariate data vectors with the same order of sample ID sequence of geno data are updated.\n";
	cout << "****************************************************************************\n";


	new_samSize = samSize;
	new_covdata = covdata;
	new_phenodata = phenodata;

}



void Bed::getBedVariantPos(Bed bed, CommandLine cmd) {

	uint32_t nSNPS = bed.n_variants;
	int count = 0;
	std::set<std::string> includeVariant;
	threads = cmd.threads;


	if (cmd.doFilters) {

		filterVariants = true;
		if (!cmd.includeVariantFile.empty()) {
			std::ifstream fInclude;
			string IDline;
			fInclude.open(cmd.includeVariantFile);
			if (!fInclude.is_open()) {
				cerr << "\nERROR: The file (" << cmd.includeVariantFile << ") could not be opened." << endl << endl;
				exit(1);
			}

			string vars;
			getline(fInclude, IDline);
			std::transform(IDline.begin(), IDline.end(), IDline.begin(), ::tolower);
			IDline.erase(std::remove(IDline.begin(), IDline.end(), '\r'), IDline.end());

			if (IDline == "snpid") {
				cout << "An include snp file was detected... \nIncluding SNPs for analysis based on their snpid... \n";
			}
			else {
				cerr << "\nERROR: Header name of " << cmd.includeVariantFile << " must be 'snpid' for PGEN files." << endl << endl;
				exit(1);
			}

			while (fInclude >> vars) {
				if (includeVariant.find(vars) != includeVariant.end()) {
					cout << "\nERROR: " << vars << " is a duplicate variant in " << cmd.includeVariantFile << ".\n\n";
					exit(1);
				}
				includeVariant.insert(vars);
				count++;
			}
			nSNPS = count;
			cout << "Detected " << nSNPS << " variants to be used for analysis... \nAll other variants will be excluded.\n\n" << endl;
			//count = ceil(count / Mbgen_begin.size());


			cout << "Detected " << boost::thread::hardware_concurrency() << " available thread(s)...\n";
			if (nSNPS < threads) {
				cout << "Number of variants (" << nSNPS << ") is less than the number of specified threads (" << threads << ")...\n";
				threads = nSNPS;
				cout << "Using " << threads << " for multithreading... \n\n";
			}
			else {
				cout << "Using " << threads << " for multithreading... \n\n";
			}

		}

		begin.resize(threads);
		end.resize(threads);

		for (uint t = 0; t < threads; t++) {
			begin[t] = floor((nSNPS / threads) * t);
			end[t] = ((t + 1) == threads) ? nSNPS - 1 : floor(((nSNPS / threads) * (t + 1)) - 1);
		}

		std::ifstream fIDMat;
		fIDMat.open(cmd.bimFile);

		string IDline;

		uint32_t k = 0;
		long long unsigned int pvalIndex = 0;
		while (getline(fIDMat, IDline)) {
			std::istringstream iss(IDline);
			string value;
			vector <string> values;

			while (getline(iss, value, bed.bimDelim)) {
				value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
				values.push_back(value);
			}
			if (includeVariant.find(values[1]) != includeVariant.end()) {
				bedVariantPos.push_back(pvalIndex);
				k++;
			}
			pvalIndex++;
		}

		if (k != nSNPS) {
			cerr << "\nERROR: There are one or more SNPs in second column of .bim file not in " << cmd.includeVariantFile << ".\n\n";
			exit(1);
		}
	}
	else {
		threads = cmd.threads;
		if (bed.n_variants < threads) {
			cout << "Number of variants (" << bed.n_variants << ") is less than the number of specified threads (" << threads << ")...\n";
			threads = bed.n_variants;
			cout << "Using " << threads << " thread(s) instead... \n\n";
		}

		begin.resize(threads);
		end.resize(threads);

		for (uint32_t t = 0; t < threads; t++) {
			begin[t] = floor((bed.n_variants / threads) * t);
			if ((t + 1) == (threads)) {
				end[t] = bed.n_variants - 1;
			}
			else {
				end[t] = floor(((bed.n_variants / threads) * (t + 1)) - 1);
			}
		}
	}
}


void gemBED(uint32_t begin, uint32_t end, string bedFile, string bimFile, int thread_num, bool filterVariants, std::vector<long long unsigned int> bedPos,
			int Sq, int numSelCol, int numIntSelCol, int numExpSelCol, int phenoType, int robust, int samSize, int stream_snps, double MAF, double missGenoCutoff,
			char bimDelim, int bimLast, uint n_samples, double sigma2, double* resid, double* XinvXTX, double* covX, vector<double> miu,
			vector<long int> include_idx, string outFile) {

	auto start_time = std::chrono::high_resolution_clock::now();
	std::string output = outFile + "_bin_" + std::to_string(thread_num) + ".tmp";
	std::ofstream results(output, std::ofstream::binary);
	std::ostringstream oss;


	int Sq1 = Sq + 1;
	int intSq = numIntSelCol;
	int intSq1 = intSq + 1;
	int expSq = numExpSelCol;
	double maxMAF = 1.0 - MAF;
	vector <string> geno_snpid(stream_snps);
	vector <double> ZGSvec(samSize * (1 + Sq) * stream_snps);
	vector <double> ZGSR2vec(samSize * (1 + Sq) * stream_snps);
	vector <double> WZGSvec(samSize * (1 + Sq) * stream_snps);
	double* WZGS = &WZGSvec[0];
	std::vector<uint> missingIndex;
	vector <double> AF(stream_snps);
	int ZGS_col = Sq1 * stream_snps;

	std::ifstream fIDMat;
	fIDMat.open(bimFile);
	string IDline;
	uint32_t skipIndex = 0;
	if (!filterVariants) {
		while (skipIndex != begin) {
			getline(fIDMat, IDline);
			skipIndex++;
		}
	}
	else {
		while (skipIndex != bedPos[begin]) {
			getline(fIDMat, IDline);
			skipIndex++;
		}
	}


	std::ifstream readbedfile(bedFile.c_str(), std::ios::binary);
	uint nblocks = (n_samples + 3) / 4, pos;
	unsigned char temp[2];
	unsigned char* buffer = new unsigned char[nblocks];


	uint32_t snploop = begin;
	int variant_index = 0;
	int keepIndex = 0;
	double geno;
	while (snploop <= end) {

		int stream_i = 0;
		while (stream_i < stream_snps) {

			if ((snploop == (end + 1)) && stream_i == 0) {
				break;
			}

			if (snploop == end + 1 && stream_i != 0) {
				stream_snps = stream_i;
				ZGS_col = Sq1 * stream_snps;
				break;
			}

			string value;
			vector <string> values;
			if (!filterVariants) {
				readbedfile.seekg((std::streamoff)snploop * nblocks + 3, readbedfile.beg);
				getline(fIDMat, IDline);
				std::istringstream iss(IDline);
				while (getline(iss, value, bimDelim)) {
					values.push_back(value);
				}
			}
			else {
				while (skipIndex != bedPos[snploop]) {
					getline(fIDMat, IDline);
					skipIndex++;
				}
				readbedfile.seekg((std::streamoff)bedPos[snploop] * nblocks + 3, readbedfile.beg);
				getline(fIDMat, IDline);
				skipIndex++;
				std::istringstream iss(IDline);
				while (getline(iss, value, bimDelim)) {
					values.push_back(value);
				}
			}
			snploop++;


			int tmp1 = stream_i * Sq1 * samSize;
			int idx_k = 0;
			int nMissing = 0;
			uint ncount = 0;
			readbedfile.read((char*)buffer, nblocks);
			for (size_t n = 0; n < nblocks; n++) {
				pos = 0;
				for (int i = 0; i < 4; i++) {
					if ((ncount == n_samples) && (n == nblocks - 1)) {
						break;
					}

					for (size_t l = 0; l < 2; ++l) {
						temp[l] = (buffer[n] >> pos) & 1;
						pos++;
					}

					if (include_idx[idx_k] != ncount) {
						ncount++;
						continue;
					}


					if (temp[0] == 0 && temp[1] == 0) {
						geno = 2.0;
					}
					else if (temp[0] == 1 && temp[1] == 1) {
						geno = 0.0;
					}
					else if (temp[0] == 0 && temp[1] == 1) {
						geno = 1.0;
					}
					else {
						missingIndex.push_back(idx_k);
						nMissing++;
						idx_k++;
						ncount++;
						continue;
					}

					ncount++;
					int tmp2 = idx_k + tmp1;
					AF[stream_i] += geno;

					if (phenoType == 1) {
						ZGSvec[tmp2] = miu[idx_k] * (1 - miu[idx_k]) * geno;
					}
					else {
						ZGSvec[tmp2] = geno;
					}
					idx_k++;
				}
			}

	
			double gmean = AF[stream_i] / double(samSize - nMissing);
			double cur_AF = AF[stream_i] / double(samSize - nMissing) / 2.0;
			double percMissing = nMissing / (samSize * 1.0);

			if ((cur_AF < MAF || cur_AF > maxMAF) || (percMissing > missGenoCutoff)) {
				AF[stream_i] = 0;
				continue;
			}
			else {
				AF[stream_i] = cur_AF;
			}

			if (nMissing > 0) {
				if (phenoType == 0) {
					for (long unsigned int nm = 0; nm < missingIndex.size(); nm++) {
						int tmp5 = tmp1 + missingIndex[nm];
						ZGSvec[tmp5] = gmean;
					}
				}
				else {
					for (long unsigned int nm = 0; nm < missingIndex.size(); nm++) {
						int tmp5 = tmp1 + missingIndex[nm];
						ZGSvec[tmp5] = miu[missingIndex[nm]] * (1 - miu[missingIndex[nm]]) * gmean;
					}
				}
				missingIndex.clear();
			}


			values[bimLast].erase(std::remove(values[bimLast].begin(), values[bimLast].end(), '\r'), values[bimLast].end());
			geno_snpid[stream_i] = values[1] + "\t" + values[0] + "\t" + values[bimLast - 2] + "\t" + values[bimLast] + "\t" + values[bimLast - 1] + "\t" + std::to_string(samSize - nMissing);

			for (int j = 0; j < Sq; j++) {
				int tmp3 = samSize * (j + 1) + tmp1;
				for (int i = 0; i < samSize; i++) {
					int tmp4 = i * (numSelCol + 1);
					ZGSvec[tmp3 + i] = covX[tmp4 + j + 1] * ZGSvec[tmp1 + i]; // here we save ZGS in column wise
				}
			}
			variant_index++;
			stream_i++;
			keepIndex++;
		}

		if ((snploop == (end + 1)) & (stream_i == 0)) {
			break;
		}

		/***************************************************************/
		//	genodata and envirment data
		double* ZGS = &ZGSvec[0];
		double* ZGSR2 = &ZGSR2vec[0];
		// transpose(X) * ZGS
		// it is non-squred matrix, attention that continuous memory is column-major due to Fortran in BLAS.
		// important!!!!
		double* XtransZGS = new double[(numSelCol + 1) * ZGS_col];
		matNmatNprod(covX, ZGS, XtransZGS, numSelCol + 1, samSize, ZGS_col);

		if (phenoType == 0) {
			matNmatNprod(XinvXTX, XtransZGS, ZGSR2, samSize, (numSelCol + 1), ZGS_col);
			matAdd(ZGS, ZGSR2, samSize * ZGS_col, -1);
			if (robust == 1) {
				for (int j = 0; j < ZGS_col; j++) {
					for (int i = 0; i < samSize; i++) {
						ZGSR2vec[j * samSize + i] = ZGS[j * samSize + i] * resid[i] * resid[i];
					}
				}
			}
		}
		else if (phenoType == 1) {
			WZGSvec = ZGSvec;
			matNmatNprod(XinvXTX, XtransZGS, ZGSR2, samSize, (numSelCol + 1), ZGS_col);
			matAdd(WZGS, ZGSR2, samSize * ZGS_col, -1);

			for (int j = 0; j < ZGS_col; j++) {
				for (int i = 0; i < samSize; i++) {
					ZGS[j * samSize + i] = WZGS[j * samSize + i] / miu[i] / (1.0 - miu[i]);
					if (robust == 1) ZGSR2vec[j * samSize + i] = ZGS[j * samSize + i] * resid[i] * resid[i];
				}
			}
		}
		delete[] XtransZGS;


		// transpose(ZGS) * resid
		double* ZGStR = new double[ZGS_col];
		matvecprod(ZGS, resid, ZGStR, ZGS_col, samSize);
		// transpose(ZGS) * ZGS
		double* ZGStZGS = new double[ZGS_col * ZGS_col];
		if (phenoType == 0) {
			matmatTprod(ZGS, ZGS, ZGStZGS, ZGS_col, samSize, ZGS_col);
		}
		if (phenoType == 1) {
			matmatTprod(ZGS, WZGS, ZGStZGS, ZGS_col, samSize, ZGS_col);
		}

		// transpose(ZGSR2) * ZGS
		double* ZGSR2tZGS = new double[ZGS_col * ZGS_col];
		if (robust == 1) matmatTprod(ZGSR2, ZGS, ZGSR2tZGS, ZGS_col, samSize, ZGS_col);


		double* betaM = new double[stream_snps];
		double* VarbetaM = new double[stream_snps];
		double** betaInt = new double* [stream_snps];
		double** VarbetaInt = new double* [stream_snps];
		double* PvalM = new double[stream_snps];
		double* PvalInt = new double[stream_snps];
		double* PvalJoint = new double[stream_snps];
		boost::math::chi_squared chisq_dist_M(1);
		boost::math::chi_squared chisq_dist_Int(expSq);
		boost::math::chi_squared chisq_dist_Joint(1 + expSq);


		if (robust == 0) {
			for (int i = 0; i < stream_snps; i++) {
				// initialize dynamic 2D array
				betaInt[i] = new double[expSq];
				VarbetaInt[i] = new double[expSq * expSq];

				// betamain
				int tmp1 = i * ZGS_col * Sq1 + i * Sq1;
				betaM[i] = ZGStR[i * Sq1] / ZGStZGS[tmp1];
				VarbetaM[i] = sigma2 / ZGStZGS[tmp1];

				// ZGStR
				double* subZGStR = new double[Sq1];
				subMatrix(ZGStR, subZGStR, Sq1, 1, Sq1, Sq1, i * Sq1);

				// inv(ZGStZGS[)
				double* invZGStZGS = new double[Sq1 * Sq1];
				subMatrix(ZGStZGS, invZGStZGS, Sq1, Sq1, ZGS_col, Sq1, tmp1);
				matInv(invZGStZGS, Sq1);


				double* betaAll = new double[Sq1 * 1];
				matvecprod(invZGStZGS, subZGStR, betaAll, Sq1, Sq1);
				double* VarBetaAll = new double[Sq1 * Sq1];
				subMatrix(invZGStZGS, VarBetaAll, Sq1, Sq1, Sq1, Sq1, 0);
				for (int k = 0; k < Sq1 * Sq1; k++) {
					VarBetaAll[k] *= sigma2;
				}

				//calculating Marginal P values
				double statM = betaM[i] * betaM[i] / VarbetaM[i];
				if (isnan(statM) || statM <= 0.0) {
					PvalM[i] = NAN;
				}
				else {
					PvalM[i] = boost::math::cdf(complement(chisq_dist_M, statM));
				}


				// VarBetaInt
				subMatrix(VarBetaAll, VarbetaInt[i], expSq, expSq, Sq1, expSq, (intSq1 * Sq1 + intSq1));

				double* invVarbetaint = new double[expSq * expSq];
				subMatrix(VarBetaAll, invVarbetaint, expSq, expSq, Sq1, expSq, (intSq1 * Sq1 + intSq1));
				matInv(invVarbetaint, expSq);

				// Beta Int
				subMatrix(betaAll, betaInt[i], expSq, 1, expSq, 1, intSq1);

				// StatInt
				double* Stemp3 = new double[expSq];
				matvecprod(invVarbetaint, betaInt[i], Stemp3, expSq, expSq);
				double statInt = 0.0;
				for (int j = 0; j < expSq; j++) {
					statInt += betaInt[i][j] * Stemp3[j];
				}

				//calculating Interaction P values
				if (isnan(statInt) || statInt <= 0.0) {
					PvalInt[i] = NAN;
				}
				else {
					PvalInt[i] = boost::math::cdf(complement(chisq_dist_Int, statInt));
				}

				vector<double> invAvec((1 + expSq) * (1 + expSq));
				vector<double> betaMIntvec(1 + expSq);
				int betaIndex = 0;
				int tmpIndex = 0;
				for (int k = 0; k < Sq1; k++) {
					if (k > 0 && k <= intSq) { continue; }
					else {
						betaMIntvec[betaIndex] = betaAll[k];
						betaIndex++;
						for (int j = 0; j < Sq1; j++) {
							if (j > 0 && j <= intSq) { continue; }
							else {
								invAvec[tmpIndex] = VarBetaAll[(k * Sq1) + j];
								tmpIndex++;
							}
						}
					}
				}
				double* invA = &invAvec[0];
				double* betaMInt = &betaMIntvec[0];
				matInv(invA, 1 + expSq);
				double* Stemp4 = new double[1 + expSq];
				matvecprod(invA, betaMInt, Stemp4, 1 + expSq, 1 + expSq);
				double statJoint = 0.0;
				for (int k = 0; k < 1 + expSq; k++) {
					statJoint += betaMInt[k] * Stemp4[k];
				}
				if (isnan(statJoint) || statJoint <= 0.0) {
					PvalJoint[i] = NAN;
				}
				else {
					PvalJoint[i] = boost::math::cdf(complement(chisq_dist_Joint, statJoint));
				}



				delete[] subZGStR;
				delete[] invZGStZGS;
				delete[] betaAll;
				delete[] VarBetaAll;
				delete[] invVarbetaint;
				delete[] Stemp3;
				delete[] Stemp4;
			}
		}
		else if (robust == 1) {
			for (int i = 0; i < stream_snps; i++) {

				betaInt[i] = new double[expSq];
				VarbetaInt[i] = new double[expSq * expSq];

				//BetaMain
				int tmp1 = i * ZGS_col * Sq1 + i * Sq1;
				betaM[i] = ZGStR[i * Sq1] / ZGStZGS[tmp1];
				VarbetaM[i] = ZGSR2tZGS[tmp1] / (ZGStZGS[tmp1] * ZGStZGS[tmp1]);

				// ZGStR
				double* subZGStR = new double[Sq1];
				subMatrix(ZGStR, subZGStR, Sq1, 1, Sq1, Sq1, i * Sq1);

				// ZGSR2tZGS
				double* subZGSR2tZGS = new double[Sq1 * Sq1];
				subMatrix(ZGSR2tZGS, subZGSR2tZGS, Sq1, Sq1, ZGS_col, Sq1, tmp1);

				// inv(ZGStZGS[)
				double* invZGStZGS = new double[Sq1 * Sq1];
				subMatrix(ZGStZGS, invZGStZGS, Sq1, Sq1, ZGS_col, Sq1, tmp1);
				matInv(invZGStZGS, Sq1);


				double* betaAll = new double[Sq1 * 1];
				matvecprod(invZGStZGS, subZGStR, betaAll, Sq1, Sq1);
				double* ZGSR2tZGSxinvZGStZGS = new double[Sq1 * Sq1];
				matNmatNprod(subZGSR2tZGS, invZGStZGS, ZGSR2tZGSxinvZGStZGS, Sq1, Sq1, Sq1);
				double* VarBetaAll = new double[Sq1 * Sq1];
				matNmatNprod(invZGStZGS, ZGSR2tZGSxinvZGStZGS, VarBetaAll, Sq1, Sq1, Sq1);


				//calculating Marginal P values
				double statM = betaM[i] * betaM[i] / VarbetaM[i];
				if (isnan(statM) || statM <= 0.0) {
					PvalM[i] = NAN;
				}
				else {
					PvalM[i] = boost::math::cdf(complement(chisq_dist_M, statM));
				}


				// VarBetaInt
				subMatrix(VarBetaAll, VarbetaInt[i], expSq, expSq, Sq1, expSq, (intSq1 * Sq1 + intSq1));

				double* invVarbetaint = new double[expSq * expSq];
				subMatrix(VarBetaAll, invVarbetaint, expSq, expSq, Sq1, expSq, (intSq1 * Sq1 + intSq1));
				matInv(invVarbetaint, expSq);

				// Beta Int
				subMatrix(betaAll, betaInt[i], expSq, 1, expSq, 1, intSq1);

				// StatInt
				double* Stemp3 = new double[expSq];
				matvecprod(invVarbetaint, betaInt[i], Stemp3, expSq, expSq);
				double statInt = 0.0;
				for (int j = 0; j < expSq; j++) {
					statInt += betaInt[i][j] * Stemp3[j];
				}

				//calculating Interaction P values
				if (isnan(statInt) || statInt <= 0.0) {
					PvalInt[i] = NAN;
				}
				else {
					PvalInt[i] = boost::math::cdf(complement(chisq_dist_Int, statInt));
				}

				vector<double> invAvec((1 + expSq) * (1 + expSq));
				vector<double> betaMIntvec(1 + expSq);
				int betaIndex = 0;
				int tmpIndex = 0;
				for (int k = 0; k < Sq1; k++) {
					if (k > 0 && k <= intSq) { continue; }
					else {
						betaMIntvec[betaIndex] = betaAll[k];
						betaIndex++;
						for (int j = 0; j < Sq1; j++) {
							if (j > 0 && j <= intSq) { continue; }
							else {
								invAvec[tmpIndex] = VarBetaAll[(k * Sq1) + j];
								tmpIndex++;
							}
						}
					}
				}
				double* invA = &invAvec[0];
				double* betaMInt = &betaMIntvec[0];
				matInv(invA, 1 + expSq);
				double* Stemp4 = new double[1 + expSq];
				matvecprod(invA, betaMInt, Stemp4, 1 + expSq, 1 + expSq);
				double statJoint = 0.0;
				for (int k = 0; k < 1 + expSq; k++) {
					statJoint += betaMInt[k] * Stemp4[k];
				}
				if (isnan(statJoint) || statJoint <= 0.0) {
					PvalJoint[i] = NAN;
				}
				else {
					PvalJoint[i] = boost::math::cdf(complement(chisq_dist_Joint, statJoint));
				}



				delete[] subZGStR;
				delete[] subZGSR2tZGS;
				delete[] invZGStZGS;
				delete[] betaAll;
				delete[] ZGSR2tZGSxinvZGStZGS;
				delete[] VarBetaAll;
				delete[] invVarbetaint;
				delete[] Stemp3;
				delete[] Stemp4;

		}
		} // end of if robust == 1


		for (int i = 0; i < stream_snps; i++) {
			oss << geno_snpid[i] << "\t" << AF[i] << "\t" << betaM[i] << "\t" << VarbetaM[i] << "\t";
			for (int ii = 0; ii < expSq; ii++) {
				oss << betaInt[i][ii] << "\t";
			}

			for (int ii = 0; ii < expSq; ii++) {
				for (int jj = 0; jj < expSq; jj++) {
					oss << VarbetaInt[i][ii * expSq + jj] << "\t";
				}
			}
			oss << PvalM[i] << "\t" << PvalInt[i] << "\t" << PvalJoint[i] << '\n';
			AF[i] = 0.0;
		}

		delete[] ZGStR;
		delete[] ZGStZGS;
		delete[] ZGSR2tZGS;

		delete[] betaM;
		delete[] VarbetaM;
		for (int i = 0; i < stream_snps; i++) {
			delete[] betaInt[i];
			delete[] VarbetaInt[i];
		}
		delete[] betaInt;
		delete[] VarbetaInt;
		delete[] PvalM;
		delete[] PvalInt;
		delete[] PvalJoint;

		if (variant_index % 10000 == 0) {
			results << oss.str();
			oss.str(std::string());
			oss.clear();
		}
	}

	if (variant_index % 10000 != 0) {
		results << oss.str();
		oss.str(std::string());
		oss.clear();
	}

	// Close files
	results.close();
	fIDMat.close();


	auto end_time = std::chrono::high_resolution_clock::now();
	cout << "Thread " << thread_num << " finished in ";
	printExecutionTime1(start_time, end_time);
}



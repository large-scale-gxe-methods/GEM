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


	std::ifstream readFam;
	string IDline;
	readFam.open(famFile);
	if (!readFam) { 
		cerr << "\nERROR: Cannot open .fam file " << famFile << ".\n\n";
		exit(1);
	}
	uint nsamples = 0;
	while (getline(readFam, IDline)) {
		   nsamples++;
	}
	readFam.close();


	std::ifstream readBim;
	string var;
	readBim.open(bimFile);
	if (!readBim) {
		cerr << "\nERROR: Cannot open .bim file " << bimFile << ".\n\n";
		exit(1);
	}
	uint nvars = 0;
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



void Bed::processFam(Bed bed, string famFile, unordered_map<string, vector<string>> phenomap, string phenoMissingKey, vector<double> phenodata, vector<double> covdata, int numSelCol, int samSize) {

	unordered_set<int> genoUnMatchID;

	std::ifstream fIDMat;
	string IDline, FID, strtmp;

	fIDMat.open(famFile);
	if (!fIDMat.is_open()) {
		cerr << "\nERROR: .fam file could not be opened.\n\n";
		exit(1);
	}


	string iid;
	vector <string> values;

	int k = 0;
	for (uint m = 0; m < bed.n_samples; m++) {
		getline(fIDMat, IDline);
		std::istringstream iss(IDline);
		iss >> FID >> strtmp;

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

	new_samSize = samSize;
	new_covdata = covdata;
	new_phenodata = phenodata;


	if (samSize == 0) {
		cout << "\nerror: sample size changed from " << samSize + genoUnMatchID.size() << " to " << samSize << ".\n";
		cout << "\ncheck if sample ids are consistent between the phenotype file and .psam file, or check if (--sampleid-name) is specified correctly. \n\n";
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
}




void gemBED(uint32_t begin, uint32_t end, string bedFile, string bimFile, int thread_num, Bed test) {

	auto start_time = std::chrono::high_resolution_clock::now();
	std::string output = test.outFile + "_bin_" + std::to_string(thread_num) + ".tmp";
	std::ofstream results(output, std::ofstream::binary);
	std::ostringstream oss;


	int Sq = test.Sq;
	int Sq1 = Sq + 1;
	int intSq = test.numIntSelCol;
	int intSq1 = intSq + 1;
	int expSq = test.numExpSelCol;
	int phenoType = test.phenoTyp;
	int robust = test.robust;
	int samSize = test.new_samSize;
	int numSelCol = test.numSelCol;
	int stream_snps = test.stream_snps;
	double MAF = test.maf;
	double maxMAF = 1.0 - MAF;
	double missGenoCutoff = test.missGeno;
	double sigma2 = test.sigma2;
	double* resid = test.resid;
	double* XinvXTX = test.XinvXTX;
	double* covX = test.covX;
	vector <double> miu = test.miu;
	vector<long int> include_idx = test.include_idx;
	vector <string> geno_snpid(stream_snps);
	vector <double> ZGSvec(samSize * (1 + Sq) * stream_snps);
	vector <double> ZGSR2vec(samSize * (1 + Sq) * stream_snps);
	vector <double> WZGSvec(samSize * (1 + Sq) * stream_snps);
	double* WZGS = &WZGSvec[0];
	std::vector<uint> missingIndex;
	vector <double> AF(stream_snps);
	int ZGS_col = Sq1 * stream_snps;
	uint n_samples = test.n_samples;
	uint n_variants = test.n_variants;

	std::ifstream fIDMat;
	fIDMat.open(bimFile);

	string IDline;
	uint32_t skipIndex = 0;
	while (skipIndex != begin) {
		getline(fIDMat, IDline);
		skipIndex++;
	}


	std::ifstream readbedfile(bedFile.c_str(), std::ios::binary);
	uint nblocks = (n_samples + 3) / 4, pos;
	unsigned char temp[2];
	unsigned char* buffer = new unsigned char[nblocks];


	uint32_t snploop = begin;
	int variant_index = 0;
	int keepIndex = 0;
	double geno;
	while (snploop <= n_variants) {

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
			readbedfile.seekg(snploop * nblocks + 3);
			getline(fIDMat, IDline);
			std::istringstream iss(IDline);
			string value;
			vector <string> values;
			while (getline(iss, value, '\t')) {
				values.push_back(value);
			}
			values[5].erase(std::remove(values[5].begin(), values[5].end(), '\r'), values[5].end());
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
						geno = 0.0;
					}
					else if (temp[0] == 1 && temp[1] == 1) {
						geno = 2.0;
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
			geno_snpid[stream_i] = values[1] + "\t" + values[0] + "\t" + values[3] + "\t" + values[4] + "\t" + values[5] + "\t" + std::to_string(samSize - nMissing);

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

				//double* S2TransS2     = new double[Sq * Sq];
				//double* S2TransR      = new double[Sq];
				double* InvVarbetaint = new double[expSq * expSq];

				// inv(ZGStZGS[tmp1])
				double* invZGStZGStmp1 = new double[intSq1 * intSq1];
				subMatrix(ZGStZGS, invZGStZGStmp1, intSq1, intSq1, ZGS_col, intSq1, tmp1);
				matInv(invZGStZGStmp1, intSq1);
				// ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col]
				double* ZGStZGSsGRow = new double[expSq * intSq1];
				subMatrix(ZGStZGS, ZGStZGSsGRow, expSq, intSq1, ZGS_col, expSq, tmp1 + intSq1);

				/* For S2TransR*/
				// ZGStR[i * Sq1 + ind1 + 1]
				double* expZGStR = new double[expSq];
				subMatrix(ZGStR, expZGStR, expSq, 1, expSq, expSq, (i * Sq1) + intSq1);
				// ZGStR[i * Sq1]
				double* int1ZGStR = new double[intSq1];
				subMatrix(ZGStR, int1ZGStR, intSq1, 1, expSq + intSq, intSq1, (i * Sq1));
				// ZGStR[i * Sq1] / ZGStZGS[tmp1]
				double* invZGStmpInt1ZGStR = new double[intSq1];
				matvecprod(invZGStZGStmp1, int1ZGStR, invZGStmpInt1ZGStR, intSq1, intSq1);
				// ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col] * ZGStZGS[tmp1 + ind2 + 1] / ZGStZGS[tmp1];
				double* S2TransRright = new double[expSq];
				matNmatNprod(ZGStZGSsGRow, invZGStmpInt1ZGStR, S2TransRright, expSq, intSq1, 1);
				matAdd(expZGStR, S2TransRright, expSq, -1.0);


				/* For S2TransS2S2*/
				double* expS2TransS2 = new double[expSq * expSq];
				subMatrix(ZGStZGS, expS2TransS2, expSq, expSq, ZGS_col, expSq, tmp1 + (intSq1 * ZGS_col + intSq1));
				// ZGStZGS[tmp1 + ind2 + 1] / ZGStZGS[tmp1]
				double* invZGStmp1IntExpZGStZGS = new double[intSq1 * expSq];
				matTmatprod(invZGStZGStmp1, ZGStZGSsGRow, invZGStmp1IntExpZGStZGS, intSq1, intSq1, expSq);
				// ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col] * ZGStZGS[tmp1 + ind2 + 1] / ZGStZGS[tmp1]
				double* S2TransS2right = new double[expSq * expSq];
				matNmatNprod(ZGStZGSsGRow, invZGStmp1IntExpZGStZGS, S2TransS2right, expSq, intSq1, expSq);
				matAdd(expS2TransS2, S2TransS2right, expSq * expSq, -1.0);


				// invert (S2TransS2)
				matInv(expS2TransS2, expSq);
				// betaInt = invert(S2TransS2) * S2TransR
				matvecprod(expS2TransS2, expZGStR, betaInt[i], expSq, expSq);

				// Inv(S2TransS2) * S2DS2
				double* Stemp2 = new double[expSq * expSq];

				for (int j = 0; j < expSq; j++) {
					for (int k = 0; k < expSq; k++) {
						VarbetaInt[i][j * expSq + k] = sigma2 * expS2TransS2[j * expSq + k];
						InvVarbetaint[j * expSq + k] = VarbetaInt[i][j * expSq + k];
					}
				}

				// calculating P values
				double statM = betaM[i] * betaM[i] / VarbetaM[i];
				if (isnan(statM) || statM <= 0.0) {
					PvalM[i] = NAN;
				}
				else {
					PvalM[i] = boost::math::cdf(complement(chisq_dist_M, statM));
				}

				// invert VarbetaInt[i]
				matInv(InvVarbetaint, expSq);
				double* Stemp3 = new double[expSq];
				matvecprod(InvVarbetaint, betaInt[i], Stemp3, expSq, expSq);

				double statInt = 0.0;
				for (int j = 0; j < expSq; j++) {
					statInt += betaInt[i][j] * Stemp3[j];
				}

				if (isnan(statInt) || statInt <= 0.0) {
					PvalInt[i] = NAN;
				}
				else {
					PvalInt[i] = boost::math::cdf(complement(chisq_dist_Int, statInt));
				}

				double statJoint = statM + statInt;
				if (isnan(statJoint) || statJoint <= 0.0) {
					PvalJoint[i] = NAN;
				}
				else {
					PvalJoint[i] = boost::math::cdf(complement(chisq_dist_Joint, statJoint));
				}

				//delete[] S2TransS2;
				//delete[] S2TransR;
				///delete[] ZGStZGSsGsInt;
				delete[]Stemp2;
				delete[]Stemp3;
				delete[]InvVarbetaint;
				delete[]invZGStZGStmp1;
				delete[]ZGStZGSsGRow;

				/* For S2TransR*/
				delete[]expZGStR;
				delete[]int1ZGStR;
				delete[]invZGStmpInt1ZGStR;
				delete[]S2TransRright;

				/* For S2TransS2S2*/
				delete[]expS2TransS2;
				delete[]invZGStmp1IntExpZGStZGS;
				delete[]S2TransS2right;
			}
		}
		else if (robust == 1) {

			for (int i = 0; i < stream_snps; i++) {
				// initialize dynamic 2D array
				betaInt[i] = new double[expSq];
				VarbetaInt[i] = new double[expSq * expSq];

				//betamain
				int tmp1 = i * ZGS_col * Sq1 + i * Sq1;
				betaM[i] = ZGStR[i * Sq1] / ZGStZGS[tmp1];
				VarbetaM[i] = ZGSR2tZGS[tmp1] / (ZGStZGS[tmp1] * ZGStZGS[tmp1]);

				//double* S2TransS2 = new double[Sq * Sq];
				//double* S2TransR  = new double[Sq];
				//double* S2DS2     = new double[Sq * Sq];
				double* InvVarbetaint = new double[expSq * expSq];


				// inv(ZGStZGS[tmp1])
				double* invZGStZGStmp1 = new double[intSq1 * intSq1];
				subMatrix(ZGStZGS, invZGStZGStmp1, intSq1, intSq1, ZGS_col, intSq1, tmp1);
				matInv(invZGStZGStmp1, intSq1);
				// ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col]
				double* ZGStZGSsGRow = new double[expSq * intSq1];
				subMatrix(ZGStZGS, ZGStZGSsGRow, expSq, intSq1, ZGS_col, expSq, tmp1 + intSq1);


				/* For S2TransR*/
				// ZGStR[i * Sq1 + ind1 + 1]
				double* expZGStR = new double[expSq];
				subMatrix(ZGStR, expZGStR, expSq, 1, expSq, expSq, (i * Sq1) + intSq1);
				// ZGStR[i * Sq1]
				double* int1ZGStR = new double[intSq1];
				subMatrix(ZGStR, int1ZGStR, intSq1, 1, expSq + intSq, intSq1, (i * Sq1));
				// ZGStR[i * Sq1] / ZGStZGS[tmp1]
				double* invZGStmpInt1ZGStR = new double[intSq1];
				matvecprod(invZGStZGStmp1, int1ZGStR, invZGStmpInt1ZGStR, intSq1, intSq1);
				// ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col] * ZGStZGS[tmp1 + ind2 + 1] / ZGStZGS[tmp1];
				double* S2TransRright = new double[expSq];
				matNmatNprod(ZGStZGSsGRow, invZGStmpInt1ZGStR, S2TransRright, expSq, intSq1, 1);
				matAdd(expZGStR, S2TransRright, expSq, -1.0);


				/* For S2TransS2S2*/
				double* expS2TransS2 = new double[expSq * expSq];
				subMatrix(ZGStZGS, expS2TransS2, expSq, expSq, ZGS_col, expSq, tmp1 + (intSq1 * ZGS_col + intSq1));
				// ZGStZGS[tmp1 + ind2 + 1] / ZGStZGS[tmp1]
				double* invZGStmp1IntExpZGStZGS = new double[intSq1 * expSq];
				matTmatprod(invZGStZGStmp1, ZGStZGSsGRow, invZGStmp1IntExpZGStZGS, intSq1, intSq1, expSq);
				// ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col] * ZGStZGS[tmp1 + ind2 + 1] / ZGStZGS[tmp1]
				double* S2TransS2right = new double[expSq * expSq];
				matNmatNprod(ZGStZGSsGRow, invZGStmp1IntExpZGStZGS, S2TransS2right, expSq, intSq1, expSq);
				matAdd(expS2TransS2, S2TransS2right, expSq * expSq, -1.0);


				/* For S2DS2*/
				double* expS2DS2 = new double[expSq * expSq];
				subMatrix(ZGSR2tZGS, expS2DS2, expSq, expSq, ZGS_col, expSq, tmp1 + (intSq1 * ZGS_col + intSq1));
				// ZGSR2tZGS[tmp1 + ind2 + 1]
				double* intExpZGSR2tZGS = new double[expSq * intSq1];
				subMatrix(ZGSR2tZGS, intExpZGSR2tZGS, expSq, intSq1, ZGS_col, expSq, tmp1 + intSq1);
				// ZGSR2tZGS[tmp1 + (ind1 + 1) * ZGS_col]
				double* ZGSR2tZGSsGRow = new double[expSq * intSq1];
				subMatrix(ZGSR2tZGS, ZGSR2tZGSsGRow, expSq, intSq1, ZGS_col, expSq, tmp1 + intSq1);
				// ZGSR2tZGS[tmp1]
				double* ZGSR2tZGStmp1 = new double[intSq1 * intSq1];
				subMatrix(ZGSR2tZGS, ZGSR2tZGStmp1, intSq1, intSq1, ZGS_col, intSq1, tmp1);


				// ZGSR2tZGS[tmp1 + ind2 + 1] / ZGStZGS[tmp1]
				double* invZGStmp1IntExpZGSR2tZGS = new double[intSq1 * expSq];
				matTmatprod(invZGStZGStmp1, intExpZGSR2tZGS, invZGStmp1IntExpZGSR2tZGS, intSq1, intSq1, expSq);
				// ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col] * ZGSR2tZGS[tmp1 + ind2 + 1] / ZGStZGS[tmp1]
				double* S2DS2second = new double[expSq * expSq];
				matNmatNprod(ZGStZGSsGRow, invZGStmp1IntExpZGSR2tZGS, S2DS2second, expSq, intSq1, expSq);
				// ZGSR2tZGS[tmp1 + (ind1 + 1) * ZGS_col] * ZGStZGS[tmp1 + ind2 + 1] / ZGStZGS[tmp1]
				double* S2DS2third = new double[expSq * expSq];
				matNmatNprod(ZGSR2tZGSsGRow, invZGStmp1IntExpZGStZGS, S2DS2third, expSq, intSq1, expSq);
				// ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col] * ZGSR2tZGS[tmp1] * ZGStZGS[tmp1 + ind2 + 1] / ZGStZGS[tmp1] / ZGStZGS[tmp1]
				double* invZGStZGStmp1_intExpZGStZGS = new double[intSq1 * expSq];
				matTmatprod(invZGStZGStmp1, ZGStZGSsGRow, invZGStZGStmp1_intExpZGStZGS, intSq1, intSq1, expSq);
				double* ZGSR2tZGS_invZGStZGStmp1_intExpZGStZGS = new double[intSq1 * expSq];
				matNmatNprod(ZGSR2tZGStmp1, invZGStZGStmp1_intExpZGStZGS, ZGSR2tZGS_invZGStZGStmp1_intExpZGStZGS, intSq1, intSq1, expSq);
				double* S2DS2Forthtemp = new double[intSq1 * expSq];
				matNmatNprod(invZGStZGStmp1, ZGSR2tZGS_invZGStZGStmp1_intExpZGStZGS, S2DS2Forthtemp, intSq1, intSq1, expSq);
				double* S2DS2forth = new double[expSq * expSq];
				matNmatNprod(ZGStZGSsGRow, S2DS2Forthtemp, S2DS2forth, expSq, intSq1, expSq);
				matAdd(expS2DS2, S2DS2second, expSq * expSq, -1.0);
				matAdd(expS2DS2, S2DS2third, expSq * expSq, -1.0);
				matAdd(expS2DS2, S2DS2forth, expSq * expSq, 1.0);


				// invert (S2TransS2)
				matInv(expS2TransS2, expSq);

				// betaInt = invert(S2TransS2) * S2TransR
				matvecprod(expS2TransS2, expZGStR, betaInt[i], expSq, expSq);

				// Inv(S2TransS2) * S2DS2
				double* Stemp2 = new double[expSq * expSq];
				matmatprod(expS2TransS2, expS2DS2, Stemp2, expSq, expSq, expSq);

				// Stemp2 * Inv(S2TransS2)
				matNmatNprod(Stemp2, expS2TransS2, VarbetaInt[i], expSq, expSq, expSq);



				for (int j = 0; j < expSq; j++) {
					for (int k = 0; k < expSq; k++) {
						InvVarbetaint[j * expSq + k] = VarbetaInt[i][j * expSq + k];
					}
				}


				//calculating P values
				double statM = betaM[i] * betaM[i] / VarbetaM[i];
				if (isnan(statM) || statM <= 0.0) {
					PvalM[i] = NAN;
				}
				else {
					PvalM[i] = boost::math::cdf(complement(chisq_dist_M, statM));
				}


				// invert VarbetaInt[i]
				matInv(InvVarbetaint, expSq);
				double* Stemp3 = new double[expSq];
				matvecprod(InvVarbetaint, betaInt[i], Stemp3, expSq, expSq);

				double statInt = 0.0;
				for (int j = 0; j < expSq; j++) {
					statInt += betaInt[i][j] * Stemp3[j];
				}

				if (isnan(statInt) || statInt <= 0.0) {
					PvalInt[i] = NAN;
				}
				else {
					PvalInt[i] = boost::math::cdf(complement(chisq_dist_Int, statInt));
				}

				double statJoint = statM + statInt;
				if (isnan(statJoint) || statJoint <= 0.0) {
					PvalJoint[i] = NAN;
				}
				else {
					PvalJoint[i] = boost::math::cdf(complement(chisq_dist_Joint, statJoint));
				}

				//delete[] S2TransS2;
				//delete[] S2TransR;
				///delete[] ZGStZGSsGsInt;
				delete[]Stemp2;
				delete[]Stemp3;
				delete[]InvVarbetaint;
				delete[]invZGStZGStmp1;
				delete[]ZGStZGSsGRow;

				/* For S2TransR*/
				delete[]expZGStR;
				delete[]int1ZGStR;
				delete[]invZGStmpInt1ZGStR;
				delete[]S2TransRright;

				/* For S2TransS2S2*/
				delete[]expS2TransS2;
				delete[]invZGStmp1IntExpZGStZGS;
				delete[]S2TransS2right;

				/* For S2DS2*/
				delete[]expS2DS2;
				delete[]intExpZGSR2tZGS;
				delete[]ZGSR2tZGSsGRow;
				delete[]ZGSR2tZGStmp1;
				delete[]invZGStmp1IntExpZGSR2tZGS;
				delete[]S2DS2second;
				delete[]S2DS2third;
				delete[]invZGStZGStmp1_intExpZGStZGS;
				delete[]ZGSR2tZGS_invZGStZGStmp1_intExpZGStZGS;
				delete[]S2DS2Forthtemp;
				delete[]S2DS2forth;
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
	printExecutionTime(start_time, end_time);
}



#include "declars.h"
#include "ReadBGEN.h"




// Read BGEN v1.1, v1.2, and v1.3 header block and flags.
void Bgen::processBgenHeaderBlock(char genofile[300]) {


	// Ensure genotype file contains BGEN extension
	string genopath(genofile);
	if (genopath.substr(genopath.length() - 5, 5) != ".bgen") {
		cout << "\nError: " << genopath << " is not a .bgen file. Currently only supporting genotype files in BGEN format. \n\n";
		exit(1);
	}


	/***************************************************************
	  Read General information of bgen data.
	  Conduct Sample IDMatching process if necessary.
	***************************************************************/
	fin = fopen(genofile, "rb");
	if (fin == 0) {
		cerr << "\nError: BGEN file could not be opened.\n\n";
		exit(1);
	}



	cout << "General information of BGEN file. \n";


	// First four bytes (offset)
	fread(&offset, 4, 1, fin); cout << "offset: " << offset << '\n';
	

	// The header block
	uint L_H;      fread(&L_H, 4, 1, fin);   cout << "L_H: \n";
	fread(&Mbgen, 4, 1, fin); cout << "BGEN snpBlocks (Mbgen): " << Mbgen << '\n'; assert(Mbgen != 0);
	fread(&Nbgen, 4, 1, fin); cout << "BGEN samples (Nbgen): " << Nbgen << '\n';
	char magic[5]; fread(magic, 1, 4, fin); magic[4] = '\0';  //cout << "magic bytes: " << string(magic) << endl;
	fseek(fin, L_H - 20, SEEK_CUR);                           //cout << "skipping L_H-20 = " << L_H-20 << " bytes (free data area)" << endl;
	uint flags;    fread(&flags, 4, 1, fin);                  //cout << "flags: " << flags << endl;


	// The header block - flag definitions
	CompressedSNPBlocks = flags & 3; cout << "CompressedSNPBlocks: " << CompressedSNPBlocks << '\n';
	assert(CompressedSNPBlocks == 1); // REQUIRE CompressedSNPBlocks==1
	Layout = (flags >> 2) & 0xf; cout << "Layout: " << Layout << '\n';
	assert(Layout == 1 || Layout == 2); // REQUIRE Layout==1 or Layout==2
	SampleIdentifiers = flags >> 31; cout << "SampleIdentifiers: " << SampleIdentifiers << '\n';


}







// This functions reads the sample block of BGEN v1.1, v1.2, and v1.3. Also finds which samples to remove if they have missing values in the pheno file.
void Bgen::processBgenSampleBlock(Bgen bgen, char samplefile[300], unordered_map<string, vector<string>> phenomap, string phenoMissingKey, vector<double> phenodata, vector<double> covdata, int numSelCol, int samSize) {


	int k = 0;
	uint maxLA = 65536;
	char* samID = new char[maxLA + 1];
	unordered_set<int> genoUnMatchID;

	
	if (bgen.SampleIdentifiers == 1) {


		uint LS1;  fread(&LS1,  4, 1, bgen.fin); // std::cout << "LS1: " << LS1 << std::endl; // LS1 + L_H <= offset
		uint Nrow; fread(&Nrow, 4, 1, bgen.fin); // cout << "Nrow: " << Nrow << " " << std::flush;
		if (Nrow != bgen.Nbgen) {
			cerr << "\nERROR: Nrow = " << Nrow << " does not match Nbgen = " << bgen.Nbgen << "\n\n";
			exit(1);
		}


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
					for (int c = 0; c < numSelCol; c++) {
						sscanf(tmp_valvec[c + 1].c_str(), "%lf", &covdata[k * (numSelCol + 1) + c + 1]);
					}

					k++;
				}

				// erase the used element in phenomap
				phenomap.erase(strtmp);
			 }

			 // save the index with unmatched ID into genoUnMatchID.
			 if (itmp == k) {
				 genoUnMatchID.insert(m);
			 }
		}
		
	} // end SampleIdentifiers == 1
	else {

		if (bgen.SampleIdentifiers == 0 && samplefile[0] == '\0') {
			cerr << "\nError: BGEN file does not contain sample identifiers (SampleIdentifiers == 0). A .sample file is required. \n"
				<< "          See https://www.well.ox.ac.uk/~gav/qctool/documentation/sample_file_formats.html for .sample file format. \n\n";
			exit(1);
		}


		std::ifstream fIDMat;
		fIDMat.open(samplefile);

		if (!fIDMat.is_open()) {
			cerr << "\nError: Sample file could not be opened." << endl << endl;
			exit(1);
		}


		string IDline;
		getline(fIDMat, IDline);
		getline(fIDMat, IDline);


		for (uint m = 0; m < bgen.Nbgen; m++) {
			 // IDMatching
			 getline(fIDMat, IDline);
			 std::istringstream iss(IDline);
			 string strtmp;
			 iss >> strtmp;

			 int itmp = k;
			 if (phenomap.find(strtmp) != phenomap.end()) {
				 auto tmp_valvec = phenomap[strtmp];
				 if (find(tmp_valvec.begin(), tmp_valvec.end(), phenoMissingKey) == tmp_valvec.end()) {
					 sscanf(tmp_valvec[0].c_str(), "%lf", &phenodata[k]);
					 // covdata[k*(numSelCol+1) + 0] = 1.0;
					 for (int c = 0; c < numSelCol; c++) {
						  sscanf(tmp_valvec[c + 1].c_str(), "%lf", &covdata[k * (numSelCol + 1) + c + 1]);
					 }
					k++;
				 }

				// erase the used element in phenomap
				phenomap.erase(strtmp);
			 }

			// save the index with unmatched ID into genoUnMatchID.
			if (itmp == k) genoUnMatchID.insert(m);
		}
		fIDMat.close();
	}

	// After IDMatching, resizing phenodata and covdata, and updating samSize;
	phenodata.resize(k);
	covdata.resize(k * (numSelCol + 1));
	samSize = k;

	new_samSize = samSize;
	new_covdata = covdata;
	new_phenodata = phenodata;


	if (samSize == 0) {
		cerr << "\nError: Sample size changed from" << samSize + genoUnMatchID.size() << " to " << samSize 
			 << "         Check if sample IDs are consistent across files. \n\n";
		exit(1);
	}


	int ii = 0;
	include_idx.resize(samSize);
	for (uint i = 0; i < bgen.Nbgen; i++) {
		 if (genoUnMatchID.find(i) == genoUnMatchID.end()) {
			 include_idx[ii] = i;
		 	 ii++;
		}
	}



	cout << "****************************************************************************\n";
	if (genoUnMatchID.empty()) {
		cout << "After processes of sample IDMatching and checking missing values, the sample size does not change.\n\n";
	} else {
		cout << "After processes of sample IDMatching and checking missing values, the sample size changes from "
			<< samSize + genoUnMatchID.size() << " to " << samSize << ".\n\n";
	}
	cout << "Sample IDMatching and checking missing values processes have been completed.\n";
	cout << "New pheno and covariate data vectors with the same order of sample ID sequence of geno data are updated.\n";
	cout << "****************************************************************************\n";
}








// This function reads just the variant block for BGEN files version v1.1, v1.2, and v1.3 and is used to grab the byte where the variant begins.
std::vector<long int> getPositionOfBgenVariant(FILE* fin, uint offset, uint Mbgen, uint Nbgen, uint CompressedSNPBlocks, uint Layout, vector<int> Mbgen_begin) {

	int t = 0;
	uint maxLA = 65536;
	uint maxLB = 65536;
	char* snpID   = new char[maxLA + 1];
	char* rsID    = new char[maxLA + 1];
	char* chrStr  = new char[maxLA + 1];
	char* allele1 = new char[maxLA + 1];
	char* allele0 = new char[maxLB + 1];

	vector <uchar> zBuf;
	vector <long int> variant_pos(Mbgen_begin.size());

	
	fseek(fin, offset + 4, SEEK_SET);

	for (int snploop = 0; snploop < Mbgen; snploop++) {

		 if (snploop == Mbgen_begin[t]) {

			 variant_pos[t] = ftell(fin);
			 t++;

			 if (t == (Mbgen_begin.size())) {
				 break;
			 }

		 }
		 
		
		 /**** Variant Data Block ********/

		 // Number of individuals. Only present when Layout == 1
		 if (Layout == 1) {
			 uint Nrow; fread(&Nrow, 4, 1, fin); // cout << "Nrow: " << Nrow << " " << std::flush;  
			 if (Nrow != Nbgen) {
				 cerr << "ERROR: Nrow = " << Nrow << " does not match Nbgen = " << Nbgen << '\n';
				 exit(1);
			 }
		 }

		 // The length of the variant identifier
		 ushort LS; fread(&LS, 2, 1, fin);  // cout << "LS: " << LS << " " << std::flush;
		 if (LS > maxLA) {
			 maxLA = 2 * LS;
			 delete[] snpID;
			 char* snpID = new char[maxLA + 1];
		 }
		 // The variant identifier
		 fread(snpID, 1, LS, fin); snpID[LS] = '\0'; // cout << "snpID: " << string(snpID) << " " << std::flush;
		


		 // The length of the rsid
		 ushort LR; fread(&LR, 2, 1, fin); // cout << "LR: " << LR << " " << std::flush;
		 if (LR > maxLA) {
			 maxLA = 2 * LR;
			 delete[] rsID;
			 char* rsID = new char[maxLA + 1];
		 }
		 // The rsid
		 fread(rsID, 1, LR, fin); rsID[LR] = '\0'; // cout << "rsID: " << string(rsID) << " " << std::flush;




		 // The length of the chromosome
		 ushort LC; fread(&LC, 2, 1, fin); // cout << "LC: " << LC << " " << std::flush;
		 // The chromosome
		 fread(chrStr, 1, LC, fin); chrStr[LC] = '\0';



		 // The variant position
		 uint physpos; fread(&physpos, 4, 1, fin); // cout << "physpos: " << physpos << " " << std::flush;





		 // The number of alleles if Layout == 2. If Layout == 1, this value is assumed to be 2
		 if (Layout == 2) {
			 ushort LKnum; fread(&LKnum, 2, 1, fin); // this is for Layout = 2, Lnum = 2 is Layout = 1

			 if (LKnum != 2) {
				 cerr << "\nError: Non-bi-allelic variant found: " << LKnum << " alleles\n\n";
				 exit(1);
			 }
		 }


		 // Length of the first allele
		 uint LA; fread(&LA, 4, 1, fin); // cout << "LA: " << LA << " " << std::flush;
		 if (LA > maxLA) {
			 maxLA = 2 * LA;
			 delete[] allele1;
			 char* allele1 = new char[maxLA + 1];
		 }
		 // The first allele
		 fread(allele1, 1, LA, fin); allele1[LA] = '\0';


		 // The length of the second allele
		 uint LB; fread(&LB, 4, 1, fin); // cout << "LB: " << LB << " " << std::flush;
		 if (LB > maxLB) {
			 maxLB = 2 * LB;
			 delete[] allele0;
			 char* allele0 = new char[maxLB + 1];
		 }
		 // The second allele
		 fread(allele0, 1, LB, fin); allele0[LB] = '\0';




		 // Seeks past the uncompressed genotype.
		 uint zLen;  fread(&zLen, 4, 1, fin); // cout << "zLen: " << zLen << endl;
		 uint DLen;
		 if (CompressedSNPBlocks == 0) {
			 fseek(fin, zLen, SEEK_CUR);

		 } else {
			 fseek(fin, 4 + zLen - 4, SEEK_CUR);
		 }

	}


	return variant_pos;
}
















void BgenParallelGWAS(int begin, int end, long int byte, char genobgen[300], int thread_num, Bgen test) {


	auto start_time = std::chrono::high_resolution_clock::now();
	std::string output = "gem_bin_" + std::to_string(thread_num) + ".tmp";
	std::ofstream results(output, std::ofstream::binary);
	std::ostringstream oss;

	vector <uchar> shortBuf;
	uint maxLA = 65536;
	uint maxLB = 65536;
	char* snpID   = new char[maxLA + 1];
	char* rsID    = new char[maxLA + 1];
	char* chrStr  = new char[maxLA + 1];
	char* allele1 = new char[maxLA + 1];
	char* allele0 = new char[maxLB + 1];
	vector <uchar> zBuf;
	string physpos_tmp;

	int Sq = test.Sq;
	int phenoType   = test.phenoTyp;
	int robust      = test.robust;
	int samSize     = test.new_samSize;
	int numSelCol   = test.numSelCol;
	int stream_snps = test.stream_snps;
	double MAF      = test.maf;
	double sigma2   = test.sigma2;
	double* resid   = test.resid;
	double* XinvXTX = test.XinvXTX;

	uint Nbgen  = test.Nbgen;
	uint Layout = test.Layout;
	uint CompressedSNPBlocks = test.CompressedSNPBlocks;
	vector <double> miu = test.miu;
	vector<long int> include_idx = test.include_idx;

	double* covX = test.covX;

	vector <double> ZGSvec(samSize * (1 + Sq) * stream_snps);
	vector <double> ZGSR2vec(samSize * (1 + Sq) * stream_snps);
	vector <double> WZGSvec(samSize * (1 + Sq) * stream_snps);


	double* WZGS = &WZGSvec[0];
	int Sq1 = Sq + 1;
	
	
	vector <string> geno_snpid(stream_snps);

	FILE* fin3;
	fin3 = fopen(genobgen, "rb");


	fseek(fin3, byte, SEEK_SET);

	

	int snploop = begin;
	int variant_index = 0;
	while (snploop <= end) {

		int Sq1 = Sq + 1;
		int ZGS_col = Sq1 * stream_snps;
		vector <double> AF(stream_snps);

		int stream_i = 0;
		while (stream_i < stream_snps) {

			if (snploop == (end + 1) && stream_i == 0) {
				break;
			}

			if (snploop == end + 1 && stream_i != 0) {
				stream_snps = stream_i;
				Sq1 = Sq + 1;
				ZGS_col = Sq1 * stream_snps;
				vector <double> AF(stream_snps);

				break;
			}
			snploop++;


			// Number of individuals. Only present when Layout == 1
			if (Layout == 1) {
				uint Nrow;  fread(&Nrow, 4, 1, fin3); // cout << "Nrow: " << Nrow << " " << std::flush;  
				if (Nrow != Nbgen) {
					cerr << "\nError: Nrow = " << Nrow << " does not match Nbgen = " << Nbgen << "\n\n";
					exit(1);
				}
			}


			// The length of the variant identifier
			ushort LS;  fread(&LS, 2, 1, fin3);  // cout << "LS: " << LS << " " << std::flush;
			if (LS > maxLA) {
				maxLA = 2 * LS;
				delete[] snpID;
				snpID = new char[maxLA + 1];
			}

			// The variant identifier
			fread(snpID, 1, LS, fin3); snpID[LS] = '\0'; // cout << "snpID: " << string(snpID) << " " << std::flush;


			// The length of the rsid
			ushort LR;  fread(&LR, 2, 1, fin3); // cout << "LR: " << LR << " " << std::flush;
			if (LR > maxLA) {
				maxLA = 2 * LR;
				delete[] rsID;
				rsID = new char[maxLA + 1];
			}

			// The rsid
			fread(rsID, 1, LR, fin3); rsID[LR] = '\0'; // cout << "rsID: " << string(rsID) << " " << std::flush;


			// The length of the chromosome
			ushort LC;  fread(&LC, 2, 1, fin3); // cout << "LC: " << LC << " " << std::flush;

			// The chromosome
			fread(chrStr, 1, LC, fin3); chrStr[LC] = '\0';

			// The variant position
			uint physpos;  fread(&physpos, 4, 1, fin3); // cout << "physpos: " << physpos << " " << std::flush;
			physpos_tmp = std::to_string(physpos);


			// The number of alleles if Layout == 2. If Layout == 1, this value is assumed to be 2 and not stored
			if (Layout == 2) {
				ushort LKnum;  fread(&LKnum, 2, 1, fin3); // this is for Layout = 2, Lnum = 2 is Layout = 1
				if (LKnum != 2) {
					cerr << "\nError: Non-bi-allelic variant found: " << LKnum << " alleles \n\n";
					exit(1);
				}
			}

			// Length of the first allele
			ushort LA;  fread(&LA, 4, 1, fin3); // cout << "LA: " << LA << " " << std::flush;
			if (LA > maxLA) {
				maxLA = 2 * LA;
				delete[] allele1;
				allele1 = new char[maxLA + 1];
			}
			// The first allele
			fread(allele1, 1, LA, fin3); allele1[LA] = '\0';


			// The length of the second allele
			ushort LB; fread(&LB, 4, 1, fin3); // cout << "LB: " << LB << " " << std::flush;
			if (LB > maxLB) {
				maxLB = 2 * LB;
				delete[] allele0;
				allele0 = new char[maxLB + 1];
			}
			// The second allele
			fread(allele0, 1, LB, fin3); allele0[LB] = '\0';



			//cout << string((*tss).snpID) + "\t" + string((*tss).rsID) + "\t" + string((*tss).chrStr) + "\t" + std::to_string((*tss).physpos) + "\t" + string((*tss).allele1) + "\t" + string((*tss).allele0);
			if (Layout == 1) {

				uint zLen; fread(&zLen, 4, 1, fin3); // cout << "zLen: " << zLen << endl;
				fread(&zBuf[0], 1, zLen, fin3);

				uLongf destLen = 6 * Nbgen;

				if (uncompress(&shortBuf[0], &destLen, &zBuf[0], zLen) != Z_OK || destLen != 6 * Nbgen) {
					cerr << "\nError: uncompress() failed\n\n";
					exit(1);
				}


				// read genotype probabilities
				const double scale = 1.0 / 32768;
				int tmp1 = stream_i * Sq1 * samSize;

				//if (IDMatching == 1) {
				int idx_k = 0;
				for (uint i = 0; i < Nbgen; i++) {

					 if (include_idx[idx_k] == i) {
						 double p11 = shortBuf[3 * i] * scale;
						 double p10 = shortBuf[3 * i + 1] * scale;
						 double p00 = shortBuf[3 * i + 2] * scale;

						 double pTot = p11 + p10 + p00;
						 double dosage = (2 * p00 + p10) / pTot;

						 int tmp2 = idx_k + tmp1;
						 AF[stream_i] += dosage;

						 if (phenoType == 1) {
							 ZGSvec[tmp1] = miu[idx_k] * (1 - miu[idx_k]) * dosage;
						 }
						 else {
							 ZGSvec[tmp1] = dosage;
						 }

						 idx_k++;
					 }
				}

				if ((AF[stream_i] / 2 / samSize) < MAF || (AF[stream_i] / 2 / samSize) > (1 - MAF)) {
					AF[stream_i] = 0;
					continue;
				} else {
					AF[stream_i] = AF[stream_i] / 2 / samSize;
				}

				for (int j = 0; j < Sq; j++) {
					 int tmp3 = samSize * (j + 1) + tmp1;

					 for (uint i = 0; i < samSize; i++) {
						  int tmp4 = i * (numSelCol + 1);
						  ZGSvec[tmp3 + i] = covX[tmp4 + j + 1] * ZGSvec[tmp1 + i]; // here we save ZGS in column wise
					 }

				}

			} // end of reading genotype data when Layout = 1




			if (Layout = 2) {

				uint zLen; fread(&zLen, 4, 1, fin3); // cout << "zLen: " << zLen << endl;
				zBuf.resize(zLen - 4);
				uint DLen;
				if (CompressedSNPBlocks == 0) {
					DLen = zLen;
					fread(&zBuf[0], 1, zLen, fin3);
				}
				else {
					fread(&DLen, 4, 1, fin3);
					fread(&zBuf[0], 1, zLen - 4, fin3);
				}



				uLongf destLen = DLen; //6*Nbgen;
				shortBuf.resize(DLen);


				if (uncompress(&shortBuf[0], &destLen, &zBuf[0], zLen - 4) != Z_OK || destLen != DLen) {
					cout << "destLen: " << destLen << " " << zLen - 4 << "\n\n";
					cerr << "\nError: uncompress() failed\n";
					exit(1);
				}


				// read genotype probabilities
				uchar* bufAt = &shortBuf[0];
				uint N = bufAt[0] | (bufAt[1] << 8) | (bufAt[2] << 16) | (bufAt[3] << 24); bufAt += 4;
				if (N != Nbgen) {
					cerr << "\nError: " << "snpName " << " has N = " << N << " (mismatch with header block)\n\n";
					exit(1);
				}



				uint K = bufAt[0] | (bufAt[1] << 8); bufAt += 2;
				if (K != 2U) {
					cout << "\nError: There are SNP(s) with more than 2 alleles (non-bi-allelic). . Currently unsupported. \n\n";
					exit(1);
				}


				uint Pmin = *bufAt; bufAt++;
				if (Pmin > 2U) {
					cerr << "\nError: " << snpID << " has minimum ploidy = " << Pmin << ". Currently unsupported. \n\n";
					exit(1);
				}


				uint Pmax = *bufAt; bufAt++;
				if (Pmax > 2U) {
					cerr << "\nError: " << snpID << " has maximum ploidy = " << Pmax << ". Currently unsupported. \n\n";
					exit(1);
				}




				vector<unsigned int> ploidy_and_missing_info(N);
				int ploidy_sum = 0;
				int idx_to_sum = 0;
				for (uint i = 0; i < N; i++) {
					 uint ploidyMiss = *bufAt;
					 ploidy_and_missing_info[i] = ploidyMiss;

					if (ploidyMiss > 2U) {
						std::cerr << "\nError: " << snpID << " has ploidy/missingness byte = " << ploidyMiss
							<< " (not 2) \n\n";
						exit(1);
					}

					if (include_idx[idx_to_sum] == i) {
						ploidy_sum += ploidyMiss;
						idx_to_sum++;
					}

					bufAt++;
				}


				// Phased information  indicating what is stored in the row. 
				//    Phased = 1; row stores one probability per allele
				//    Phased = 0; row stores one probability per possible genotype
				//    Everything else is error.
				uint Phased = *bufAt; bufAt++;
				if (Phased != 0U && Phased != 1U) {
					//if(Phased == 1U){ cerr << "\nERROR: Phased data is not currently supported by GEM.\n"; exit(1);}
					cerr << "\nError: " << snpID << " has Phased = " << Phased << ". Must be 0 or 1. \n"
						<< "       See https://www.well.ox.ac.uk/~gav/bgen_format/spec/latest.html for more details. \n\n";
					exit(1);
				}

				uint B = *bufAt; bufAt++;
				uint Bbits = std::pow(2, B);
				if ((B != 8U) && (B != 16U) && (B != 24U) && (B != 32U)) {
					std::cerr << "\nError: " << "snpName " << " has B = " << B << " (not divisible by 8)\n\n";
					exit(1);
				}



				int tmp1 = stream_i * Sq1 * samSize;
				int idx_k = 0;
				for (uint i = 0; i < N; i++) {
					 uint chartem;
					 uint chartem1;

					if (B == 8U)
						chartem = bufAt[0];
					else if (B == 16U)
						chartem = bufAt[0] | (bufAt[1] << 8);
					else if (B == 24U)
						chartem = bufAt[0] | (bufAt[1] << 8) | (bufAt[2] << 16);
					else if (B == 32U)
						chartem = bufAt[0] | (bufAt[1] << 8) | (bufAt[2] << 16) | (bufAt[3] << 24);
					bufAt += B / 8;


					if (ploidy_and_missing_info[i] == 2U) {
						if (B == 8U)
							chartem1 = bufAt[0];
						else if (B == 16U)
							chartem1 = bufAt[0] | (bufAt[1] << 8);
						else if (B == 24U)
							chartem1 = bufAt[0] | (bufAt[1] << 8) | (bufAt[2] << 16);
						else if (B == 32U)
							chartem1 = bufAt[0] | (bufAt[1] << 8) | (bufAt[2] << 16) | (bufAt[3] << 24);
						bufAt += B / 8;
					}




					if (include_idx[idx_k] == i) {
						double p11 = chartem / double(1.0 * (Bbits - 1));
						double p10 = chartem1 / double(1.0 * (Bbits - 1));
						double dosage;

						if (Phased == 1U) {
							if (ploidy_and_missing_info[i] == 1U) {
								dosage = 1 - p11;
							}
							else {
								dosage = 2 - (p11 + p10);
							}

						}
						else {
							dosage = 2 * (1 - p11 - p10) + p10;
						}

						int tmp2 = idx_k + tmp1;
						AF[stream_i] += dosage;

						if (phenoType == 1) {
							ZGSvec[tmp2] = miu[idx_k] * (1 - miu[idx_k]) * dosage;
						}
						else {
							ZGSvec[tmp2] = dosage; // replace your new data from other genotype files here.
						}

						idx_k++;
					}

				}

				//cout << string(snpID) + "\t" + string(rsID) + "\t" + string(chrStr) + "\t" + std::to_string(physpos) + "\t" + string(allele1) + "\t" + string(allele0) + "\t AF: " + std::to_string(cur_AF / ploidy_sum) << "\n";
				if ((AF[stream_i] / ploidy_sum) < MAF || (AF[stream_i] / ploidy_sum) > (1 - MAF)) {
					AF[stream_i] = 0;
					continue;
				}
				else {
					AF[stream_i] = AF[stream_i] / ploidy_sum;
				}
		
				geno_snpid[stream_i] = string(snpID) + "\t" + string(rsID) + "\t" + string(chrStr) + "\t" + physpos_tmp + "\t" + string(allele1) + "\t" + string(allele0);



				for (int j = 0; j < Sq; j++) {
					 int tmp3 = samSize * (j + 1) + tmp1;
					 for (uint i = 0; i < samSize; i++) {
						int tmp4 = i * (numSelCol + 1);
						ZGSvec[tmp3 + i] = covX[tmp4 + j + 1] * ZGSvec[tmp1 + i]; // here we save ZGS in column wise
					 }
				}


			} // end of layout 2
	
			variant_index++;
			stream_i++;
		} // end of stream_i

		if (snploop == end + 1 & stream_i == 0) {
			break;
		}


		/***************************************************************/
		//	genodata and envirment data
		double* ZGS = &ZGSvec[0];

		// transpose(X) * ZGS
		// it is non-squred matrix, attention that continuous memory is column-major due to Fortran in BLAS.
		// important!!!!
		double* XtransZGS = new double[(numSelCol + 1) * ZGS_col];
		matNmatNprod(covX, ZGS, XtransZGS, numSelCol + 1, samSize, ZGS_col);



		if (phenoType == 0) {
			for (int j = 0; j < ZGS_col; j++) {
				for (int i = 0; i < samSize; i++) {
					for (int k = 0; k <= numSelCol; k++) {
						ZGS[j * samSize + i] -= XinvXTX[k * samSize + i] * XtransZGS[j * (numSelCol + 1) + k];
					}
					if (robust == 1) { ZGSR2vec[j * samSize + i] = ZGS[j * samSize + i] * resid[i] * resid[i]; }
				}
			}
		}
		else if (phenoType == 1) {

			for (int j = 0; j < ZGS_col; j++) {
				for (int i = 0; i < samSize; i++) {
					double ZGStemp = 0.0;
					for (int k = 0; k <= numSelCol; k++) {
						if (k == 0) {
							ZGStemp += ZGS[j * samSize + i] / miu[i] / (1.0 - miu[i]) - XinvXTX[k * samSize + i] * XtransZGS[j * (numSelCol + 1) + k];
						}
						else {
							ZGStemp -= XinvXTX[k * samSize + i] * XtransZGS[j * (numSelCol + 1) + k];
						}
					}

					ZGS[j * samSize + i] = ZGStemp;
					if (robust == 1) ZGSR2vec[j * samSize + i] = ZGS[j * samSize + i] * resid[i] * resid[i];
					WZGS[j * samSize + i] = miu[i] * (1 - miu[i]) * ZGS[j * samSize + i];
				}
			}
		}



		double* ZGSR2 = &ZGSR2vec[0];
		delete[] XtransZGS;
		// transpose(ZGS) * resid
		double* ZGStR = new double[ZGS_col];
		matvecprod(ZGS, resid, ZGStR, ZGS_col, samSize);
		// transpose(ZGS) * ZGS
		double* ZGStZGS = new double[ZGS_col * ZGS_col];
		if (phenoType == 0) {
			matmatTprod(ZGS, ZGS, ZGStZGS, ZGS_col, samSize, ZGS_col);
		}
		else if (phenoType == 1) {
			matmatTprod(ZGS, WZGS, ZGStZGS, ZGS_col, samSize, ZGS_col);
		}
		else {
			cout << "phenoTyp is not equal to 0 or 1. Kill the job!! \n";
			exit(1);
		}


		// transpose(ZGSR2) * ZGS
		double* ZGSR2tZGS = new double[ZGS_col * ZGS_col];
		if (robust == 1) matmatTprod(ZGSR2, ZGS, ZGSR2tZGS, ZGS_col, samSize, ZGS_col);





		double*  betaM      = new double[stream_snps];
		double*  VarbetaM   = new double[stream_snps];
		double** betaInt    = new double* [stream_snps];
		double** VarbetaInt = new double* [stream_snps];
		double*  PvalM      = new double[stream_snps];
		double*  PvalInt    = new double[stream_snps];
		double*  PvalJoint  = new double[stream_snps];
		boost::math::chi_squared chisq_dist_M(1);
		boost::math::chi_squared chisq_dist_Int(Sq);
		boost::math::chi_squared chisq_dist_Joint(1 + Sq);




		if (robust == 0) {

			for (int i = 0; i < stream_snps; i++) {
				// initialize dynamic 2D array
				betaInt[i] = new double[Sq];
				VarbetaInt[i] = new double[Sq * Sq];

				// betamain
				int tmp1 = i * ZGS_col * Sq1 + i * Sq1;
				betaM[i] = ZGStR[i * Sq1] / ZGStZGS[tmp1];
				VarbetaM[i] = sigma2 / ZGStZGS[tmp1];

				double* S2TransS2 = new double[Sq * Sq];
				double* S2TransR = new double[Sq];
				double* S2DS2 = new double[Sq * Sq];
				double* InvVarbetaint = new double[Sq * Sq];
				for (int ind1 = 0; ind1 < Sq; ind1++) {
					for (int ind2 = 0; ind2 < Sq; ind2++) {
						// transpose(Snew2) * Snew2
						S2TransS2[ind1 * Sq + ind2] = ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col + ind2 + 1] - ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col] * ZGStZGS[tmp1 + ind2 + 1] / ZGStZGS[tmp1];
					}
					//transpose(Snew2) * resid
					S2TransR[ind1] = ZGStR[i * Sq1 + ind1 + 1] - ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col] * ZGStR[i * Sq1] / ZGStZGS[tmp1];
				}

				// invert (S2TransS2)
				matInv(S2TransS2, Sq);

				// betaInt = invert(S2TransS2) * S2TransR
				matvecprod(S2TransS2, S2TransR, betaInt[i], Sq, Sq);

				// Inv(S2TransS2) * S2DS2
				double* Stemp2 = new double[Sq * Sq];

				for (int j = 0; j < Sq; j++) {
					for (int k = 0; k < Sq; k++) {
						VarbetaInt[i][j * Sq + k] = sigma2 * S2TransS2[j * Sq + k];
						InvVarbetaint[j * Sq + k] = VarbetaInt[i][j * Sq + k];
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
				matInv(InvVarbetaint, Sq);
				double* Stemp3 = new double[Sq];
				matvecprod(InvVarbetaint, betaInt[i], Stemp3, Sq, Sq);

				double statInt = 0.0;
				for (int j = 0; j < Sq; j++) {
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

				delete[] S2TransS2;
				delete[] S2TransR;
				delete[] S2DS2;
				delete[]Stemp2;
				delete[]Stemp3;
				delete[]InvVarbetaint;
			}
		}
		else if (robust == 1) {

			for (int i = 0; i < stream_snps; i++) {
				// initialize dynamic 2D array
				betaInt[i] = new double[Sq];
				VarbetaInt[i] = new double[Sq * Sq];

				//betamain
				int tmp1 = i * ZGS_col * Sq1 + i * Sq1;
				betaM[i] = ZGStR[i * Sq1] / ZGStZGS[tmp1];
				VarbetaM[i] = ZGSR2tZGS[tmp1] / (ZGStZGS[tmp1] * ZGStZGS[tmp1]);

				double* S2TransS2 = new double[Sq * Sq];
				double* S2TransR = new double[Sq];
				double* S2DS2 = new double[Sq * Sq];
				double* InvVarbetaint = new double[Sq * Sq];

				for (int ind1 = 0; ind1 < Sq; ind1++) {
					for (int ind2 = 0; ind2 < Sq; ind2++) {
						// transpose(Snew2) * Snew2
						S2TransS2[ind1 * Sq + ind2] = ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col + ind2 + 1] - ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col] * ZGStZGS[tmp1 + ind2 + 1] / ZGStZGS[tmp1];
						// transpose(Snew2) * D * Snew2
						S2DS2[ind1 * Sq + ind2] = ZGSR2tZGS[tmp1 + (ind1 + 1) * ZGS_col + ind2 + 1] - ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col] * ZGSR2tZGS[tmp1 + ind2 + 1] / ZGStZGS[tmp1] - ZGSR2tZGS[tmp1 + (ind1 + 1) * ZGS_col] * ZGStZGS[tmp1 + ind2 + 1] / ZGStZGS[tmp1] + ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col] * ZGSR2tZGS[tmp1] * ZGStZGS[tmp1 + ind2 + 1] / ZGStZGS[tmp1] / ZGStZGS[tmp1];
					}
					//transpose(Snew2) * resid
					S2TransR[ind1] = ZGStR[i * Sq1 + ind1 + 1] - ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col] * ZGStR[i * Sq1] / ZGStZGS[tmp1];
				}

				// invert (S2TransS2)
				matInv(S2TransS2, Sq);

				// betaInt = invert(S2TransS2) * S2TransR
				matvecprod(S2TransS2, S2TransR, betaInt[i], Sq, Sq);

				// Inv(S2TransS2) * S2DS2
				double* Stemp2 = new double[Sq * Sq];
				matmatprod(S2TransS2, S2DS2, Stemp2, Sq, Sq, Sq);

				// Stemp2 * Inv(S2TransS2)
				matNmatNprod(Stemp2, S2TransS2, VarbetaInt[i], Sq, Sq, Sq);



				for (int j = 0; j < Sq; j++) {
					for (int k = 0; k < Sq; k++) {
						InvVarbetaint[j * Sq + k] = VarbetaInt[i][j * Sq + k];
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
				matInv(InvVarbetaint, Sq);
				double* Stemp3 = new double[Sq];
				matvecprod(InvVarbetaint, betaInt[i], Stemp3, Sq, Sq);

				double statInt = 0.0;
				for (int j = 0; j < Sq; j++) {
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


				delete[] S2TransS2;
				delete[] S2TransR;
				delete[] S2DS2;
				delete[] Stemp2;
				delete[] Stemp3;
				delete[] InvVarbetaint;
			}
		} // end of if robust == 1



		for (int i = 0; i < stream_snps; i++) {
			oss << geno_snpid[i] << "\t" << AF[i] << "\t" << betaM[i] << "\t" << VarbetaM[i] << "\t";
			for (int ii = 0; ii < Sq; ii++) {
				 oss << betaInt[i][ii] << "\t";
			}

			for (int ii = 0; ii < Sq; ii++) {
				for (int jj = 0; jj < Sq; jj++) {
					oss << VarbetaInt[i][ii * Sq + jj] << "\t";
				}
			}
			oss << PvalM[i] << "\t" << PvalInt[i] << "\t" << PvalJoint[i] << '\n';
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
		AF.clear();


		if ((variant_index % 10000 == 0) || snploop == end + 1) {
			results << oss.str();
			oss.str(std::string());
			oss.clear();
			variant_index = 0;
		}
	} // end of snploop 


	// Close files
	results.close();
	fclose(fin3);

	auto end_time = std::chrono::high_resolution_clock::now();
	cout << "Thread " << thread_num << " finished in ";
	printExecutionTime(start_time, end_time);
}






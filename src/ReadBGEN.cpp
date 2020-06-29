#include "declars.h"
#include "ReadBGEN.h"
#include "zstd.h"
#include "../thirdparty/libdeflate/libdeflate.h"

/**************************************
This function is revised based on the Parse function in BOLT-LMM v2.3 source code
*************************************/

// Read BGEN v1.1, v1.2, and v1.3 header block and flags.
void Bgen::processBgenHeaderBlock(char genofile[300]) {


	// Ensure genotype file contains BGEN extension
	string genopath(genofile);
	if (genopath.substr(genopath.length() - 5, 5) != ".bgen") {
		cout << "\nERROR: " << genopath << " is not a .bgen file. Currently only supporting genotype files in BGEN format. \n\n";
		exit(1);
	}


	/***************************************************************
	  Read General information of bgen data.
	  Conduct Sample IDMatching process if necessary.
	***************************************************************/
	fin = fopen(genofile, "rb");
	if (fin == 0) {
		cerr << "\nERROR: BGEN file could not be opened.\n\n";
		exit(1);
	}



	cout << "General information of BGEN file. \n";


	// First four bytes (offset)
	fread(&offset, 4, 1, fin);
	

	// The header block
	uint L_H;      fread(&L_H, 4, 1, fin);
	fread(&Mbgen, 4, 1, fin); cout << "Number of variants: " << Mbgen << '\n'; assert(Mbgen != 0);
	fread(&Nbgen, 4, 1, fin); cout << "Number of samples: " << Nbgen << '\n';
	char magic[5]; fread(magic, 1, 4, fin); magic[4] = '\0';  //cout << "magic bytes: " << string(magic) << endl;
	fseek(fin, L_H - 20, SEEK_CUR);                           //cout << "skipping L_H-20 = " << L_H-20 << " bytes (free data area)" << endl;
	uint flags;    fread(&flags, 4, 1, fin);                  //cout << "flags: " << flags << endl;


	// The header block - flag definitions
	CompressedSNPBlocks = flags & 3; cout << "Genotype Block Compression Type: ";
	switch (CompressedSNPBlocks) {

	case 0:
		cout << "Uncompressed\n";
		break;
	case 1:
		cout << "Zlib\n";
		break;
	case 2:
		cout << "Zstd\n";
		break;
	default:
		cout << "\nERROR: SNP block must be uncompressed or compressed using Zlib or Zstd.\n\n";
		exit(1);
	}
	Layout = (flags >> 2) & 0xf; cout << "Layout: " << Layout << '\n';
	if (Layout < 1U && Layout > 2U) {
		cerr << "\nERROR: BGEN layout flag must be 1 or 2.\n\n";
		exit(1);
	}
	SampleIdentifiers = flags >> 31; cout << "Sample Identifiers Present: ";  SampleIdentifiers == 0 ? cout << "False \n" : cout << "True \n";
	if (SampleIdentifiers < 0 && SampleIdentifiers > 1) {
		cerr << "\nERROR: BGEN sample identifier flag must be 0 or 1\n\n";
		exit(1);
	}

	if (Layout == 2) {
		if (CompressedSNPBlocks < 2) {
			cout << "BGEN file version: " << "v1.2\n";
		}
		else {
			cout << "BGEN file version: " << "v1.3\n";
		}
	}
	else {
		cout << "BGEN file version: " << "v1.1\n";
	}
}






/**********************************************************************************
This function is revised based on the Parse function in BOLT-LMM v2.3 source code
***********************************************************************************/

// This functions reads the sample block of BGEN v1.1, v1.2, and v1.3. Also finds which samples to remove if they have missing values in the pheno file.
void Bgen::processBgenSampleBlock(Bgen bgen, char samplefile[300], bool useSample, unordered_map<string, vector<string>> phenomap, string phenoMissingKey, vector<double> phenodata, vector<double> covdata, int numSelCol, int samSize) {


	int k = 0;
	uint maxLA = 65536;
	char* samID = new char[maxLA + 1];
	unordered_set<int> genoUnMatchID;

	
	std::vector<string> tempID;

	if ((bgen.SampleIdentifiers == 0) || useSample) {

		if (bgen.SampleIdentifiers == 0 && !useSample) {
			cerr << "\nERROR: BGEN file does not contain sample identifiers. A .sample file is required. \n"
				<< "       See https://www.well.ox.ac.uk/~gav/qctool/documentation/sample_file_formats.html for .sample file format. \n\n";
			exit(1);
		}


		std::ifstream fIDMat;
		fIDMat.open(samplefile);

		if (!fIDMat.is_open()) {
			cerr << "\nERROR: Sample file could not be opened." << endl << endl;
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

	if ((bgen.SampleIdentifiers == 1) && !useSample) {

		uint LS1;  fread(&LS1, 4, 1, bgen.fin); // std::cout << "LS1: " << LS1 << std::endl; // LS1 + L_H <= offset
		uint Nrow; fread(&Nrow, 4, 1, bgen.fin); // cout << "Nrow: " << Nrow << " " << std::flush;
		if (Nrow != bgen.Nbgen) {
			cerr << "\nERROR: Nrow = " << Nrow << " does not match Nbgen = " << bgen.Nbgen << "\n\n";
			exit(1);
		}


		int tempIndex = 0;
		for (uint m = 0; m < bgen.Nbgen; m++) {
			ushort LSID; fread(&LSID, 2, 1, bgen.fin); // std::cout << "LSID: " << LSID << " ";
			fread(samID, 1, LSID, bgen.fin); // std::cout << "samID: " << samID << " " << std::endl;

			// IDMatching
			string strtmp(samID);
			int itmp = k;
			if (tempIndex < 5) {
				tempID.push_back(strtmp);
				tempIndex++;
			}

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




	// After IDMatching, resizing phenodata and covdata, and updating samSize;
	phenodata.resize(k);
	covdata.resize(k * (numSelCol + 1));
	samSize = k;

	new_samSize = samSize;
	new_covdata = covdata;
	new_phenodata = phenodata;


	if (samSize == 0) {

		cerr << "\nERROR: Sample size changed from " << samSize + genoUnMatchID.size() << " to " << samSize << ".\n\n";
		if (bgen.SampleIdentifiers == 1 && !useSample) {
			cout << "ID matching was done using the BGEN sample identifier block. \nHere are the first 5 sample identifiers in the block: \n";
			for (int i = 0; i < 5; i++) {
				cout << " " << tempID[i] << "\n";
			}
			cout << "\nRename the sample identifiers in the BGEN file or use a .sample file (--sample) instead. \n\n";

		}

		if (bgen.SampleIdentifiers == 0 || useSample) {
			cout << "Check if sample IDs are consistent between the phenotype file and sample file, or check if (--sampleid-name) is specified correctly. \n\n";
		}

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







/***********************************************************************************
This function is revised based on the Parse function in BOLT-LMM v2.3 source code
************************************************************************************/

// This function reads just the variant block for BGEN files version v1.1, v1.2, and v1.3 and is used to grab the byte where the variant begins.
//    Necesary when there's no bgen index file.
void Bgen::getPositionOfBgenVariant(Bgen bgen, CommandLine cmd) {



	int t = 0;
	int count = 0;
	unsigned int index = 0;
	uint CompressedSNPBlocks = bgen.CompressedSNPBlocks;
	uint Layout = bgen.Layout;
	uint offset = bgen.offset;
	uint Mbgen = bgen.Mbgen;
	uint nSNPS = Mbgen;
	uint Nbgen = bgen.Nbgen;
	uint maxLA = 65536;
	uint maxLB = 65536;
	char* snpID = new char[maxLA + 1];
	char* rsID = new char[maxLA + 1];
	char* chrStr = new char[maxLA + 1];
	char* allele1 = new char[maxLA + 1];
	char* allele0 = new char[maxLB + 1];
	threads = cmd.threads;
	string IDline;

	std::set<std::string> includeVariant;
	std::vector<std::vector<uint>> includeVariantIndex;
	bool checkSNPID = false;
	bool checkRSID = false;
	bool checkInclude = false;
	bool checkExclude = false;


	if (cmd.doFilters) {

		filterVariants = true;
		if (!cmd.includeVariantFile.empty()) {
			checkInclude = true;
			std::ifstream fInclude;
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
				checkSNPID = true;
			}
			else if (IDline == "rsid") {
				cout << "An include snp file was detected... \nIncluding SNPs for analysis based on their rsid... \n";
				checkRSID = true;
			}
			else {
				cerr << "\nERROR: Header name of " << cmd.includeVariantFile << " must be snpid or rsid." << endl << endl;
				exit(1);
			}


			while (fInclude >> vars) {
				includeVariant.insert(vars);
				count++;
			}
			nSNPS = count;
			cout << "Detected " << nSNPS << " variants to be used for analysis... \nAll other variants will be excluded.\n\n" << endl;
			//count = ceil(count / Mbgen_begin.size());


			cout << "Detected " << boost::thread::hardware_concurrency() << " available thread(s)...\n";
			if (nSNPS < threads) {
				threads = nSNPS;
				cout << "Number of variants (" << nSNPS << ") is less than the number of specified threads (" << cmd.threads << ")...\n";
				cout << "Using " << threads << " for multithreading... \n\n";
			}
			else {
				cout << "Using " << threads << " for multithreading... \n\n";
			}
		}



		cout << "Dividing BGEN file into " << threads << " block(s)...\n";
		cout << "Identifying start position of each block...\n";
		vector<uint> endIndex(threads);
		int nBlocks = ceil(nSNPS / threads);
		int lastSNP = nSNPS - 1;
		int index = 0;
		int k = 0;
		Mbgen_begin.resize(threads);
		Mbgen_end.resize(threads);
		bgenVariantPos.resize(threads);
		keepVariants.resize(threads);


		for (int t = 0; t < threads; t++) {
			endIndex[t] = ((t + 1) == threads) ? nSNPS - 1 : floor(((nSNPS / threads) * (t + 1)) - 1);
		}



		FILE* fin = bgen.fin;
		fseek(fin, offset + 4, SEEK_SET);

		// Identifying the start position of each BGEN variant block for multithreading when no bgen index file present.
		for (int snploop = 0; snploop < Mbgen; snploop++) {

			long unsigned int prev = ftell(fin);

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
			if (checkSNPID) {
				if ((checkInclude) && (includeVariant.find(snpID) != includeVariant.end())) {
					keepVariants[k].push_back(snploop);
					if (index == (nBlocks * k)) {
						Mbgen_begin[k] = snploop;
						long int curr = ftell(fin);
						bgenVariantPos[k] = curr - (curr - (prev));
					}
					if (index == endIndex[k]) {
						Mbgen_end[k] = snploop;
						k++;
						if (k == threads) {
							break;
						}
					}
					index++;
				}
			}


			// The length of the rsid
			ushort LR; fread(&LR, 2, 1, fin); // cout << "LR: " << LR << " " << std::flush;
			if (LR > maxLA) {
				maxLA = 2 * LR;
				delete[] rsID;
				char* rsID = new char[maxLA + 1];
			}



			// The rsid
			fread(rsID, 1, LR, fin); rsID[LR] = '\0'; // cout << "rsID: " << string(rsID) << " " << std::flush;
			if (checkRSID) {
				if (checkInclude && includeVariant.find(rsID) != includeVariant.end()) {
					keepVariants[k].push_back(snploop);
					if (index == (nBlocks * k)) {
						Mbgen_begin[k] = snploop;
						long unsigned int curr = ftell(fin);
						bgenVariantPos[k] = curr - (curr - (prev));
					}
					if (index == endIndex[k]) {
						Mbgen_end[k] = snploop;
						k++;
						if (k == threads) {
							break;
						}
					}
					index++;
				}
			}



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
					cerr << "\nERROR: Non-bi-allelic variant found: " << LKnum << " alleles\n\n";
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
			if (Layout == 1) {
				if (CompressedSNPBlocks == 0) {
					fseek(fin, 6 * Nbgen, SEEK_CUR);

				}
				else {
					fseek(fin, zLen, SEEK_CUR);
				}
			}
			if (Layout == 2) {
				if (CompressedSNPBlocks == 0) {
					fseek(fin, zLen, SEEK_CUR);

				}
				else {
					fseek(fin, 4 + zLen - 4, SEEK_CUR);
				}
			}
		}

		if (index != (nSNPS - 1)) {
			cerr << "\nError: There are one or more SNPs in BGEN file with " << IDline << " not in " << cmd.includeVariantFile << "\n\n";
			exit(1);
		}
	}
	else {

		filterVariants = false;

		cout << "Detected " << boost::thread::hardware_concurrency() << " available thread(s)...\n";
		if (Mbgen < threads) {
			threads = Mbgen;
			nSNPS = Mbgen;
			cout << "Number of variants (" << Mbgen << ") is less than the number of specified threads (" << cmd.threads << ")...\n";
			cout << "Using " << threads << " for multithreading... \n\n";
		}
		else {
			nSNPS = Mbgen;
			cout << "Using " << threads << " for multithreading... \n\n";
		}


		cout << "Dividing BGEN file into " << threads << " block(s)..." << endl;
		Mbgen_begin.resize(threads);
		Mbgen_end.resize(threads);
		bgenVariantPos.resize(threads);
		keepVariants.resize(threads);

		for (int t = 0; t < threads; t++) {
			Mbgen_begin[t] = floor((nSNPS / threads) * t);

			if ((t + 1) == (threads)) {
				Mbgen_end[t] = nSNPS - 1;
			}
			else {
				Mbgen_end[t] = floor(((nSNPS / threads) * (t + 1)) - 1);
			}
		}


		FILE* fin = bgen.fin;
		fseek(fin, offset + 4, SEEK_SET);


		// Identifying the start position of each BGEN variant block for multithreading when no bgen index file present.
		for (int snploop = 0; snploop < Mbgen; snploop++) {


			if (snploop == Mbgen_begin[t]) {

				bgenVariantPos[t] = ftell(fin);
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
					cerr << "\nERROR: Non-bi-allelic variant found: " << LKnum << " alleles\n\n";
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
			if (Layout == 2) {
				if (CompressedSNPBlocks > 0) {
					uint zLen;  fread(&zLen, 4, 1, fin);
					fseek(fin, 4 + zLen - 4, SEEK_CUR);

				}
				else {
					uint zLen;  fread(&zLen, 4, 1, fin);
					fseek(fin, zLen, SEEK_CUR);
				}
			}
			else {
				if (CompressedSNPBlocks == 1) {
					uint zLen;  fread(&zLen, 4, 1, fin);
					fseek(fin, zLen, SEEK_CUR);

				}
				else {
					fseek(fin, 6 * Nbgen, SEEK_CUR);
				}
			}
		}
	}
}





/**************************************************************************************************************************************************************************************
The two functions Bgen13GetOneVal and Bgen13GetTwoVals ar being redistributed as part of Plink 2.0 source code
***************************************************************************************************************************************************************************************/

uintptr_t Bgen13GetOneVal(const unsigned char* prob_start, uint64_t prob_offset, uint32_t bit_precision, uintptr_t numer_mask) {
	const uint64_t bit_offset = prob_offset * bit_precision;
	uint64_t relevant_bits;
	// This can read slightly past the end of the buffer.
	memcpy(&relevant_bits, &(prob_start[bit_offset / CHAR_BIT]), sizeof(int64_t));
	return (relevant_bits >> (bit_offset % CHAR_BIT)) & numer_mask;
}


void Bgen13GetTwoVals(const unsigned char* prob_start, uint64_t prob_offset, uint32_t bit_precision, uintptr_t numer_mask, uintptr_t* first_val_ptr, uintptr_t* second_val_ptr) {
	const uint64_t bit_offset = prob_offset * bit_precision;
	uint64_t relevant_bits;
	// This can read slightly past the end of the buffer.
	// Note that with bit_precision=29 and variable ploidy,
	// (bit_offset % CHAR_BIT) == 7 is possible, so we may only get 57 bits when
	// we need 58; thus we don't support 29-31 bits for now.
	memcpy(&relevant_bits, &(prob_start[bit_offset / CHAR_BIT]), sizeof(int64_t));
	relevant_bits = relevant_bits >> (bit_offset % CHAR_BIT);
	*first_val_ptr = relevant_bits & numer_mask;
	*second_val_ptr = (relevant_bits >> bit_precision) & numer_mask;
}




/*************************************************************************************************************************
This function contains code that has been revised based on the Parse function in BOLT-LMM v2.3 source code and Plink 2.00
**************************************************************************************************************************/

void BgenParallelGWAS(int begin, int end, long int byte, vector<uint> keepVariants, char genobgen[300], bool filterVariants, int thread_num, Bgen test) {

	auto start_time = std::chrono::high_resolution_clock::now();
	std::string output = test.outFile + "_bin_" + std::to_string(thread_num) + ".tmp";
	std::ofstream results(output, std::ofstream::binary);
	std::ostringstream oss;


	uint Nbgen  = test.Nbgen;
	uint Layout = test.Layout;
	uint CompressedSNPBlocks = test.CompressedSNPBlocks;

	uint maxLA = 65536;
	uint maxLB = 65536;
	char* snpID   = new char[maxLA + 1];
	char* rsID    = new char[maxLA + 1];
	char* chrStr  = new char[maxLA + 1];
	char* allele1 = new char[maxLA + 1];
	char* allele0 = new char[maxLB + 1];
	string physpos_tmp;
	vector <uchar> zBuf;
	vector <uchar> shortBuf;


	int Sq     = test.Sq;
	int Sq1    = Sq + 1;
	int intSq  = test.numIntSelCol;
	int intSq1 = intSq + 1;
	int expSq        = test.numExpSelCol;
	int phenoType    = test.phenoTyp;
	int robust       = test.robust;
	int samSize      = test.new_samSize;
	int numSelCol    = test.numSelCol;
	int stream_snps  = test.stream_snps;
	double MAF       = test.maf;
	double sigma2    = test.sigma2;
	double* resid    = test.resid;
	double* XinvXTX  = test.XinvXTX;
	double* covX = test.covX;
	vector <double> miu = test.miu;
	vector<long int> include_idx = test.include_idx;

	uint16_t* shortBuf1; 
	uchar* zBuf1;
	if (Layout == 1) {
		shortBuf1 = (uint16_t*)malloc(6 * Nbgen);
		zBuf1 = (unsigned char*)malloc(6 * Nbgen);
	}

	vector <double> ZGSvec(samSize   * (1 + Sq) * stream_snps);
	vector <double> ZGSR2vec(samSize * (1 + Sq) * stream_snps);
	vector <double> WZGSvec(samSize  * (1 + Sq) * stream_snps);
	double* WZGS = &WZGSvec[0];
	
	

	vector <string> geno_snpid(stream_snps);
	struct libdeflate_decompressor* decompressor = libdeflate_alloc_decompressor();

	FILE* fin3;
	fin3 = fopen(genobgen, "rb");
	fseek(fin3, byte, SEEK_SET);

	
	
	int snploop = begin;
	int variant_index = 0;
	int keepIndex = 0;
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

			maxLA = 65536;
			maxLB = 65536;

			// Number of individuals. Only present when Layout == 1
			if (Layout == 1) {
				uint Nrow; fread(&Nrow, 4, 1, fin3); 
				if (Nrow != Nbgen) {
					cerr << "\nERROR: Nrow = " << Nrow << " does not match Nbgen = " << Nbgen << "\n\n";
					exit(1);
				}
			}

			// The variant identifier
			ushort LS; fread(&LS, 2, 1, fin3); 
			if (LS > maxLA) {
				maxLA = 2 * LS;
				delete[] snpID;
				snpID = new char[maxLA + 1];
			}
			fread(snpID, 1, LS, fin3); snpID[LS] = '\0';

			// The rsid
			ushort LR; fread(&LR, 2, 1, fin3);
			if (LR > maxLA) {
				maxLA = 2 * LR;
				delete[] rsID;
				rsID = new char[maxLA + 1];
			}
			fread(rsID, 1, LR, fin3); rsID[LR] = '\0'; 
		
			// The chromosome
			ushort LC; fread(&LC, 2, 1, fin3);
			fread(chrStr, 1, LC, fin3); chrStr[LC] = '\0';

			// The variant position
			uint physpos; fread(&physpos, 4, 1, fin3); 
			physpos_tmp = std::to_string(physpos);

			// The number of alleles if Layout == 2. If Layout == 1, this value is assumed to be 2 and not stored
			if (Layout == 2) {
				ushort LKnum; fread(&LKnum, 2, 1, fin3); // this is for Layout = 2, Lnum = 2 is Layout = 1
				if (LKnum != 2) {
					cerr << "\nERROR: " << snpID << " is a non-bi-allelic variant with " << LKnum << " alleles \n\n";
					exit(1);
				}
			}

			// The first allele
			ushort LA;  fread(&LA, 4, 1, fin3);
			if (LA > maxLA) {
				maxLA = 2 * LA;
				delete[] allele1;
				allele1 = new char[maxLA + 1];
			}
			fread(allele1, 1, LA, fin3); allele1[LA] = '\0';

			// The second allele
			ushort LB; fread(&LB, 4, 1, fin3); 
			if (LB > maxLB) {
				maxLB = 2 * LB;
				delete[] allele0;
				allele0 = new char[maxLB + 1];
			}
			fread(allele0, 1, LB, fin3); allele0[LB] = '\0';





			if (Layout == 1) {

				uLongf destLen1 = 6 * Nbgen;

				if (CompressedSNPBlocks == 1) {
					uint zLen; fread(&zLen, 4, 1, fin3);
					fread(zBuf1, 1, zLen, fin3);
					if (libdeflate_zlib_decompress(decompressor, &zBuf1[0], zLen, &shortBuf1[0], destLen1, NULL) != LIBDEFLATE_SUCCESS) {
						cerr << "\nERROR: Decompressing " << snpID << " block failed with libdeflate.\n\n";
						exit(1);
					}
				}
				else {
					fread(zBuf1, 1, destLen1, fin3);
					shortBuf1 = reinterpret_cast<uint16_t*>(zBuf1);
				}

				// read genotype probabilities
				const double scale = 1.0 / 32768;
				int tmp1 = stream_i * Sq1 * samSize;

				int idx_k = 0;
				for (uint i = 0; i < Nbgen; i++) {
					if (include_idx[idx_k] == i) {
						double p11 = shortBuf1[3 * i] * scale;
						double p10 = shortBuf1[3 * i + 1] * scale;
						double p00 = shortBuf1[3 * i + 2] * scale;

						double pTot = p11 + p10 + p00;
						double dosage = (2 * p00 + p10) / pTot;
						
						int tmp2 = idx_k + tmp1;
						AF[stream_i] += dosage;

						if (phenoType == 1) {
							ZGSvec[tmp2] = miu[idx_k] * (1 - miu[idx_k]) * dosage;
						}
						else {
							ZGSvec[tmp2] = dosage;
						}

						idx_k++;
					}
				}

				if ((AF[stream_i] / 2 / samSize) < MAF || (AF[stream_i] / 2 / samSize) > (1 - MAF)) {
					AF[stream_i] = 0;
					continue;
				}
				else {
					AF[stream_i] = AF[stream_i] / 2 / samSize;
				}

				geno_snpid[stream_i] = string(snpID) + "\t" + string(rsID) + "\t" + string(chrStr) + "\t" + physpos_tmp + "\t" + string(allele1) + "\t" + string(allele0);
				
				for (int j = 0; j < Sq; j++) {
					 int tmp3 = samSize * (j + 1) + tmp1;

					 for (uint i = 0; i < samSize; i++) {
						  int tmp4 = i * (numSelCol + 1);
						  ZGSvec[tmp3 + i] = covX[tmp4 + j + 1] * ZGSvec[tmp1 + i]; // here we save ZGS in column wise
					 }

				}

			} // end of reading genotype data when Layout = 1

			if (Layout == 2) {
				uint zLen; fread(&zLen, 4, 1, fin3);


				if (filterVariants && keepVariants[keepIndex]+1 != snploop) {
					if (CompressedSNPBlocks == 0) {
						fseek(fin3, zLen, SEEK_CUR);
					}
					else {
						fseek(fin3, 4 + zLen - 4, SEEK_CUR);
					}
					continue;
				}

				uint DLen;
				uchar* bufAt;
				if (CompressedSNPBlocks == 1) {
					zBuf.resize(zLen - 4);
					fread(&DLen, 4, 1, fin3);
					fread(&zBuf[0], 1, zLen - 4, fin3);
					shortBuf.resize(DLen);
					uLongf destLen = DLen;

					if (libdeflate_zlib_decompress(decompressor, &zBuf[0], zLen - 4, &shortBuf[0], destLen, NULL) != LIBDEFLATE_SUCCESS) {
						cerr << "\nERROR: Decompressing " << snpID << " block failed\n\n";
						exit(1);
					}
					bufAt = &shortBuf[0];
				}
				else if (CompressedSNPBlocks == 2) {
					zBuf.resize(zLen - 4);
					fread(&DLen, 4, 1, fin3);
					fread(&zBuf[0], 1, zLen - 4, fin3);
					shortBuf.resize(DLen);
					uLongf destLen = DLen;

					size_t ret = ZSTD_decompress(&shortBuf[0], destLen, &zBuf[0], zLen - 4);
					if (ret > destLen) {
						if (ZSTD_isError(ret)) {
							cout << "ZSTD ERROR: " << ZSTD_getErrorName(ret);
						}
					}
					bufAt = &shortBuf[0];
				}
				else {
					zBuf.resize(zLen);
					fread(&zBuf[0], 1, zLen, fin3);
					bufAt = &zBuf[0];
				}


				uint32_t N;
				memcpy(&N, bufAt, sizeof(int32_t));
				if (N != Nbgen) {
					cerr << "\nERROR: " << snpID << " number of samples (" << N << ") is different than header block N\n\n";
					exit(1);
				}

				uint16_t K;
				memcpy(&K, &(bufAt[4]), sizeof(int16_t));
				if (K != 2U) {
					cout << "\nERROR: There are SNP(s) with more than 2 alleles (non-bi-allelic). . Currently unsupported. \n\n";
					exit(1);
				}
			
				const uint32_t min_ploidy = bufAt[6];
				if (min_ploidy > 2U) {
					cerr << "\nERROR: " << snpID << " has minimum ploidy " << min_ploidy << ". Currently unsupported. \n\n";
					exit(1);
				}

				const uint32_t max_ploidy = bufAt[7];
				if (max_ploidy > 2U) {
					cerr << "\nERROR: " << snpID << " has maximum ploidy " << max_ploidy << ". Currently unsupported. \n\n";
					exit(1);
				}

				const unsigned char* missing_and_ploidy_iter = &(bufAt[8]);
				const unsigned char* probs_start = &(bufAt[10 + N]);
				const uint32_t is_phased = probs_start[-2];
				const uint32_t bit_precision = probs_start[-1];
				const uintptr_t numer_mask = (1U << bit_precision) - 1;

				if (is_phased > 1U) {
					cerr << "\nERROR: " << snpID << " has phased value of " << is_phased << ". This must be 0 or 1. \n"
						<< "       See https://www.well.ox.ac.uk/~gav/bgen_format/spec/latest.html for more details. \n\n";
					exit(1);
				}

				int ploidy_sum = 0;
				int idx_to_sum = 0;
				int tmp1 = stream_i * Sq1 * samSize;
				int idx_k = 0;
				uintptr_t prob_offset = 0;
				if (!is_phased) {
					for (int i = 0; i < N; i++) {
						const uint32_t missing_and_ploidy = *missing_and_ploidy_iter++;

						uintptr_t numer_aa;
						uintptr_t numer_ab;

						switch (missing_and_ploidy) {
						case 1:
							Bgen13GetOneVal(probs_start, prob_offset, bit_precision, numer_mask);
							prob_offset++;
							break;
						case 2:
							Bgen13GetTwoVals(probs_start, prob_offset, bit_precision, numer_mask, &numer_aa, &numer_ab);
							prob_offset += 2;
							break;
						default:
							cerr << "\nERROR: " << snpID << " contains ploidy value " << missing_and_ploidy << ". Currently unsupported.\n\n";
							exit(1);
						}

						if (include_idx[idx_k] == i) {
							ploidy_sum += missing_and_ploidy;
							double p11 = numer_aa / double(1.0 * (numer_mask));
							double p10 = numer_ab / double(1.0 * (numer_mask));
							double dosage;

							dosage = 2 * (1 - p11 - p10) + p10;


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

				} else {
					for (int i = 0; i < N; i++) {

						const uint32_t missing_and_ploidy = *missing_and_ploidy_iter++;

					    uintptr_t numer_aa;
						uintptr_t numer_ab;

						switch (missing_and_ploidy) {
						case 1:
							Bgen13GetOneVal(probs_start, prob_offset, bit_precision, numer_mask);
							prob_offset++;
							break;
						case 2:
							Bgen13GetTwoVals(probs_start, prob_offset, bit_precision, numer_mask, &numer_aa, &numer_ab);
							prob_offset += 2;
							break;
						default:
							cerr << "\nERROR: " << snpID << " contains ploidy " << missing_and_ploidy << ". Currently unsupported.\n\n";
							exit(1);
						}

						if (include_idx[idx_k] == i) {
							ploidy_sum += missing_and_ploidy;
							double p11 = numer_aa / double(1.0 * (numer_mask));
							double p10 = numer_ab / double(1.0 * (numer_mask));
							double dosage;

								
							if (missing_and_ploidy == 1U) {
								dosage = 1 - p11;
							}
							else {
								dosage = 2 - (p11 + p10);
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
				}
				
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
			keepIndex++;
		} // end of stream_i

		if (snploop == end + 1 & stream_i == 0) {
			break;
		}

		/***************************************************************/
		//	genodata and envirment data
		double* ZGS   = &ZGSvec[0];
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
			matAdd(WZGS, ZGSR2, samSize* ZGS_col, -1);

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


		double*  betaM      = new double[stream_snps];
		double*  VarbetaM   = new double[stream_snps];
		double** betaInt    = new double* [stream_snps];
		double** VarbetaInt = new double* [stream_snps];
		double*  PvalM      = new double[stream_snps];
		double*  PvalInt    = new double[stream_snps];
		double*  PvalJoint  = new double[stream_snps];
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

			
				//Old code
				//for (int ind1 = 0; ind1 < Sq; ind1++) {
				//	 for (int ind2 = 0; ind2 < Sq; ind2++) {
				//		  // transpose(Snew2) * Snew2
				//		  S2TransS2[ind1 * Sq + ind2] = ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col + ind2 + 1] - ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col] * ZGStZGS[tmp1 + ind2 + 1] / ZGStZGS[tmp1];
	
 			    //   }
				//	 //transpose(Snew2) * resid
				//	 S2TransR[ind1] = ZGStR[i * Sq1 + ind1 + 1] - ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col] * ZGStR[i * Sq1] / ZGStZGS[tmp1];
				//}
				

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
				subMatrix(ZGStR, expZGStR, expSq, 1, expSq, expSq, (i* Sq1) + intSq1);
				// ZGStR[i * Sq1]
				double* int1ZGStR = new double[intSq1];
				subMatrix(ZGStR, int1ZGStR, intSq1, 1, expSq + intSq, intSq1, (i* Sq1));
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
				matAdd(expS2TransS2, S2TransS2right, expSq* expSq, -1.0);


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
				matAdd(expS2DS2, S2DS2third,  expSq * expSq, -1.0);
				matAdd(expS2DS2, S2DS2forth,  expSq * expSq, 1.0);

				// Old code
				//for (int ind1 = 0; ind1 < Sq; ind1++) {
				//	for (int ind2 = 0; ind2 < Sq; ind2++) {
				//		// transpose(Snew2) * Snew2
				//		S2TransS2[ind1 * Sq + ind2] = ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col + ind2 + 1] - ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col] * ZGStZGS[tmp1 + ind2 + 1] / ZGStZGS[tmp1];
				//		// transpose(Snew2) * D * Snew2
				//	S2DS2[ind1 * Sq + ind2] = ZGSR2tZGS[tmp1 + (ind1 + 1) * ZGS_col + ind2 + 1] - ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col] * ZGSR2tZGS[tmp1 + ind2 + 1] / ZGStZGS[tmp1] - ZGSR2tZGS[tmp1 + (ind1 + 1) * ZGS_col] * ZGStZGS[tmp1 + ind2 + 1] / ZGStZGS[tmp1] + ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col] * ZGSR2tZGS[tmp1] * ZGStZGS[tmp1 + ind2 + 1] / ZGStZGS[tmp1] / ZGStZGS[tmp1];
				//	}
				//	//transpose(Snew2) * resid
				//	S2TransR[ind1] = ZGStR[i * Sq1 + ind1 + 1] - ZGStZGS[tmp1 + (ind1 + 1) * ZGS_col] * ZGStR[i * Sq1] / ZGStZGS[tmp1];
				//}
				
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

	if (Layout == 1 && CompressedSNPBlocks == 1) {
		free(shortBuf1);
		free(zBuf1);
	}
	else {
		free(zBuf1);
	}

	// Close files
	results.close();
	fclose(fin3);

	auto end_time = std::chrono::high_resolution_clock::now();
	cout << "Thread " << thread_num << " finished in ";
	printExecutionTime(start_time, end_time);
}






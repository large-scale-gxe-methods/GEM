#include "declars.h"
#include "ReadPGEN.h"
#include "../thirdparty/plink-2.0/plink2_bits.h"
#include "../thirdparty/plink-2.0/plink2_base.h"
#include "../thirdparty/plink-2.0/pgenlib_misc.h"
#include "../thirdparty/plink-2.0/pgenlib_read.h"
#include "../thirdparty/plink-2.0/pgenlib_ffi_support.h"




void Pgen::processPgenHeader(string pgenFile) {

	cout << "General information of PGEN file. \n";

	const char* geno_filename = pgenFile.c_str();

	plink2::PgenFileInfo _info_ptr;
	plink2::PreinitPgfi(&_info_ptr);
	plink2::PgenHeaderCtrl header_ctrl;

	uintptr_t pgfi_alloc_cacheline_ct;
	char errstr_buf[plink2::kPglErrstrBufBlen];
	if (PgfiInitPhase1(geno_filename, UINT32_MAX, UINT32_MAX, 0, &header_ctrl, &_info_ptr, &pgfi_alloc_cacheline_ct, errstr_buf) != plink2::kPglRetSuccess) {
		throw std::runtime_error(errstr_buf);
	}
	assert((header_ctrl & 0x30) == 0); // no alt allele counts
	assert((header_ctrl & 0xc0) != 0xc0); // no explicit nonref_flags

	const uint32_t raw_variant_tmp = _info_ptr.raw_variant_ct;
	const uint32_t raw_sample_tmp  = _info_ptr.raw_sample_ct;
	const uint32_t file_sample_ct = _info_ptr.raw_sample_ct;

	if (raw_variant_tmp == 0) {
		cerr << "\nERROR: No variants in the .pgen file.\n\n";
		exit(1);
	}
	if (raw_sample_tmp == 0) {
		cerr << "\nERROR: No samples in the .pgen file.\n\n";
	}

	raw_sample_ct = raw_sample_tmp;
	raw_variant_ct = raw_variant_tmp;


	unsigned char* pgfi_alloc = nullptr;
	if (plink2::cachealigned_malloc(pgfi_alloc_cacheline_ct * plink2::kCacheline, &pgfi_alloc)) {
		cerr << "Out of memory" << endl;
		exit(1);
	}

	uint32_t max_vrec_width;
	uintptr_t pgr_alloc_cacheline_ct;
	if (PgfiInitPhase2(header_ctrl, 1, 0, 0, 0, raw_variant_ct, &max_vrec_width, &_info_ptr, pgfi_alloc, &pgr_alloc_cacheline_ct, errstr_buf)) {
		if (pgfi_alloc && (!_info_ptr.vrtypes)) {
			plink2::aligned_free(pgfi_alloc);
		}
		throw std::runtime_error(errstr_buf);
	}

	plink2::PgenReader _state_ptr;
	plink2::PreinitPgr(&_state_ptr);
	plink2::PgrSetFreadBuf(nullptr, &_state_ptr);
	const uintptr_t pgr_alloc_main_byte_ct = pgr_alloc_cacheline_ct * plink2::kCacheline;
	const uintptr_t sample_subset_byte_ct = plink2::DivUp(file_sample_ct, plink2::kBitsPerVec) * plink2::kBytesPerVec;
	const uintptr_t cumulative_popcounts_byte_ct = plink2::DivUp(file_sample_ct, plink2::kBitsPerWord * plink2::kInt32PerVec) * plink2::kBytesPerVec;
	const uintptr_t genovec_byte_ct = plink2::DivUp(file_sample_ct, plink2::kNypsPerVec) * plink2::kBytesPerVec;
	const uintptr_t dosage_main_byte_ct = plink2::DivUp(file_sample_ct, (2 * plink2::kInt32PerVec)) * plink2::kBytesPerVec;
	uintptr_t multiallelic_hc_byte_ct = 0;

	unsigned char* pgr_alloc;
	if (plink2::cachealigned_malloc(pgr_alloc_main_byte_ct + (2 * plink2::kPglNypTransposeBatch + 5) * sample_subset_byte_ct + cumulative_popcounts_byte_ct + (1 + plink2::kPglNypTransposeBatch) * genovec_byte_ct + multiallelic_hc_byte_ct + dosage_main_byte_ct + plink2::kPglBitTransposeBufbytes + 4 * (plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 8), &pgr_alloc)) {
		cerr << "Out of memory" << endl;
		exit(1);
	}

	plink2::PglErr reterr = PgrInit(geno_filename, max_vrec_width, &_info_ptr, &_state_ptr, pgr_alloc);
	if (reterr) {
		throw std::runtime_error("Out of memory.");
	}

	uint32_t max_allele_ct = _info_ptr.max_allele_ct;
	if (max_allele_ct != 2) {
		cerr << "\nERROR: There are non-biallelic variants in the .pgen file.\n\n";
		exit(1);
	}

	cout << "Number of variants: " << raw_variant_ct << '\n'; 
	cout << "Number of samples: " << raw_sample_ct << '\n';
}



// This functions reads the sample block of BGEN v1.1, v1.2, and v1.3. Also finds which samples to remove if they have missing values in the pheno file.
void Pgen::processPsam(Pgen pgen, string psamFile, unordered_map<string, vector<string>> phenomap, string phenoMissingKey, vector<double> phenodata, vector<double> covdata, int numSelCol, int samSize) {

	 unordered_set<int> genoUnMatchID;

	 std::ifstream fIDMat;
	 fIDMat.open(psamFile);

	 if (!fIDMat.is_open()) {
		 cerr << "\nERROR: .psam file could not be opened.\n\n";
	 	 exit(1);
	 }

	 string IDline;
	 getline(fIDMat, IDline);

	uint nSamples = 0;
	while (getline(fIDMat, IDline)) {
		nSamples++;
	}
	if (nSamples != pgen.raw_sample_ct) {
		cerr << "\nERROR: Number of sample identifiers in .psam file (" << nSamples << ") does not match the number of samples specified in pgen file (" << pgen.raw_sample_ct << ").\n\n";
		exit(1);
	}
	else {
		fIDMat.clear();
		fIDMat.seekg(0, fIDMat.beg);
		getline(fIDMat, IDline);
	}

	int k = 0;
	for (uint m = 0; m < pgen.raw_sample_ct; m++) {
	     getline(fIDMat, IDline);
	     std::istringstream iss(IDline);
	     string strtmp;
	     iss >> strtmp;

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
	

	// After IDMatching, resizing phenodata and covdata, and updating samSize;
	phenodata.resize(k);
	covdata.resize(k * (numSelCol + 1));
	samSize = k;

	new_samSize = samSize;
	new_covdata = covdata;
	new_phenodata = phenodata;


	if (samSize == 0) {
		cout << "\nERROR: Sample size changed from " << samSize + genoUnMatchID.size() << " to " << samSize << ".\n";
		cout << "\nCheck if sample IDs are consistent between the phenotype file and .psam file, or check if (--sampleid-name) is specified correctly. \n\n";
		exit(1);
	}


	int ii = 0;
	include_idx.resize(samSize);
	for (uint i = 0; i < pgen.raw_sample_ct; i++) {
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


void Pgen::processPvar(Pgen pgen, string pvarFile) {


	std::ifstream fIDMat;
	fIDMat.open(pvarFile);

	if (!fIDMat.is_open()) {
		cerr << "\nERROR: .pvar file could not be opened.\n\n";
		exit(1);
	}

	string IDline;
	getline(fIDMat, IDline);
	std::istringstream iss(IDline);
	string value;
	vector <string> values;

	while (getline(iss, value, '\t')) {
		value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
		values.push_back(value);
	}
	if (values[0] != "#CHROM") {
		cerr << "\nERROR: First line in .pvar file should begin with #CHROM.\n\n";
		exit(1);
	}
	if (values[1] != "POS") {
		cerr << "\nERROR: Second column should be POS in .pvar file.\n\n";
		exit(1);
	}
	if (values[2] != "ID") {
		cerr << "\nERROR: Third column should be ID in .pvar file.\n\n";
		exit(1);
	}
	if (values[3] != "REF") {
		cerr << "\nERROR: Forth column should be REF in .pvar file.\n\n";
		exit(1);
	}
	if (values[4] != "ALT") {
		cerr << "\nERROR: Fifth column should be ALT in .pvar file.\n\n";
		exit(1);
	}

	uint nVariants = 0;
	while (getline(fIDMat, IDline)) {
		nVariants++;
	}
	if (nVariants != pgen.raw_variant_ct) {
		cout << "\nERROR: Number of variants in .pvar file (" << nVariants << ") does not match the number of variants in pgen file (" << pgen.raw_variant_ct << ").\n\n";
		exit(1);
	}

	fIDMat.close();

}






//void Pgen::getPgenVariantPos(Pgen pgen, CommandLine cmd) {
//
//	uint32_t nSNPS = pgen.raw_variant_ct;
//	bool checkSNPID = false;
//	bool checkInclude = false;
//	int count = 0;
//	std::set<std::string> includeVariant;
//
//
//	if (cmd.doFilters) {
//
//		filterVariants = true;
//		if (!cmd.includeVariantFile.empty()) {
//			checkInclude = true;
//			std::ifstream fInclude;
//			string IDline;
//			fInclude.open(cmd.includeVariantFile);
//			if (!fInclude.is_open()) {
//				cerr << "\nERROR: The file (" << cmd.includeVariantFile << ") could not be opened." << endl << endl;
//				exit(1);
//			}
//
//			string vars;
//			getline(fInclude, IDline);
//			std::transform(IDline.begin(), IDline.end(), IDline.begin(), ::tolower);
//			IDline.erase(std::remove(IDline.begin(), IDline.end(), '\r'), IDline.end());
//
//			if (IDline == "id") {
//				cout << "An include snp file was detected... \nIncluding SNPs for analysis based on their id... \n";
//				checkSNPID = true;
//			}
//			else {
//				cerr << "\nERROR: Header name of " << cmd.includeVariantFile << " must be 'id' for PGEN files." << endl << endl;
//				exit(1);
//			}
//
//			while (fInclude >> vars) {
//				includeVariant.insert(vars);
//				count++;
//			}
//			nSNPS = count;
//			cout << "Detected " << nSNPS << " variants to be used for analysis... \nAll other variants will be excluded.\n\n" << endl;
//			//count = ceil(count / Mbgen_begin.size());
//
//
//			cout << "Detected " << boost::thread::hardware_concurrency() << " available thread(s)...\n";
//			if (nSNPS < threads) {
//				threads = nSNPS;
//				cout << "Number of variants (" << nSNPS << ") is less than the number of specified threads (" << threads << ")...\n";
//				cout << "Using " << threads << " for multithreading... \n\n";
//			}
//			else {
//				cout << "Using " << threads << " for multithreading... \n\n";
//			}
//
//		}
//
//		begin.resize(threads);
//		end.resize(threads);
//
//
//		std::ifstream fIDMat;
//		fIDMat.open(cmd.pvarFile);
//
//		string IDline;
//		getline(fIDMat, IDline);
//		std::istringstream iss(IDline);
//		string value;
//		vector <string> values;
//
//		getline(iss, value, '\t');
//
//		int k = 0;
//		for (uint32_t snploop = 0; snploop < raw_variant_ct; snploop++) {
//			while (getline(iss, value, '\t')) {
//				value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
//				values.push_back(value);
//			}
//
//			if (includeVariant.find(rsID) != includeVariant.end()) {
//				pgenVariantPos[k] = snploop;
//				k++;
//			}
//			
//		}
//
//	}
//}




void gemPGEN(uint32_t begin, uint32_t end, string pgenFile, string pvarFile, int thread_num, Pgen test) {

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
		vector <double> ZGSvec(samSize   * (1 + Sq) * stream_snps);
		vector <double> ZGSR2vec(samSize * (1 + Sq) * stream_snps);
		vector <double> WZGSvec(samSize  * (1 + Sq) * stream_snps);
		double* WZGS = &WZGSvec[0];
		std::vector<uint> missingIndex;
		vector <double> AF(stream_snps);
		int ZGS_col = Sq1 * stream_snps;
		std::ifstream fIDMat;
		fIDMat.open(pvarFile);

		string IDline;
		getline(fIDMat, IDline);
		uint32_t skipIndex = 0;
		while (skipIndex != begin) {
			getline(fIDMat, IDline);
			skipIndex++;
		}

		const char* geno_filename = pgenFile.c_str();

		plink2::PgenFileInfo _info_ptr;
		plink2::PreinitPgfi(&_info_ptr);
		plink2::PgenHeaderCtrl header_ctrl;

		uintptr_t pgfi_alloc_cacheline_ct;
		char errstr_buf[plink2::kPglErrstrBufBlen];
		if (PgfiInitPhase1(geno_filename, UINT32_MAX, UINT32_MAX, 0, &header_ctrl, &_info_ptr, &pgfi_alloc_cacheline_ct, errstr_buf) != plink2::kPglRetSuccess) {
			throw std::runtime_error(errstr_buf);
		}
		assert((header_ctrl & 0x30) == 0); // no alt allele counts
		assert((header_ctrl & 0xc0) != 0xc0); // no explicit nonref_flags

		const uint32_t raw_variant_ct = _info_ptr.raw_variant_ct;
		const uint32_t file_sample_ct = _info_ptr.raw_sample_ct;

		unsigned char* pgfi_alloc = nullptr;
		if (plink2::cachealigned_malloc(pgfi_alloc_cacheline_ct * plink2::kCacheline, &pgfi_alloc)) {
			cerr << "Out of memory" << endl;
		}

		uint32_t max_vrec_width;
		uintptr_t pgr_alloc_cacheline_ct;
		if (PgfiInitPhase2(header_ctrl, 1, 0, 0, 0, raw_variant_ct, &max_vrec_width, &_info_ptr, pgfi_alloc, &pgr_alloc_cacheline_ct, errstr_buf)) {
			if (pgfi_alloc && (!_info_ptr.vrtypes)) {
				plink2::aligned_free(pgfi_alloc);
			}
			throw std::runtime_error(errstr_buf);
		}
	
		plink2::PgenVariant _pgv;
		plink2::PgenReader _state_ptr;
		plink2::PreinitPgr(&_state_ptr);
		plink2::PgrSetFreadBuf(nullptr, &_state_ptr);
		const uintptr_t pgr_alloc_main_byte_ct = pgr_alloc_cacheline_ct * plink2::kCacheline;
		const uintptr_t sample_subset_byte_ct = plink2::DivUp(file_sample_ct, plink2::kBitsPerVec) * plink2::kBytesPerVec;
		const uintptr_t cumulative_popcounts_byte_ct = plink2::DivUp(file_sample_ct, plink2::kBitsPerWord * plink2::kInt32PerVec) * plink2::kBytesPerVec;
		const uintptr_t genovec_byte_ct = plink2::DivUp(file_sample_ct, plink2::kNypsPerVec) * plink2::kBytesPerVec;
		const uintptr_t dosage_main_byte_ct = plink2::DivUp(file_sample_ct, (2 * plink2::kInt32PerVec)) * plink2::kBytesPerVec;
		uintptr_t multiallelic_hc_byte_ct = 0;

		unsigned char* pgr_alloc;
		if (plink2::cachealigned_malloc(pgr_alloc_main_byte_ct + (2 * plink2::kPglNypTransposeBatch + 5) * sample_subset_byte_ct + cumulative_popcounts_byte_ct + (1 + plink2::kPglNypTransposeBatch) * genovec_byte_ct + multiallelic_hc_byte_ct + dosage_main_byte_ct + plink2::kPglBitTransposeBufbytes + 4 * (plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 8), &pgr_alloc)) {
			cerr << "Out of memory" << endl;
		}

		plink2::PglErr reterr = PgrInit(geno_filename, max_vrec_width, &_info_ptr, &_state_ptr, pgr_alloc);
		if (reterr) {
			throw std::runtime_error("Out of memory.");
		}

		unsigned char* pgr_alloc_iter = &(pgr_alloc[pgr_alloc_main_byte_ct]);
		uintptr_t* _subset_include_vec = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
		pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
		uintptr_t* _subset_include_interleaved_vec = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
		pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
		_subset_include_interleaved_vec[-1] = 0;

		//uint32_t* _subset_cumulative_popcounts = reinterpret_cast<uint32_t*>(pgr_alloc_iter);
		pgr_alloc_iter = &(pgr_alloc_iter[cumulative_popcounts_byte_ct]);
		_pgv.genovec = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
		pgr_alloc_iter = &(pgr_alloc_iter[genovec_byte_ct]);
		_pgv.patch_01_set = nullptr;
		_pgv.patch_01_vals = nullptr;
		_pgv.patch_10_set = nullptr;
		_pgv.patch_10_vals = nullptr;

		_pgv.phasepresent = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
		pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
		_pgv.phaseinfo = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
		pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
		_pgv.dosage_present = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
		pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
		_pgv.dosage_main = reinterpret_cast<uint16_t*>(pgr_alloc_iter);
		pgr_alloc_iter = &(pgr_alloc_iter[dosage_main_byte_ct]);

		uint32_t _subset_size = file_sample_ct;

		plink2::PgrSampleSubsetIndex _subset_index;
		pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglBitTransposeBufbytes]);

		uint32_t snploop = begin;
		int variant_index = 0;
		int keepIndex = 0;
	    vector<double> buf(file_sample_ct);
		while (snploop <= raw_variant_ct) {

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

				if (snploop >= _info_ptr.raw_variant_ct) {
					char errstr_buf[256];
					sprintf(errstr_buf, "variant_num out of range (%d; must be 1..%u)", snploop + 1, _info_ptr.raw_variant_ct);
					cerr << errstr_buf << "\n";
				}

				uint32_t dosage_ct;
				reterr = plink2::PgrGet1D(_subset_include_vec, _subset_index, _subset_size, snploop, 1, &_state_ptr, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, &dosage_ct);
				plink2::Dosage16ToDoubles(plink2::kGenoDoublePairs, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, _subset_size, dosage_ct, &buf[0]);

				getline(fIDMat, IDline);
				std::istringstream iss(IDline);
				string value;
				vector <string> values;
				while (getline(iss, value, '\t')) {
					values.push_back(value);
				}
				values[4].erase(std::remove(values[4].begin(), values[4].end(), '\r'), values[4].end());
				snploop++;

				int tmp1 = stream_i * Sq1 * samSize;
				int idx_k = 0;
				int nMissing = 0;
				for (uint32_t n = 0; n < file_sample_ct; n++) {
					if (buf[n] == -9.0) {
						if (include_idx[idx_k] == n) {
							nMissing++;
							missingIndex.push_back(idx_k);
							idx_k++;
						}
						continue;
					 }

				     if (include_idx[idx_k] == n) {
						 int tmp2 = idx_k + tmp1;
						 AF[stream_i] += buf[n];

					     if (phenoType == 1) {
				     	     ZGSvec[tmp2] = miu[idx_k] * (1 - miu[idx_k]) * buf[n];
					     }
					     else {
					         ZGSvec[tmp2] = buf[n];
					     }
						 idx_k++;
				     }

				}
				double gmean = AF[stream_i] / double(samSize - nMissing);
				double cur_AF = AF[stream_i] / double(samSize - nMissing) / 2.0;;
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
				geno_snpid[stream_i] = values[2] + "\t" + values[0] + "\t" + values[1] + "\t" + values[3] + "\t" + values[4] + "\t" + std::to_string(samSize - nMissing);

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



		if (_info_ptr.vrtypes) {
			plink2::aligned_free(_info_ptr.vrtypes);
		}
		plink2::PglErr reterr2 = plink2::kPglRetSuccess;
		plink2::CleanupPgfi(&_info_ptr, &reterr2);
		
		reterr2 = plink2::kPglRetSuccess;
		plink2::CleanupPgr(&_state_ptr, &reterr2);
		if (PgrGetFreadBuf(&_state_ptr)) {
			plink2::aligned_free(PgrGetFreadBuf(&_state_ptr));
		}
		_subset_size = 0;

		auto end_time = std::chrono::high_resolution_clock::now();
		cout << "Thread " << thread_num << " finished in ";
		printExecutionTime(start_time, end_time);
}
	


#include "declars.h"
#include "ReadPGEN.h"
#include "../thirdparty/plink-2.0/plink2_bits.h"
#include "../thirdparty/plink-2.0/plink2_base.h"
#include "../thirdparty/plink-2.0/pgenlib_misc.h"
#include "../thirdparty/plink-2.0/pgenlib_read.h"
#include "../thirdparty/plink-2.0/pgenlib_ffi_support.h"




void Pgen::processPgenHeader(string pgenFile) 
{
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



// This functions reads the .psam file
void Pgen::processPsam(Pgen pgen, string psamFile, unordered_map<string, vector<string>> phenomap, string phenoMissingKey, int numSelCol, int samSize) 
{
     unordered_set<int> genoUnMatchID;
     new_phenodata.resize(samSize);
     new_covdata.resize(samSize * (numSelCol+1));
     vector<double> new_covdata_orig(samSize * (numSelCol+1));
     std::ifstream fIDMat;
     fIDMat.open(psamFile);
     if (!fIDMat.is_open()) {
         cerr << "\nERROR: .psam file could not be opened.\n\n";
          exit(1);
     }

     string IDline;
     string value;
     vector <string> values;
     int prev = fIDMat.tellg();
     while (getline(fIDMat, IDline)) {
         std::istringstream iss(IDline);
         while (getline(iss, value, '\t')) {
             value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
             values.push_back(value);
         }
         if (values[0].rfind("##", 0) != 0) {
             break;
         }
         prev = fIDMat.tellg();
         values.clear();
     }

     int iidIndex;
     if ((values[0].compare("#FID") == 0)) {
         std::vector<std::string>::iterator it;
         it = find(values.begin(), values.end(), "IID");
         iidIndex = std::distance(values.begin(), it);

         uint nSamples = 0;
         while (getline(fIDMat, IDline)) {
             nSamples++;
         }

         if (nSamples != pgen.raw_sample_ct) {
             cerr << "\nERROR: Number of sample identifiers in .psam file (" << nSamples << ") does not match the number of samples specified in pgen file (" << pgen.raw_sample_ct << ").\n\n";
             exit(1);
         }
         fIDMat.clear();
         fIDMat.seekg(prev);
         getline(fIDMat, IDline);
         values.clear();

     }
     else if ((values[0].compare("#IID") == 0)) {
         iidIndex = 0;

         uint nSamples = 0;
         while (getline(fIDMat, IDline)) {
             nSamples++;
         }

         if (nSamples != pgen.raw_sample_ct) {
             cerr << "\nERROR: Number of sample identifiers in .psam file (" << nSamples << ") does not match the number of samples specified in pgen file (" << pgen.raw_sample_ct << ").\n\n";
             exit(1);
         }
         fIDMat.clear();
         fIDMat.seekg(prev);
         getline(fIDMat, IDline);
         values.clear();

     }
     else {
         cout << "\nWARNING: No header line with #FID or #IID present in .psam file.\n";
         cout << "Assuming the .psam file column order is FID, IID, PAT, MAT, SEX, PHENO1, as shown here https://www.cog-genomics.org/plink/2.0/formats#psam. \n";
         iidIndex = 1;

         uint nSamples = 1;
         while (getline(fIDMat, IDline)) {
             nSamples++;
         }
         if (nSamples != pgen.raw_sample_ct) {
             cerr << "\nERROR: Number of sample identifiers in .psam file (" << nSamples << ") does not match the number of samples specified in pgen file (" << pgen.raw_sample_ct << ").\n\n";
             exit(1);
         }
         fIDMat.clear();
         fIDMat.seekg(prev);
         values.clear();
     }


    int k = 0;
    for (uint m = 0; m < pgen.raw_sample_ct; m++) {
         getline(fIDMat, IDline);
         std::istringstream iss(IDline);

         string value;
         vector <string> values;
         while (getline(iss, value, '\t')) {
             value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
             values.push_back(value);
         }

         string strtmp = values[iidIndex];
         int itmp = k;
         if (phenomap.find(strtmp) != phenomap.end()) {
             auto tmp_valvec = phenomap[strtmp];
             if (find(tmp_valvec.begin(), tmp_valvec.end(), phenoMissingKey) == tmp_valvec.end() && find(tmp_valvec.begin(), tmp_valvec.end(), "") == tmp_valvec.end()) {
                 sscanf(tmp_valvec[0].c_str(), "%lf", &new_phenodata[k]);
                 new_covdata_orig[k * (numSelCol+1)] = 1.0;
                for (int c = 0; c < numSelCol; c++) {
                    sscanf(tmp_valvec[c + 1].c_str(), "%lf", &new_covdata_orig[k * (numSelCol + 1) + c + 1]);
                }
                sampleID.push_back(strtmp);
                k++;
             }
         }

         if (itmp == k) genoUnMatchID.insert(m);
    }
    fIDMat.close();
    

    // After IDMatching, resizing phenodata and covdata, and updating samSize;
    new_phenodata.resize(k);
    new_covdata_orig.resize(k * (numSelCol+1));
    samSize = k;

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
     

    new_samSize = samSize;
    if (new_samSize<(numSelCol+1) || new_samSize == (numSelCol+1)){
        cout << "\nERROR: The sample size should be greater than the number of predictors!" <<endl;
        exit(1);
    }
    MatrixXd matcovX (samSize,(numSelCol+1));
    for (int i=0; i<samSize; i++){    
        for (int j=0; j<(numSelCol+1); j++) {
          matcovX(i,j) =new_covdata_orig [i * (numSelCol+1) +j];
        }
    }
    Eigen::HouseholderQR<MatrixXd> qr;
    qr.compute(matcovX);
    Eigen::MatrixXd R = qr.matrixQR();
    int colR=R.cols();
    int rowR=R.rows();
    int diagR_size;
    if (colR > rowR){
        diagR_size = rowR;
    }
    else{
        diagR_size = colR;
    }
    VectorXd diagR (diagR_size);
    for (int i=0; i<diagR_size; i++){
        diagR(i)=R(i,i);
    }
    double sqrtEps =sqrt(std::numeric_limits<double>::epsilon());
    double maxdiag = *std::max_element( diagR.begin(), diagR.end() ) ;
    double colinear_cut = abs(maxdiag * sqrtEps);
    for (int i=0; i<diagR_size; i++){
        if (abs(diagR(i)) < colinear_cut){
            excludeCol.push_back(i);
        }
    }
    matcovX.resize(0,0);
    R.resize(0,0);

    int NumExcludeCol = excludeCol.size();

    if (excludeCol.size()>0){    
        vector <int> remove_colinear;
        for (int i=0; i<excludeCol.size(); i++){
            for (int j=0; j<samSize; j++) {
                remove_colinear.push_back(j * (numSelCol+1) + excludeCol[i]);
            }
        }
        numSelCol=numSelCol- excludeCol.size();
        new_covdata.resize(samSize * (numSelCol+1));
        vector<double> temp;
        for (int i=0; i<new_covdata_orig.size(); i++)
        {
            if (std::find(remove_colinear.begin(), remove_colinear.end(), i) == remove_colinear.end())
            {
                new_covdata.push_back(new_covdata_orig[i]);
                temp.push_back(new_covdata_orig[i]);
                
            }
        }
        new_covdata = temp;
    } 
    else {
            new_covdata.resize(samSize * (numSelCol+1));
            new_covdata = new_covdata_orig;
    }
}


void Pgen::processPvar(Pgen pgen, string pvarFile) 
{
    std::ifstream fIDMat;
    fIDMat.open(pvarFile);

    if (!fIDMat.is_open()) {
        cerr << "\nERROR: .pvar file could not be opened.\n\n";
        exit(1);
    }

    string IDline;
    string value;
    vector <string> values;
    while (getline(fIDMat, IDline)) {
        std::istringstream iss(IDline);
        while (getline(iss, value, '\t')) {
            value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
            values.push_back(value);
        }
        if (values[0].rfind("##", 0) != 0) {
            break;
        }
        values.clear();
    }

    if (values[0].compare("#CHROM") == 0) {
        std::vector<std::string>::iterator it;
        
        it = find(values.begin(), values.end(), "ID");
        if (it != values.end()) {
            int idx = std::distance(values.begin(), it);
            pvarIndex.push_back(idx);
        }
        else {
            cerr << "\nERROR: Cannot find ID column header name in .pvar file when #CHROM header line is present.\n\n";
        }

        pvarIndex.push_back(0);

        it = find(values.begin(), values.end(), "POS");
        if (it != values.end()) {
            int idx = std::distance(values.begin(), it);
            pvarIndex.push_back(idx);
        }
        else {
            cerr << "\nERROR: Cannot find POS column header name in .pvar file when #CHROM header line is present.\n\n";
        }

        it = find(values.begin(), values.end(), "REF");
        if (it != values.end()) {
            int idx = std::distance(values.begin(), it);
            pvarIndex.push_back(idx);
        }
        else {
            cerr << "\nERROR: Cannot find REF column header name in .pvar file when #CHROM header line is present.\n\n";
        }

        it = find(values.begin(), values.end(), "ALT");
        if (it != values.end()) {
            int idx = std::distance(values.begin(), it);
            pvarIndex.push_back(idx);
        }
        else {
            cerr << "\nERROR: Cannot find ALT column header name in .pvar file when #CHROM header line is present.\n\n";
        }

        uint nVariants = 0;
        while (getline(fIDMat, IDline)) {
            nVariants++;
        }
        if (nVariants != pgen.raw_variant_ct) {
            cout << "\nERROR: Number of variants in .pvar file (" << nVariants << ") does not match the number of variants in pgen file (" << pgen.raw_variant_ct << ").\n\n";
            exit(1);
        }
        pvarLast = 4;

    }
    else {
        if (values.size() == 5) {
            cout << "\nWARNING: No header line with #CHROM present in .pvar file. Number of columns in .pvar file is 5.\n";
            cout << "Assuming the .pvar file column order is CHROM, ID, POS, ALT, REF as shown here https://www.cog-genomics.org/plink/2.0/formats#pvar. \n";
            pvarIndex.push_back(1);
            pvarIndex.push_back(0);
            pvarIndex.push_back(2);
            pvarIndex.push_back(4);
            pvarIndex.push_back(3);
            pvarLast = 4;
        }
        else if (values.size() == 6) {
            cout << "\nWARNING: No header line with #CHROM present in .pvar file. Number of columns in .pvar file is 6.\n";
            cout << "Assuming the .pvar file column order is CHROM, ID, CM, POS, ALT, REF as shown here https://www.cog-genomics.org/plink/2.0/formats#pvar. \n";
            pvarIndex.push_back(1);
            pvarIndex.push_back(0);
            pvarIndex.push_back(3);
            pvarIndex.push_back(5);
            pvarIndex.push_back(4);
            pvarLast = 5;
        }
        else {
            cerr << "\nERROR: Number of columns is not 5 or 6  in .pvar file when the #CHROM header line is not present ( https://www.cog-genomics.org/plink/2.0/formats#pvar ).\n\n";
            exit(1);
        }

        uint nVariants = 1;
        while (getline(fIDMat, IDline)) {
            nVariants++;
        }
        if (nVariants != pgen.raw_variant_ct) {
            cout << "\nERROR: Number of variants in .pvar file (" << nVariants << ") does not match the number of variants in pgen file (" << pgen.raw_variant_ct << ").\n\n";
            exit(1);
        }
    }

    fIDMat.close();
}



void Pgen::getPgenVariantPos(Pgen pgen, CommandLine cmd) 
{
    uint32_t nSNPS = pgen.raw_variant_ct;
    int count = 0;
    std::set<std::string> includeVariant;
    threads = cmd.threads;
    vector<int> pvarIndex = pgen.pvarIndex;

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
                cerr << "\nERROR: Header name of " << cmd.includeVariantFile << " must be 'snpid' for PGEN files.\n\n";
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
            cout << "Detected " << nSNPS << " variants to be used for analysis... \nAll other variants will be excluded.\n\n\n";

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
        fIDMat.open(cmd.pvarFile);

        string IDline;
        string value;
        vector <string> values;
        int prev = fIDMat.tellg();
        int pi;
        while (getline(fIDMat, IDline)) {
            std::istringstream iss(IDline);
            while (getline(iss, value, '\t')) {
                value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                values.push_back(value);
            }
            if (values[0].rfind("##", 0) != 0) {
                break;
            }
            prev = fIDMat.tellg();
            values.clear();
        }
        if (values[0].compare("#CHROM") == 0) {
            std::vector<std::string>::iterator it;
            it = find(values.begin(), values.end(), "ID");
            pi = std::distance(values.begin(), it);
        }
        else {
            fIDMat.seekg(prev);
            pi = 1;
        }


        uint32_t k = 0;
        long long unsigned int pvalIndex = 0;
        while (getline(fIDMat, IDline)) {
            std::istringstream iss(IDline);
            string value;
            vector <string> values;

            while (getline(iss, value, '\t')) {
                value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                values.push_back(value);
            }
            if (includeVariant.find(values[pi]) != includeVariant.end()) {
                pgenVariantPos.push_back(pvalIndex);
                k++;

                values.clear();
            }
            pvalIndex++;
        }

        if (k != nSNPS) {
            cerr << "\nERROR: There are one or more SNPs in " << cmd.includeVariantFile << " file that is not in the .pvar file.\n\n";
            exit(1);
        }

        fIDMat.close();
    }
    else {
         threads = cmd.threads;
         if (pgen.raw_variant_ct < threads) {
             cout << "Number of variants (" << pgen.raw_variant_ct << ") is less than the number of specified threads (" << threads << ")...\n";
             threads = pgen.raw_variant_ct;
              cout << "Using " << threads << " thread(s) instead... \n\n";
         }
         begin.resize(threads);
         end.resize(threads);

         for (uint32_t t = 0; t < threads; t++) {
              begin[t] = floor((pgen.raw_variant_ct / threads) * t);
              if ((t + 1) == (threads)) {
                  end[t] = pgen.raw_variant_ct - 1;
              }
              else {
                  end[t] = floor(((pgen.raw_variant_ct / threads) * (t + 1)) - 1);
              }
         }
    }

}




void gemPGEN(int thread_num, double sigma2, double* resid, double* XinvXTX, vector<double> miu, BinE binE, Pgen pgen, CommandLine cmd) 
{
    auto start_time = std::chrono::high_resolution_clock::now();
    std::string output = cmd.outFile + "_bin_" + std::to_string(thread_num) + ".tmp";
    std::ofstream results(output, std::ofstream::binary);
    std::ostringstream oss;


    bool filterVariants = pgen.filterVariants;
    string outStyle = cmd.outStyle;
    int phenoType   = pgen.phenoType;
    int stream_snps = cmd.stream_snps;
    int samSize     = pgen.new_samSize;
    int robust      = cmd.robust;
    int intSq1      = pgen.numIntSelCol_new + 1;
    int expSq       = pgen.numExpSelCol_new;
    int expSq1      = expSq+1;
    int Sq1         = intSq1 + expSq;
    int Sq          = Sq1-1;
    int numSelCol1  = pgen.numSelCol_new + Sq1; 
    double MAF      = cmd.MAF;
    double maxMAF   = 1 - MAF;
    double missGenoCutoff = cmd.missGenoRate;
    vector<long int> include_idx = pgen.include_idx;


    int numBinE   = binE.nBinE;
    bool strata   = (numBinE > 0 ) ? true : false;
    int strataLen = binE.strataLen;
    vector<int> stratum_idx = binE.stratum_idx;
    vector<double> binE_AF(stream_snps * strataLen, 0.0), binE_N(stream_snps * strataLen, 0.0);

    int pvarLast = pgen.pvarLast;
    vector<int> pvarIndex = pgen.pvarIndex;
    uint32_t snploop = pgen.begin[thread_num], end = pgen.end[thread_num];
    std::vector<long long unsigned int> pgenPos = pgen.pgenVariantPos;

    int ZGS_col = Sq1 * stream_snps;
    //vector <uint>   missingIndex;
    vector <double> ZGSvec(samSize   * (Sq1) * stream_snps);
    vector <double> ZGSR2vec(samSize * (Sq1) * stream_snps);
    vector <double> WZGSvec(samSize  * (Sq1) * stream_snps);
    vector <double> AF(stream_snps);
    vector <string> geno_snpid(stream_snps);
    double* WZGS = &WZGSvec[0];
    double* covX = &pgen.new_covdata[0];   

    int pvarLength = pgen.pvarIndex.size();
    std::ifstream fIDMat;
    fIDMat.open(cmd.pvarFile);
    string IDline;
    string tmpvalue;
    vector <string> tmpvalues;
    int prev = fIDMat.tellg();
    while (getline(fIDMat, IDline)) {
        std::istringstream iss(IDline);
        while (getline(iss, tmpvalue, '\t')) {
            tmpvalue.erase(std::remove(tmpvalue.begin(), tmpvalue.end(), '\r'), tmpvalue.end());
            tmpvalues.push_back(tmpvalue);
        }
        if (tmpvalues[0].rfind("##", 0) != 0) {
            break;
        }
        prev = fIDMat.tellg();
        tmpvalues.clear();
    }
    if (tmpvalues[0].compare("#CHROM") != 0) {
        fIDMat.seekg(prev);
    }

    uint32_t skipIndex = 0;
    if (!filterVariants) {
        while (skipIndex != snploop) {
            getline(fIDMat, IDline);
            skipIndex++;
        }
    }
    else {
        while (skipIndex != pgenPos[snploop]) {
            getline(fIDMat, IDline);
            skipIndex++;
        }
    }

    const char* geno_filename = cmd.pgenFile.c_str();

    plink2::PgenFileInfo _info_ptr;
    plink2::PreinitPgfi(&_info_ptr);
    plink2::PgenHeaderCtrl header_ctrl;

    uintptr_t pgfi_alloc_cacheline_ct;
    char errstr_buf[plink2::kPglErrstrBufBlen];
    if (PgfiInitPhase1(geno_filename, UINT32_MAX, UINT32_MAX, 0, &header_ctrl, &_info_ptr, &pgfi_alloc_cacheline_ct, errstr_buf) != plink2::kPglRetSuccess) {
        throw std::runtime_error(errstr_buf);
    }

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

    pgr_alloc_iter = &(pgr_alloc_iter[cumulative_popcounts_byte_ct]);
    _pgv.genovec = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[genovec_byte_ct]);

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


    int printStart = 1; int printEnd = expSq1;
    bool printFull = false;
    bool printMeta = false;
    if (expSq == 0) {
        printStart = 0; printEnd = 0;
    }
    else if (cmd.outStyle.compare("meta") == 0) {
        printStart = 0; printEnd = Sq1; printMeta = true;
    } else if (cmd.outStyle.compare("full") == 0) {
        printStart = 0; printEnd = Sq1; printFull = true;
    }
        
    boost::math::chi_squared chisq_dist_M(1);

    int variant_index = 0;
    int keepIndex = 0;
    vector<double> buf(file_sample_ct);
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

            if (snploop >= _info_ptr.raw_variant_ct) {
                char errstr_buf[256];
                sprintf(errstr_buf, "variant_num out of range (%d; must be 1..%u)", snploop + 1, _info_ptr.raw_variant_ct);
                cerr << errstr_buf << "\n";
            }

            uint32_t dosage_ct;
            string value;
            vector <string> values;
            if (!filterVariants) {
                reterr = plink2::PgrGet1D(_subset_include_vec, _subset_index, _subset_size, snploop, 1, &_state_ptr, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, &dosage_ct);
                getline(fIDMat, IDline);
                std::istringstream iss(IDline);
                while (getline(iss, value, '\t')) {
                    values.push_back(value);
                }
            }
            else {
                while (skipIndex != pgenPos[snploop]) {
                    getline(fIDMat, IDline);
                    skipIndex++;
                }
                getline(fIDMat, IDline);
                skipIndex++;
                std::istringstream iss(IDline);
                while (getline(iss, value, '\t')) {
                    values.push_back(value);
                }
                reterr = plink2::PgrGet1D(_subset_include_vec, _subset_index, _subset_size, pgenPos[snploop], 1, &_state_ptr, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, &dosage_ct);
            }
            plink2::Dosage16ToDoubles(plink2::kGenoDoublePairs, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, _subset_size, dosage_ct, &buf[0]);
            snploop++;

            int tmp1 = stream_i * Sq1 * samSize;
            int strata_i = stream_i * strataLen;
            int idx_k = 0;
            int nMissing = 0;
            vector <uint>   missingIndex;
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
                        
                    if (strata) {
                        binE_N[strata_i + stratum_idx[idx_k]]+=1.0;
                        binE_AF[strata_i + stratum_idx[idx_k]]+=buf[n];
                    }

                    idx_k++;
                }
            }

            double gmean = AF[stream_i] / double(samSize - nMissing);
            double cur_AF = AF[stream_i] / double(samSize - nMissing) / 2.0;
            double percMissing = nMissing / (samSize * 1.0);
            if ((cur_AF < MAF || cur_AF > maxMAF) || (percMissing > missGenoCutoff)) {
                AF[stream_i] = 0.0;
                if (strata) {
                    for (int i = 0; i < strataLen; i++) {
                        binE_N[strata_i + i] = 0.0;
                        binE_AF[strata_i + i] = 0.0;
                    }
                }
                variant_index++;
                keepIndex++;
                continue;
            }
            else {
                AF[stream_i] = cur_AF;
            }

            if (strata) { 
                for (int i = 0; i < strataLen; i++) {
                    binE_AF[strata_i + i] = binE_AF[strata_i + i] / binE_N[strata_i + i] / 2.0;
                }
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

            string tmpString = "";
            values[pvarLast].erase(std::remove(values[pvarLast].begin(), values[pvarLast].end(), '\r'), values[pvarLast].end());
            for (int p = 0; p < pvarLength; p++) {
                tmpString = tmpString + values[pvarIndex[p]] + "\t";
            }
            geno_snpid[stream_i] = tmpString + std::to_string(samSize - nMissing);

            for (int j = 0; j < Sq; j++) {
                int tmp3 = samSize * (j + 1) + tmp1;
                for (int i = 0; i < samSize; i++) {
                    int tmp4 = i * numSelCol1;
                    ZGSvec[tmp3 + i] = covX[tmp4 + j + 1] * ZGSvec[tmp1 + i]; // here we save ZGS in column wise
                }
            }
            variant_index++;
            stream_i++;
            keepIndex++;
        }

        if ((snploop == (end + 1)) & (stream_i == 0)) { break; }

        /***************************************************************/
        double* ZGS = &ZGSvec[0];
        double* ZGSR2 = &ZGSR2vec[0];

        // transpose(X) * ZGS
        // it is non-squred matrix, attention that continuous memory is column-major due to Fortran in BLAS.
        // important!!!!
        double* XtransZGS = new double[numSelCol1 * ZGS_col];
        matNmatNprod(covX, ZGS, XtransZGS, numSelCol1, samSize, ZGS_col);

        if (phenoType == 0) {
            matNmatNprod(XinvXTX, XtransZGS, ZGSR2, samSize, numSelCol1, ZGS_col);
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
            matNmatNprod(XinvXTX, XtransZGS, ZGSR2, samSize, numSelCol1, ZGS_col);
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


        double*  betaM      = new double[stream_snps];
        double*  VarbetaM   = new double[stream_snps];
        double*  PvalM      = new double[stream_snps];
        double*  PvalInt    = new double[stream_snps];
        double*  PvalJoint  = new double[stream_snps];
        double** betaAll    = new double* [stream_snps];
        double** VarBetaAll = new double* [stream_snps];
        double*  mbVarbetaM   = nullptr;
        double*  mbPvalM      = nullptr;
        double*  mbPvalInt    = nullptr;
        double*  mbPvalJoint  = nullptr;
        double** invZGStZGS   = nullptr;
        if (robust == 1) {
            invZGStZGS =  new double* [stream_snps];
        }
        if (printMeta || printFull) {
            mbVarbetaM  = new double[stream_snps];          
            mbPvalM     = new double[stream_snps];            
            mbPvalInt   = new double[stream_snps];
            mbPvalJoint = new double[stream_snps] ;
        }

        if (robust == 0) {
            for (int i = 0; i < stream_snps; i++) {

                // betamain
                int tmp1 = i * ZGS_col * Sq1 + i * Sq1;
                betaM[i] = ZGStR[i * Sq1] / ZGStZGS[tmp1];
                VarbetaM[i] = sigma2 / ZGStZGS[tmp1];

                //calculating Marginal P values
                double statM = betaM[i] * betaM[i] / VarbetaM[i];
                PvalM[i] = (isnan(statM) || statM <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_M, statM));

                if (expSq != 0) {
                    boost::math::chi_squared chisq_dist_Int(expSq);
                    boost::math::chi_squared chisq_dist_Joint(expSq1);
                    
                    // ZGStR
                    double* subZGStR = new double[Sq1];
                    subMatrix(ZGStR, subZGStR, Sq1, 1, Sq1, Sq1, i * Sq1);

                    // inv(ZGStZGS)
                    VarBetaAll[i] = new double[Sq1 * Sq1];
                    subMatrix(ZGStZGS, VarBetaAll[i], Sq1, Sq1, ZGS_col, Sq1, tmp1);
                    matInv(VarBetaAll[i], Sq1);

                    betaAll[i] = new double[Sq1];
                    matvecprod(VarBetaAll[i], subZGStR, betaAll[i], Sq1, Sq1);
                    for (int k = 0; k < Sq1 * Sq1; k++) {
                        VarBetaAll[i][k] *= sigma2;
                    }
                    //invVarBetaInt
                    double* invVarbetaint = new double[expSq * expSq];
                    subMatrix(VarBetaAll[i], invVarbetaint, expSq, expSq, Sq1, expSq, Sq1+1);
                    matInv(invVarbetaint, expSq);
                    
                    // StatInt
                    double* Stemp3 = new double[expSq];
                    matvecSprod(invVarbetaint, betaAll[i], Stemp3, expSq, expSq, 1);

                    double statInt = 0.0;
                    for (int j = 1; j < expSq1; j++) {
                        statInt += betaAll[i][j] * Stemp3[j-1];
                    }
                    
                    //calculating Interaction P values
                    PvalInt[i] = (isnan(statInt) || statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int, statInt));

                    double* invA = new double[expSq1 * expSq1];
                    subMatrix(VarBetaAll[i], invA, expSq1, expSq1, Sq1, expSq1, 0);
                    matInv(invA, expSq1);

                    double* Stemp4 = new double[expSq1];
                    matvecprod(invA, betaAll[i], Stemp4, expSq1, expSq1);
                    
                    double statJoint = 0.0;
                    for (int k = 0; k < expSq1; k++) {
                        statJoint += betaAll[i][k] * Stemp4[k];
                    }

                    PvalJoint[i] = (isnan(statJoint) || statJoint <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Joint, statJoint));

                    delete[] subZGStR;
                    delete[] invVarbetaint;
                    delete[] Stemp3;
                    delete[] Stemp4;
                    delete[] invA;
                }
            }
        }
        else if (robust == 1) {
            for (int i = 0; i < stream_snps; i++) {

                //BetaMain
                int tmp1 = i * ZGS_col * Sq1 + i * Sq1;
                betaM[i] = ZGStR[i * Sq1] / ZGStZGS[tmp1];
                VarbetaM[i] = ZGSR2tZGS[tmp1] / (ZGStZGS[tmp1] * ZGStZGS[tmp1]);

                //calculating Marginal P values
                double statM = betaM[i] * betaM[i] / VarbetaM[i];
                if (isnan(statM) || statM <= 0.0) {
                    PvalM[i] = NAN;
                }
                else {
                    PvalM[i] = boost::math::cdf(complement(chisq_dist_M, statM));
                }

                if (expSq != 0) {
                    boost::math::chi_squared chisq_dist_Int(expSq);
                    boost::math::chi_squared chisq_dist_Joint(expSq1);
                    
                    // ZGStR
                    double* subZGStR = new double[Sq1];
                    subMatrix(ZGStR, subZGStR, Sq1, 1, Sq1, Sq1, i * Sq1);

                    // ZGSR2tZGS
                    double* subZGSR2tZGS = new double[Sq1 * Sq1];
                    subMatrix(ZGSR2tZGS, subZGSR2tZGS, Sq1, Sq1, ZGS_col, Sq1, tmp1);

                    // inv(ZGStZGS)
                    invZGStZGS[i] = new double[Sq1 * Sq1];
                    subMatrix(ZGStZGS, invZGStZGS[i], Sq1, Sq1, ZGS_col, Sq1, tmp1);   
                    matInv(invZGStZGS[i], Sq1);

                    betaAll[i]= new double[Sq1];
                    matvecprod(invZGStZGS[i], subZGStR, betaAll[i], Sq1, Sq1);
                    double* ZGSR2tZGSxinvZGStZGS = new double[Sq1 * Sq1];
                    matNmatNprod(subZGSR2tZGS, invZGStZGS[i], ZGSR2tZGSxinvZGStZGS, Sq1, Sq1, Sq1);
                    VarBetaAll[i] = new double[Sq1 * Sq1];
                    matNmatNprod(invZGStZGS[i], ZGSR2tZGSxinvZGStZGS, VarBetaAll[i], Sq1, Sq1, Sq1);


                    // invVarBetaInt
                    double* invVarbetaint = new double[expSq * expSq];
                    subMatrix(VarBetaAll[i], invVarbetaint, expSq, expSq, Sq1, expSq, 1+Sq1);
                    matInv(invVarbetaint, expSq);

                    // StatInt
                    double* Stemp3 = new double[expSq];
                    matvecSprod(invVarbetaint, betaAll[i], Stemp3, expSq, expSq, 1);
                    double statInt = 0.0;
                    for (int j = 1; j < expSq1; j++) {
                        statInt += betaAll[i][j] * Stemp3[j-1];
                    }
            
                    //calculating Interaction P values
                    PvalInt[i] = (isnan(statInt) || statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int, statInt));

                    double* invA = new double[expSq1 * expSq1];
                    subMatrix(VarBetaAll[i], invA, expSq1, expSq1, Sq1, expSq1, 0);
                    matInv(invA, expSq1);

                    double* Stemp4 = new double[expSq1];
                    matvecprod(invA, betaAll[i], Stemp4, expSq1, expSq1);

                    double statJoint = 0.0;
                    for (int k = 0; k < expSq1; k++) {
                        statJoint += betaAll[i][k] * Stemp4[k];
                    }

                    PvalJoint[i] = (isnan(statJoint) || statJoint <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Joint, statJoint));

                    if (printMeta || printFull) {

                        mbVarbetaM[i] = sigma2 / ZGStZGS[tmp1];

                        //calculating model-based Marginal P values
                        double statM = betaM[i] * betaM[i] / mbVarbetaM[i];
                        mbPvalM[i] = (isnan(statM) || statM <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_M, statM));

                        for (int k = 0; k < Sq1 * Sq1; k++) {
                            invZGStZGS[i][k] *= sigma2;
                        }
      
                        double* mbInvVarbetaint = new double[expSq * expSq];
                        subMatrix(invZGStZGS[i], mbInvVarbetaint, expSq, expSq, Sq1, expSq, 1+Sq1);
                        matInv(mbInvVarbetaint, expSq);

                        double* Stemp3 = new double[expSq];
                        matvecSprod(mbInvVarbetaint, betaAll[i], Stemp3, expSq, expSq, 1);

                        double statInt = 0.0;
                        for (int j = 1; j < expSq1; j++) {
                            statInt += betaAll[i][j] * Stemp3[j-1];
                        }
                    
                        //calculating model-based interaction P values
                        mbPvalInt[i] = (isnan(statInt) || statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int, statInt));
 
                        subMatrix(invZGStZGS[i], invA, expSq1, expSq1, Sq1, expSq1, 0);
                        matInv(invA, expSq1);
                        matvecprod(invA, betaAll[i], Stemp4, expSq1, expSq1);
                    
                        double statJoint = 0.0;
                        for (int k = 0; k < expSq1; k++) {
                            statJoint += betaAll[i][k] * Stemp4[k];
                        }
                        
                        //calculating model-based joint P values
                        mbPvalJoint[i] = (isnan(statJoint) || statJoint <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Joint, statJoint));
                    }

                    delete[] subZGStR;
                    delete[] subZGSR2tZGS;
                    delete[] ZGSR2tZGSxinvZGStZGS;
                    delete[] invVarbetaint;
                    delete[] Stemp3;
                    delete[] Stemp4;
                    delete[] invA;
                }

            }
        } // end of if robust == 1


        for (int i = 0; i < stream_snps; i++) {
            oss << geno_snpid[i] << "\t" << AF[i] << "\t";

            int strata_ii = i * strataLen;
            for (int k = 0; k < strataLen; k++) {
                oss << binE_N[strata_ii + k] << "\t" << binE_AF[strata_ii + k] << "\t";
            }
            
            oss << betaM[i] << "\t" << sqrt(VarbetaM[i]) << "\t";

            if ((robust == 1) && (printFull || printMeta)) {
                oss << sqrt(mbVarbetaM[i]) << "\t";
            }

            for (int ii = printStart; ii < printEnd; ii++) {
                 oss << betaAll[i][ii] << "\t";
            }
            for (int ii = printStart; ii < printEnd; ii++) {
                oss << sqrt(VarBetaAll[i][ii * Sq1 + ii]) << "\t";
            }

            for (int ii = printStart; ii < printEnd; ii++) {
                for (int jj = printStart; jj < printEnd; jj++) {
                    if (ii < jj) {
                        oss << VarBetaAll[i][ii * Sq1 + jj] << "\t";
                    }
                }
            }

            if ((robust == 1) && (printFull || printMeta)) {
                for (int ii = printStart; ii < printEnd; ii++) {
                    for (int jj = printStart; jj < printEnd; jj++) {
                        if (ii == jj) {
                            oss << sqrt(invZGStZGS[i][ii*Sq1 + jj]) << "\t";
                        }
                    }
                }

                for (int ii = printStart; ii < printEnd; ii++) {
                    for (int jj = printStart; jj < printEnd; jj++) {
                        if (ii < jj) {
                            oss << invZGStZGS[i][ii*Sq1 + jj] << "\t";
                        }
                    }
                }

            }

            if (expSq != 0) {
                oss << PvalM[i] << "\t" << PvalInt[i] << "\t" << PvalJoint[i];
                if ((robust == 1) && (printMeta || printFull)) {
                    oss << "\t" << mbPvalM[i] << "\t" << mbPvalInt[i] << "\t" << mbPvalJoint[i];
                }
                oss << "\n";
            }
            else {
                oss << PvalM[i] << "\n";
            }
            AF[i] = 0.0;
        }

        if (strata) {       
            std::fill(binE_N.begin(), binE_N.end(), 0.0);
            std::fill(binE_AF.begin(), binE_AF.end(), 0.0);
        }
        delete[] ZGStR;
        delete[] ZGStZGS;
        delete[] ZGSR2tZGS;
        delete[] betaM;
        delete[] VarbetaM;

        if (expSq != 0 ) {
            for (int i = 0; i < stream_snps; i++) {
                delete[] betaAll[i];
                delete[] VarBetaAll[i];
            }
            if (robust == 1) {
                for (int i = 0; i < stream_snps; i++) {
                    delete[] invZGStZGS[i];
                }
                if (printMeta || printFull) {
                    delete[] mbVarbetaM;
                    delete[] mbPvalM;
                    delete[] mbPvalInt;
                    delete[] mbPvalJoint;
                }
            }
        }

        delete[] betaAll;
        delete[] VarBetaAll;
        delete[] invZGStZGS;
        delete[] PvalInt;
        delete[] PvalJoint;
        delete[] PvalM;
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
        printExecutionTime1(start_time, end_time);
}
    


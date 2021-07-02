
#include "declars.h"
#include "ReadBGEN.h"
#include "ReadParameters.h"
#include "../thirdparty/zstd-1.4.5/lib/zstd.h"
#include "../thirdparty/libdeflate-1.6/libdeflate.h"

/**************************************
This function is revised based on the Parse function in BOLT-LMM v2.3 source code
*************************************/
typedef std::numeric_limits< double > dbl;
// Read BGEN v1.1, v1.2, and v1.3 header block and flags.
void Bgen::processBgenHeaderBlock(string bgenfile) {

    char genofile[300];
    strcpy(genofile, bgenfile.c_str());
    string genopath(genofile);
    if (genopath.substr(genopath.length() - 5, 5) != ".bgen") {
        cout << "\nERROR: " << genopath << " does not have a .bgen extension. \n\n";
        exit(1);
    }

    fin = fopen(genofile, "rb");
    if (fin == 0) {
        cerr << "\nERROR: BGEN file could not be opened.\n\n";
        exit(1);
    }

    cout << "General information of BGEN file: \n";
    if (!fread(&offset, 4, 1, fin)) {
        cerr << "\nERROR: Cannot read BGEN header (offset).\n\n";
        exit(1);
    }

    uint L_H; 
    if (!fread(&L_H, 4, 1, fin)) {
        cerr << "\nERROR: Cannot read BGEN header (LH).\n\n";
        exit(1);
    }

    if (fread(&Mbgen, 4, 1, fin)) { 
        if (Mbgen <= 0) {
            cerr << "\nERROR: The number of variants in the BGEN file is 0.\n\n";
            exit(1);
        }
        cout << "Number of variants: " << Mbgen << '\n';
    } 
    else { 
        cerr << "\nERROR: Cannot read BGEN header (M).\n\n"; 
        exit(1);
    }

    if (fread(&Nbgen, 4, 1, fin)) {
        if (Nbgen <= 0) {
            cerr << "\nERROR: The number of samples in the BGEN file is 0.\n\n";
            exit(1);
        }
        cout << "Number of samples: " << Nbgen << '\n';
    }
    else {
        cerr << "\nERROR: Cannot read BGEN header (N). \n\n";
        exit(1);
    }


    char magic[5]; 
    if (!fread(magic, 1, 4, fin)) { 
        cerr << "\nERROR: Cannot read BGEN header (magic bytes). \n\n"; 
        exit(1);
    }
    magic[4] = '\0';
    if (!(magic[0] == 'b' && magic[1] == 'g' && magic[2] == 'e' && magic[3] == 'n')) {
        cerr << "\nERROR: BGEN file's four magic number bytes does not match 'b' 'g' 'e' 'n'.\n\n";
        exit(1);
    }

    fseek(fin, L_H - 20, SEEK_CUR);         

    uint flags; 
    if (!fread(&flags, 4, 1, fin)) {
        cerr << "\nERROR: Cannot read BGEN header (flags). \n\n";
        exit(1); 
    }


    // The header block - flag definitions
    CompressedSNPBlocks = flags & 3;
    
    switch (CompressedSNPBlocks) {
    case 0:
        cout << "Genotype Block Compression Type: Uncompressed\n";
        break;
    case 1:
        cout << "Genotype Block Compression Type: Zlib\n";
        break;
    case 2:
        cout << "Genotype Block Compression Type: Zstd\n";
        break;
    default:
        cout << "\nERROR: BGEN compression flag must be 0 (uncompressed), 1 (zlib compression), or 2 (zstd compression). Value in file: " << CompressedSNPBlocks << ".\n\n";
        exit(1);
    }

    Layout = (flags >> 2) & 0xf; 
    cout << "Layout: " << Layout << '\n';
    if (Layout != 1U && Layout != 2U) {
        cerr << "\nERROR: BGEN layout flag must be 1 or 2.\n\n";
        exit(1);
    }

    SampleIdentifiers = flags >> 31; 
    cout << "Sample Identifiers Present: ";  
    SampleIdentifiers == 0 ? cout << "False \n" : cout << "True \n";
    if (SampleIdentifiers != 0 && SampleIdentifiers != 1) {
        cerr << "\nERROR: BGEN sample identifier flag must be 0 or 1.\n\n";
        exit(1);
    }
    
}






/**********************************************************************************
This function is revised based on the Parse function in BOLT-LMM v2.3 source code
***********************************************************************************/

// This functions reads the sample block of BGEN v1.1, v1.2, and v1.3. Also finds which samples to remove if they have missing values in the pheno file.
void Bgen::processBgenSampleBlock(Bgen bgen, char samplefile[300], bool useSample, unordered_map<string, vector<string>> phenomap, string phenoMissingKey, int numSelCol, int samSize) {


    int k = 0;
    unordered_set<int> genoUnMatchID;
    std::vector<string> tempID;
    new_phenodata.resize(samSize);
    new_covdata.resize(samSize * (numSelCol+1));
    if ((bgen.SampleIdentifiers == 0) || useSample) {

        if (bgen.SampleIdentifiers == 0 && !useSample) {
            cerr << "\nERROR: BGEN file does not contain sample identifiers. A .sample file is required. \n"
                << "       See https://www.well.ox.ac.uk/~gav/qctool/documentation/sample_file_formats.html for .sample file format. \n\n";
            exit(1);
        }

        std::ifstream fIDMat;
        fIDMat.open(samplefile);
        if (!fIDMat.is_open()) {
            cerr << "\nERROR: Sample file could not be opened.\n\n";
            exit(1);
        }

        string IDline;
        getline(fIDMat, IDline);
        getline(fIDMat, IDline);
        uint nSamples = 0;
        while (getline(fIDMat, IDline)) {
            nSamples++;
        }
        if (nSamples != bgen.Nbgen) {
            cout << "\nERROR: Number of sample identifiers in .sample file (" << nSamples << ") does not match the number of samples specified in BGEN file (" << bgen.Nbgen << ").\n\n";
            exit(1);
        }
        else {
            fIDMat.clear();
            fIDMat.seekg(0, fIDMat.beg);
            getline(fIDMat, IDline);
            getline(fIDMat, IDline);
        }

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
                    sscanf(tmp_valvec[0].c_str(), "%lf", &new_phenodata[k]);
                    new_covdata[k * (numSelCol+1)] = 1.0;
                    for (int c = 0; c < numSelCol; c++) {
                        sscanf(tmp_valvec[c + 1].c_str(), "%lf", &new_covdata[k * (numSelCol + 1) + c + 1]);
                    }
                    sampleID.push_back(strtmp);
                    k++;
                }
            }

            // save the index with unmatched ID into genoUnMatchID.
            if (itmp == k) genoUnMatchID.insert(m);
        }
        fIDMat.close();
    }

    if ((bgen.SampleIdentifiers == 1) && !useSample) {

        uint maxLA = 65536;
        char* samID = new char[maxLA + 1];

        uint LS1;  
        if (!fread(&LS1, 4, 1, bgen.fin)) {
            cerr << "\nERROR: Cannot read BGEN sample block (LS).\n\n";
            exit(1);
        }

        uint Nrow; 
        if (!fread(&Nrow, 4, 1, bgen.fin)) {
            cerr << "\nERROR: Cannot read BGEN sample block (N).\n\n";
            exit(1);
        }
        if (Nrow != bgen.Nbgen) {
            cerr << "\nERROR: Number of sample identifiers (" << Nrow << ") does not match number of samples specified in BGEN file (" << bgen.Nbgen << ").\n\n";
            exit(1);
        }


        for (uint m = 0; m < bgen.Nbgen; m++) {
            ushort LSID; 
            if (!fread(&LSID, 2, 1, bgen.fin)) { 
                cerr << "\nERROR: Cannot read BGEN sample block (LSID).\n\n"; 
                exit(1); 
            }
            if (!fread(samID, 1, LSID, bgen.fin)) {
                cerr << "\nERROR: Cannot read BGEN sample block (sample id).\n\n";
                exit(1);
            }
            samID[LSID] = '\0';

            string strtmp(samID);
            int itmp = k;
            
            if (m < 5) {
                tempID.push_back(strtmp);
            }

            if (phenomap.find(strtmp) != phenomap.end()) {
                auto tmp_valvec = phenomap[strtmp];
                if (find(tmp_valvec.begin(), tmp_valvec.end(), phenoMissingKey) == tmp_valvec.end()) {
                    sscanf(tmp_valvec[0].c_str(), "%lf", &new_phenodata[k]);
                    new_covdata[k * (numSelCol+1)] = 1.0;
                    for (int c = 0; c < numSelCol; c++) {
                        sscanf(tmp_valvec[c + 1].c_str(), "%lf", &new_covdata[k * (numSelCol + 1) + c + 1]);
                    }
                    sampleID.push_back(strtmp);
                    k++;
                }
            }

            if (itmp == k) {
                genoUnMatchID.insert(m);
            }
        }

        delete[] samID;
    } // end SampleIdentifiers == 1


    // After IDMatching, resizing phenodata and covdata, and updating samSize;
    new_phenodata.resize(k);
    new_covdata.resize(k * (numSelCol + 1));
    samSize = k;


    if (samSize == 0) {
        cerr << "\nERROR: Sample size changed from " << samSize + genoUnMatchID.size() << " to " << samSize << ".\n\n";
        if (bgen.SampleIdentifiers == 1 && !useSample) {
            int print_i = 5;
            if (bgen.Nbgen < 5) { 
                print_i = bgen.Nbgen; 
            }
            cout << "ID matching was done using the BGEN sample identifier block. \nHere are the first " << print_i << " sample identifiers in BGEN file: \n";
            for (int i = 0; i < print_i; i++) {
                cout << " " << tempID[i] << "\n";
            }
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
    }
    else {
        cout << "After processes of sample IDMatching and checking missing values, the sample size changes from "
            << samSize + genoUnMatchID.size() << " to " << samSize << ".\n\n";
    }
    cout << "Sample IDMatching and checking missing values processes have been completed.\n";
    cout << "New pheno and covariate data vectors with the same order of sample ID sequence of geno data are updated.\n";
    cout << "****************************************************************************\n";


    new_samSize = samSize;
}







/***********************************************************************************
This function contains code that is revised based on BOLT-LMM v2.3 source code
************************************************************************************/

// This function reads just the variant block for BGEN files version v1.1, v1.2, and v1.3 and is used to grab the byte where the variant begins.
//    Necesary when there's no bgen index file.
void Bgen::getPositionOfBgenVariant(Bgen bgen, CommandLine cmd) {


    int count = 0;
    uint CompressedSNPBlocks = bgen.CompressedSNPBlocks;
    uint Layout = bgen.Layout;
    uint offset = bgen.offset;
    uint Mbgen = bgen.Mbgen;
    uint nSNPS = Mbgen;
    uint Nbgen = bgen.Nbgen;
    uint maxLA = 65536;
    char* snpID   = new char[maxLA + 1];
    char* rsID    = new char[maxLA + 1];
    char* chrStr  = new char[maxLA + 1];
    char* allele1 = new char[maxLA + 1];
    char* allele0 = new char[maxLA + 1];
    threads = cmd.threads;
    string IDline;

    std::set<std::string> includeVariant;
    std::vector<std::vector<uint>> includeVariantIndex;
    bool checkSNPID = false;
    bool checkRSID = false;
    bool checkInclude = false;
    int ret;

    if (cmd.doFilters) {

        filterVariants = true;
        if (!cmd.includeVariantFile.empty()) {
            checkInclude = true;

            std::ifstream fInclude;
            fInclude.open(cmd.includeVariantFile);
            if (!fInclude.is_open()) {
                cerr << "\nERROR: The file (" << cmd.includeVariantFile << ") could not be opened.\n\n";
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
                cerr << "\nERROR: Header name of " << cmd.includeVariantFile << " must be snpid or rsid.\n\n";
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


        cout << "Dividing BGEN file into " << threads << " block(s)...\n";
        cout << "Identifying start position of each block...\n";
        vector<uint> endIndex(threads);
        int nBlocks = ceil(nSNPS / threads);
        uint index = 0;
        uint k = 0;
        uint sucessCount = 0;
        Mbgen_begin.resize(threads);
        Mbgen_end.resize(threads);
        bgenVariantPos.resize(threads);
        keepVariants.resize(threads);

        for (uint t = 0; t < threads; t++) {
            endIndex[t] = ((t + 1) == threads) ? nSNPS - 1 : floor(((nSNPS / threads) * (t + 1)) - 1);
        }



        FILE* fin = bgen.fin;
        fseek(fin, offset + 4, SEEK_SET);
        for (uint snploop = 0; snploop < Mbgen; snploop++) {
            long long unsigned int prev = ftell(fin);

            uint Nrow;
            if (Layout == 1) {
                ret = fread(&Nrow, 4, 1, fin); 
                if (Nrow != Nbgen) {
                    cerr << "\nERROR: Number of samples (" << Nrow << ") with genotype probabilities does not match number of samples specified in BGEN file (" << Nbgen << ").\n\n";
                    exit(1);
                }
            }

            ushort LS; 
            ret = fread(&LS, 2, 1, fin);

            ret = fread(snpID, 1, LS, fin);
            snpID[LS] = '\0';
            if (checkSNPID) {
                if ((checkInclude) && (includeVariant.find(snpID) != includeVariant.end())) {
                    sucessCount++;
                    keepVariants[k].push_back(snploop);
                    if (index == (nBlocks * k)) {
                        Mbgen_begin[k] = snploop;
                        long long int curr = ftell(fin);
                        bgenVariantPos[k] = curr - (curr - (prev));
                    }
                    if (index == endIndex[k]) {
                        Mbgen_end[k] = snploop;
                        k++;
                        if (k == threads) {break;}
                    }
                    index++;
                }
            }

            ushort LR; 
            ret = fread(&LR, 2, 1, fin);

            ret = fread(rsID, 1, LR, fin); 
            rsID[LR] = '\0';
            if (checkRSID) {
                if (checkInclude && (includeVariant.find(rsID) != includeVariant.end())) {
                    sucessCount++;
                    keepVariants[k].push_back(snploop);
                    if (index == (nBlocks * k)) {
                        Mbgen_begin[k] = snploop;
                        long long unsigned int curr = ftell(fin);
                        bgenVariantPos[k] = curr - (curr - (prev));
                    }
                    if (index == endIndex[k]) {
                        Mbgen_end[k] = snploop;
                        k++;
                        if (k == threads) {break;}
                    }
                    index++;
                }
            }

            ushort LC; 
            ret = fread(&LC, 2, 1, fin);

            ret = fread(chrStr, 1, LC, fin); 
            chrStr[LC] = '\0';

            uint32_t physpos; 
            ret = fread(&physpos, 4, 1, fin);

            uint16_t LKnum;
            if (Layout == 2) {
                ret = fread(&LKnum, 2, 1, fin);
                if (LKnum != 2) {
                    cerr << "\nERROR: " << string(snpID) << " is a non-bi-allelic variant with " << LKnum << " alleles. Please filter these variants for now.\n\n";
                    exit(1);
                }
            }

            uint32_t LA; 
            ret = fread(&LA, 4, 1, fin);
            ret = fread(allele1, 1, LA, fin); 
            allele1[LA] = '\0';

            uint32_t LB; 
            ret = fread(&LB, 4, 1, fin);
            ret = fread(allele0, 1, LB, fin); 
            allele0[LB] = '\0';


            if (Layout == 2) {
                if (CompressedSNPBlocks > 0) {
                    uint zLen; 
                    ret = fread(&zLen, 4, 1, fin);
                    fseek(fin, 4 + zLen - 4, SEEK_CUR);

                }
                else {
                    uint zLen; 
                    ret = fread(&zLen, 4, 1, fin);
                    fseek(fin, zLen, SEEK_CUR);
                }
            }
            else {
                if (CompressedSNPBlocks == 1) {
                    uint zLen; 
                    ret = fread(&zLen, 4, 1, fin);
                    fseek(fin, zLen, SEEK_CUR);

                }
                else {
                    fseek(fin, 6 * Nbgen, SEEK_CUR);
                }
            }
        }

        if (sucessCount != nSNPS) {
            cerr << "\nERROR: There are one or more SNPs in BGEN file with " << IDline << " not in " << cmd.includeVariantFile << ".\n\n";
            exit(1);
        }
    }
    else {

        filterVariants = false;

        cout << "Detected " << boost::thread::hardware_concurrency() << " available thread(s)...\n";
        if (Mbgen < threads) {
            threads = Mbgen;
            cout << "Number of variants (" << Mbgen << ") is less than the number of specified threads (" << cmd.threads << ")...\n";
            cout << "Using " << threads << " for multithreading... \n\n";
        }
        else {
            cout << "Using " << threads << " for multithreading... \n\n";
        }

        cout << "Dividing BGEN file into " << threads << " block(s)..." << endl;
        Mbgen_begin.resize(threads);
        Mbgen_end.resize(threads);
        bgenVariantPos.resize(threads);
        keepVariants.resize(threads);

        for (uint t = 0; t < threads-1; t++) {
            Mbgen_begin[t] = floor((Mbgen / threads) * t);
            Mbgen_end[t] = floor(((Mbgen / threads) * (t + 1)) - 1);
        }
        Mbgen_begin[threads-1] = floor((Mbgen / threads) * (threads - 1));
        Mbgen_end[threads-1] = Mbgen - 1;


        uint t = 0;
        FILE* fin = bgen.fin;
        fseek(fin, offset + 4, SEEK_SET);
        for (uint snploop = 0; snploop < Mbgen; snploop++) {

            if (snploop == Mbgen_begin[t]) {
                bgenVariantPos[t] = ftell(fin);
                t++;
                if (t == (Mbgen_begin.size())) {
					break;
				}
            }

            uint Nrow;
            if (Layout == 1) {
                ret = fread(&Nrow, 4, 1, fin);
                if (Nrow != Nbgen) {
                    cerr << "\nERROR: Number of samples (" << Nrow << ") with genotype probabilities does not match number of samples specified in BGEN file (" << Nbgen << ").\n\n";
                    exit(1);
                }
            }

            ushort LS; 
            ret = fread(&LS, 2, 1, fin);
            ret = fread(snpID, 1, LS, fin); 
            snpID[LS] = '\0';

            ushort LR; 
            ret = fread(&LR, 2, 1, fin);
            ret = fread(rsID, 1, LR, fin); 
            rsID[LR] = '\0';

            ushort LC; 
            ret = fread(&LC, 2, 1, fin);
            ret = fread(chrStr, 1, LC, fin); 
            chrStr[LC] = '\0';

            uint32_t physpos; 
            ret = fread(&physpos, 4, 1, fin);

            uint16_t LKnum;
            if (Layout == 2) {
                ret = fread(&LKnum, 2, 1, fin);
                if (LKnum != 2) {
                    cerr << "\nERROR: " << string(snpID) << " is a non-bi-allelic variant with " << LKnum << " alleles. Please filter these variants for now.\n\n";
                    exit(1);
                }
            }

            uint32_t LA; 
            ret = fread(&LA, 4, 1, fin);
            ret = fread(allele1, 1, LA, fin); 
            allele1[LA] = '\0';

            uint32_t LB; 
            ret = fread(&LB, 4, 1, fin);
            ret = fread(allele0, 1, LB, fin); 
            allele0[LB] = '\0';


            // Seeks past the uncompressed genotype.
            if (Layout == 2) {
                if (CompressedSNPBlocks > 0) {
                    uint zLen;  
                    ret = fread(&zLen, 4, 1, fin);
                    ret = fseek(fin, 4 + zLen - 4, SEEK_CUR);

                }
                else {
                    uint zLen; 
                    ret = fread(&zLen, 4, 1, fin);
                    ret = fseek(fin, zLen, SEEK_CUR);
                }
            }
            else {
                if (CompressedSNPBlocks == 1) {
                    uint zLen;  
                    ret = fread(&zLen, 4, 1, fin);
                    ret = fseek(fin, zLen, SEEK_CUR);

                }
                else {
                    ret = fseek(fin, 6 * Nbgen, SEEK_CUR);
                }
            }
        }
    }

    (void)ret;
    delete[] snpID;
    delete[] rsID;
    delete[] chrStr;
    delete[] allele1;
    delete[] allele0;
}




void Bgen13GetTwoVals(const unsigned char* prob_start, uint32_t bit_precision, uintptr_t offset, uintptr_t* first_val_ptr, uintptr_t* second_val_ptr) {

    switch (bit_precision) {
    case 8:
        *first_val_ptr = prob_start[0];
        prob_start += offset;
        *second_val_ptr = prob_start[0];
        break;
    case 16:
        *first_val_ptr = prob_start[0] | (prob_start[1] << 8);
        prob_start += offset;
        *second_val_ptr = prob_start[0] | (prob_start[1] << 8);
        break;
    case 24:
        *first_val_ptr = prob_start[0] | (prob_start[1] << 8) | (prob_start[2] << 16);
        prob_start += offset;
        *second_val_ptr = prob_start[0] | (prob_start[1] << 8) | (prob_start[2] << 16);
        break;
    case 32:
        *first_val_ptr = prob_start[0] | (prob_start[1] << 8) | (prob_start[2] << 16) | (prob_start[3] << 24);
        prob_start += offset;
        *second_val_ptr = prob_start[0] | (prob_start[1] << 8) | (prob_start[2] << 16) | (prob_start[3] << 24);
        break;
    }

}





/*************************************************************************************************************************
This function contains code that has been revised based on BOLT-LMM v2.3 source code
**************************************************************************************************************************/

void gemBGEN(int thread_num, double sigma2, double* resid, double* XinvXTX, vector<double> miu,  BinE binE, Bgen bgen, CommandLine cmd) {

    auto start_time = std::chrono::high_resolution_clock::now();
    std::string output = cmd.outFile + "_bin_" + std::to_string(thread_num) + ".tmp";
    std::ofstream results(output, std::ofstream::binary);
    std::ostringstream oss;


    uint maxLA = 65536;
    char* snpID   = new char[maxLA + 1];
    char* rsID    = new char[maxLA + 1];
    char* chrStr  = new char[maxLA + 1];
    char* allele1 = new char[maxLA + 1];
    char* allele0 = new char[maxLA + 1];
    uint Nbgen  = bgen.Nbgen;
    uint Layout = bgen.Layout;
    uint CompressedSNPBlocks = bgen.CompressedSNPBlocks;

    string physpos_tmp;
    vector <uchar> zBuf;
    vector <uchar> shortBuf;

    vector <uchar> zBuf1;
    vector <uint16_t> shortBuf1;
    uLongf destLen1 = 6 * Nbgen;
    if (Layout == 1) {
        if (CompressedSNPBlocks == 0) {
            zBuf1.resize(destLen1);
        }
        else {
            shortBuf1.resize(destLen1);
        }
    }

    bool filterVariants = bgen.filterVariants;
    string outStyle = cmd.outStyle;
    int phenoType   = bgen.phenoType;
    int stream_snps = cmd.stream_snps;
    int samSize     = bgen.new_samSize;
    int robust      = cmd.robust;
    int intSq1      = cmd.numIntSelCol + 1;
    int expSq       = cmd.numExpSelCol;
    int expSq1      = expSq+1;
    int Sq1         = intSq1 + expSq;
    int Sq          = Sq1-1;
    int numSelCol1  = cmd.numSelCol + Sq1;
    double MAF      = cmd.MAF;
    double maxMAF   = 1 - MAF;
    double missGenoCutoff = cmd.missGenoRate;

    vector<long int> include_idx = bgen.include_idx;
    vector<uint> keepVariants = bgen.keepVariants[thread_num];
    uint snploop = bgen.Mbgen_begin[thread_num], end = bgen.Mbgen_end[thread_num];
   
    int numBinE   = binE.nBinE;
    int strataLen = binE.strataLen;
    bool strata   = (numBinE > 0 ) ? true : false;
    vector<int> stratum_idx = binE.stratum_idx;
    vector<double> binE_AF(stream_snps * strataLen, 0.0), binE_N(stream_snps * strataLen, 0.0);
   
    int ZGS_col  = Sq1 * stream_snps;
    vector <double> ZGSvec(samSize   * (Sq1) * stream_snps);
    vector <double> ZGSR2vec(samSize * (Sq1) * stream_snps);
    vector <double> WZGSvec(samSize  * (Sq1) * stream_snps);
    vector <double> AF(stream_snps);
    vector<uint> missingIndex;
    vector <string> geno_snpid(stream_snps);
    double* WZGS = &WZGSvec[0];
    double* covX = &bgen.new_covdata[0];
    boost::math::chi_squared chisq_dist_M(1);

    struct libdeflate_decompressor* decompressor = libdeflate_alloc_decompressor();

    int printStart = 1; 
    int printEnd   = expSq1;
    bool printMeta = false;
    bool printFull = false;
    if (expSq == 0) {
        printStart = 0; printEnd = 0;
    }
    else if (outStyle.compare("meta") == 0) {
        printStart = 0; printEnd = Sq1; printMeta = true;
    } else if (outStyle.compare("full") == 0) {
        printStart = 0; printEnd = Sq1; printFull = true;
    }


    FILE* fin3;
    fin3 = fopen(cmd.bgenFile.c_str(), "rb");
    long long unsigned int byte = bgen.bgenVariantPos[thread_num];
    fseek(fin3, byte, SEEK_SET);

    int variant_index = 0;
    int keepIndex = 0;
    int ret;
    while (snploop <= end) {

        int stream_i = 0;
        while (stream_i < stream_snps) {

            if (snploop == (end + 1) && stream_i == 0) {
                break;
            }
            if (snploop == end + 1 && stream_i != 0) {
                stream_snps = stream_i;
                ZGS_col = Sq1 * stream_snps;
                break;
            }
            snploop++;


            uint Nrow;
            if (Layout == 1) {
                ret = fread(&Nrow, 4, 1, fin3);
                if (Nrow != Nbgen) {
                    cerr << "\nERROR: Number of samples (" << Nrow << ") with genotype probabilities does not match number of samples specified in BGEN file (" << Nbgen << ").\n\n";
                    exit(1);
                }
            }

            ushort LS; 
            ret = fread(&LS, 2, 1, fin3);
            ret = fread(snpID, 1, LS, fin3); 
            snpID[LS] = '\0';

            ushort LR; 
            ret = fread(&LR, 2, 1, fin3);
            ret = fread(rsID, 1, LR, fin3); 
            rsID[LR] = '\0';

            ushort LC; 
            ret = fread(&LC, 2, 1, fin3);
            ret = fread(chrStr, 1, LC, fin3); 
            chrStr[LC] = '\0';

            uint32_t physpos; 
            ret = fread(&physpos, 4, 1, fin3);
            physpos_tmp = std::to_string(physpos);

            uint16_t LKnum;
            if (Layout == 2) {
                ret = fread(&LKnum, 2, 1, fin3);
                if (LKnum != 2) {
                    cout << "\nERROR: " << snpID << " is a non-bi-allelic variant with " << LKnum << " alleles. Please filter these variant for now. \n\n";
                    exit(1);
                }
            }

            uint32_t LA;  
            ret = fread(&LA, 4, 1, fin3);
            ret = fread(allele1, 1, LA, fin3); 
            allele1[LA] = '\0';

            uint32_t LB; 
            ret = fread(&LB, 4, 1, fin3);
            ret = fread(allele0, 1, LB, fin3); 
            allele0[LB] = '\0';

            uint nMissing = 0;
            int tmp1 = stream_i * Sq1 * samSize;
            int strata_i = stream_i * strataLen;
            if (Layout == 1) {
                uint16_t* probs_start;
                if (CompressedSNPBlocks == 1) {
                    uint zLen; 
                    ret = fread(&zLen, 4, 1, fin3);
                    zBuf1.resize(zLen);
                    ret = fread(&zBuf1[0], 1, zLen, fin3);
                    if (libdeflate_zlib_decompress(decompressor, &zBuf1[0], zLen, &shortBuf1[0], destLen1, NULL) != LIBDEFLATE_SUCCESS) {
                        cerr << "\nERROR: Decompressing " << snpID << " block failed with libdeflate.\n\n";
                        exit(1);
                    }
                    probs_start = &shortBuf1[0];
                }
                else {
                    ret = fread(&zBuf1[0], 1, destLen1, fin3);
                    probs_start = reinterpret_cast<uint16_t*>(&zBuf1[0]);
                }


                // read genotype probabilities
                const double scale = 1.0 / 32768;

                int idx_k = 0;
                for (uint i = 0; i < Nbgen; i++) {
                    if (include_idx[idx_k] == i) {
                        double p11 = probs_start[3 * i] * scale;
                        double p10 = probs_start[3 * i + 1] * scale;
                        double p00 = probs_start[3 * i + 2] * scale;

                        if (p11 == 0 && p10 == 0 && p00 == 0) {
                            missingIndex.push_back(idx_k);
                            nMissing++;
                        }
                        else {
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

                            if (strata) {
                                binE_N[strata_i + stratum_idx[idx_k]]+=1.0;
                                binE_AF[strata_i + stratum_idx[idx_k]]+=dosage;
                            }
                        }

                        idx_k++;
                    }
                }

            } // end of reading genotype data when Layout = 1


            if (Layout == 2) {
                uint zLen; 
                ret = fread(&zLen, 4, 1, fin3);

                if (filterVariants && keepVariants[keepIndex] + 1 != snploop) {
                    CompressedSNPBlocks > 0 ? fseek(fin3, 4 + zLen - 4, SEEK_CUR) : fseek(fin3, zLen, SEEK_CUR);
                    continue;
                }

                uint DLen;
                uchar* bufAt;
                if (CompressedSNPBlocks == 1) {
                    zBuf.resize(zLen - 4);
                    ret = fread(&DLen, 4, 1, fin3);
                    ret = fread(&zBuf[0], 1, zLen - 4, fin3);
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
                    ret = fread(&DLen, 4, 1, fin3);
                    ret = fread(&zBuf[0], 1, zLen - 4, fin3);
                    shortBuf.resize(DLen);
                    uLongf destLen = DLen;

                    size_t ret = ZSTD_decompress(&shortBuf[0], destLen, &zBuf[0], zLen - 4);
                    if (ZSTD_isError(ret)) {
                        cout << "ZSTD ERROR: " << ZSTD_getErrorName(ret);
                    }
                    bufAt = &shortBuf[0];
                }
                else {
                    zBuf.resize(zLen);
                    ret = fread(&zBuf[0], 1, zLen, fin3);
                    bufAt = &zBuf[0];
                }


                uint32_t N; 
                memcpy(&N, bufAt, sizeof(int32_t));
                if (N != Nbgen) {
                    cerr << "\nERROR: " << snpID << " number of samples (" << N << ") with genotype probabilties does not match number of samples specified in BGEN file (" << Nbgen << ").\n\n";
                    exit(1);
                }

                uint16_t K; 
                memcpy(&K, &(bufAt[4]), sizeof(int16_t));
                if (K != 2U) {
                    cout << "\nERROR: There are SNP(s) with more than 2 alleles (non-bi-allelic). Currently unsupported. \n\n";
                    exit(1);
                }
            
                const uint32_t min_ploidy = bufAt[6];
                if (min_ploidy != 2U) {
                    cerr << "\nERROR: " << snpID << " has minimum ploidy " << min_ploidy << ". Currently unsupported. \n\n";
                    exit(1);
                }

                const uint32_t max_ploidy = bufAt[7];
                if (max_ploidy != 2U) {
                    cerr << "\nERROR: " << snpID << " has maximum ploidy " << max_ploidy << ". Currently unsupported. \n\n";
                    exit(1);
                }

                const unsigned char* missing_and_ploidy_info = &(bufAt[8]);
                const unsigned char* probs_start = &(bufAt[10 + N]);
                
                const uint32_t is_phased = probs_start[-2];
                if (is_phased != 1 && is_phased != 0) {
                    cerr << "\nERROR: " << snpID << " has phased value of " << is_phased << ". This must be 0 or 1. \n\n";
                    exit(1);
                }
                
                const uint32_t bit_precision = probs_start[-1];
                if (bit_precision != 8 && bit_precision != 16 && bit_precision != 24 && bit_precision != 32) {
                    cerr << "\nERROR: Bits to store probabilities must be 8, 16, 24, or 32. \n\n";
                    exit(1);
                }
                const uintptr_t numer_mask = (1U << bit_precision) - 1;
                const uintptr_t probs_offset = bit_precision / 8;


                int idx_k = 0;
                if (!is_phased) {
                    for (uint32_t i = 0; i < N; i++) {
                        const uint32_t missing_and_ploidy = missing_and_ploidy_info[i];
                        uintptr_t numer_aa;
                        uintptr_t numer_ab;

                        if (missing_and_ploidy == 2) {
                            Bgen13GetTwoVals(probs_start, bit_precision, probs_offset, &numer_aa, &numer_ab);
                            probs_start += (probs_offset * 2);
                        }
                        else if (missing_and_ploidy == 130) {
                            if (include_idx[idx_k] == i) {
                                nMissing++;
                                missingIndex.push_back(idx_k);
                                idx_k++;
                            }
                            probs_start += (probs_offset * 2);
                            continue;
                        }
                        else {
                            cerr << "\nERROR: Ploidy value " << missing_and_ploidy << " is unsupported. Must be 2 or 130 (missing). \n\n";
                            exit(1);
                        }

                        if (include_idx[idx_k] == i) {
                            double p11 = numer_aa / double(1.0 * (numer_mask));
                            double p10 = numer_ab / double(1.0 * (numer_mask));
                            double dosage = 2 * (1 - p11 - p10) + p10;

                            int tmp2 = idx_k + tmp1;
                            AF[stream_i] += dosage;

                            if (phenoType == 1) {
                                ZGSvec[tmp2] = miu[idx_k] * (1 - miu[idx_k]) * dosage;
                            }
                            else {
                                ZGSvec[tmp2] = dosage;
                            }

                            if (strata) {
                                binE_N[strata_i + stratum_idx[idx_k]]+=1.0;
                                binE_AF[strata_i + stratum_idx[idx_k]]+=dosage;
                            }
                            idx_k++;
                        }
                    }

                } else {
                    for (uint32_t i = 0; i < N; i++) {
                        const uint32_t missing_and_ploidy = missing_and_ploidy_info[i];
                        uintptr_t numer_aa;
                        uintptr_t numer_ab;

                        if (missing_and_ploidy == 2) {
                            Bgen13GetTwoVals(probs_start, bit_precision, probs_offset, &numer_aa, &numer_ab);
                            probs_start += (probs_offset * 2);
                        }
                        else if (missing_and_ploidy == 130) {
                            if (include_idx[idx_k] == i) {
                                nMissing++;
                                missingIndex.push_back(idx_k);
                                idx_k++;
                            }
                            probs_start += (probs_offset * 2);
                            continue;
                        }
                        else {
                            cerr << "\nERROR: Ploidy value " << missing_and_ploidy << " is unsupported. Must be 2 or 130. \n\n";
                            exit(1);
                        }

                        if (include_idx[idx_k] == i) {
                            double p11 = numer_aa / double(1.0 * (numer_mask));
                            double p10 = numer_ab / double(1.0 * (numer_mask));
                            double dosage = 2 - (p11 + p10);

                            int tmp2 = idx_k + tmp1;
                            AF[stream_i] += dosage;

                            if (phenoType == 1) {
                                ZGSvec[tmp2] = miu[idx_k] * (1 - miu[idx_k]) * dosage;
                            }
                            else {
                                ZGSvec[tmp2] = dosage;
                            }

                            if (strata) {
                                binE_N[strata_i + stratum_idx[idx_k]]+=1.0;
                                binE_AF[strata_i + stratum_idx[idx_k]]+=dosage;
                            }

                            idx_k++;
                        }
                    }
                }
            } // end of layout 2
    
            double gmean  = AF[stream_i] / double(samSize - nMissing);
            double cur_AF = AF[stream_i] / 2.0 / double(samSize - nMissing);
            double percMissing = nMissing / (samSize * 1.0);
            if ((cur_AF < MAF || cur_AF > maxMAF) || (percMissing > missGenoCutoff)) {
                AF[stream_i] = 0.0;
                if (strata) {
                    for (int i = 0; i < strataLen; i++) {
                        binE_N[strata_i + i] = 0.0;
                        binE_AF[strata_i + i] = 0.0;
                    }
                }
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
                    for (size_t nm = 0; nm < missingIndex.size(); nm++) {
                        int tmp5 = tmp1 + missingIndex[nm];
                        ZGSvec[tmp5] = gmean;
                    }
                }
                else {
                    for (size_t nm = 0; nm < missingIndex.size(); nm++) {
                        int tmp5 = tmp1 + missingIndex[nm];
                        ZGSvec[tmp5] = miu[missingIndex[nm]] * (1 - miu[missingIndex[nm]]) * gmean;
                    }
                }
                missingIndex.clear();
            }

            geno_snpid[stream_i] = ((LS > 0) ? string(snpID) : "NA") + "\t" + ((LR > 0) ? string(rsID) : "NA") + "\t" + ((LC > 0) ? string(chrStr) : "NA") + "\t" + physpos_tmp + "\t" + string(allele1) + "\t" + string(allele0) + "\t" + std::to_string(samSize - nMissing);
                
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
        } // end of stream_i

        if ((snploop == (end + 1)) & (stream_i == 0)) {
            break;
        }

        /***************************************************************/
        //	genodata and envirment data
        double* ZGS   = &ZGSvec[0];
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

            int tmp_strata = i * strataLen;
            for (int k = 0; k < strataLen; k++) {
                oss << binE_N[tmp_strata + k] << "\t" << binE_AF[tmp_strata + k] << "\t";
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

    } // end of snploop 

    if (variant_index % 10000 != 0) {
        results << oss.str();
        oss.str(std::string());
        oss.clear();
    }


    libdeflate_free_decompressor(decompressor);
    delete[] snpID;
    delete[] rsID;
    delete[] chrStr;
    delete[] allele1;
    delete[] allele0;
    (void)ret;

    // Close files
    results.close();
    fclose(fin3);

    auto end_time = std::chrono::high_resolution_clock::now();
    cout << "Thread " << thread_num << " finished in ";
    printExecutionTime1(start_time, end_time);
}






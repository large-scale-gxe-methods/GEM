/*  GEM : Gene-Environment interaction analysis for Millions of samples
 *  Copyright (C) 2018-2021  Liang Hong, Han Chen, Duy Pham
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


int  checkBinary(unordered_map<string, vector<string>> phenoMap, vector<string> sampleID);
void center(int center, int scale, int samSize, int numSelCol, vector<double> covdata, vector<double>* covdata_ret);
void fitNullModel(int samSize, int numSelCol, int phenoType, double epsilon, int robust, std::vector<string> covSelHeadersName, std::vector<double> phenodata, std::vector<double> covdata, std::vector<double>* XinvXTX_ret, vector<double>* miu_ret, vector<double>* resid_ret, double* sigma2_ret);
void printCovVarMat(int numCovs, vector<string> covNames, double* covVarMat, double* beta, int phenoType, int samSize);
void printOutputHeader(bool useBgen, int numExpSelCol, int Sq1, vector<string> covNames, string output, string outStyle, int robust, double sigma2, BinE binE);

int main(int argc, char* argv[]) {

    // Process command line
    CommandLine cmd;
    cmd.processCommandLine(argc, argv);


    // Parameters
    int samSize;
    int phenoCol;
    int samIDCol;
    int robust      = cmd.robust;
    char delim      = cmd.pheno_delim;
    double epsilon  = cmd.tol;
    string phenoHeaderName = cmd.phenoName;
    string samIDHeaderName = cmd.sampleID;
    string phenoMissingKey = cmd.missing;

    string output = cmd.outFile;
    int numSelCol    = cmd.numSelCol;
    int numExpSelCol = cmd.numExpSelCol;
    int numIntSelCol = cmd.numIntSelCol;
    int Sq = numExpSelCol + numIntSelCol;
    vector<string> covSelHeadersName    = cmd.cov;
    vector<string> expCovSelHeadersName = cmd.exp;
    vector<string> intCovSelHeadersName = cmd.icov;


    // Rearranging exposures, interaction covariates, and covariates for matrix operations
    numSelCol = numSelCol + numIntSelCol + numExpSelCol;
    vector<int> colSelVec(numSelCol);
    if (numIntSelCol != 0) {
        for (int i = numIntSelCol - 1; i >= 0; i--) { covSelHeadersName.insert(covSelHeadersName.begin(), intCovSelHeadersName[i]); }
    }
    for (int i = numExpSelCol - 1; i >= 0; i--) { covSelHeadersName.insert(covSelHeadersName.begin(), expCovSelHeadersName[i]); }


    // Start clock
    auto wall0 = std::chrono::system_clock::now();
    std::clock_t cpu0 = std::clock();



    //Reading phenotype file headers
    std::unordered_map<string, int> colNames;
    long unsigned int  phenoncols;

    string phenopath(cmd.phenoFile);
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

    for (int i = 0; i < numExpSelCol; i++) {
        if (colNames.find(expCovSelHeadersName[i]) == colNames.end()) {
            cerr << "\nERROR: Cannot find exposure column " << expCovSelHeadersName[i] << " in phenotype file. \n\n";
            exit(1);
        }
    }
    for (int i = 0; i < numIntSelCol; i++) {
        if (colNames.find(intCovSelHeadersName[i]) == colNames.end()) {
            cerr << "\nERROR: Cannot find interaction covariate column " << intCovSelHeadersName[i] << " in phenotype file. \n\n";
            exit(1);
        }
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
    colNames.clear();



    // Count sample size
    int nrows = 0;
    while (getline(finph, line)) nrows++;
    samSize = nrows;
    finph.clear();
    finph.seekg(0, finph.beg);
    getline(finph, line);
    cout << "Before ID Matching and checking missing values... \n";
    cout << "Size of the phenotype vector is: " << samSize << " X 1\n";
    cout << "Size of the selected covariate matrix (including first column for intercept values) is: " << samSize << " X " << numSelCol + 1 << '\n';


    // A Hashmap phenodata for IDMatching process.
    // key is smapleID in phenotype file,
    // value is a vector of pheno data as string for the sampleID
    unordered_map<string, vector<string>> phenomap;
    for (int r = 0; r < samSize; r++) {
        getline(finph, line);
        std::istringstream iss(line);
        string value;
        string temvalue;
        vector <string> values;
        while (getline(iss, value, delim)) values.push_back(value);
        if (values.size() != phenoncols) {
            cerr << "ERROR: Wrong number of entries in " << r  << " row.";
            cerr << "Expected " << phenoncols << " fields; parsed " << values.size() << '\n';
            exit(1);
        }
        phenomap[values[samIDCol]] = { values[phenoCol] };
        for (int c = 0; c < numSelCol; c++) {
            phenomap[values[samIDCol]].push_back(values[colSelVec[c]]);
        }
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


    double sigma2;
    BinE binE;
    if (cmd.usePgenFile) {
        Pgen pgen;

        pgen.processPgenHeader(cmd.pgenFile);
        pgen.processPvar(pgen, cmd.pvarFile);
        pgen.processPsam(pgen, cmd.psamFile, phenomap, phenoMissingKey, numSelCol, samSize);

        samSize = pgen.new_samSize;
        
        pgen.phenoType = checkBinary(phenomap, pgen.sampleID);
        binE.checkBinaryCovariates(binE, phenomap, pgen.sampleID, pgen.include_idx, samSize, numExpSelCol, numSelCol, covSelHeadersName);
        phenomap.clear();

        if ((cmd.center == 1) || (cmd.scale == 1)) {
            center(cmd.center, cmd.scale, samSize, numSelCol, pgen.new_covdata, &pgen.new_covdata);
        }
        
        cout << "Starting GWAS... \n\n";
        vector <double> miuvec(samSize), residvec(samSize);
        vector <double> XinvXTXvec(samSize* (numSelCol + 1));
        fitNullModel(samSize, numSelCol, pgen.phenoType, epsilon, robust, covSelHeadersName, pgen.new_phenodata, pgen.new_covdata, &XinvXTXvec, &miuvec, &residvec, &sigma2);
        pgen.new_phenodata.clear();
        
        pgen.getPgenVariantPos(pgen, cmd);
        cout << "The ALT allele in the .pvar file will be used for association testing.\n";
        auto start_time = std::chrono::high_resolution_clock::now();
        if (pgen.threads > 1) {
            cout << "Running multithreading...\n";
            boost::thread_group thread_grp;
		    for (uint i = 0; i < pgen.threads; i++) {
				thread_grp.create_thread(boost::bind(&gemPGEN, i, sigma2, &residvec[0], &XinvXTXvec[0], miuvec, boost::ref(binE), boost::ref(pgen), boost::ref(cmd)));
			}
            cout << "Joining threads... \n";
            thread_grp.join_all();
        }
        else {
            cout << "Running with single thread...\n";
            gemPGEN(0, sigma2, &residvec[0], &XinvXTXvec[0], miuvec, binE, pgen, cmd);

        }
        cmd.threads = pgen.threads;
        auto end_time = std::chrono::high_resolution_clock::now();
        printExecutionTime(start_time, end_time);
    }

    
    if (cmd.useBedFile) {
        Bed bed;
        bed.processBed(cmd.bedFile, cmd.bimFile, cmd.famFile);
        bed.processFam(bed, cmd.famFile, phenomap, phenoMissingKey, numSelCol, samSize);

        samSize = bed.new_samSize;

        bed.phenoType = checkBinary(phenomap, bed.sampleID);
        binE.checkBinaryCovariates(binE, phenomap, bed.sampleID, bed.include_idx, samSize, numExpSelCol, numSelCol, covSelHeadersName);
        phenomap.clear();

        if ((cmd.center == 1) || (cmd.scale == 1)) {
            center(cmd.center, cmd.scale, samSize, numSelCol, bed.new_covdata, &bed.new_covdata);
        }
        
        cout << "Starting GWAS... \n\n";
        vector <double> miuvec(samSize), residvec(samSize);
        vector <double> XinvXTXvec(samSize* (numSelCol + 1));
        fitNullModel(samSize, numSelCol, bed.phenoType, epsilon, robust, covSelHeadersName, bed.new_phenodata, bed.new_covdata, &XinvXTXvec, &miuvec, &residvec, &sigma2);
        bed.new_phenodata.clear();

        bed.getBedVariantPos(bed, cmd);
        cout << "The ALT allele in the .bim file will be used for association testing.\n";
        auto start_time = std::chrono::high_resolution_clock::now();
        if (bed.threads > 1) {
            cout << "Running multithreading...\n";
            boost::thread_group thread_grp;
			for (uint i = 0; i < bed.threads; i++) {
				thread_grp.create_thread(boost::bind(&gemBED, i, sigma2, &residvec[0], &XinvXTXvec[0], miuvec, boost::ref(binE), boost::ref(bed), boost::ref(cmd)));
			}
            cout << "Joining threads... \n";
			thread_grp.join_all();
        }
        else {
            cout << "Running with single thread...\n";
            gemBED(0, sigma2, &residvec[0], &XinvXTXvec[0], miuvec, binE, bed, cmd);
        }
        cmd.threads = bed.threads;
        auto end_time = std::chrono::high_resolution_clock::now();
        printExecutionTime(start_time, end_time);
    }

    
    if (cmd.useBgenFile) {
        Bgen bgen;
        bgen.processBgenHeaderBlock(cmd.bgenFile);
        bgen.processBgenSampleBlock(bgen, cmd.samplefile, cmd.useSampleFile, phenomap, phenoMissingKey, numSelCol, samSize);

        samSize = bgen.new_samSize;

        bgen.phenoType = checkBinary(phenomap, bgen.sampleID);
        binE.checkBinaryCovariates(binE, phenomap, bgen.sampleID, bgen.include_idx, samSize, numExpSelCol, numSelCol, covSelHeadersName);
        phenomap.clear();

        if ((cmd.center == 1) || (cmd.scale == 1)) {
            center(cmd.center, cmd.scale, samSize, numSelCol, bgen.new_covdata, &bgen.new_covdata);
        }
        
        cout << "Starting GWAS... \n\n";
        vector <double> miuvec(samSize), residvec(samSize);
        vector <double> XinvXTXvec(samSize * (numSelCol + 1));
        fitNullModel(samSize, numSelCol, bgen.phenoType, epsilon, robust, covSelHeadersName, bgen.new_phenodata, bgen.new_covdata, &XinvXTXvec, &miuvec, &residvec, &sigma2);
        bgen.new_phenodata.clear();

        auto start_time = std::chrono::high_resolution_clock::now();
        bgen.getPositionOfBgenVariant(bgen, cmd);
        auto end_time = std::chrono::high_resolution_clock::now();
        printExecutionTime(start_time, end_time);

        //Preparing for parallelizing of BGEN file
        cout << "The second allele in the BGEN file will be used for association testing.\n";
        start_time = std::chrono::high_resolution_clock::now();
        if (bgen.threads > 1) {
            cout << "Running multithreading...\n";
            boost::thread_group thread_grp;
			for (uint i = 0; i < bgen.threads; i++) {
				thread_grp.create_thread(boost::bind(&gemBGEN, i, sigma2, &residvec[0], &XinvXTXvec[0], miuvec, boost::ref(binE), boost::ref(bgen), boost::ref(cmd)));
			}
            cout << "Joining threads... \n";
			thread_grp.join_all();
        }
        else {
            cout << "Running with single thread...\n";
            gemBGEN(0, sigma2, &residvec[0], &XinvXTXvec[0], miuvec, binE, bgen, cmd);
        }
        cmd.threads = bgen.threads;
        end_time = std::chrono::high_resolution_clock::now();
        printExecutionTime(start_time, end_time);
    }


    // Write all results from each thread to 1 file
    cout << "Combining results... \n";
    auto start_time = std::chrono::high_resolution_clock::now();
    printOutputHeader(cmd.useBgenFile, numExpSelCol, Sq+1, covSelHeadersName, output, cmd.outStyle, cmd.robust, sigma2, binE);
    std::ofstream results(output, std::ios_base::app);
    for (int i = 0; i < cmd.threads; i++) {
         std::string threadOutputFile = cmd.outFile + "_bin_" + std::to_string(i) + ".tmp";
         std::ifstream thread_output(threadOutputFile);
         results << thread_output.rdbuf();

         thread_output.close();
         boost::filesystem::remove(threadOutputFile.c_str());
    }
    results.close();
    auto end_time = std::chrono::high_resolution_clock::now();
    printExecutionTime(start_time, end_time);


    // Finished
    std::chrono::duration<double> wallduration = (std::chrono::system_clock::now() - wall0);
    double cpuduration = (std::clock() - cpu0) / (double)CLOCKS_PER_SEC;
    cout << "Total Wall Time = " << wallduration.count() << "  Seconds\n";
    cout << "Total CPU Time  = " << cpuduration << "  Seconds\n";
    cout << "*********************************************************\n";
    
    
    return 0;
}

int checkBinary(unordered_map<string, vector<string>> phenoMap, vector<string> sampleID) {

    int is_bin = 1;
    std::unordered_map<string, int> map;
    size_t sampleSize = sampleID.size();
    for (size_t i = 0; i < sampleSize; i++) {
        auto tmp = phenoMap[sampleID[i]];

        if (!map.count(tmp[0])) {
            map[tmp[0]] = 1;
        }
        if ((map.size()/sampleSize) > 0.05) {
            is_bin = 0;
            break;
        }
    }
    if (map.size() < 2) {
        cerr << "ERROR: All values of the phenotype are the same. \n\n";
    }
    
    return is_bin;
}



void center(int center, int scale, int samSize, int numSelCol, vector<double> covdata, vector<double>* covdata_ret) 
{
    vector<double> tmp1(samSize, 1);
    double* tmpMean = new double[numSelCol + 1];
    vector<double> tmpSD(numSelCol + 1);
    if (center) 
    {
        matmatprod(&tmp1[0], &covdata[0], tmpMean, 1, samSize, numSelCol + 1);
        if (!scale) {
            cout << "Centering ALL exposures and covariates..." << endl;
            for (int i = 1; i < numSelCol + 1; i++) {
                tmpMean[i] /= double(samSize * 1.0);
                tmpSD[i] = 1.0;
            }
        }
        else {
            cout << "Centering and scaling ALL exposures and covariates..." << endl;
            for (int i = 1; i < numSelCol + 1; i++) {
                tmpMean[i] /= double(samSize * 1.0);
            }
            for (int i = 0; i < samSize; i++) {
                for (int j = 1; j < numSelCol + 1; j++) {
                    tmpSD[j] += pow(covdata[i * (numSelCol + 1) + j] - tmpMean[j], 2.0);
                }
            }
            for (int i = 1; i < numSelCol + 1; i++) {
                tmpSD[i] = sqrt(tmpSD[i] / double(samSize * 1.0 - 1.0));
            }
        }

        for (int i = 0; i < samSize; i++) {
            for (int j = 1; j < numSelCol + 1; j++) {
                covdata[i * (numSelCol + 1) + j] = (covdata[i * (numSelCol + 1) + j] - tmpMean[j]) / tmpSD[j];
            }
        }

    }
    else {
        if (scale) {
            cout << "Scaling ALL exposures and covariates..." << endl;
            matmatprod(&tmp1[0], &covdata[0], tmpMean, 1, samSize, numSelCol + 1);
            for (int i = 1; i < numSelCol + 1; i++) {
                tmpMean[i] /= double(samSize * 1.0);
            }
            for (int i = 0; i < samSize; i++) {
                for (int j = 1; j < numSelCol + 1; j++) {
                    tmpSD[j] += pow(covdata[i * (numSelCol + 1) + j] - tmpMean[j], 2.0);
                }
            }
            for (int i = 1; i < numSelCol + 1; i++) {
                tmpSD[i] = sqrt(tmpSD[i] / double(samSize * 1.0 - 1.0));
            }
            for (int i = 0; i < samSize; i++) {
                for (int j = 1; j < numSelCol + 1; j++) {
                    covdata[i * (numSelCol + 1) + j] /= tmpSD[j];
                }
            }
        }
    }
    delete[] tmpMean;
    *covdata_ret = covdata;
}


void fitNullModel(int samSize, int numSelCol, int phenoType, double epsilon, int robust, std::vector<string> covSelHeadersName, std::vector<double> phenodata, std::vector<double> covdata, std::vector<double>* XinvXTX_ret, vector<double>* miu_ret, vector<double>* resid_ret, double* sigma2_ret) 
{
    double* phenoY = &phenodata[0];
    double* covX = &covdata[0];
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

    // logistic regression
    while ((phenoType == 1) && (Check != (numSelCol + 1))) 
    {
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
    vector<double> XinvXTXvec(samSize * (numSelCol + 1));
    double* XinvXTX = &XinvXTXvec[0];
    if (phenoType == 1) {
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
    double* resid = &residvec[0];

    // sqr(sigma) = transpose(resid)*resid/[samSize-(numSelCol+1)]
    sigma2 = sigma2 / (samSize - (numSelCol + 1));
    if (phenoType == 1) sigma2 = 1.0;


    vector<double> XR2vec;
    if (!robust) 
    {
        for (int i = 0; i < (numSelCol + 1) * (numSelCol + 1); i++) {
            XTransX[i] = XTransX[i] * sigma2;
        }
        printCovVarMat(numSelCol + 1, covSelHeadersName, XTransX, beta, phenoType, samSize);
    }
    else {
        vector<double> XR2vec = covdata;
        for (int i = 0; i < samSize; i++) {
            for (int j = 0; j < numSelCol + 1; j++) {
                XR2vec[i * (numSelCol + 1) + j] = XR2vec[i * (numSelCol + 1) + j] * resid[i] * resid[i];
            }
        }

        double* XR2 = &XR2vec[0];
        double* XR2tX = new double[(numSelCol + 1) * (numSelCol + 1)];
        matTmatprod(XR2, covX, XR2tX, samSize, numSelCol + 1, numSelCol + 1);
        double* XTransXtXR2tX = new double[(numSelCol + 1) * (numSelCol + 1)];
        matmatTprod(XR2tX, XTransX, XTransXtXR2tX, numSelCol + 1, numSelCol + 1, numSelCol + 1);
        double* XTransXR2 = new double[(numSelCol + 1) * (numSelCol + 1)];
        matmatTprod(XTransX, XTransXtXR2tX, XTransXR2, numSelCol + 1, numSelCol + 1, numSelCol + 1);

        printCovVarMat(numSelCol + 1, covSelHeadersName, XTransXR2, beta, phenoType, samSize);
        delete[] XR2tX;
        delete[] XTransXtXR2tX;
        delete[] XTransXR2;
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    printExecutionTime(start_time, end_time);

    delete[] XTransX;
    delete[] XTransY;
    delete[] beta;
    delete[] Xbeta;

    *miu_ret = miu;
    *sigma2_ret = sigma2;
    *resid_ret = residvec;
    *XinvXTX_ret = XinvXTXvec;
}



void printCovVarMat(int numCovs, vector<string> covNames, double* covVarMat, double* beta, int phenoType, int samSize) 
{
    covNames.insert(covNames.begin(), "Intercept");
    boost::math::chi_squared chisq_dist_M(1);

    cout << "\nCoefficients: \n";
    cout << boost::format("%-26s %-17s %-22s %-19s %-15s\n") % "" % "Estimate" % "Std. Error" % "Z-value" % "P-value";
    for (int i = 0; i < numCovs; i++) 
    {
        double stdError = sqrt(covVarMat[i * numCovs + i]);
        double zvalue = beta[i] / stdError;
        double pr = (isnan(zvalue)) ? NAN : boost::math::cdf(complement(chisq_dist_M, (beta[i] * beta[i]) / covVarMat[i * numCovs + i]));
        cout << boost::format("%+15s %19.6e %19.6e %19.6e %19.6e\n") % covNames[i] % beta[i] % stdError % zvalue % pr;
    }

    cout << "\nVariance-Covariance Matrix: \n";
    cout << boost::format("%+35s") % covNames[0];
    for (int i = 1; i < numCovs; i++) {
        cout << boost::format("%+20s") % covNames[i];
    }
    cout << "\n";
    for (int i = 0; i < numCovs; i++) {
        cout << boost::format("%+15s") % covNames[i];
        for (int j = 0; j < numCovs; j++) {
            cout << boost::format("%20.6e") % covVarMat[j * numCovs + i];
        }
        cout << "\n";
    }
    cout << "\n";
}


void printOutputHeader(bool useBgen, int numExpSelCol, int Sq1, vector<string> covNames, string output, string outStyle, int robust, double sigma2, BinE binE) 
{
    std::ofstream results(output, std::ofstream::binary);

    bool printFull = false;
    bool printMeta = false;
    int printStart = 1; 
    int printEnd   = numExpSelCol+1; 
    if (outStyle.compare("meta") == 0) {
        printStart = 0; 
        printEnd   = Sq1;
        printMeta  = true;
    } else if (outStyle.compare("full") == 0) {
        printStart = 0; 
        printEnd   = Sq1; 
        printFull  = true;
        results << "#dispersion: " << sigma2 << "\n";
    }

    results << "SNPID" << ((useBgen) ? "\tRSID\t" : "\t") << "CHR" << "\t" << "POS" << "\t" << "Non_Effect_Allele" << "\t" << "Effect_Allele" << "\t" << "N_Samples" << "\t" << "AF" << "\t";
    int nBinE = binE.nBinE;
    if (nBinE > 0) {
        cout << "Here" << endl;
        vector<string> bn = binE.bn_header;
        for (size_t i = 0; i < bn.size(); i++) {
            results << "N_" << bn[i] << "\t";
            results << "AF_" << bn[i] << "\t";
        }
    }

    for (int i = 0; i < Sq1-1; i++) {
        covNames[i] = "G-" + covNames[i];
    }
    covNames.insert(covNames.begin(), "G");


    string seMHeader = "SE_Beta_Marginal";
    string seHeader  = "SE_Beta_";
    string covHeader = "Cov_Beta_";
    if (robust == 1) {
        seMHeader = "robust_" + seMHeader;
        seHeader  = "robust_" + seHeader;
        covHeader = "robust_" + covHeader;
    }

    results << "Beta_Marginal" << "\t" << seMHeader << "\t";
    if ((robust == 1) && (printMeta || printFull)) {
        results << "SE_Beta_Marginal" << "\t";
    }
    if (numExpSelCol != 0) {
        for (int i = printStart; i < printEnd; i++) {
            results << "Beta_" << covNames[i] << "\t";
        }
        for (int i = printStart; i < printEnd; i++) {
             results << seHeader << covNames[i] << "\t";  
        }
        for (int i = printStart; i < printEnd; i++) {
            for (int j = printStart; j < printEnd; j++) {
                if (i < j) {
                   results << covHeader << covNames[i] << "_" << covNames[j] << "\t";  
                } 
            }
        }
        if (robust == 1) {
            if (printMeta || printFull) {
                for (int i = printStart; i < printEnd; i++) {
                    for (int j = printStart; j < printEnd; j++) {
                        if (i == j) {
                            results << "SE_Beta_" << covNames[j] << "\t"; 
                        }
                    }
                }
                for (int i = printStart; i < printEnd; i++) {
                    for (int j = printStart; j < printEnd; j++) {
                        if (i < j) {
                            results << "Cov_Beta_" << covNames[i] << "_" << covNames[j] << "\t"; 
                        }
                    }
                }

                results << "robust_P_Value_Marginal" << "\t" << "robust_P_Value_Interaction" << "\t" << "robust_P_Value_Joint" << "\t";
                results << "P_Value_Marginal" << "\t" << "P_Value_Interaction" << "\t" << "P_Value_Joint\n";
            } else {
                results << "robust_P_Value_Marginal" << "\t" << "robust_P_Value_Interaction" << "\t" << "robust_P_Value_Joint\n";
            }
        } else {
                results << "P_Value_Marginal" << "\t" << "P_Value_Interaction" << "\t" << "P_Value_Joint\n";
        }
    }
    else {
        results << "P_Value_Marginal\n";
    }

    results.close();
}

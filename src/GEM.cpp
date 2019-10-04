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
  version Log: logistic regression for binary phenodata
  7/27/18: input column header name in pheno file not index number
  8/01/18: sample ID matching
  8/30/18: initial sample Size is not necessary equal
  8/30/18: hash table of genoUnMatchID for significant unmatching numbers;
  12/10/18: aritrary stream_snps implemented, omp parallel with > 20 covariates (HC)
  2/7/19: speed up reading genotype data without looking up in genoUnMatchID

  To-Do List:
  1. OOP
*/

#include "declars.h"

int main(int argc, char *argv[]) {
  char paramfile[300], genofile[300], phenofile[300], samplefile[300];
  PARAMETERS prm;

  /********************************
    Parse command line
  ********************************/
  {
    cout << "*********************************************************\n";
    cout << "Number of command-line arguments: " << argc << '\n';
    if (argc != 2) {
      cout << "Only one argument for the parameter input file name is a must. \n";
      exit(1);
    }
    if (argv[1][0] == '-') {
      cout << "Usage options are not provided. All settings can be found in the parameter input file. \n";
      exit(1);
    }
    sscanf(argv[1], "%s", paramfile);
    cout << "Parameter input file is: " << paramfile << '\n';
  }

  /****************************************************
    Call subroutine to read parameters from a file
  ****************************************************/
  ReadParameters(paramfile, genofile, phenofile, samplefile, &prm);

  string genopath(genofile); 
  if (genopath.substr(genopath.length()-5,5) != ".bgen") {
    cout << genopath << " is not a bgen format. Only support bgen format. \n";
    exit(1);
  }

  /****************************************************
    Parameters
  ****************************************************/ 
//  int samSize = prm.samSize;
  int samSize;
  int phenoTyp = prm.phenoTyp;
  int phenoCol;
  string phenoHeaderName(prm.phenoHeader);
  int samIDCol;
  string samIDHeaderName(prm.samIDHeader);
  string output(prm.outputfile);
  int numSelCol = prm.covSelHeaders.size();
  vector<int> colSelVec(numSelCol);
  vector<string> covSelHeadersName;
  for (int i = 0; i < numSelCol; i++) covSelHeadersName.push_back(prm.covSelHeaders[i]);
  int robust = prm.robust;
//  int IDMatching = prm.IDMatching;

  int stream_snps = prm.stream_snps;
  int Sq = prm.Sq;

  double epsilon = prm.epsilon;

  double exetime = 0;

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
    colNames[headerName] = header_i;
    ++header_i;
  }
  phenoncols = colNames.size();
  if (colNames.find(phenoHeaderName) == colNames.end()) {
    cerr << "Pheno header name is wrong.\n";
    exit(1);
  } else 
    phenoCol = colNames[phenoHeaderName];
  if (colNames.find(samIDHeaderName) == colNames.end()) {
    cerr << "Sample ID header name is wrong.\n";
    exit(1);
  } else
    samIDCol = colNames[samIDHeaderName];
  for (int i = 0; i < numSelCol; i++) {
    if (colNames.find(covSelHeadersName[i]) == colNames.end()) {
      cerr << "Covariate header name is wrong.\n";
      exit(i);
    } else
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
  cout << "Size of the selected covariate matrix (including first column for interception values) is: " << samSize << " X " << numSelCol+1 << '\n';

  // initialize data matrix
  vector <double> phenodata(samSize);
  vector <string> sampleIds(samSize);
  vector <double> covdata(samSize * (numSelCol+1));
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
    covdata[r*(numSelCol+1)] = 1.0;
    for (int c = 0; c < numSelCol; c++) 
      sscanf(values[colSelVec[c]].c_str(), "%lf", &covdata[r*(numSelCol+1) + c + 1]);

    // IDMatching
//    if (IDMatching == 1) {
    phenomap[values[samIDCol]] = {values[phenoCol]};
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
  FILE *fin = fopen(genofile, "rb");
  cout << "General information of bgen file. \n";
  uint offset; fread(&offset, 4, 1, fin); cout << "offset: " << offset << '\n';
  uint L_H; fread(&L_H, 4, 1, fin); cout << "L_H: \n";
  uint Mbgen; fread(&Mbgen, 4, 1, fin); cout << "BGEN snpBlocks (Mbgen): " << Mbgen << '\n';
  assert(Mbgen != 0);
  uint Nbgen; fread(&Nbgen, 4, 1, fin); cout << "BGEN samples (Nbgen): " << Nbgen << '\n';
//  if (Nbgen != samSize) {
//    cerr << "ERROR: Number of samples in BGEN header does not match sample file" << '\n';
//    exit(1);
//  }

  char magic[5]; fread(magic, 1, 4, fin); magic[4] = '\0'; //cout << "magic bytes: " << string(magic) << endl;
  fseek(fin, L_H-20, SEEK_CUR); //cout << "skipping L_H-20 = " << L_H-20 << " bytes (free data area)" << endl;
  uint flags; fread(&flags, 4, 1, fin); //cout << "flags: " << flags << endl;
  uint CompressedSNPBlocks = flags&3; cout << "CompressedSNPBlocks: " << CompressedSNPBlocks << '\n';
  assert(CompressedSNPBlocks==1); // REQUIRE CompressedSNPBlocks==1
  uint Layout = (flags>>2)&0xf; cout << "Layout: " << Layout << '\n';
  assert(Layout==1 || Layout==2); // REQUIRE Layout==1 or Layout==2
  uint SampleIdentifiers = flags>>31; cout << "SampleIdentifiers: " << SampleIdentifiers << '\n';
  //  char snpID[65536], rsID[65536], chrStr[65536];
  uint maxLA = 65536, maxLB = 65536;
  char* snpID = new char[maxLA+1];
  char* rsID = new char[maxLA+1];
  char* chrStr = new char[maxLA+1];
  char* allele1 = new char[maxLA+1];
  char* allele0 = new char[maxLB+1];
  char* samID = new char[maxLA+1];
  vector <uchar> zBuf;
  vector <uchar> shortBuf;

  /**** sample identifier block ********/
  /* exists when SampleIdentifiers = 1 */
  if (SampleIdentifiers == 1) {
    uint LS1; fread(&LS1, 4, 1, fin); // std::cout << "LS1: " << LS1 << std::endl; // LS1 + L_H <= offset
    uint Nrow; fread(&Nrow, 4, 1, fin); // cout << "Nrow: " << Nrow << " " << std::flush;
    if (Nrow != Nbgen) {
      cerr << "ERROR: Nrow = " << Nrow << " does not match Nbgen = " << Nbgen << '\n';
        exit(1);
    }
 
//    if (IDMatching == 0) {
//      for (uint m = 0; m < Nbgen; m++) {
//        ushort LSID; fread(&LSID, 2, 1, fin); // std::cout << "LSID: " << LSID << " ";
//        fread(samID, 1, LSID, fin); // std::cout << "samID: " << samID << " " << std::endl;
//      }
//    } 
//    else if (IDMatching == 1) {
      // checking sample ID matching
      int k = 0;
      for (uint m = 0; m < Nbgen; m++) {
        ushort LSID; fread(&LSID, 2, 1, fin); // std::cout << "LSID: " << LSID << " ";
        fread(samID, 1, LSID, fin); // std::cout << "samID: " << samID << " " << std::endl;

        // IDMatching
        string strtmp(samID);
        int itmp = k;
        if (phenomap.find(strtmp) != phenomap.end()) {
          auto tmp_valvec = phenomap[strtmp];
          if (find(tmp_valvec.begin(), tmp_valvec.end(), phenoMissingKey) == tmp_valvec.end()) {
            sscanf(tmp_valvec[0].c_str(), "%lf", &phenodata[k]);
            // covdata[k*(numSelCol+1) + 0] = 1.0;
	    for (int c = 0; c < numSelCol; c++)
              sscanf(tmp_valvec[c + 1].c_str(), "%lf", &covdata[k*(numSelCol+1) + c + 1]);
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
      covdata.resize(k*(numSelCol + 1));
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
    for (uint m = 0; m < Nbgen; m++) {
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
            sscanf(tmp_valvec[c + 1].c_str(), "%lf", &covdata[k*(numSelCol+1) + c + 1]);
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
      covdata.resize(k*(numSelCol + 1));
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
  double* covX = &covdata[0];
  vector <double> residvec(samSize);

  // for logistic regression
  vector <double> miu(samSize);
  int Check = 1; // convergence condition of beta^(i+1) - beta^(i)
  int iter = 1;

  // write to file
  std::ofstream results(output, std::ofstream::binary);
  results << "SNPID" << "\t" << "rsID" << "\t" << "CHR" << "\t" << "POS" << "\t" << "Allele1" << "\t" << "Allele2" << "\t" << "AF" << "\t" << "Beta_Main" << "\t" << "Var_Beta_Main" << "\t";
  for (int i = 1; i <= Sq; i++) 
    results << "Beta_Interaction" << "_" << i << "\t"; 
  for (int i = 1; i <= Sq; i++)
    for (int j = 1; j <= Sq; j++)
      results << "Var_Beta_Interaction" << "_" << i << "_" << j << "\t";
  results << "P_Value_Main" << "\t" << "P_Value_Interaction" << "\t" << "P_Value_Joint\n";

  // transpose(X) * X
  double* XTransX = new double[(numSelCol+1)*(numSelCol+1)];  
  matTmatprod(covX, covX, XTransX, samSize, numSelCol+1, numSelCol+1);
  // invert (XTransX)
  matInv(XTransX, numSelCol+1);
  // transpose(X) * Y
  double* XTransY = new double[(numSelCol+1)];
  matTvecprod(covX, phenoY, XTransY, samSize, numSelCol+1);
  // beta = invert(XTransX) * XTransY
  double* beta = new double[(numSelCol+1)];
  matvecprod(XTransX, XTransY, beta, numSelCol+1, numSelCol+1);

  while ((phenoTyp == 1) && (Check != (numSelCol+1))) { // logistic regression
    iter++;
    // X * beta
    double* XbetaFL = new double[samSize];
    matvecprod(covX, beta, XbetaFL, samSize, numSelCol+1);
    double* Yip1 = new double[samSize];	
    // W * X and W * Y
    double* WX = new double[samSize*(numSelCol+1)];
    double* WYip1 = new double[samSize];
    for (int i = 0; i < samSize; i++) {
      miu[i] = exp(XbetaFL[i])/(1.0+exp(XbetaFL[i]));
      Yip1[i] = XbetaFL[i] + (phenoY[i] - miu[i])/(miu[i]*(1-miu[i]));
      WYip1[i] = miu[i]*(1-miu[i])*Yip1[i]; 
      for (int j = 0; j < numSelCol+1; j++) {
	WX[i*(numSelCol+1)+j] = miu[i]*(1-miu[i])*covX[i*(numSelCol+1)+j];
      }
    } 
    // transpose(X) * WX
    matTmatprod(covX, WX, XTransX, samSize, numSelCol+1, numSelCol+1);
    // invert (XTransX)
    matInv(XTransX, numSelCol+1);
    // transpose(X) * WYip1
    matTvecprod(covX, WYip1, XTransY, samSize, numSelCol+1);
    // beta = invert(XTransX) * XTransY
    double* betaT = new double[(numSelCol+1)];
    matvecprod(XTransX, XTransY, betaT, numSelCol+1, numSelCol+1);
    Check = 0;
    for (int i = 0; i < numSelCol+1; i++) {
      if (std::abs(betaT[i] - beta[i]) <= epsilon) Check++;
      beta[i] = betaT[i];
    }

    delete [] Yip1;
    delete [] WYip1;
    delete [] WX;
    delete [] XbetaFL;
    delete [] betaT;
  }

  // X * beta
  double* Xbeta = new double[samSize];
  matvecprod(covX, beta, Xbeta, samSize, numSelCol+1);
  if (phenoTyp == 1) {
    double* WX = new double[samSize*(numSelCol+1)];
    for (int i = 0; i < samSize; i++) {
      miu[i] = exp(Xbeta[i])/(1.0+exp(Xbeta[i]));
      Xbeta[i] = miu[i];
      for (int j = 0; j < numSelCol+1; j++) {
        WX[i*(numSelCol+1)+j] = miu[i]*(1.0-miu[i])*covX[i*(numSelCol+1)+j];
      } 
    }
    // transpose(X) * WX
    matTmatprod(covX, WX, XTransX, samSize, numSelCol+1, numSelCol+1);
    // invert (XTransX)
    matInv(XTransX, numSelCol+1);
    delete [] WX;

    cout << "Logistic regression reaches convergence after " << iter << " steps.\n";
    cout << "*********************************************************\n";
  }

  // X*[invert (XTransX)]
  double* XinvXTX = new double[samSize*(numSelCol+1)];
  matmatprod(covX, XTransX, XinvXTX, samSize, numSelCol+1, numSelCol+1);
  // residual = Y - X * beta
  double sigma2 = 0;
  for (int i = 0; i < samSize; i++) {
    residvec[i] = phenoY[i] - Xbeta[i];
    sigma2 += residvec[i]*residvec[i];
  }
  // sqr(sigma) = transpose(resid)*resid/[samSize-(numSelCol+1)]
  sigma2 = sigma2/(samSize-(numSelCol+1));
  if (phenoTyp == 1) sigma2 = 1.0;

  double* resid = &residvec[0];

  delete [] XTransY;
  delete [] beta;
  delete [] Xbeta;

  cout << "Streaming SNPs for speeding up GWAS analysis in parallel. \n";
  cout << "Number of SNPs in each batch is: " << stream_snps << '\n';
  cout << "*********************************************************\n";
  vector <double> ZGSvec(samSize * (1+Sq)*stream_snps);
  vector <double> ZGSR2vec(samSize * (1+Sq)*stream_snps);
  // for logistic regression
  vector <double> WZGSvec(samSize * (1+Sq)*stream_snps);
  double* WZGS = &WZGSvec[0];

  // for output unique snpids in genotype file
  vector <string> geno_snpid(stream_snps);

  fseek(fin, offset + 4, SEEK_SET);
  int stream_snps1 = stream_snps;

  uint* include_idx = new uint[samSize];
  int ii = 0;
  for (uint i = 0; i < Nbgen; i++) {
    if (genoUnMatchID.find(i) == genoUnMatchID.end()) {
      include_idx[ii] = i;
      ii++;
    }
  }
  genoUnMatchID.clear();

  for (int snploop = 0; snploop*stream_snps1 < Mbgen; snploop++) {
    /*    if ((Mbgen % stream_snps) != 0) {
      cerr << "ERROR: total snps is not divisible by streaming snp numbers ( " << stream_snps << " ). \n";
      exit(1);
      } */
    if((snploop + 1) * stream_snps1 >= Mbgen) {
      stream_snps = Mbgen - snploop * stream_snps1;
    }

    int Sq1 = Sq+1;
    int ZGS_col = Sq1*stream_snps;
    vector <double> AF(stream_snps);

    for (int stream_i = 0; stream_i < stream_snps; stream_i++) {
    /**** variant data block ********/
      if (Layout == 1) {
        uint Nrow; fread(&Nrow, 4, 1, fin); // cout << "Nrow: " << Nrow << " " << std::flush;
        if (Nrow != Nbgen) {
          cerr << "ERROR: Nrow = " << Nrow << " does not match Nbgen = " << Nbgen << '\n';
          exit(1);
        }
      }
      ushort LS; fread(&LS, 2, 1, fin);  // cout << "LS: " << LS << " " << std::flush;
      if (LS > maxLA) {
        maxLA = 2*LS;
        delete [] snpID;
        char* snpID = new char[maxLA+1];
      }
      fread(snpID, 1, LS, fin); snpID[LS] = '\0'; // cout << "snpID: " << string(snpID) << " " << std::flush;
      ushort LR; fread(&LR, 2, 1, fin); // cout << "LR: " << LR << " " << std::flush;
      if (LR > maxLA) {
        maxLA = 2*LR;
        delete [] rsID;
        char* rsID = new char[maxLA+1];
      }
      fread(rsID, 1, LR, fin); rsID[LR] = '\0'; // cout << "rsID: " << string(rsID) << " " << std::flush;
      ushort LC; fread(&LC, 2, 1, fin); // cout << "LC: " << LC << " " << std::flush;
      fread(chrStr, 1, LC, fin); chrStr[LC] = '\0';
      uint physpos; fread(&physpos, 4, 1, fin); // cout << "physpos: " << physpos << " " << std::flush;

      if (Layout == 2) {
        ushort LKnum; fread(&LKnum, 2, 1, fin); // this is for Layout = 2, Lnum = 2 is Layout = 1
        if (LKnum != 2) {
          cerr << "ERROR: Non-bi-allelic variant found: " << LKnum << " alleles\n";
          exit(1);
        }
      }
      uint LA; fread(&LA, 4, 1, fin); // cout << "LA: " << LA << " " << std::flush;
      if (LA > maxLA) {
        maxLA = 2*LA;
        delete [] allele1;
        char* allele1 = new char[maxLA+1];
      }
      fread(allele1, 1, LA, fin); allele1[LA] = '\0';
      uint LB; fread(&LB, 4, 1, fin); // cout << "LB: " << LB << " " << std::flush;
      if (LB > maxLB) {
        maxLB = 2*LB;
        delete [] allele0;
        char* allele0 = new char[maxLB+1];
      }
      fread(allele0, 1, LB, fin); allele0[LB] = '\0';

      geno_snpid[stream_i] = string(snpID) + "\t" + string(rsID) + "\t" + string(chrStr) + "\t" + std::to_string(physpos) + "\t" + string(allele1) + "\t" + string(allele0);

      if (Layout == 1) {
        uint zLen; fread(&zLen, 4, 1, fin); // cout << "zLen: " << zLen << endl;
	fread(&zBuf[0], 1, zLen, fin);
	uLongf destLen = 6*Nbgen;
        if (uncompress(&shortBuf[0], &destLen, &zBuf[0], zLen) != Z_OK || destLen != 6*Nbgen) {
          cerr << "ERROR: uncompress() failed\n";
          exit(1);
        }
        // read genotype probabilities
        double sum_eij = 0, sum_fij_minus_eij2 = 0;
        const double scale = 1.0/32768;
        int tmp1 = stream_i*Sq1*samSize;
//        if (IDMatching == 1) {
          int k = 0;
          for (uint i = 0; i < Nbgen; i++) {

//	    if (find(genoUnMatchID.begin(), genoUnMatchID.end(), i) == genoUnMatchID.end()) {
//          if (genoUnMatchID.find(i) == genoUnMatchID.end()) {
	    if (include_idx[k] == i) {
	      double p11 = shortBuf[3*i] * scale;
	      double p10 = shortBuf[3*i+1] * scale;
	      double p00 = shortBuf[3*i+2] * scale;

	      double pTot = p11 + p10 + p00;
	      double dosage = (2*p00 + p10) / pTot;

              int tmp2 = k+tmp1;
	      AF[stream_i] += dosage;
              if (phenoTyp == 1) 
		ZGSvec[tmp2] = miu[k]*(1-miu[k])*dosage;
	      else
		ZGSvec[tmp2] = dosage;
              k++;
            }
          }
//        }
/*        else if (IDMatching == 0) {
          for (uint i = 0; i < Nbgen; i++) {
            double p11 = shortBuf[3*i] * scale;
            double p10 = shortBuf[3*i+1] * scale;
            double p00 = shortBuf[3*i+2] * scale;

            double pTot = p11 + p10 + p00;
            double dosage = (2*p00 + p10) / pTot;

            int tmp2 = i+tmp1;
            ZGSvec[tmp2] = dosage;
            if (phenoTyp == 1) ZGSvec[tmp2] = miu[i]*(1-miu[i])*dosage;
          }
        }
*/
        for (int j = 0; j < Sq; j++) {
          int tmp3 = samSize*(j+1)+tmp1;
          for (uint i = 0; i < samSize; i++) {
            int tmp4 = i*(numSelCol+1);
            ZGSvec[tmp3+i] = covX[tmp4+j+1]*ZGSvec[tmp1+i]; // here we save ZGS in column wise
          }
        }
      } // end of reading genotype data when Layout = 1

      if (Layout == 2) {
        uint zLen; fread(&zLen, 4, 1, fin); // cout << "zLen: " << zLen << endl;
        uint DLen;
        zBuf.resize(zLen-4);
        if (CompressedSNPBlocks == 0) {
          DLen = zLen;
          fread(&zBuf[0], 1, zLen, fin);
        } else {
          fread(&DLen, 4, 1, fin);
          fread(&zBuf[0], 1, zLen-4, fin);
        }
        uLongf destLen = DLen; //6*Nbgen;
        shortBuf.resize(DLen);

        if (uncompress(&shortBuf[0], &destLen, &zBuf[0], zLen-4) != Z_OK || destLen != DLen) {
          cout << "destLen: " << destLen << " " << zLen-4 << '\n';
          cerr << "ERROR: uncompress() failed\n";
          exit(1);
        }
        // read genotype probabilities
        uchar *bufAt = &shortBuf[0];
        uint N = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16)|(bufAt[3]<<24); bufAt += 4;
        if (N != Nbgen) {
          cerr << "ERROR: " << "snpName " << " has N = " << N << " (mismatch with header block)\n";
          exit(1);
        }
        uint K = bufAt[0]|(bufAt[1]<<8); bufAt += 2;
        if (K != 2U) {
          cerr << "ERROR: " << "snpName " << " has K = " << K << " (non-bi-allelic)\n";
          exit(1);
        }
        uint Pmin = *bufAt; bufAt++;
        if (Pmin != 2U) {
          cerr << "ERROR: " << "snpName " << " has minimum ploidy = " << Pmin << " (not 2)\n";
          exit(1);
        }
        uint Pmax = *bufAt; bufAt++;
        if (Pmax != 2U) {
          cerr << "ERROR: " << "snpName " << " has maximum ploidy = " << Pmax << " (not 2)\n";
          exit(1);
        }
        for (uint i = 0; i < N; i++) {
          uint ploidyMiss = *bufAt; bufAt++;
          if (ploidyMiss != 2U) {
          //  std::cerr << "ERROR: " << "snpName " << " has ploidy/missingness byte = " << ploidyMiss
          //       << " (not 2)" << endl;
          //  exit(1);
          }
        }
        uint Phased = *bufAt; bufAt++;
        if (Phased != 0U) {
          cerr << "ERROR: " << "snpName " << " has Phased = " << Phased << " (not 0)\n";
          exit(1);
        }
        uint B = *bufAt; bufAt++;
        uint Bbits = std::pow(2,B);
        if ((B != 8U) && (B != 16U) && (B != 24U) && (B != 32U)) {
          std::cerr << "ERROR: " << "snpName " << " has B = " << B << " (not divisible by 8)\n";
          exit(1);
        }

        int tmp1 = stream_i*Sq1*samSize;

//	if (IDMatching == 1) {
          int k = 0;
          for (uint i = 0; i < N; i++) {
    	    uint chartem;
	    if (B == 8U)
              chartem = bufAt[0];
            else if (B == 16U)
              chartem = bufAt[0]|(bufAt[1]<<8);
            else if (B == 24U)
              chartem = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16);
            else if (B == 32U)
              chartem = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16)|(bufAt[3]<<24);
	    bufAt += B/8;
	    uint chartem1;
            if (B == 8U)
              chartem1 = bufAt[0];
            else if (B == 16U)
              chartem1 = bufAt[0]|(bufAt[1]<<8);
            else if (B == 24U)
              chartem1 = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16);
            else if (B == 32U)
              chartem1 = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16)|(bufAt[3]<<24);
            bufAt += B/8;

//            if (find(genoUnMatchID.begin(), genoUnMatchID.end(), i) == genoUnMatchID.end()) {
//          if (genoUnMatchID.find(i) == genoUnMatchID.end()) {
	    if (include_idx[k] == i) {
	      double p11 = chartem/double(1.0*(Bbits-1));
	      double p10 = chartem1/double(1.0*(Bbits-1));
	      double dosage = 2*(1-p11-p10) + p10;

              int tmp2 = k+tmp1;
              AF[stream_i] += dosage;
              if (phenoTyp == 1) 
		ZGSvec[tmp2] = miu[k]*(1-miu[k])*dosage;
	      else
		ZGSvec[tmp2] = dosage; // replace your new data from other genotype files here.
	      k++;
	    }
          }
//	}
/*        else if (IDMatching == 0) {
          for (uint i = 0; i < N; i++) {
            uint chartem;
            if (B == 8U)
              chartem = bufAt[0];
            else if (B == 16U)
              chartem = bufAt[0]|(bufAt[1]<<8);
            else if (B == 24U)
              chartem = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16);
            else if (B == 32U)
              chartem = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16)|(bufAt[3]<<24);
            double p11 = chartem/double(1.0*(Bbits-1));
            bufAt += B/8;
            uint chartem1;
            if (B == 8U)
              chartem1 = bufAt[0];
            else if (B == 16U)
              chartem1 = bufAt[0]|(bufAt[1]<<8);
            else if (B == 24U)
              chartem1 = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16);
            else if (B == 32U)
              chartem1 = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16)|(bufAt[3]<<24);
            double p10 = chartem1/double(1.0*(Bbits-1));
            bufAt += B/8;
            double dosage = 2*(1-p11-p10) + p10;

            int tmp2 = i+tmp1;
            ZGSvec[tmp2] = dosage;
            if (phenoTyp == 1) ZGSvec[tmp2] = miu[i]*(1-miu[i])*dosage;
          }
        }
*/	
        for (int j = 0; j < Sq; j++) {
          int tmp3 = samSize*(j+1)+tmp1;
          for (uint i = 0; i < samSize; i++) {
            int tmp4 = i*(numSelCol+1);
            ZGSvec[tmp3+i] = covX[tmp4+j+1]*ZGSvec[tmp1+i]; // here we save ZGS in column wise
          }
        }
      }	// end of reading genotype data when Layout = 2
    } // end of stream_i

    /***************************************************************/
    // genodata and envirment data
    auto wallexe0 = std::chrono::system_clock::now();

    double* ZGS = &ZGSvec[0];

    // transpose(X) * ZGS
    // it is non-squred matrix, attention that continuous memory is column-major due to Fortran in BLAS.
    // important!!!!
    double* XtransZGS = new double[(numSelCol+1)*ZGS_col];
    matNmatNprod(covX, ZGS, XtransZGS, numSelCol+1, samSize, ZGS_col);

    if (phenoTyp == 0) {
      #pragma omp parallel // used when multi-thread is called
      for (int j = 0; j < ZGS_col; j++) {
        #pragma omp for nowait
        for (int i = 0; i < samSize; i++) {
          // notice that column number of X is 4
          // we do not write for loop is because it is within the openmp parallel for
          ZGS[j*samSize+i] -= XinvXTX[0*samSize+i]*XtransZGS[j*(numSelCol+1)+0];
          ZGS[j*samSize+i] -= XinvXTX[1*samSize+i]*XtransZGS[j*(numSelCol+1)+1];
	  if (numSelCol > 1) ZGS[j*samSize+i] -= XinvXTX[2*samSize+i]*XtransZGS[j*(numSelCol+1)+2];
          if (numSelCol > 2) ZGS[j*samSize+i] -= XinvXTX[3*samSize+i]*XtransZGS[j*(numSelCol+1)+3];
          if (numSelCol > 3) ZGS[j*samSize+i] -= XinvXTX[4*samSize+i]*XtransZGS[j*(numSelCol+1)+4];
          if (numSelCol > 4) ZGS[j*samSize+i] -= XinvXTX[5*samSize+i]*XtransZGS[j*(numSelCol+1)+5];
          if (numSelCol > 5) ZGS[j*samSize+i] -= XinvXTX[6*samSize+i]*XtransZGS[j*(numSelCol+1)+6];
          if (numSelCol > 6) ZGS[j*samSize+i] -= XinvXTX[7*samSize+i]*XtransZGS[j*(numSelCol+1)+7];
          if (numSelCol > 7) ZGS[j*samSize+i] -= XinvXTX[8*samSize+i]*XtransZGS[j*(numSelCol+1)+8];
          if (numSelCol > 8) ZGS[j*samSize+i] -= XinvXTX[9*samSize+i]*XtransZGS[j*(numSelCol+1)+9];
          if (numSelCol > 9) ZGS[j*samSize+i] -= XinvXTX[10*samSize+i]*XtransZGS[j*(numSelCol+1)+10];
          if (numSelCol > 10) ZGS[j*samSize+i] -= XinvXTX[11*samSize+i]*XtransZGS[j*(numSelCol+1)+11];
          if (numSelCol > 11) ZGS[j*samSize+i] -= XinvXTX[12*samSize+i]*XtransZGS[j*(numSelCol+1)+12];
          if (numSelCol > 12) ZGS[j*samSize+i] -= XinvXTX[13*samSize+i]*XtransZGS[j*(numSelCol+1)+13];
          if (numSelCol > 13) ZGS[j*samSize+i] -= XinvXTX[14*samSize+i]*XtransZGS[j*(numSelCol+1)+14];
          if (numSelCol > 14) ZGS[j*samSize+i] -= XinvXTX[15*samSize+i]*XtransZGS[j*(numSelCol+1)+15];
          if (numSelCol > 15) ZGS[j*samSize+i] -= XinvXTX[16*samSize+i]*XtransZGS[j*(numSelCol+1)+16];
          if (numSelCol > 16) ZGS[j*samSize+i] -= XinvXTX[17*samSize+i]*XtransZGS[j*(numSelCol+1)+17];
          if (numSelCol > 17) ZGS[j*samSize+i] -= XinvXTX[18*samSize+i]*XtransZGS[j*(numSelCol+1)+18];
          if (numSelCol > 18) ZGS[j*samSize+i] -= XinvXTX[19*samSize+i]*XtransZGS[j*(numSelCol+1)+19];
          if (numSelCol > 19) ZGS[j*samSize+i] -= XinvXTX[20*samSize+i]*XtransZGS[j*(numSelCol+1)+20];
	  if (numSelCol > 20)
	    for (int k = 21; k <= numSelCol; k++) ZGS[j*samSize+i] -= XinvXTX[k*samSize+i]*XtransZGS[j*(numSelCol+1)+k]; 
          if (robust == 1) ZGSR2vec[j*samSize+i] = ZGS[j*samSize+i]*resid[i]*resid[i];
        }
      }
    } 
    else if (phenoTyp == 1) {
      #pragma omp parallel // used when multi-thread is called
      for (int j = 0; j < ZGS_col; j++) {
        #pragma omp for nowait
        for (int i = 0; i < samSize; i++) {
	  double ZGStemp = 0.0;
          ZGStemp += ZGS[j*samSize+i]/miu[i]/(1.0-miu[i]) - XinvXTX[0*samSize+i]*XtransZGS[j*(numSelCol+1)+0];
          ZGStemp -= XinvXTX[1*samSize+i]*XtransZGS[j*(numSelCol+1)+1];
          if (numSelCol > 1) ZGStemp -= XinvXTX[2*samSize+i]*XtransZGS[j*(numSelCol+1)+2];
          if (numSelCol > 2) ZGStemp -= XinvXTX[3*samSize+i]*XtransZGS[j*(numSelCol+1)+3];
          if (numSelCol > 3) ZGStemp -= XinvXTX[4*samSize+i]*XtransZGS[j*(numSelCol+1)+4];
          if (numSelCol > 4) ZGStemp -= XinvXTX[5*samSize+i]*XtransZGS[j*(numSelCol+1)+5];
          if (numSelCol > 5) ZGStemp -= XinvXTX[6*samSize+i]*XtransZGS[j*(numSelCol+1)+6];
          if (numSelCol > 6) ZGStemp -= XinvXTX[7*samSize+i]*XtransZGS[j*(numSelCol+1)+7];
          if (numSelCol > 7) ZGStemp -= XinvXTX[8*samSize+i]*XtransZGS[j*(numSelCol+1)+8];
          if (numSelCol > 8) ZGStemp -= XinvXTX[9*samSize+i]*XtransZGS[j*(numSelCol+1)+9];
          if (numSelCol > 9) ZGStemp -= XinvXTX[10*samSize+i]*XtransZGS[j*(numSelCol+1)+10];
          if (numSelCol > 10) ZGStemp -= XinvXTX[11*samSize+i]*XtransZGS[j*(numSelCol+1)+11];
          if (numSelCol > 11) ZGStemp -= XinvXTX[12*samSize+i]*XtransZGS[j*(numSelCol+1)+12];
          if (numSelCol > 12) ZGStemp -= XinvXTX[13*samSize+i]*XtransZGS[j*(numSelCol+1)+13];
          if (numSelCol > 13) ZGStemp -= XinvXTX[14*samSize+i]*XtransZGS[j*(numSelCol+1)+14];
          if (numSelCol > 14) ZGStemp -= XinvXTX[15*samSize+i]*XtransZGS[j*(numSelCol+1)+15];
          if (numSelCol > 15) ZGStemp -= XinvXTX[16*samSize+i]*XtransZGS[j*(numSelCol+1)+16];
          if (numSelCol > 16) ZGStemp -= XinvXTX[17*samSize+i]*XtransZGS[j*(numSelCol+1)+17];
          if (numSelCol > 17) ZGStemp -= XinvXTX[18*samSize+i]*XtransZGS[j*(numSelCol+1)+18];
          if (numSelCol > 18) ZGStemp -= XinvXTX[19*samSize+i]*XtransZGS[j*(numSelCol+1)+19];
          if (numSelCol > 19) ZGStemp -= XinvXTX[20*samSize+i]*XtransZGS[j*(numSelCol+1)+20];
	  if (numSelCol > 20)
	    for (int k = 21; k <= numSelCol; k++) ZGStemp -= XinvXTX[k*samSize+i]*XtransZGS[j*(numSelCol+1)+k]; 
	  ZGS[j*samSize+i] = ZGStemp; 
          if (robust == 1) ZGSR2vec[j*samSize+i] = ZGS[j*samSize+i]*resid[i]*resid[i];
          WZGS[j*samSize+i] = miu[i]*(1-miu[i])*ZGS[j*samSize+i];
        }
      }
    }

    double* ZGSR2 = &ZGSR2vec[0];

    delete [] XtransZGS;

    // transpose(ZGS) * resid
    double* ZGStR = new double[ZGS_col];
    matvecprod(ZGS, resid, ZGStR, ZGS_col, samSize);      
    // transpose(ZGS) * ZGS
    double* ZGStZGS = new double[ZGS_col*ZGS_col];
    if (phenoTyp == 0) {
      matmatTprod(ZGS, ZGS, ZGStZGS, ZGS_col, samSize, ZGS_col);
    } 
    else if (phenoTyp == 1) {
      matmatTprod(ZGS, WZGS, ZGStZGS, ZGS_col, samSize, ZGS_col);
    } 
    else {
      cout << "phenoTyp is not equal to 0 or 1. Kill the job!! \n";
      exit(1);
    }
    // transpose(ZGSR2) * ZGS
    double* ZGSR2tZGS = new double[ZGS_col*ZGS_col];
    if (robust == 1) matmatTprod(ZGSR2, ZGS, ZGSR2tZGS, ZGS_col, samSize, ZGS_col);

    double* betaM = new double[stream_snps];
    double* VarbetaM = new double[stream_snps];
    double** betaInt = new double*[stream_snps];
    double** VarbetaInt = new double*[stream_snps];
    double* PvalM = new double[stream_snps];
    double* PvalInt = new double[stream_snps];
    double* PvalJoint = new double[stream_snps];
    boost::math::chi_squared chisq_dist_M(1);
    boost::math::chi_squared chisq_dist_Int(Sq);
    boost::math::chi_squared chisq_dist_Joint(1+Sq);

    if (robust == 0) {
      for (int i = 0; i < stream_snps; i++) {
        // initialize dynamic 2D array
        betaInt[i] = new double[Sq];
        VarbetaInt[i] = new double[Sq*Sq];

        // betamain
        int tmp1 = i*ZGS_col*Sq1+i*Sq1;
        betaM[i] = ZGStR[i*Sq1]/ZGStZGS[tmp1];
	VarbetaM[i] = sigma2/ZGStZGS[tmp1];

        double* S2TransS2 = new double[Sq*Sq];
        double* S2TransR = new double[Sq];
        double* S2DS2 = new double[Sq*Sq];
	double* InvVarbetaint = new double[Sq*Sq];
        for (int ind1 = 0; ind1 < Sq; ind1++) {
          for (int ind2 = 0; ind2 < Sq; ind2++) {
            // transpose(Snew2) * Snew2
            S2TransS2[ind1*Sq+ind2] = ZGStZGS[tmp1+(ind1+1)*ZGS_col+ind2+1]-ZGStZGS[tmp1+(ind1+1)*ZGS_col]*ZGStZGS[tmp1+ind2+1]/ZGStZGS[tmp1];
	  }
          // transpose(Snew2) * resid
          S2TransR[ind1] = ZGStR[i*Sq1+ind1+1]-ZGStZGS[tmp1+(ind1+1)*ZGS_col]*ZGStR[i*Sq1]/ZGStZGS[tmp1];
        }

        // invert (S2TransS2)
        matInv(S2TransS2, Sq);

        // betaInt = invert(S2TransS2) * S2TransR
        matvecprod(S2TransS2, S2TransR, betaInt[i], Sq, Sq);

        // Inv(S2TransS2) * S2DS2
        double* Stemp2 = new double[Sq*Sq];

        for (int j = 0; j < Sq; j++) {
          for (int k = 0; k < Sq; k++) {
            VarbetaInt[i][j*Sq+k] = sigma2*S2TransS2[j*Sq+k];
	    InvVarbetaint[j*Sq+k] = VarbetaInt[i][j*Sq+k]; 
          }
        }

        // calculating P values
        double statM = betaM[i]*betaM[i]/VarbetaM[i];
        if (isnan(statM) || statM <= 0.0) PvalM[i] = NAN;
        else PvalM[i] = boost::math::cdf(complement(chisq_dist_M, statM));
        // invert VarbetaInt[i]
        matInv(InvVarbetaint, Sq);
        double* Stemp3 = new double[Sq];
        matvecprod(InvVarbetaint, betaInt[i], Stemp3, Sq, Sq);
        double statInt = 0.0;
        for (int j = 0; j < Sq; j++) statInt += betaInt[i][j]*Stemp3[j];
        if (isnan(statInt) || statInt <= 0.0) PvalInt[i] = NAN;
        else PvalInt[i] = boost::math::cdf(complement(chisq_dist_Int, statInt));
        double statJoint = statM + statInt;
        if (isnan(statJoint) || statJoint <= 0.0) PvalJoint[i] = NAN;
        else PvalJoint[i] = boost::math::cdf(complement(chisq_dist_Joint, statJoint));

        delete [] S2TransS2;
        delete [] S2TransR;
        delete [] S2DS2;
        delete [] Stemp2;
        delete [] Stemp3;
	delete [] InvVarbetaint;
      }
    } 
    else if (robust == 1) {
      for (int i = 0; i < stream_snps; i++) {
        // initialize dynamic 2D array
        betaInt[i] = new double[Sq];
        VarbetaInt[i] = new double[Sq*Sq];

        // betamain
        int tmp1 = i*ZGS_col*Sq1+i*Sq1;
        betaM[i] = ZGStR[i*Sq1]/ZGStZGS[tmp1];
        VarbetaM[i] = ZGSR2tZGS[tmp1]/(ZGStZGS[tmp1]*ZGStZGS[tmp1]);

        double* S2TransS2 = new double[Sq*Sq];
        double* S2TransR = new double[Sq];
        double* S2DS2 = new double[Sq*Sq];
        double* InvVarbetaint = new double[Sq*Sq];
        for (int ind1 = 0; ind1 < Sq; ind1++) {
          for (int ind2 = 0; ind2 < Sq; ind2++) {
            // transpose(Snew2) * Snew2
            S2TransS2[ind1*Sq+ind2] = ZGStZGS[tmp1+(ind1+1)*ZGS_col+ind2+1]-ZGStZGS[tmp1+(ind1+1)*ZGS_col]*ZGStZGS[tmp1+ind2+1]/ZGStZGS[tmp1];
            // transpose(Snew2) * D * Snew2
            S2DS2[ind1*Sq+ind2] = ZGSR2tZGS[tmp1+(ind1+1)*ZGS_col+ind2+1]-ZGStZGS[tmp1+(ind1+1)*ZGS_col]*ZGSR2tZGS[tmp1+ind2+1]/ZGStZGS[tmp1]-ZGSR2tZGS[tmp1+(ind1+1)*ZGS_col]*ZGStZGS[tmp1+ind2+1]/ZGStZGS[tmp1]+ZGStZGS[tmp1+(ind1+1)*ZGS_col]*ZGSR2tZGS[tmp1]*ZGStZGS[tmp1+ind2+1]/ZGStZGS[tmp1]/ZGStZGS[tmp1];
          }
          // transpose(Snew2) * resid
          S2TransR[ind1] = ZGStR[i*Sq1+ind1+1]-ZGStZGS[tmp1+(ind1+1)*ZGS_col]*ZGStR[i*Sq1]/ZGStZGS[tmp1];
        }

        // invert (S2TransS2)
        matInv(S2TransS2, Sq);

        // betaInt = invert(S2TransS2) * S2TransR
        matvecprod(S2TransS2, S2TransR, betaInt[i], Sq, Sq);

        // Inv(S2TransS2) * S2DS2
        double* Stemp2 = new double[Sq*Sq];
        matmatprod(S2TransS2, S2DS2, Stemp2, Sq, Sq, Sq);
        // Stemp2 * Inv(S2TransS2)
        matNmatNprod(Stemp2, S2TransS2, VarbetaInt[i], Sq, Sq, Sq);

        for (int j = 0; j < Sq; j++) {
          for (int k = 0; k < Sq; k++) {
            InvVarbetaint[j*Sq+k] = VarbetaInt[i][j*Sq+k];
          }
        }

        // calculating P values
        double statM = betaM[i]*betaM[i]/VarbetaM[i];
        if (isnan(statM) || statM <= 0.0)  PvalM[i] = NAN;
        else
          PvalM[i] = boost::math::cdf(complement(chisq_dist_M, statM));
        // invert VarbetaInt[i]
        matInv(InvVarbetaint, Sq);
        double* Stemp3 = new double[Sq];
        matvecprod(InvVarbetaint, betaInt[i], Stemp3, Sq, Sq);
        double statInt = 0.0;
        for (int j = 0; j < Sq; j++) statInt += betaInt[i][j]*Stemp3[j];
        if (isnan(statInt) || statInt <= 0.0)  PvalInt[i] = NAN;
        else
          PvalInt[i] = boost::math::cdf(complement(chisq_dist_Int, statInt));
        double statJoint = statM + statInt;
        if (isnan(statJoint) || statJoint <= 0.0)  PvalJoint[i] = NAN;
        else
          PvalJoint[i] = boost::math::cdf(complement(chisq_dist_Joint, statJoint));

        delete [] S2TransS2;
        delete [] S2TransR;
        delete [] S2DS2;
        delete [] Stemp2;
        delete [] Stemp3;
        delete [] InvVarbetaint;
      }
    } // end of if robust == 1

    delete [] ZGStR;
    delete [] ZGStZGS;
    delete [] ZGSR2tZGS;

    auto wallexe1 = std::chrono::system_clock::now();
    std::chrono::duration<double> wallexeduration = (wallexe1 - wallexe0);
    exetime += wallexeduration.count();

    for (int i = 0; i < stream_snps; i++) {
      results << geno_snpid[i] << "\t" << AF[i]/2/samSize << "\t" << betaM[i] << "\t" << VarbetaM[i] << "\t";
      for (int ii = 0; ii < Sq; ii++)
	results << betaInt[i][ii] << "\t";
      for (int ii = 0; ii < Sq; ii++)
	for (int jj = 0; jj < Sq; jj++)
	  results << VarbetaInt[i][ii*Sq+jj] << "\t";
      results << PvalM[i] << "\t" << PvalInt[i] << "\t" << PvalJoint[i] << '\n';
    }

    delete [] betaM;
    delete [] VarbetaM;
    for(int i = 0; i < stream_snps; i++) {
      delete [] betaInt[i];
      delete [] VarbetaInt[i];
    }
    delete [] betaInt;
    delete [] VarbetaInt;
    delete [] PvalM;
    delete [] PvalInt;
    delete [] PvalJoint;
    AF.clear();
  } // end of snploop 
  phenodata.clear();
  covdata.clear();
  delete [] XTransX;
  delete [] XinvXTX;

  fclose(fin);

  delete [] snpID;
  delete [] rsID;
  delete [] chrStr;
  delete [] allele1;
  delete [] allele0;
  delete [] samID;
  delete [] include_idx;
  
  results.close();
  std::chrono::duration<double> wallduration = (std::chrono::system_clock::now() - wall0);
  double cpuduration = (std::clock() - cpu0)/(double)CLOCKS_PER_SEC;
  cout << "Total Wall Time = " << wallduration.count() << "  Seconds\n";
  cout << "Total CPU Time = " << cpuduration << "  Seconds\n";
  cout << "Execution Wall Time = " << exetime << "  Seconds\n";
  cout << "*********************************************************\n";
  return 0;
}
/* end of main*/


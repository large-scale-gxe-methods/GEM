#pragma once
#include "GMMAT.h"
#include "MAGEE.h"
#include "declars.h"

int  checkBinary(unordered_map<string, vector<vector<string>>> phenoMap, vector<string> sampleID, double epsilon);
void center(int center, int scale, int samSize, int numSelCol, vector<double> covdata, vector<double>* covdata_ret);
void fitNullModel2(int samSize, int numSelCol, int phenoType, double epsilon, 
                    int robust, std::vector<string> covSelHeadersName, std::vector<double> phenodata, 
                    std::vector<double> covdata, std::vector<double>* XinvXTX_ret, vector<double>* miu_ret, 
                    vector<double>* resid_ret, double* sigma2_ret, std::vector<double>& beta_ret,
                    std::vector<double>& Xbeta_ret);
void fitNullModel(int samSize, int numSelCol, int phenoType, double epsilon, 
                    int robust, std::vector<string> covSelHeadersName, std::vector<double> phenodata, 
                    std::vector<double> covdata, std::vector<double>* XinvXTX_ret, vector<double>* miu_ret, 
                    vector<double>* resid_ret, double* sigma2_ret);
void printCovVarMat(int numCovs, vector<string> covNames, double* covVarMat, double* beta, int phenoType, int samSize);
void printOutputHeader(bool useBgen, int numExpSelCol_new, int Sq1, vector<string> covNames, string output, string outStyle, 
                       int robust, double sigma2, BinE binE);
            


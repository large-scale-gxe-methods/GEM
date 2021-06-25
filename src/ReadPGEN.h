#include <string>
#include <stdio.h>
#include <stdlib.h>

#include "ReadParameters.h"
#include "TimeUtils.h"

#ifndef READPGEN_H
#define READPGEN_H

class Pgen {

public:

    uint32_t raw_variant_ct;
    uint32_t raw_sample_ct;

    int new_samSize;
    std::vector<string>   sampleID;
    std::vector<double>   new_covdata;
    std::vector<double>   new_phenodata;
    std::vector<long int> include_idx;
    
    int phenoType;
    vector<uint32_t> begin;
    vector<uint32_t> end;
    uint32_t threads;
    bool filterVariants = false;
    std::vector<long long unsigned int> pgenVariantPos;
    vector<int> pvarIndex;
    int pvarLast;
    vector<int> binE_idx;

    void processPgenHeader(string pgenFile);
    void processPvar(Pgen pgen, string pvarFile);
    void processPsam(Pgen pgen, string psamFile, unordered_map<string, vector<string>> phenomap, string phenoMissingKey, int numSelCol, int samSize);
    void getPgenVariantPos(Pgen pgen, CommandLine cmd);
};

void gemPGEN(int thread_num, double sigma2, double* resid, double* XinvXTX, vector<double> miu, BinE binE, Pgen pgen, CommandLine cmd);

#endif

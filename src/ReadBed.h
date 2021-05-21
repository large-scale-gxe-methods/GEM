#pragma once
#include <string>
#include <stdio.h>
#include <stdlib.h>

#include "ReadParameters.h"
#include "TimeUtils.h"

#ifndef READBED_H
#define READBED_H

class Bed {

public:

    uint32_t n_samples = 0;
    uint32_t n_variants = 0;
    char famDelim;
    char bimDelim;
    int bimLast; 
    int new_samSize;
    std::vector<double>   new_covdata;
    std::vector<double>   new_phenodata;
    std::vector<long int> include_idx;

    vector<uint32_t> begin;
    vector<uint32_t> end;
    uint32_t threads;
    bool filterVariants = false;
    std::vector<long long unsigned int> bedVariantPos;


    void processBed(string bedFile, string bimFile, string famFile);
    void processFam(Bed bed, string famFile, unordered_map<string, vector<string>> phenomap, string phenoMissingKey, vector<double> phenodata, vector<double> covdata, int numSelCol, int samSize, double center, double scale);
    void getBedVariantPos(Bed bed, CommandLine cmd);
};

void gemBED(int thread_num, int phenoType, double sigma2, double* resid, double* XinvXTX, double* covX, vector<double> miu, Bed bed, CommandLine cmd);

#endif


#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "ReadParameters.h"
#include "TimeUtils.h"
#ifndef READBGEN_H
#define READBGEN_H
class Bgen {

    public:

        // For file
        FILE* fin;

        // For BGEN offset
        uint offset;

        // For BGEN header block
        uint Mbgen;
        uint Nbgen;
        uint CompressedSNPBlocks;
        uint Layout;

        // For BGEN header-flag block;
        uint SampleIdentifiers;


        // For ID matching
        int new_samSize;
        std::vector<double>   new_covdata;
        std::vector<double>   new_phenodata;
        std::vector<long int> include_idx;
        vector <long int> variant_pos;
        std::vector<unsigned int> includeVariantIndex;

        // For multithreading BGEN file
        uint threads;
        bool filterVariants;
        std::vector<uint> Mbgen_begin;
        std::vector<uint> Mbgen_end;
        std::vector<long long unsigned int> bgenVariantPos;
        std::vector<vector<uint>> keepVariants;


        void processBgenHeaderBlock(string bgenfile);
        void processBgenSampleBlock(Bgen bgen, char samplefile[300], bool useSample, unordered_map<string, vector<string>> phenomap, string phenoMissingKey, vector<double> phenodata, vector<double> covdata, int numSelCol, int samSize, double center, double scale);
        void getPositionOfBgenVariant(Bgen bgen, CommandLine cmd);
};

void gemBGEN(int thread_num, double sigma2, double* resid, double* XinvXTX, double* covX, vector<double> miu, Bgen bgen, CommandLine cmd);
#endif

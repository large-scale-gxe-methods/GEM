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
	std::vector<double>   new_covdata;
	std::vector<double>   new_phenodata;
	std::vector<long int> include_idx;
	
	vector<uint32_t> begin;
	vector<uint32_t> end;
	uint32_t threads;
	bool filterVariants;
	// Temporary
	int stream_snps;
	double  maf;
	double missGeno;
	vector < double > miu;
	int     phenoTyp;
	double  sigma2;
	double* covX;
	double* XinvXTX;
	double* resid;
	int numSelCol;
	int numIntSelCol;
	int numExpSelCol;
	int Sq;
	int samSize;
	int robust;
	std::string outFile;
	std::vector<long long unsigned int> pgenVariantPos;

	void processPgenHeader(string pgenFile);
	void processPvar(Pgen pgen, string pvarFile);
	void processPsam(Pgen pgen, string psamFile, unordered_map<string, vector<string>> phenomap, string phenoMissingKey, vector<double> phenodata, vector<double> covdata, int numSelCol, int samSize);
	void getPgenVariantPos(Pgen pgen, CommandLine cmd);
};

void gemPGEN(uint32_t begin, uint32_t end, string pgenFile, string pvarFile, int thread_num, Pgen test);

#endif

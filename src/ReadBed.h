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

	void processBed(string bedFile, string bimFile, string famFile);
	//void processBim(Pgen pgen, string pvarFile);
	void processFam(Bed bed, string famFile, unordered_map<string, vector<string>> phenomap, string phenoMissingKey, vector<double> phenodata, vector<double> covdata, int numSelCol, int samSize);
	//void getPgenVariantPos(Pgen pgen, CommandLine cmd);
};

void gemBED(uint32_t begin, uint32_t end, string bedFile, string bimFile, int thread_num, Bed test);

#endif

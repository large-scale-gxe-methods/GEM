#ifndef READBGEN_H
#define READBGEN_H


#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "TimeUtils.h"

class Bgen: public Time {

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
	std::vector<string>   geno_snpid;
	std::vector<double>   new_covdata;
	std::vector<double>   new_phenodata;
	std::vector<long int> include_idx;



	// Temporary
	int stream_snps;
	double  maf;
	vector < double > miu;
	int     phenoTyp;
	double  sigma2;
	double* covX;
	double* XinvXTX;
	double* resid;
	int numSelCol;
	int Sq;
	int samSize;
	int robust;
	void processBgenHeaderBlock(char genofile[300]);
	void processBgenSampleBlock(Bgen bgen, char samplefile[300], unordered_map<string, vector<string>> phenomap, string phenoMissingKey, vector<double> phenodata, vector<double> covdata, int numSelCol, int samSize);
};


std::vector<long int> getPositionOfBgenVariant(FILE* fin, uint offset, uint Mbgen, uint Nbgen, uint CompressedSNPBlocks, uint Layout, vector<int> Mbgen_begin);
void BgenParallelGWAS(int begin, int end, long int byte, char genobgen[300], int thread_num, Bgen test);
void GEMforBGEN13(vector<uint> DLens, vector<vector <uchar>> zBufs, vector<uint> zLens, vector<string> geno_snpid, Bgen test);
#endif

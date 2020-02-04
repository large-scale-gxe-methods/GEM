#ifndef READBGEN_H
#define READBGEN_H


#include <string>
#include <stdio.h>
#include <stdlib.h>


class Bgen {

public:

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
	vector<long int> include_idx;
	std::vector<double> covdata;
	std::vector<double> phenodata;
	unordered_set<int>  genoUnMatchID;

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
	//void processBgenSamples(BgenHeader bgen, char samplefile[300], unordered_map<string, vector<string>> phenomap, string phenoMissingKey, vector<double> phenodata, vector<double> covdata, int numSelCol);
};


std::vector<long int> getPositionOfBgenVariant(FILE* fin, uint offset, uint Mbgen, uint Nbgen, uint CompressedSNPBlocks, uint Layout, vector<int> Mbgen_begin);
void BgenParallelGWAS(int begin, int end, long int byte, char genobgen[300], int thread_num,  Bgen test);

#endif

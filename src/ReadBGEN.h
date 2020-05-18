
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
		int threads;
		bool filterVariants;
		vector<uint> Mbgen_begin;
		vector<uint> Mbgen_end;
		vector<uint> bgenVariantPos;
		vector<vector<uint>> keepVariants;


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
		int numIntSelCol;
		int numExpSelCol;
		int Sq;
		int samSize;
		int robust;
		std::string outFile;
		void processBgenHeaderBlock(char genofile[300]);
		void processBgenSampleBlock(Bgen bgen, char samplefile[300], unordered_map<string, vector<string>> phenomap, string phenoMissingKey, vector<double> phenodata, vector<double> covdata, int numSelCol, int samSize);
		void getPositionOfBgenVariant(Bgen bgen, CommandLine cmd);
};

void BgenParallelGWAS(int begin, int end, long int byte, vector<uint> keepVariants, char genobgen[300], bool filterVariants, int thread_num, Bgen test);
#endif

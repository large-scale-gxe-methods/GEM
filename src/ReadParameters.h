#ifndef READPARAMETERS_H
#define READPARAMETERS_H

#include <string>
#include <stdio.h>
#include <stdlib.h>

class CommandLine {

	public:

		char genofile[300];
		char samplefile[300];


		// Input files
		std::string pFile;
		std::string bgenFile;
		std::string sampleFile;
		std::string phenoFile;
		std::string outFile;
		bool useSampleFile = false;

		// Inputs
		std::vector<std::string> cov;
		std::vector<std::string> icov;
		std::vector<std::string> exp;
		std::string phenoName;
		std::string sampleID;
		int phenoType;
	
		// Filtering options
		double MAF;
		double missGenoRate;
		std::string includeVariantFile;
		bool doFilters;

		std::string delim;
		char pheno_delim[300];
		std::string missing;


		// Performance options
		int threads;
		int stream_snps;
		int robust;
		double tol;
		void processCommandLine(int argc, char* argv[]);
};


#endif

#ifndef READPARAMETERS_H
#define READPARAMETERS_H

#include <string>
#include <stdio.h>
#include <stdlib.h>

class CommandLine {

public:
	char samplefile[300];

	// Genotype files
	std::string bgenFile;
	std::string sampleFile;
	std::string pgenFile;
	std::string pvarFile;
	std::string psamFile;
	std::string bedFile;
	std::string bimFile;
	std::string famFile;
	bool useSampleFile = false;
	bool useBgenFile = false;
	bool usePgenFile = false;
	bool useBedFile = false;

	// Phenotype file
	std::string phenoFile;
	std::string delim;
	std::string missing;
	char pheno_delim;

	// Out file
	std::string outFile;
	std::string outStyle;

	// Inputs
	std::string phenoName;
	std::string sampleID;
	int phenoType;

	int numSelCol    = 0;
	int numExpSelCol = 0;
	int numIntSelCol = 0;
	std::vector<std::string> cov;
	std::vector<std::string> icov;
	std::vector<std::string> exp;
	std::unordered_map<string, int> covHM;
	std::unordered_map<string, int> intHM;
	std::unordered_map<string, int> expHM;

	// Filtering options
	double MAF;
	double missGenoRate;
	std::string includeVariantFile;
	bool doFilters;

	// Performance options
	int center;
	int scale;
	double tol;
	int robust;
	int threads;
	int stream_snps;
	void processCommandLine(int argc, char* argv[]);
};


#endif

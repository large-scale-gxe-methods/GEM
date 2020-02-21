#ifndef READPARAMETERS_H
#define READPARAMETERS_H

#include <string>
#include <stdio.h>
#include <stdlib.h>

class CommandLine {

public:

	// Input files
	std::string pFile;


	// Filtering options
	double MAF;


    // Performance options
	int threads;

	void processCommandLine(int argc, char* argv[]);
};



void ReadParameters(char* paramfile, char* genofile, char* phenofile, char* samplefile, PARAMETERS* prm);

#endif

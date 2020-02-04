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



class GwasThreadSpecificVariables {
public:
	vector <double> ZGSvec;
	vector <double> ZGSR2vec;
	vector <double> WZGSvec;
	double* WZGS;
	double* ZGS;
	double* XtransZGS;
	double   ZGStemp;
	double* ZGSR2;
	double* ZGStR;
	double* ZGStZGS;
	double* ZGSR2tZGS;

	double* betaM;
	double* VarbetaM;
	double** betaInt;
	double** VarbetaInt;
	double* PvalM;
	double* PvalInt;
	double* PvalJoint;

	double* S2TransS2;
	double* S2TransR;
	double* S2DS2;
	double* InvVarbetaint;

	double* Stemp2;
	double   statM;
	double* Stemp3;
	double   statInt;
	double   statJoint;
};



class BgenThreadSpecificVariables: public GwasThreadSpecificVariables {
public:

	// For idexing. May change within code
	int stream_snps;
	int stream_i;
	int snploop;
	int ZGS_col;
	int Sq1;

	FILE* fin3;

	// For reading variant block
	uint maxLA = 65536;
	uint maxLB = 65536;
	uint Nrow;
	ushort LS;
	char* snpID   = new char[maxLA + 1];
	ushort LR;
	char* rsID    = new char[maxLA + 1];
	ushort LC;
	char* chrStr  = new char[maxLA + 1];
	uint physpos;
	ushort LKnum;
	ushort LA;
	char* allele1 = new char[maxLA + 1];
	ushort LB;
	char* allele0 = new char[maxLB + 1];
	vector <std::string> geno_snpid;
	uint zLen;
	uint DLen;
	vector <uchar> zBuf;



	// For uncompressing genotype block
	uLongf destLen;
	vector <uchar> shortBuf;
	uchar* bufAt;


	// For genotype block
	uint N;
	uint K;
	uint Pmin;
	uint Pmax;
	uint ploidyMiss;
	vector<unsigned int> ploidy_and_missing_info;
	int ploidy_sum;
	int idx_to_sum;
	uint Phased;
	uint B;
	uint Bbits;
	uint chartem;
	uint chartem1;


	// For computing dosages
	int tmp1;
	int tmp2;
	int tmp3;
	int tmp4;
	uint idx_k;
	vector <double> AF;
	double p10;
	double p11;
	double p00;
	double dosage;
	double pTot;
};




void ReadParameters(char* paramfile, char* genofile, char* phenofile, char* samplefile, PARAMETERS* prm);
#endif


#include <vector>
#include <string>
#include "ReadParameters.h"

#ifndef BinaryEUtils_H
#define BinaryEUtils_H

class BinE {

    public:
        int nBinE = 0;
        int strataLen = 0;

        std::vector<int> stratum_idx;
        
        std::vector<string> bin_headers;

       void checkBinaryCovariates(BinE binE, CommandLine cmd, unordered_map<string, vector<string>> phenoMap, vector<string> sampleID, vector<long int> include_idx, int samSize, int phenoType, string phenoHeaderName, std::vector<string> covSelHeadersName, std::vector<string> covSelHeadersName_new, int Sq_new);
};

std::unordered_map <int, vector<int>> GetPowerSet(std::string v);
std::vector<std::string> vector_binstrings(int width);

#endif


#include <vector>
#include <string>

#ifndef BinaryEUtils_H
#define BinaryEUtils_H

class BinE {

    public:
        int nBinE = 0;
        int nss = 0;
        int strataLen = 0;
        std::vector<int> binE_idx;

        std::vector<int> stratum_idx;
        std::vector<int> sub_stratum_idx;
        std::vector<int> sub_stratum_size;
        
        std::vector<string> bn_header;

        void checkBinaryCovariates(BinE binE, unordered_map<string, vector<string>> phenoMap, vector<string> sampleID, vector<long int> include_idx, int samSize, int numExpSelCol, int numSelCol, std::vector<string> covSelHeadersName);
};

std::unordered_map <int, vector<int>> GetPowerSet(std::string v);
std::vector<std::string> vector_binstrings(int width);

#endif

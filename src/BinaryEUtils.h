
#include <vector>
#include <string>

#ifndef BinaryEUtils_H
#define BinaryEUtils_H

class BinE {

    public:
        size_t numBinE = 0;
        size_t numSubStrata = 0;
        std::vector<int> binE_idx;

        std::vector<int> stratum_idx;
        std::vector<int> sub_stratum_idx;
        std::vector<int> sub_stratum_size;
        
        std::vector<std::string> strata_names;

        void checkBinaryCovariates(BinE binE, unordered_map<string, vector<string>> phenoMap, vector<string> sampleID, vector<long int> include_idx, int numExpSelCol, int numSelCol);
};

std::unordered_map <int, vector<int>> GetPowerSet(std::string v);
std::vector<std::string> vector_binstrings(int width);

#endif

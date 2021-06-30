#include "declars.h"
#include "BinaryEUtils.h"

std::unordered_map <int, vector<int>> GetPowerSet(std::string v) 
{
    std::string emptyString;
    std::vector<std::string> powerSet;
    int n = (int) std::pow(2.0, (double) v.size());
    powerSet.reserve(n);
    powerSet.push_back(emptyString);

    for (std::string::iterator it = v.begin(); it < v.end(); it++) {
        unsigned int tempSize = powerSet.size();
        for (std::size_t j = 0; j < tempSize; j++)
            powerSet.push_back(powerSet[j] + *it);
    }

    powerSet.erase(powerSet.begin());

    std::vector<int> Data;
    std::transform(powerSet.begin(), powerSet.end(), std::back_inserter(Data),
               [](const std::string& str) { return std::stoi(str); });
    sort(Data.begin(), Data.end());
    std::vector<std::string> Data2;
    std::transform(Data.begin(), Data.end(), std::back_inserter(Data2),
        [](const int& str) { return std::to_string(str); });

    std::unordered_map <int, vector<int>> map;
    for (size_t i = 0; i < Data2.size(); i++) {
        std::vector<int> vec;
        for (size_t j = 0; j < Data2[i].size(); ++j) {                                 
            vec.push_back((Data2[i][j] - '0') - 1);
        }
        map[i] = vec;
    }

    return(map);
}

std::string to_binstring(int x, int width)
{
    std::string s(width, '0');
    for (int i = 0; i < width; i++)
    {
        s[i] = '0' + x % 2;
        x /= 2;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

std::unordered_map<std::string, unsigned int> map_binstrings(int width)
{
    const unsigned int limit = 1 << width;
    std::unordered_map<std::string, unsigned int> bins;
    for (unsigned int i = 0; i < limit; i++) {
        string binString = to_binstring(i, width);
        bins[binString] = i;
    }
    return bins;
}

std::vector<std::string> vector_binstrings(int width)
{
    const unsigned int limit = 1 << width;
    std::vector<std::string> bins;
    for (unsigned int i = 0; i < limit; i++) {
        bins.push_back(to_binstring(i, width));
    }
    return bins;
}


void BinE::checkBinaryCovariates(BinE binE, int samSize, unordered_map<string, vector<string>> phenoMap, vector<string> sampleID, vector<long int> include_idx, int numExpSelCol, int numSelCol) 
{

    stratum_idx.resize(include_idx.size());
    std::unordered_map<string, int> map;

    for (int i = 0; i < numExpSelCol; i++) 
    {
        for (int j = 0; j < samSize; j++ ) {
            auto tmp = phenoMap[sampleID[j]];
            if (!map.count(tmp[1 + i])) { map[tmp[1 + i]] = 1; }
            if (map.size() > 2) { break; }
        }

        if (map.size() < 2) { cerr << "\nERROR: All values of covariate columns are the same.\n\n"; }

        if (map.size() == 2) {
            binE_idx.push_back(i + 1);
            for (auto kv : map) { strata_names.push_back(kv.first); }
            numBinE++;
        }
        map.clear();
    }


    if (numBinE > 0) 
    {
        std::unordered_map<std::string, unsigned int> stratum_map = map_binstrings(numBinE);
        for (int j = 0; j < samSize; j++ ) {
            auto tmp = phenoMap[sampleID[j]];
            string stratum = "";
            for (int i = 0; i < numBinE; i++) {
                int tmp1 = i*2;
                (strata_names[tmp1].compare(tmp[binE_idx[i]])) ? stratum+="1" : stratum+="0";
            }
            stratum_idx[j] = stratum_map[stratum];          
        }
        
        if (numBinE > 1) { 
                       
            std::string tmp = "";
            for(int i = 1; i <= numBinE; i++) {
                tmp+=std::to_string(i);
            }

            std::unordered_map<int, vector<int>> map_combo = GetPowerSet(tmp);

            std::vector<std::string> stratum_vector = vector_binstrings(numBinE);
            for (size_t i = 0; i < (map_combo.size()-1); i++) {
                vector<int> vec = map_combo[i];
                std::vector<std::string> vec_combo = vector_binstrings(vec.size());
                numSubStrata+=vec_combo.size();
                for (size_t j = 0; j < vec_combo.size(); j++) {
                    int cnt = 0;
                    for (size_t k = 0; k < stratum_vector.size(); k++) {
                        string tmp_main = stratum_vector[k];
                        string tmp_s = "";
                        for (size_t l = 0; l < vec.size(); l++) {
                            tmp_s += tmp_main[vec[l]];
                        }
                        if (tmp_s.compare(vec_combo[j]) == 0) {
                            sub_stratum_idx.push_back(k);
                            cnt++;
                        }
                    }
                    sub_stratum_size.push_back(cnt);
                }
            }
        }
    }
}
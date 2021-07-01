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

std::unordered_map<std::string, unsigned int> cartesian_map( vector<vector<string> >& v ) {
  auto product = []( long long a, vector<string>& b ) { return a*b.size(); };
  const long long N = accumulate( v.begin(), v.end(), 1LL, product );
  vector<string> u(v.size());
  std::unordered_map<std::string, unsigned int> map;
  for( long long n=0 ; n<N ; ++n ) {
    lldiv_t q { n, 0 };
    for( long long i=v.size()-1 ; 0<=i ; --i ) {
      q = div( q.quot, v[i].size() );
      u[i] = v[i][q.rem];
    }

    string tmp = "";
    for( size_t i = 0; i < u.size(); i++) { 
        tmp+=u[i]; 
    }
    map[tmp] = n;
  }
  return map;
}

std::vector<std::string> cartesian_vec( vector<vector<string> >& v ) {
  auto product = []( long long a, vector<string>& b ) { return a*b.size(); };
  const long long N = accumulate( v.begin(), v.end(), 1LL, product );
  vector<string> u(v.size());
  std::vector<std::string> vec(N);
  for( long long n=0 ; n<N ; ++n ) {
    lldiv_t q { n, 0 };
    for( long long i=v.size()-1 ; 0<=i ; --i ) {
      q = div( q.quot, v[i].size() );
      u[i] = v[i][q.rem];
    }

    string tmp = "";
    for( size_t i = 0; i < u.size(); i++) { 
        tmp+=u[i]; 
    }
    vec[n] = tmp;
  }
  return vec;
}

std::vector<std::string> cartesian_vec_sep( vector<vector<string> >& v) {
  auto product = []( long long a, vector<string>& b ) { return a*b.size(); };
  const long long N = accumulate( v.begin(), v.end(), 1LL, product );
  vector<string> u(v.size());
  std::vector<std::string> vec(N);
  for( long long n=0 ; n<N ; ++n ) {
    lldiv_t q { n, 0 };
    for( long long i=v.size()-1 ; 0<=i ; --i ) {
      q = div( q.quot, v[i].size() );
      u[i] = v[i][q.rem];
    }

    string tmp = "";
    for( size_t i = 0; i < u.size(); i++) { 
        tmp = (i == 0) ? u[i] : tmp + "_" + u[i]; 
    }
    vec[n] = tmp;
  }
  return vec;
}

void BinE::checkBinaryCovariates(BinE binE, unordered_map<string, vector<string>> phenoMap, vector<string> sampleID, vector<long int> include_idx, int samSize, int numExpSelCol, int numSelCol, std::vector<string> covNames) 
{

    size_t t_sampleSize = include_idx.size(); 
    stratum_idx.resize(t_sampleSize);
    std::unordered_map<string, int> map;
    std::vector<std::vector<std::string>> stratum_names;

    for (int i = 0; i < numExpSelCol; i++) {
        for (int j = 0; j < samSize; j++ ) {
            auto tmp = phenoMap[sampleID[j]];
            if (!map.count(tmp[1 + i])) { map[tmp[1 + i]] = 1; }
            if ((map.size()/t_sampleSize) > 0.05) { break; }
        }

        if (map.size() < 2) { cerr << "\nERROR: All values of covariate columns are the same.\n\n"; }

        if ((map.size()/t_sampleSize) < 0.05) {
            binE_idx.push_back(i + 1);
            vector<string> v;
            for (auto kv : map) {  
                v.push_back(kv.first); 
            }
            stratum_names.push_back(v);
            nBinE++;
        }
        map.clear();
    }

    if (nBinE > 0) 
    {
        //Assign each sample an index from stratum_map
        std::unordered_map<string, unsigned int> stratum_map = cartesian_map(stratum_names);
        for (int j = 0; j < samSize; j++ ) {
            auto sv = phenoMap[sampleID[j]];
            std::string stratum = "";
            for (int i = 0; i < nBinE; i++) {
                stratum += sv[binE_idx[i]];
            }
            stratum_idx[j] = stratum_map[stratum]; 
        }
        for (auto kv : stratum_map) { 
            strataLen++;
        }
        std::vector<string> bn;
        for(int i = 0; i < nBinE; i++) {
            bn.push_back(covNames[binE_idx[i] - 1]);
        }

        if (nBinE > 1) {     
            std::string powerString = "";
            for(int i = 1; i <= nBinE; i++) {
                powerString+=std::to_string(i);
            }
            std::unordered_map<int, vector<int>> powerSet = GetPowerSet(powerString);
            
            std::vector<string> stratum_vec = cartesian_vec(stratum_names);
            vector<vector<string>> sn;
            for (size_t i = 0; i < (powerSet.size()-1); i++) {
                vector<int> ps = powerSet[i];
                string bns = "";
                for (size_t j = 0; j < ps.size(); j++) {
                    sn.push_back(stratum_names[ps[j]]);
                    bns = bns + bn[ps[j]] + "_";
                }
 
                std::vector<string> cart = cartesian_vec(sn);
                std::vector<string> cart_sep = cartesian_vec_sep(sn);
                for (size_t k = 0; k < cart.size(); k++) {
                    int cnt = 0;
                    for (size_t l = 0; l < stratum_vec.size(); l++) {
                        string svs = "";
                        for (size_t j = 0; j < ps.size(); j++) {
                            svs += stratum_vec[l][ps[j]];
                        }
                        if (svs.compare(cart[k]) == 0) {
                            sub_stratum_idx.push_back(l);
                            cnt++;
                        }
                    }
                    bn_header.push_back(bns + cart_sep[k]);
                    sub_stratum_size.push_back(cnt);
                    nss++;
                }
                sn.clear();
            }

            // Final combination for header output
            string bns = "";
            for (int i = 0; i < nBinE; i++) {
                bns = bns + bn[i] + "_";
            }
            std::vector<string> cart_sep = cartesian_vec_sep(stratum_names);
            for (size_t j = 0; j < cart_sep.size(); j++) {
                bn_header.push_back(bns + cart_sep[j]);
            }

        } else {
            vector<vector<string>> sn;
            for (int i = 0; i < nBinE; i++) {
                string bns = bn[i] + "_";
                sn.push_back(stratum_names[i]);
                std::vector<string> cart_sep = cartesian_vec_sep(sn);
                for (size_t j = 0; j < cart_sep.size(); j++) {
                    bn_header.push_back(bns+cart_sep[j]);
                }
                sn.clear();
            }
        }
    }
}
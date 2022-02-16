#include "declars.h"
#include "BinaryEUtils.h"

std::unordered_map<std::string, unsigned int> cartesian_map( vector<vector<string> >& v ) 
{
  auto product = []( long long a, vector<string>& b ) { return a*b.size(); };
  const long long N = accumulate( v.begin(), v.end(), 1LL, product );
  vector<string> u(v.size());
  std::unordered_map<std::string, unsigned int> map;
  for( long long n=0 ; n<N ; ++n ) {
    lldiv_t q { n, 0 };
    for( long long i=v.size()-1 ; 0<=i ; --i ) {
      q = lldiv( q.quot, v[i].size() );
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

std::vector<std::string> cartesian_vec_sep( vector<vector<string> >& v) 
{
  auto product = []( long long a, vector<string>& b ) { return a*b.size(); };
  const long long N = accumulate( v.begin(), v.end(), 1LL, product );
  vector<string> u(v.size());
  std::vector<std::string> vec(N);
  for( long long n=0 ; n<N ; ++n ) {
    lldiv_t q { n, 0 };
    for( long long i=v.size()-1 ; 0<=i ; --i ) {
      q = lldiv( q.quot, v[i].size() );
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

void BinE::checkBinaryCovariates(BinE binE, CommandLine cmd, unordered_map<string, vector<string>> phenoMap, vector<string> sampleID, vector<long int> include_idx, int samSize, std::vector<string> covNames) 
{

    stratum_idx.resize(samSize);

    std::vector<int> binE_idx;
    std::unordered_map<string, int> map;
    std::vector<std::vector<std::string>> stratum_names;

    int Sq = cmd.numExpSelCol + cmd.numIntSelCol;
    int cat_threshold = cmd.cat_threshold;
    std::vector<std::string> cat_names = cmd.cat_names;
    for (int i = 0; i < Sq; i++) {
        int cnt = 0;
        for (int j = 0; j < samSize; j++ ) {
            auto tmp = phenoMap[sampleID[j]];
            if (!map.count(tmp[1 + i])) { 
                map[tmp[1 + i]] = 1;
                cnt++;
            }
        }

        if (cnt == 1) { 
            cerr << "\nERROR: All values of " << covNames[i] << " column are the same.\n\n"; 
            exit(1);	
	    }

        if ((cnt <= cat_threshold) || (std::find(cat_names.begin(), cat_names.end(), covNames[i]) != cat_names.end())) {
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
        
        for (auto kv : stratum_map) { 
            strataLen++;
        }

        // Output header prefix
        string bin_names = "";
        for(int i = 0; i < nBinE; i++) {
            bin_names = bin_names + covNames[binE_idx[i] - 1] + "_";
        }
  
        std::vector<string> cart_sep = cartesian_vec_sep(stratum_names);
        for (size_t i = 0; i < cart_sep.size(); i++) {
            bin_headers.push_back(bin_names+cart_sep[i]);
        }  
	    
	std::vector<string> bin_headers_tmp=bin_headers;  
        vector<int> headerMap(cart_sep.size(), 0);
	    
        for (int i = 0 ; i < cart_sep.size() ; i++) {
         headerMap[i] = i;
        }
	    
        sort(headerMap.begin(), headerMap.end(),
       [&](const int& a, const int& b) {
         return (bin_headers[a] < bin_headers[b]);
        }
            );


        sort(bin_headers.begin(),bin_headers.end());
 
        for (int i = 0 ; i < cart_sep.size() ; i++){
          auto it=find(bin_headers.begin(),bin_headers.end(),bin_headers_tmp[i]);
          if (it != bin_headers.end()){
             headerMap[i] = it- bin_headers.begin();
            
          }
        }


        for (int j = 0; j < samSize; j++ ) {
            auto sv = phenoMap[sampleID[j]];
            std::string stratum = "";
            for (int i = 0; i < nBinE; i++) {
                stratum += sv[binE_idx[i]];

            }

            stratum_idx[j] = headerMap[stratum_map[stratum]];
    
        }
    }

    cout << "Number of categorical variables: " << nBinE << endl;
}

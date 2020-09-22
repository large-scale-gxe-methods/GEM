#include <vector>
#include <cstring>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <utility>
#include <cctype>
#include <sstream>
#include <cstdlib>
#include <cassert>
#include <ctime>
#include <chrono>
#include <omp.h>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <exception>


#include <boost/thread.hpp>
#include <boost/thread/tss.hpp>
#include <boost/asio.hpp>
#include <boost/make_shared.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string.hpp> 
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/program_options.hpp>
#include <boost/make_shared.hpp>
#include <boost/thread.hpp>
#include <boost/thread/tss.hpp>
#include <boost/asio.hpp>
#include <boost/filesystem.hpp> 
#include "zlib.h"


using std::string;
using std::vector;
using std::cin;
using std::cout;
using std::endl;
using std::cerr;
using std::getline;
using std::find;
using std::unordered_map;
using std::unordered_set;
namespace po = boost::program_options;


typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long uint64;
typedef unsigned short ushort;


#include "MatrixUtils.h"
#include "ReadParameters.h"
#include "ReadBGEN.h"
#include "ReadPGEN.h"
#include "ReadBed.h"

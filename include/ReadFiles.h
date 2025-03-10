#pragma once

/**
 * @file ReadFiles.h
 * @author In this file we define a class to read files with different delimiter 
 * @brief 
 * @version 0.1
 * @date 2024-03-01
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <optional>
#include <string>
#include <string_view>
#include <map>
#include <unordered_map>
#include <utility>
#include <set>
#include <vector>
#include <variant>
#include "fmt/core.h"
#include "fmt/format.h"



namespace std
{
    namespace ext
    {
        using V_string = std::vector<std::string>; 
        using VV_string = std::vector<V_string>; 
        using Data = unordered_map<string, V_string>;
        using V_double = std::vector<double>;
        using V_float = std::vector<float>;
        using V_int = std::vector<int>;
        using V_size_t = std::vector<std::size_t>;
        using VV_int = std::vector<V_int>;
        using V_bool = std::vector<bool>;
        using V_lluint = std::vector<long long unsigned int>;
        using map_str_int = std::map<std::string, int>;
        using V_opt_string = std::vector<std::optional<std::string>>;
        using Var_bool_int = std::variant<V_bool, V_int>;
        using Map_str_Vint = std::map<string, V_int>;
    }
}

#pragma once

/**
 * @brief A class for storing infromation in the tabular files
 * 
 */
class DataFrame
{
    public: 
        std::ext::Data m_data;
        int m_nrows = 0;
        int m_ncols = 0;
        std::ext::V_string m_headers;
        std::string m_missing_key;
        std::ext::V_string m_geno_ids;//bgenIDs

        bool isDataFrameCreated() const;
        int n_rows() const;
        int n_cols() const;

        /**
         * @brief show the n row of the dataframe
         * 
         * @param n : number of rows to show (by default 5)
         */
        void head(int n = 5);
        /**
         * @brief copy a data from one data frame based on a given headers to new one
         * 
         * @param v_hdrs : vector of requested headers to copy
         * @return DataFrame 
         */
        DataFrame copy_by_hdrs(std::ext::V_string const& v_hdrs);
        /**
         * @brief Read input file and store its data
         * 
         * @param path : path to the input file
         * @param delim : delimeter to separate columns
         */
        void read_file(std::string_view path, char delim = ',');
        /**
         * @brief Get the whole columns of a given hdr
         * 
         * @param hdr : a given header to return its values
         * @return std::ext::V_string 
         */
        std::ext::V_string get_header(std::string const& hdr) const;
        /**
         * @brief remove duplication based on the values in  given headers
         * 
         * @param v_hdrs 
         */
        void remove_duplicates(std::ext::V_string const& v_hdrs);
        size_t size_wo_duplicates(std::string const& hdr);
        DataFrame remove_duplicates(std::string const& hdr);
        bool any_duplicated(std::string const& hdr);
        std::set<std::string> list_duplicates(std::string const& hdr);
        void match_genoids(std::string hdr_id, std::ext::V_string const& v_hdrs);
        
    private:
        std::ext::V_string read_lines(std::string_view path);
        void fill_data(std::ext::V_string const& v_strs, char delim = ',');
};








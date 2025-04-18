// #include "../include/ReadFiles.h"

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include "ReadFiles.h"
#include <unordered_set>

bool DataFrame::isDataFrameCreated() const 
{
    return (m_nrows != 0 && m_ncols != 0);
}


int DataFrame::n_rows() const
{ 
    if (!isDataFrameCreated()) 
    {
        std::cerr << "Error: DataFrame has not been created.\n";
        exit(EXIT_FAILURE); 
    }
    return m_nrows; 
}


int DataFrame::n_cols() const
{
    if (!isDataFrameCreated()) 
    {
        std::cerr << "Error: DataFrame has not been created.\n";
        exit(EXIT_FAILURE); 
    }
    return m_ncols; 
}


std::ext::V_string  DataFrame::read_lines(std::string_view path)
{
    std::ifstream ifs(path.data());

    if(ifs.fail())
    {
        fmt::print("Error in reading file {}.\nPlease check your inputs.\n", path);
        exit(EXIT_FAILURE);
    }

    std::string line;
    std::ext::V_string v_strs;

    while(std::getline(ifs, line))
    {
        line.erase(std::remove(line.begin(), line.end(), '\r'),line.end());
        v_strs.emplace_back(line);
    }

    return v_strs;
}

void DataFrame::fill_data(std::ext::V_string const& lines, char delim)
{
    std::ext::VV_string vv_strs;
    for(auto const& line : lines)
    {
        std::istringstream iss(line);
        std::string cell;
        std::ext::V_string v_str_tmp;
        while(std::getline(iss, cell, delim))
        {
            cell.erase(std::remove(cell.begin(), cell.end(), '\"'), cell.end());
            // v_str_tmp.emplace_back(cell);
            if (cell.empty())
            {   
                v_str_tmp.emplace_back(m_missing_key);
            }
                else
            {
                v_str_tmp.emplace_back(cell);
            }
        }
        vv_strs.emplace_back(v_str_tmp);
    }
    m_nrows = vv_strs.size();
    m_ncols = vv_strs[0].size();

    m_headers.resize(m_ncols);

    for(int i{0}; i < m_ncols; ++i)
    {
        m_headers[i] = vv_strs[0][i];
        std::ext::V_string values;

        for(int j{1}; j < m_nrows; ++j)
        {
            if(vv_strs[j].size() != m_ncols)
            {
                values.emplace_back(m_missing_key);
                continue;
            }
            
            values.emplace_back(vv_strs[j][i]);
        }
        m_data[m_headers[i]] = values;
    }
    m_nrows = m_data[m_headers[0]].size();
}

void DataFrame::head(int n)
{
        for(auto const& curr_hdr : m_headers)
        {
            std::cout << curr_hdr << " ";
        }

        std::cout << "\n";

        for(int i{0}; i <= n; ++i)
        {
            for(auto const& curr_hdr : m_headers)
            {
                std::cout << m_data[curr_hdr][i] << " ";
            }

            std::cout << "\n";
        }

}

void DataFrame::read_file(std::string_view path, char delim)
{
    auto v_strs = read_lines(path);
    fill_data(v_strs, delim);
}

std::ext::V_string DataFrame::get_header(std::string const& hdr) const
{
    auto it = m_data.find(hdr);

    if (it != m_data.end()) 
    {
        return it->second;
    } 
    else 
    {
        fmt::print("Error: Header {} not found. It must be one of the variables in the submitted file\n", hdr);
        exit(EXIT_FAILURE); // Return an empty vector
    }
}

//remove rows with missing data and match phenoIDs with genofile IDs
void DataFrame::match_genoids(std::string hdr_id, const std::ext::V_string& v_hdrs)
{
    // Convert m_geno_IDs(bgenIDs) to an unordered set for fast lookup
    std::unordered_set<std::string> geno_id_set(m_geno_ids.begin(), m_geno_ids.end());
    const std::ext::V_string& sam_ids = m_data[hdr_id];
    std::vector<int> keep_indices;

    for (int i = 0; i < sam_ids.size(); ++i)
    {
        if (geno_id_set.count(sam_ids[i]) > 0)
        {
            keep_indices.push_back(i);
        } 
    }

   std::vector<int> valid_indices;

    for (int idx : keep_indices)
    {
        bool is_valid_row = true;

        for (const auto& hdr : v_hdrs)
        {
            auto it = m_data.find(hdr);
            if (it != m_data.end())
            {
                if (idx < it->second.size() && it->second[idx] == m_missing_key)
                {
                    // std::cerr << "Warning: missing value at row: " << idx + 1 << " for header: " << hdr << "\n";
                    is_valid_row = false; // If there's any missing value, mark the row as invalid
                    break; // No need to check further headers for this row
                }
            }
        }

        if (is_valid_row)
        {
            valid_indices.push_back(idx);
        }
    }

    // Update each header in v_hdrs to contain only the valid rows
    for (const auto& hdr : v_hdrs)
    {
        auto it = m_data.find(hdr);
        if (it != m_data.end())
        {
            std::ext::V_string filtered_data;
            filtered_data.reserve(valid_indices.size()); // Reserve space for efficiency

            for (int idx : valid_indices)
            {
                if (idx < it->second.size())
                {
                    filtered_data.push_back(it->second[idx]);
                }
            }

            // Swap the filtered data back into the original data structure
            it->second.swap(filtered_data);
        }
    }

    // Update row count after filtering
    m_nrows = m_data[m_headers[0]].size();
}


DataFrame DataFrame::copy_by_hdrs(std::ext::V_string const& v_hdrs)
{
    DataFrame new_DataFrame;
    new_DataFrame.m_ncols = v_hdrs.size();
    new_DataFrame.m_nrows = m_data[v_hdrs[0]].size();
    new_DataFrame.m_headers = v_hdrs;
    new_DataFrame.m_missing_key = m_missing_key;
    new_DataFrame.m_geno_ids = m_geno_ids;

    if (v_hdrs.empty()) 
    {
        throw std::invalid_argument("Header vector is empty");
    }
    
    for(unsigned int i{0}; i < v_hdrs.size(); ++i)
    {
        new_DataFrame.m_data[v_hdrs[i]] = m_data[v_hdrs[i]];
    }
    return new_DataFrame;
}


void DataFrame::remove_duplicates(std::ext::V_string const& v_hdrs)
{
    int n_rm_row = 0;
    std::set<std::ext::V_string> st_rows;
    for(size_t i{0}; i < m_data[m_headers[0]].size(); ++i)
    {
        std::ext::V_string tmp_v_strs;

        for(auto hdr : v_hdrs)
        {
            tmp_v_strs.push_back(m_data[hdr][i]);
        }

        [[maybe_unused]]auto [it, success] = st_rows.insert(tmp_v_strs);

        if(!success)
        {
            for(auto hdr : m_headers)
            {
                m_data[hdr].erase(m_data[hdr].begin() + i - n_rm_row);
                fmt::print("Remove duplicate at row {}\n", i);
            }
            -- m_nrows;
            ++ n_rm_row;
        }
    }   
}


DataFrame DataFrame::remove_duplicates(std::string const& hdr)
{
    DataFrame df;
    std::unordered_set<std::string> tmp_st;
    for(auto elm : m_data[hdr])  
    {
        auto [it, success] = tmp_st.insert(elm);
        if(success)
        {
            df.m_data[hdr].push_back(elm);
        }
    } 
    return df;
}


size_t DataFrame::size_wo_duplicates(std::string const& hdr)
{
    std::unordered_set<std::string> tmp_st;

    for(auto elm : m_data[hdr])  
    {
        tmp_st.insert(elm);
        
    } 

    return tmp_st.size();
}


bool DataFrame::any_duplicated(std::string const& hdr)
{
    std::unordered_set<std::string> tmp_st;
    for(auto elm : m_data[hdr])  
    {
        auto [it, success] = tmp_st.insert(elm);
        if(!success)
        {
            return true;
        }
    } 
    return false;
}


std::set<std::string> DataFrame::list_duplicates(std::string const& hdr)
{
    std::unordered_set<std::string> unique_id;
    std::ext::V_string list_dup;

    for(auto elm : m_data[hdr])
    {
        auto [it, success] = unique_id.insert(elm);
        if(!success)
        {
            list_dup.push_back(elm);
        }
    }
    std::set<std::string> unique_dup(list_dup.begin(), list_dup.end());
    return unique_dup;
}


#include "../include/Pheno.h"
#include <iostream>
#include <iterator>
#include <unordered_set>


void Pheno::read_file(std::string_view path, char delim)
{
    m_data_frame.read_file(path, m_v_hdrs, delim);
}

std::pair<std::string, std::string> Pheno::check_binary(std::string pheno_name )
{
    auto v_phenos = m_data_frame.m_data[pheno_name];
    std::unordered_set<std::string> s_tmp;
    for(auto ph : v_phenos)
    {
        s_tmp.insert(ph);
        if(s_tmp.size() > 2)
        {
            return {"gaussian", "identity"};
        }
    }
    return {"binomial", "logit"};
}

unsigned int Pheno::size()
{
    return m_data_frame.m_data[m_sam_id].size();
}


void Pheno::set_path(std::string_view path)
{
    m_path = path.data();
}


std::string Pheno::get_path() const
{
    return m_path;
}

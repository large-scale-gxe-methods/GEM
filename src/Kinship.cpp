#include "../include/Kinship.h"
namespace details
{
    void add_dquot(std::ext::V_string& v_strs)
    {
        for(auto& str : v_strs)
        {
            str = fmt::format("\"{}\"", str);
        }  
    }
}


void Kinship::read_file(std::string_view path, char delim)
{
    m_data_frame.read_file(path, delim);

    if(m_data_frame.m_headers.size() != 3)
    {
        spdlog::error("Error in kinship files: the number of columns in '{}' should be 3",
            path);
        exit(EXIT_FAILURE);
    }
    add_dquot();
} 


void Kinship::add_dquot()
{
    details::add_dquot(m_data_frame.m_data[m_data_frame.m_headers[0]]);
    details::add_dquot(m_data_frame.m_data[m_data_frame.m_headers[1]]);
}


void Kinship::set_path(std::string const& path)
{
    m_path = path;
}


unsigned int Kinship::size()
{
    return m_data_frame.m_data[m_data_frame.m_headers[0]].size();
}
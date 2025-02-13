#pragma once


#include "ReadFiles.h"

class Kinship
{
    public:
        DataFrame m_data_frame;
        std::string m_path;
        bool m_null_kin = false;
        double m_diag;
        /**
         * @brief Set the path object
         * 
         * @param path 
         */
        void set_path(std::string const& path);
        /**
         * @brief Read kinship files with a given delimeter (default ",") and path
         * 
         * @param delim 
         */
        void read_file(std::string_view, char delim = ',');
        /**
         * @brief A function to get the size of kinship
         * 
         * @return unsigned int 
         */
        unsigned int size();
    private:
        void add_dquot();
};


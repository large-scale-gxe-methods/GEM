#pragma once

#include "ReadFiles.h"


class Pheno
{
    public:
        DataFrame m_data_frame;
        std::string m_sam_id;
        std::ext::V_string m_v_hdrs;
        // std::string m_phenoMissingKey;
        /**
         * @brief a function to get the size of pheno data
         * 
         * @return unsigned int 
         */
        unsigned int size();
        /**
         * @brief Set the path
         * 
         * @param path 
         */
        void set_path(std::string_view path);
        /**
         * @brief Get the path
         * 
         * @return std::string 
         */
        std::string get_path() const;
        /**
         * @brief A function to read pheno file with a given delimeter (default ",") and path
         * 
         * @param delim 
         */
        void read_file(std::string_view, char delim = ',');  
        std::pair<std::string, std::string> check_binary(std::string pheno_name );
    private:
        std::string m_path;
};
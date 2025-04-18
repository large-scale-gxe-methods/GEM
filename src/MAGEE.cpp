#include "MAGEE.h"

// glmm_gei_* functions to run G*E test
extern void glmm_gei_bgen13(Magee_Arma const& null_obj, string const &bgenfile,
					string const &outfile, double minmaf, double missrate,
					size_t npb, int ei, int qi, std::ext::Map_str_Vint const &strata_list, 
					uint begin, uint end, long long unsigned int byte, uint Nbgen,
					uint compression, bool meta_output = false);

extern void glmm_gei_pgen13(Magee_Arma const& null_obj, string const &pgenfile, std::string pvarFile,
					string const &outfile, double minmaf, double missrate,
					size_t npb, int ei, int qi, std::ext::Map_str_Vint const &strata_list, 
					uint begin, uint end, std::ext::V_lluint pgenPos, bool filterVariants, int pvarLength,
					int pvarLast, std::ext::V_int pvarIndex, bool meta_output = false);

void glmm_gei_bed13(Magee_Arma const& null_obj, string const &bedfile, 
					std::string bimFile, string const &outfile, double minmaf,
					double missrate, size_t npb, int ei, int qi, 
					std::ext::Map_str_Vint const &strata_list, 
					uint begin, uint end, std::ext::V_lluint bedPos, 
					bool filterVariants, char bimDelim, int bimLast,
					uint32_t n_samples, bool meta_output = false);

void conver_eigen_to_arma(SpaMat const& eigenMat, arma::sp_mat& armaMat) 
{
    // Determine the number of non-zero elements
    arma::umat locations(2, eigenMat.nonZeros());  // Store row and column indices
    arma::vec values(eigenMat.nonZeros());  // Store the values

    size_t index = 0;

    for (int k = 0; k < eigenMat.outerSize(); ++k) 
    {
        for (SpaMat::InnerIterator it(eigenMat, k); it; ++it) 
        {
            // Collect the row, column, and value for non-zero elements
            locations(0, index) = it.row();
            locations(1, index) = it.col();
            values(index) = it.value();
            ++index;
        }
    }
    // Use the collected data to construct the Armadillo sparse matrix
    armaMat = arma::sp_mat(locations, values, eigenMat.rows(), eigenMat.cols());
}


std::ext::V_int match_indices(std::ext::V_string const& original, std::ext::V_opt_string const& filtered) 
{
    std::ext::map_str_int filtered_map;
    for (size_t i = 0; i < filtered.size(); ++i) 
    {
        if (filtered[i].has_value()) 
        {
            filtered_map[filtered[i].value()] = i;
        }
    }

    std::ext::V_int indices;

    for (const auto& id : original) 
    {
        auto it = filtered_map.find(id);

        if (it != filtered_map.end()) 
        {
            indices.push_back(it->second);
        } 
        else 
        {
            indices.push_back(-1); // Use -1 to indicate not found
        }
    }
    return indices;
}

//Return id for match indexes or -1 for not match
std::ext::V_int match_indices(std::ext::V_string const& original, std::ext::V_string const& filtered) {
    std::ext::map_str_int filtered_map;
    
    for (size_t i = 0; i < filtered.size(); ++i) 
    {
        filtered_map[filtered[i]] = i;
    }
    
    std::ext::V_int indices;
    for (const auto& id : original) 
    {
        auto it = filtered_map.find(id);

        if (it != filtered_map.end()) 
        {
            indices.push_back(it->second);
        } 
        else 
        {
            indices.push_back(-1); // Use -1 to indicate not found
        }
    }
    return indices;
}


std::ext::V_string match_id_include(std::ext::V_string const& id_include, std::ext::V_string const& sample_id) 
{
    std::unordered_set<std::string> sample_set(sample_id.begin(), sample_id.end());
    std::ext::V_string matched_ids;

    for (const auto& id : id_include) 
    {
        if (sample_set.find(id) != sample_set.end()) 
        {
            matched_ids.push_back(id);
        }
    }
    return matched_ids;
}


bool any_isna(std::ext::V_int vec)
{
    for (auto& i : vec)
    {
        if (i == -1) return true;
    }
    return false;
}

std::ext::V_int remove_minus_one(std::ext::V_int vec)
{
    std::ext::V_int ret;
    for (auto& i : vec)
    {
        if (i != -1)
        {
            ret.push_back(i);
        }
    }
    return ret;
}


std::ext::V_string create_strata(std::ext::VV_string const& Ecat) {
    std::ext::V_string strata;

    for (const auto& row : Ecat) 
    {
        std::ostringstream oss;

        for (size_t i = 0; i < row.size(); ++i) 
        {
            if (i != 0) oss << "_";
            oss << std::stoi(row[i]);
        }
        strata.push_back(oss.str());
    }
    return strata;
}



template <typename T>
std::vector<T> unique_id(const std::vector<T> &vec) 
{
    std::unordered_set<T> seen; 
    std::vector<T> result;

    for (auto const& val : vec) 
    {
        if (seen.find(val) == seen.end()) 
        {  
            result.push_back(val);          
            seen.insert(val);               
        }
    }

    return result;
}

std::ext::V_bool in_op(std::ext::V_string const& vec1, std::ext::V_string const& vec2) 
{
    std::ext::V_bool result;
    result.reserve(vec1.size());
    
    for (const auto& elem : vec1) 
    {
        result.push_back(std::find(vec2.begin(), vec2.end(), elem) != vec2.end());
    }
    
    return result;
}

std::ext::V_int in_op_indices(std::ext::V_bool const& vec) 
{
    std::ext::V_int filtered;

    for (size_t i = 0; i < vec.size(); ++i) 
    {
        if (vec[i]) 
        {
            filtered.push_back(i);
        }
    }    
    return filtered;
}

std::ext::V_string filtered_ids(std::ext::V_string const& vec1, std::ext::V_bool const& vec2) 
{
    std::ext::V_string filtered;

    for (size_t i = 0; i < vec1.size(); ++i) 
    {
        if (vec2[i]) 
        {
            filtered.push_back(vec1[i]);
        }
    }
        return filtered;
}

void check_not_empty(std::ext::V_string const& sample_id) 
{
    if (sample_id.empty()) 
    {
        throw std::runtime_error("Error: null.obj$id_include does not match sample id in geno.file!");
    }
}


template <typename T>
std::ext::V_opt_string slice(std::vector<T> const& vec, std::ext::V_int const& indices) 
{
    std::ext::V_opt_string result;

    for (const auto& idx : indices) 
    {
        if (idx != -1) 
        {
            result.push_back(vec[idx]);
        } 
        else 
        {
            result.push_back(std::nullopt); // Use nullopt to represent NA
        }
    }
    return result;
}


bool any_duplicated(std::ext::V_string const& vec)
{
    std::set<std::string> tmp_st;

    for(auto elm : vec)  
    {
        auto [it, success] = tmp_st.insert(elm);
        if(!success)
        {
            return true;
        }
    } 
    return false;
}


std::ext::V_bool list_duplicates_bool(std::ext::V_string const& vec) 
{
    std::set<std::string> seen;
    std::vector<bool> result(vec.size(), false);

    for (size_t i = 0; i < vec.size(); ++i) 
    {
        if (seen.find(vec[i]) != seen.end()) 
        {
            result[i] = true;
        } else 
        {
            seen.insert(vec[i]);
        }
    }
    return result;
}


Mat filter_unique_rows(Mat const& mat, std::ext::V_string& id_include) 
{
    std::ext::V_bool is_duplicate = list_duplicates_bool(id_include);
    std::ext::V_int non_duplicate_indices;

    for (size_t i = 0; i < is_duplicate.size(); ++i) 
    {
        if (!is_duplicate[i]) 
        {
            non_duplicate_indices.push_back(i);
        }
    }

    Mat filtered(non_duplicate_indices.size(), mat.cols());

    for (size_t i = 0; i < non_duplicate_indices.size(); ++i) 
    {
        filtered.row(i) = mat.row(non_duplicate_indices[i]);
    }
    return filtered;
}


std::ext::V_bool apply_on_columns(Mat const& mat, const std::function<bool(VectorXd const&, int)>& func, int threshold) 
{
    std::ext::V_bool result(mat.cols(), false);

    for (int col = 0; col < mat.cols(); ++col) 
    {
        result[col] = func(mat.col(col), threshold);
    }
    return result;
}

// Custom function to check if unique elements in a vector are <= 20
bool unique_less_equal(DensVec const& vec, int threshold) 
{
    std::unordered_set<double> unique_elements(vec.data(), vec.data() + vec.size());
    return unique_elements.size() <= threshold;
}

std::ext::Map_str_Vint generate_strata_list(std::ext::V_string const& vec) 
{
  std::ext::Map_str_Vint strata_list;

  for (int i = 0; i < vec.size(); ++i) 
  {
    strata_list[vec[i]].push_back(i);
  }
  return strata_list;
}

void scale(Mat& mat, bool center, bool scale) 
{    
    if (center) 
    {
        mat = mat.rowwise() - mat.colwise().mean();
    }

    if (scale) 
    {
        Eigen::ArrayXd stddev = ((mat.array().square().colwise().sum()) / (mat.rows() - 1)).sqrt();
        mat = mat.array().rowwise() / stddev.transpose().array();
    }
}


Mat cbind(Mat const& A, Mat const& B) 
{
    Mat combined(A.rows(), A.cols() + B.cols());
    combined << A, B;
    return combined;
}

void spa_mat_ones(arma::sp_mat &arma_mat)
{
    auto rows = arma_mat.n_rows;
	auto cols = arma_mat.n_cols;
	arma::umat locations(2, rows * cols); 
    arma::vec values(rows * cols, arma::fill::ones); 

    for (int i = 0, k = 0; i < rows; ++i) 
	{
        for (int j = 0; j < cols; ++j, ++k) 
        {
            locations(0, k) = i; // Row index
            locations(1, k) = j; // Column index
        }
    }

    arma_mat = arma::sp_mat(locations, values, rows, cols);
}


MAGEE::MAGEE(GMMAT &gmmat, Glmmkin &glmmkin, CommandLine &cmd, Bgen bgen,
            std::ext::V_string interaction_exp, std::ext::V_string interaction_cov, 
            int numSelCol) 
            : m_gmmat(&gmmat), m_glmmkin_fitnull(glmmkin), m_cmd(cmd), m_genotype(std::move(bgen)),
            m_active_genotype(GENOTYPE::Bgen), m_interaction_exp(interaction_exp), 
            m_interaction_cov(interaction_cov), m_numSelCol(numSelCol){}
MAGEE::MAGEE(GMMAT &gmmat, Glmmkin &glmmkin, CommandLine &cmd, Pgen pgen,
            std::ext::V_string interaction_exp, std::ext::V_string interaction_cov, 
            int numSelCol) 
            : m_gmmat(&gmmat), m_glmmkin_fitnull(glmmkin), m_cmd(cmd), m_genotype(std::move(pgen)),
            m_active_genotype(GENOTYPE::Pgen), m_interaction_exp(interaction_exp), 
            m_interaction_cov(interaction_cov), m_numSelCol(numSelCol){}
MAGEE::MAGEE(GMMAT &gmmat, Glmmkin &glmmkin, CommandLine &cmd, Bed bed,
            std::ext::V_string interaction_exp, std::ext::V_string interaction_cov, 
            int numSelCol) 
            : m_gmmat(&gmmat), m_glmmkin_fitnull(glmmkin), m_cmd(cmd), m_genotype(std::move(bed)),
            m_active_genotype(GENOTYPE::Bed), m_interaction_exp(interaction_exp), 
            m_interaction_cov(interaction_cov), m_numSelCol(numSelCol){}


void MAGEE::clear_m_magee_glmmkin() 
{
    m_magee_glmmkin.J.resize(0, 0);
    m_magee_glmmkin.E.resize(0, 0);
    m_magee_glmmkin.EC.resize(0, 0);
    m_magee_glmmkin.residuals.resize(0);
    m_magee_glmmkin.sigma_i.resize(0, 0);
    m_magee_glmmkin.sigma_ix.resize(0, 0);
    m_magee_glmmkin.cov.resize(0, 0);
    m_magee_glmmkin.select.clear();
}

void MAGEE::create_E()
{
    std::ext::V_int colIndices;
    std::string header;

    try
    {
        for (const auto& h : m_interaction) 
        {
            colIndices.push_back(m_gmmat->m_hdrsMap.at(h));
        }
    }
    catch(std::out_of_range const& e)
    {
        //std::cerr << "{header} does not exist in headers" << e.what() << '\n';
        fmt::print("The header: \"{}\" does not exist in headers. Error: {}", header, e.what());
        std::exit(EXIT_FAILURE);
    }
    m_magee_glmmkin.E = slice_mat_cols(m_gmmat->m_X, colIndices);// should we add m_x to glmmkin return object? 
}


void MAGEE::fill_sel(std::ext::V_string& sample_id)
{
    std::ext::V_int missing_id= match_indices(m_glmmkin_fitnull.id_include, sample_id);
    
    if(any_isna(missing_id))
    {
        std::cout << "Warnning: check your data... Some individuals of pheno file are missing in sample file!\n";
       m_glmmkin_fitnull.id_include =  match_id_include(m_glmmkin_fitnull.id_include, sample_id);
        std::cout << "Missing IDs were removed...\n" << "Remaind IDs are: " << m_glmmkin_fitnull.id_include.size() << "\n";
    }

    std::ext::V_string sample_id_original =  sample_id;
    std::ext::V_bool sample_id_bool = in_op(sample_id, unique_id(m_glmmkin_fitnull.id_include));
    sample_id = filtered_ids(sample_id, sample_id_bool); 
    m_magee_glmmkin.select  = match_indices(sample_id_original, sample_id);
    
    try
    {
        check_not_empty(sample_id);
    }
    catch(const std::runtime_error& e)
    {
        fmt::print("Error: ids from GMMAT does not match sample in geno.file! Error: {} \n", e.what());
        std::exit(EXIT_FAILURE);
    }
}


void MAGEE::fill_mat(const int numRows, const int numCols, std::ext::V_int& indixes_col, int value)
{
    std::ext::VecTriple_i tripletList;

    for (int i = 0; i < numRows; ++i) 
    {
        if (indixes_col[i] != -1)
        {
            tripletList.emplace_back(i, indixes_col[i], value); 
        }
    } 

    m_magee_glmmkin.J.setFromTriplets(tripletList.begin(), tripletList.end());   
}

void MAGEE::fill_J(std::ext::V_int& match_id, std::ext::V_string& sample_id)
{
    std::ext::V_bool match_id_mid = in_op(m_glmmkin_fitnull.id_include, sample_id);
    match_id = in_op_indices(match_id_mid);
    m_glmmkin_fitnull.id_include = filtered_ids(m_glmmkin_fitnull.id_include, match_id_mid);
    std::ext::V_string unique_id_include = unique_id(m_glmmkin_fitnull.id_include);
    std::ext::V_int mid_indixes =  match_indices(sample_id, unique_id_include);
    auto slic_indices = slice(unique_id_include, mid_indixes);
    std::ext::V_int indixes_col = match_indices(m_glmmkin_fitnull.id_include, slic_indices); 
    std::ext::V_int indixes_col_unique = unique_id(indixes_col);
    m_magee_glmmkin.J.resize(m_glmmkin_fitnull.id_include.size(), indixes_col_unique.size());
    fill_mat(m_glmmkin_fitnull.id_include.size(), indixes_col_unique.size(), indixes_col, 1);
    m_magee_glmmkin.J = m_magee_glmmkin.J.transpose();
}

void MAGEE::calculate_bin_header()
{
    // Remove duplicated rows
    Mat E_unique = filter_unique_rows(m_magee_glmmkin.E, m_glmmkin_fitnull.id_include);

    //Apply function to check if unique elements in each column are <= 20
    std::ext::V_bool Ebin = apply_on_columns(E_unique, unique_less_equal, m_cmd.cat_threshold);
	
    if (std::any_of(Ebin.begin(), Ebin.end(), [](bool b) { return b;}))
    {
        //  Ecat and further operations if any column has <= 20 unique elements
        Mat Ecat(E_unique.rows(), std::count(Ebin.begin(), Ebin.end(), true));
        int col_idx = 0;
        for (int col = 0; col < E_unique.cols(); ++col) 
        {
            if (Ebin[col]) 
            {
                Ecat.col(col_idx++) = E_unique.col(col);
            }
        }
        // Create strata by concatenating values in each row of Ecat
        std::ext::VV_string Ecat_str(Ecat.rows(), std::ext::V_string(Ecat.cols()));
        for (int row = 0; row < Ecat.rows(); ++row) 
        {
            for (int col = 0; col < Ecat.cols(); ++col) 
            {
                Ecat_str[row][col] = std::to_string(Ecat(row, col));
            }
        }

        // Concatinate values of each column 0-1 
        std::ext::V_string strata = create_strata(Ecat_str);
        std::ext::V_string uni_strata = unique_id(strata);
        // It help to keep the last occurance of each element
        std::reverse(uni_strata.begin(), uni_strata.end());
        std::sort(uni_strata.begin(), uni_strata.end());

        //kepp only headers that binary
        for (auto inter = 0; inter < m_interaction.size(); ++inter)
        {
            if(Ebin[inter])
            {
                m_interaction_new.emplace_back(m_interaction[inter]);
            }
        }

        // Concatenate m_interaction terms
        std::string cat_inter = accumulate(m_interaction_new.begin(), m_interaction_new.end(), std::string(),
            [](const std::string& a, const std::string& b) {
                return a.empty() ? b : a + "_" + b;
            });

        // Combine cat_inter with unique strata
        std::ext::V_string tmp1(uni_strata.size());
        for (size_t i = 0; i < uni_strata.size(); ++i) {
            tmp1[i] = cat_inter + "_" + uni_strata[i];
        }
        // Create bin header by appending "N_" and "AF_" to each unique strata combination
        
        for (const auto& str : tmp1) {
            m_bin_headers.push_back("N_" + str);
            m_bin_headers.push_back("AF_" + str);
        }

        // Create hash table labeling which observation belong to which strata using its index
        m_strata_list = generate_strata_list(strata);      
    }
	else 
	{
        // Handle case when no columns have <= 20 unique elements
		std::cout << "There is no strata\n";
        m_bin_headers.clear();  // Ensures it is an empty vector
        m_interaction_new.clear();
    }
}

void MAGEE::printOutputHeader_magee() 
{ 
    std::ofstream results(m_cmd.outFile, std::ofstream::binary); 
    bool printFull = false; 
    bool printMeta = false; 
    int printStart = 1;  //escape first index as we will add G for beta_G
    int printEnd  = m_interaction_exp.size() + 1;  
    if (m_cmd.outStyle.compare("meta") == 0) { 
        printStart = 0;  
        printEnd   = m_interaction_exp.size() + 1; 
        printMeta  = true; 
    } else if (m_cmd.outStyle.compare("full") == 0) { 
        printStart = 0;  
        printEnd   = m_interaction_exp.size() + 1;  
        printFull  = true; 
        results << "#dispersion: " << m_glmmkin_fitnull.sigma2 << "\n"; //sigma2 is comming from fitNullModel
    } 
 
    results << "SNPID" << ((m_cmd.useBgenFile) ? "\tRSID\t" : "\t") << "CHR" << "\t" << "POS" << "\t" << "Non_Effect_Allele" << "\t" << "Effect_Allele" << "\t" << "N_Samples" << "\t" << "AF" << "\t"; 
    if (m_bin_headers.size() > 0)  
    { 
        for (size_t i = 0; i < m_bin_headers.size(); i++)  
        { 
            results <<  m_bin_headers[i] << "\t"; 
        } 
    } 
 
    for (int i = 0; i < m_interaction_exp.size(); i++) 
    { 
        m_interaction_exp[i] = "G-" + m_interaction_exp[i]; 
    } 
    m_interaction_exp.insert(m_interaction_exp.begin(), "G"); //To create beta_G
 
    string seMHeader = "SE_Beta_Marginal"; 
    string seHeader  = "SE_Beta_"; 
    string covHeader = "Cov_Beta_"; 
    if (m_cmd.robust == 1) { 
        seMHeader = "robust_" + seMHeader; 
        seHeader  = "robust_" + seHeader; 
        covHeader = "robust_" + covHeader; 
    } 
 
    results << "Beta_Marginal" << "\t" << seMHeader << "\t"; 

    if ((m_cmd.robust == 1) && (printMeta || printFull)) 
    { 
        results << "SE_Beta_Marginal" << "\t"; 
    } 

    if (m_interaction_exp.size() != 0)  
    { 
        for (int i = printStart; i < printEnd; i++) 
        { 
            results << "Beta_" << m_interaction_exp[i] << "\t"; 
        } 

        for (int i = printStart; i < printEnd; i++) 
        { 
            results << seHeader << m_interaction_exp[i] << "\t";   
        } 

        for (int i = printStart; i < printEnd; i++) 
        { 
            for (int j = printStart; j <  printEnd; j++) 
            { 
                if (i < j) 
                { 
                    results << covHeader << m_interaction_exp[i] << "_" << m_interaction_exp[j] << "\t";   
                }  
            } 
        } 
  
        if (m_cmd.robust == 1) 
        { 
            if (printMeta || printFull) { 
                for (int i = printStart; i < printEnd; i++) 
                { 
                    for (int j = printStart; j < printEnd; j++) 
                    { 
                        if (i == j) 
                        { 
                            results << "SE_Beta_" << m_interaction_exp[j] << "\t";  
                        } 
                    } 
                } 

                for (int i = printStart; i < printEnd; i++) 
                { 
                    for (int j = printStart; j < printEnd; j++) 
                    { 
                        if (i < j) 
                        { 
                            results << "Cov_Beta_" << m_interaction_exp[i] << "_" << m_interaction_exp[j] << "\t";  
                        } 
                    } 
                } 
 
                results << "robust_P_Value_Marginal" << "\t" << "robust_P_Value_Interaction" << "\t" << "robust_P_Value_Joint" << "\t"; 
                results << "P_Value_Marginal" << "\t" << "P_Value_Interaction" << "\t" << "P_Value_Joint\n"; 
            } 
            else 
            { 
                results << "robust_P_Value_Marginal" << "\t" << "robust_P_Value_Interaction" << "\t" << "robust_P_Value_Joint\n"; 
            } 
        }
        else 
        { 
            results << "P_Value_Marginal" << "\t" << "P_Value_Interaction" << "\t" << "P_Value_Joint\n"; 
        } 
    } 
    else 
    { 
        if (m_cmd.robust == 1) 
        { 
            if (printMeta || printFull) 
            { 
                results << "robust_P_Value_Marginal" << "\t" << "P_Value_Marginal\n"; 
            } 
            else
            { 
                results << "robust_P_Value_Marginal\n"; 
            } 
             
        } 
        else 
        { 
            results << "P_Value_Marginal\n"; 
        }          
    } 
    results.close();
}


void MAGEE::fitglmm()
{
    // m_magee_glmmkin.glmmkin = m_glmmkin_fitnull;
    Magee_Arma magee_arma;
	int qi = m_interaction_cov.size();//onlycovinteraction
	int ei = m_interaction_exp.size();//only exposure interaction
    m_interaction = m_interaction_exp;//both cov and exp interactions
	
    for (int i = 0; i < m_interaction_cov.size(); i++)
	{ 
		m_interaction.insert(m_interaction.end(), m_interaction_cov[i]); 
	}
	
	std::string pheno_missing_key = m_cmd.missing;//NA values
    bool meta_output = (m_cmd.outStyle == "meta") ? true : false;//output style
    double miss_cutoff = m_cmd.missGenoRate;
	int covar_center = m_cmd.center;
    int nperbatch = m_cmd.stream_snps;
    double minmaf = m_cmd.MAF;
    create_E();
    std::ext::V_int match_id;  
    std::ext::V_string sample_id; 
	int samSize =  m_gmmat->m_vkins_sp[0].pheno.m_data_frame.m_nrows;

    if(m_active_genotype == GENOTYPE::Bgen)
    {
        sample_id = std::get<Bgen>(m_genotype).sampleID_all;//sampleID; 
        std::get<Bgen>(m_genotype).getPositionOfBgenVariant(std::get<Bgen>(m_genotype), m_cmd);
    }
 
    if(m_active_genotype == GENOTYPE::Pgen)
    {
        sample_id = std::get<Pgen>(m_genotype).sampleID_all;//sampleID;
        std::get<Pgen>(m_genotype).getPgenVariantPos(std::get<Pgen>(m_genotype), m_cmd);
    }
    
    if(m_active_genotype == GENOTYPE::Bed)
    {
        sample_id = std::get<Bed>(m_genotype).sampleID_all;//sampleID;
        std::get<Bed>(m_genotype).getBedVariantPos(std::get<Bed>(m_genotype), m_cmd);
    }

    //fill select to be passed to glmm_gei_bgen13
    fill_sel(sample_id);
    //fill J if there are duplicated IDs
    if (any_duplicated(m_glmmkin_fitnull.id_include))
    {
        fill_J(match_id, sample_id);
        m_magee_glmmkin.dupflag = true;
    }
    else
    {
        std::ext::V_int match_id_with_minus_one = match_indices(sample_id, m_glmmkin_fitnull.id_include);
        //keep only indexes remove minus ones
        match_id = remove_minus_one(match_id_with_minus_one);
    }

    m_magee_glmmkin.E = slice_mat(m_magee_glmmkin.E, match_id);
    calculate_bin_header();

    if (covar_center == 1)
    {
        scale(m_magee_glmmkin.E, true, false);
    }
    else if (covar_center == 2) {
        if (!m_interaction_cov.empty()) {
            Mat left = m_magee_glmmkin.E.leftCols(ei);
            Mat right = m_magee_glmmkin.E.middleCols(ei, qi);
            scale(right, true, false);
            m_magee_glmmkin.E = cbind(left, right);
        }
    }
    m_magee_glmmkin.residuals = slice_vec(m_glmmkin_fitnull.scaled_residuals, match_id);
    m_glmmkin_fitnull.fit.sigma_ix = slice_mat(m_glmmkin_fitnull.fit.sigma_ix, match_id);
    m_magee_glmmkin.sigma_ix = m_glmmkin_fitnull.fit.sigma_ix.sparseView();
    m_glmmkin_fitnull.fit.sigma_ix.resize(0, 0);// free the memory
    m_magee_glmmkin.sigma_i = slice_mat(m_glmmkin_fitnull.fit.sigma_i, match_id, match_id, false);
    m_glmmkin_fitnull.fit.sigma_i.resize(0, 0);
    m_magee_glmmkin.cov = m_glmmkin_fitnull.fit.cov.sparseView();
    m_glmmkin_fitnull.fit.cov.resize(0, 0);// free the memory

	if(!m_interaction_cov.empty())
	{
		m_magee_glmmkin.EC = m_magee_glmmkin.E.middleCols(ei, qi);
	}
  
    {
        arma::sp_mat J(m_magee_glmmkin.J.rows(), m_magee_glmmkin.J.cols());
        arma::sp_mat sigma_i;
        arma::sp_mat sigma_ix;
        arma::vec residuals;
        arma::mat E;
        conver_eigen_to_arma(m_magee_glmmkin.J, J);
        conver_eigen_to_arma(m_magee_glmmkin.sigma_i, sigma_i);
        conver_eigen_to_arma(m_magee_glmmkin.sigma_ix, sigma_ix);
        conver_eigen_to_arma(m_magee_glmmkin.cov, magee_arma.cov);
        magee_arma.EC = arma::mat(m_magee_glmmkin.EC.data(), m_magee_glmmkin.EC.rows(), m_magee_glmmkin.EC.cols(), true);
        E = arma::mat(m_magee_glmmkin.E.data(), m_magee_glmmkin.E.rows(), m_magee_glmmkin.E.cols(), true);
        residuals = arma::vec(m_magee_glmmkin.residuals.data(), m_magee_glmmkin.residuals.size(), true);
        magee_arma.dupflag = m_magee_glmmkin.dupflag;
        magee_arma.select = m_magee_glmmkin.select;
        clear_m_magee_glmmkin();
        size_t n = residuals.size();
        size_t n_obs = residuals.size();
        size_t block_size = (ei + qi + 1) * n_obs;
        arma::mat resblock = kron(arma::ones(ei + qi + 1), residuals);
        arma::sp_mat Eblock(block_size, block_size);
        Eblock.submat(0, 0, n_obs - 1, n_obs - 1) = arma::speye<arma::sp_mat>(n_obs, n_obs);
        for (auto block_n = 1; block_n < ei + qi + 1; block_n++) 
        {
            // Create a sparse diagonal matrix from column block_n - 1 of E
            arma::sp_mat block = arma::sp_mat(arma::diagmat(E.col(block_n - 1)));
            Eblock.submat(block_n * n_obs, block_n * n_obs, 
                        (block_n + 1) * n_obs - 1, (block_n + 1) * n_obs - 1) = block;
        }

        arma::sp_mat all_ones_sigma_i(ei + qi + 1, ei + qi + 1);
        arma::sp_mat all_ones_sigma_ix(ei + qi + 1, 1);
        spa_mat_ones(all_ones_sigma_i);
        spa_mat_ones(all_ones_sigma_ix);
        arma::sp_mat sigma_iblock = kron(all_ones_sigma_i, sigma_i);
        arma::sp_mat sigma_ixblock = kron(all_ones_sigma_ix, sigma_ix);
        magee_arma.n_obs = residuals.size();
        magee_arma.n = residuals.size();
        if (magee_arma.dupflag)
        {
            magee_arma.Jres = J * residuals;
            arma::sp_mat Jblock = kron(arma::speye(ei + qi + 1, ei + qi + 1), J);
            magee_arma.JEresblock = Jblock * Eblock * resblock;
            magee_arma.Psi = Jblock * Eblock * sigma_iblock * Eblock * Jblock.t();
            magee_arma.Xi = Jblock * Eblock * sigma_ixblock;
            magee_arma.sigma_iJJ = J * sigma_i * J.t(); 
            magee_arma.sigma_ixJ = J * sigma_ix; 
            magee_arma.n = J.n_rows;
        }
        else
        {
            //There is no J for crosse-sectional data
            magee_arma.Jres = residuals;
            magee_arma.JEresblock = Eblock * resblock;
            magee_arma.Psi = Eblock * sigma_iblock * Eblock;
            magee_arma.Xi = Eblock * sigma_ixblock;
            magee_arma.sigma_iJJ = sigma_i ; 
            magee_arma.sigma_ixJ = sigma_ix;
        }
    }

    if(m_active_genotype == GENOTYPE::Bgen)
    {
        std::string bgenfile = m_cmd.bgenFile;
        if (std::get<Bgen>(m_genotype).threads > 1) 
        {
            std::cout << "Running multithreading...\n";
            std::vector<std::thread> threads;
            for (uint i = 0; i < std::get<Bgen>(m_genotype).threads; i++) 
            {
                threads.emplace_back(std::thread(&glmm_gei_bgen13, std::cref(magee_arma),
                bgenfile, m_cmd.outFile, minmaf, miss_cutoff, nperbatch, ei, qi,
                std::cref(m_strata_list),std::get<Bgen>(m_genotype).Mbgen_begin[i], std::get<Bgen>(m_genotype).Mbgen_end[i],
                std::get<Bgen>(m_genotype).bgenVariantPos[i], std::get<Bgen>(m_genotype).Nbgen, std::get<Bgen>(m_genotype).CompressedSNPBlocks, meta_output));
            }
            
            std::cout << "Continuing GEI test... joining threads...\n";
            for (auto& thread : threads) 
            {
                thread.join(); 
            } 
        } 
        else 
        {
            std::cout << "Running with single thread...\n";
            glmm_gei_bgen13(magee_arma, bgenfile, m_cmd.outFile, minmaf, 
            miss_cutoff, nperbatch, ei, qi, m_strata_list,  std::get<Bgen>(m_genotype).Mbgen_begin[0],
            std::get<Bgen>(m_genotype).Mbgen_end[0], std::get<Bgen>(m_genotype).bgenVariantPos[0], std::get<Bgen>(m_genotype).Nbgen, std::get<Bgen>(m_genotype).CompressedSNPBlocks,
            meta_output);
        }     
    }

    if(m_active_genotype == GENOTYPE::Pgen)
    {
        std::string pgenfile = m_cmd.pgenFile;
        std::string pvarFile = m_cmd.pvarFile;
        if (std::get<Pgen>(m_genotype).threads > 1) 
        {
            std::cout << "Running multithreading...\n";
            std::vector<std::thread> threads;
            for (uint i = 0; i < std::get<Pgen>(m_genotype).threads; i++) 
            {
                threads.emplace_back(std::thread(&glmm_gei_pgen13, std::cref(magee_arma),
                pgenfile, pvarFile, m_cmd.outFile, minmaf, miss_cutoff, nperbatch, ei, qi,
                std::cref(m_strata_list),std::get<Pgen>(m_genotype).begin[i],
                std::get<Pgen>(m_genotype).end[i], std::get<Pgen>(m_genotype).pgenVariantPos, 
                std::get<Pgen>(m_genotype).filterVariants, std::get<Pgen>(m_genotype).pvarIndex.size(),
                std::get<Pgen>(m_genotype).pvarLast, std::get<Pgen>(m_genotype).pvarIndex, 
                meta_output));
            }
            
            std::cout << "Continuing GEI test... joining threads...\n";
            for (auto& thread : threads) 
            {
                thread.join(); 
            } 
        } 
        else 
        {
            std::cout << "Running with single thread...\n";
            glmm_gei_pgen13(magee_arma, pgenfile, pvarFile, m_cmd.outFile, minmaf, 
            miss_cutoff, nperbatch, ei, qi, m_strata_list, std::get<Pgen>(m_genotype).begin[0],
            std::get<Pgen>(m_genotype).end[0], std::get<Pgen>(m_genotype).pgenVariantPos, 
            std::get<Pgen>(m_genotype).filterVariants, std::get<Pgen>(m_genotype).pvarIndex.size(),
            std::get<Pgen>(m_genotype).pvarLast, std::get<Pgen>(m_genotype).pvarIndex, 
            meta_output);
        }
    }

    if(m_active_genotype == GENOTYPE::Bed)
    {
        std::string bedfile = m_cmd.bedFile;
        std::string bimFile = m_cmd.bimFile;
        if (std::get<Bed>(m_genotype).threads > 1) 
        {
            std::cout << "Running multithreading...\n";
            std::vector<std::thread> threads;
            for (uint i = 0; i < std::get<Bed>(m_genotype).threads; i++) 
            {
                threads.emplace_back(std::thread(&glmm_gei_bed13, std::cref(magee_arma),
                bedfile, bimFile, m_cmd.outFile, minmaf, miss_cutoff, nperbatch, ei, qi,
                std::cref(m_strata_list),std::get<Bed>(m_genotype).begin[i],
                std::get<Bed>(m_genotype).end[i], std::get<Bed>(m_genotype).bedVariantPos, 
                std::get<Bed>(m_genotype).filterVariants, std::get<Bed>(m_genotype).bimDelim, 
                std::get<Bed>(m_genotype).bimLast, std::get<Bed>(m_genotype).n_samples, meta_output));
            }
            
            std::cout << "Continuing GEI test... joining threads...\n";
            for (auto& thread : threads) 
            {
                thread.join(); 
            } 
        } 
        else 
        {
            std::cout << "Running with single thread...\n";
            glmm_gei_bed13(magee_arma, bedfile, bimFile, m_cmd.outFile, 
            minmaf, miss_cutoff, nperbatch, ei, qi,std::cref(m_strata_list),
            std::get<Bed>(m_genotype).begin[0], std::get<Bed>(m_genotype).end[0],
            std::get<Bed>(m_genotype).bedVariantPos, std::get<Bed>(m_genotype).filterVariants,
            std::get<Bed>(m_genotype).bimDelim, std::get<Bed>(m_genotype).bimLast,
            std::get<Bed>(m_genotype).n_samples, meta_output);
        }
    }

	// Write all results from each thread to 1 file
    cout << "Combining results... \n";
	printOutputHeader_magee();
	std::cout << std::flush;
	std::ofstream results(m_cmd.outFile, std::ios::binary | std::ios_base::app);
    for (int i = 0; i < m_cmd.threads; i++) 
	{
        std::string threadOutputFile;
        if(m_active_genotype == GENOTYPE::Bgen)
        {
            threadOutputFile = m_cmd.outFile + "_bin_" + std::to_string(std::get<Bgen>(m_genotype).Mbgen_begin[i]) + ".tmp";
        }
        
        if(m_active_genotype == GENOTYPE::Pgen)
        {
            threadOutputFile = m_cmd.outFile + "_bin_" + std::to_string(std::get<Pgen>(m_genotype).begin[i]) + ".tmp";
        }

        if(m_active_genotype == GENOTYPE::Bed)
        {
            threadOutputFile = m_cmd.outFile + "_bin_" + std::to_string(std::get<Bed>(m_genotype).begin[i]) + ".tmp";
        }

        std::ifstream thread_output(threadOutputFile);
        if (thread_output.peek() != std::ifstream::traits_type::eof())
	    {
           results<<thread_output.rdbuf();
        }
        thread_output.close();
        std::remove(threadOutputFile.c_str());
    }
    results.close();
}



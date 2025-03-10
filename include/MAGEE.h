#pragma once

#include "GMMAT.h"
#include "declars.h"
#include "ReadBGEN.h"
#include <numeric>
#include <thread>

// change it to define object og glmmkin in magee_glmmkin
/**
 * @brief struct to include magee data
 * 
 */
struct Magee_Glmmkin
{
    // Glmmkin &glmmkin;
	Mat E;
	Mat EC;
    SpaMat J;
    std::ext::V_int select;
    bool dupflag = false;
    //new sigma_ix which is sparse
    SpaMat sigma_i;
    SpaMat sigma_ix;
    SpaMat cov;
    DensVec residuals;
    Mat JEres; 
	SpaMat JPJ;
	SpaMat JEPJ;
	SpaMat JEPEJ;
};

struct Magee_Arma {
    arma::mat EC;
    std::ext::V_int select;
    bool dupflag = false;
    arma::sp_mat cov;
    size_t n_obs;
    size_t n;
    arma::mat Jres;
    arma::mat JEresblock;
    arma::sp_mat Psi;
    arma::sp_mat Xi;
    arma::sp_mat sigma_iJJ;
    arma::sp_mat sigma_ixJ;
};

/**
 * @brief A class to run GEI test.
 * 
 */
class MAGEE
{
    public:
        enum class GENOTYPE {Bgen, Pgen, Bed};
        MAGEE() = default;
        MAGEE(GMMAT &gmmat, Glmmkin &glmmkin,  CommandLine &cmd, Bgen bgen,
                std::ext::V_string interaction_exp, std::ext::V_string interaction_cov,  
                int numSelCol); 
        MAGEE(GMMAT &gmmat, Glmmkin &glmmkin,  CommandLine &cmd, Pgen pgen,
                std::ext::V_string interaction_exp, std::ext::V_string interaction_cov,  
                int numSelCol); 
        MAGEE(GMMAT &gmmat, Glmmkin &glmmkin,  CommandLine &cmd, Bed bed,
                std::ext::V_string interaction_exp, std::ext::V_string interaction_cov,  
                int numSelCol); 
        /**
         * @brief Run GEI based on glmm 
         * 
         */
        void fitglmm();
    private:
        GMMAT *m_gmmat;
        Magee_Glmmkin m_magee_glmmkin;
        Glmmkin &m_glmmkin_fitnull;  
		CommandLine &m_cmd;
        // Bgen &m_bgen;
        std::variant<Bgen, Pgen, Bed> m_genotype;
        GENOTYPE m_active_genotype;
        std::ext::V_string m_interaction;
        std::ext::V_string m_interaction_exp;
		std::ext::V_string m_interaction_new;
        std::ext::V_string m_interaction_cov;
        int m_numSelCol;
        std::ext::Map_str_Vint m_strata_list;
		std::ext::V_string m_bin_headers;
        /**
         * @brief free memory for m_magee_glmmkin
         * 
        */
        void clear_m_magee_glmmkin();
        /**
         * @brief Fill matrix E in the `Magee_Glmmkin` using environmental data coming from GMMAT
         * 
        */
        void create_E();

        /**
         * @brief Fill the selection vector based on the provided sample IDs.
         * 
         * This function fills the `select` vector in the `Magee_Glmmkin` struct 
         * with the appropriate indices based on the given `sample_id`.
         * 
         * @param sample_id Vector of sample IDs used to determine the selection.
         * @see Magee_Glmmkin::select
        */
        void fill_sel(std::ext::V_string& sample_id);
        
        /**
         * @brief Fills a sparse matrix with specified values based on column indices.
         * 
         * @param mat The sparse matrix to fill.
         * @param numRows The number of rows in the matrix.
         * @param numCols The number of columns in the matrix.
         * @param col_indices The vector of column indices.
         * @param value The value to fill in the matrix.
        */
        void fill_mat(const int numRows, const int numCols, std::ext::V_int& indixes_col, int value);

        /**
         * @brief  Fill the J matrix based on matched IDs and sample IDs.
         * 
         * This function fills the `J` matrix in the `Magee_Glmmkin` struct 
         * it fills only when we have duplicated IDs.
         * @param match_id 
         * @param sample_id 
        */
        void fill_J(std::ext::V_int& match_id, std::ext::V_string& sample_id);

        /**
         * @brief Generates binary data headers based on the E in the `Magee_Glmmkin` interaction and interaction_cov terms.
         * 
         * This function processes the environmental matrix `E` to identify and categorize binary data columns.
         * Specifically, it performs the following steps:
         * 
         * - **Filtering Unique Rows**: Removes duplicate rows from the environmental matrix `E` to ensure that each observation is unique.
         * - **Binary Column Detection**: Identifies columns in the matrix with a number of unique values less than or equal to a specified threshold (e.g., 20). These columns are considered binary or categorical.
         * - **Strata Creation**: For columns identified as binary, creates a strata by concatenating the values of these columns for each observation. The strata are then used to categorize observations.
         * - **Interaction Term Filtering**: Filters the interaction terms to retain only those corresponding to the binary columns. This creates new interaction terms used for further analysis.
         * - **Header Generation**: Constructs binary data headers by combining the interaction terms with unique strata values. The headers include prefixes such as "N_" and "AF_" to indicate different types of data.
         * - **Strata List Generation**: Creates a mapping of observations to their corresponding strata, enabling efficient categorization and lookup.
         * 
         * The function handles cases where no columns are identified as binary by clearing the relevant data structures.
         * 
         * @note This function is crucial for setting up the binary data headers required for subsequent analysis steps in the MAGEE framework.
        */
        void calculate_bin_header();

        /**
         * @brief This function write the headers in the output file
         * 
         */
		void printOutputHeader_magee();
};


void conver_eigen_to_arma( SpaMat const& eigenMat, arma::sp_mat& armaMat);

/**
 * @brief Matches indices from the original vector to the filtered vector.
 * 
 * This function returns a vector of indices that match the elements in the original vector
 * with those in the filtered vector. If an element in the original vector does not have a match,
 * the corresponding index is set to -1.
 * 
 * @param original The original vector of strings.
 * @param filtered The filtered vector of optional strings.
 * @return A vector of matched indices or -1 for non-matches.
*/
std::ext::V_int match_indices(std::ext::V_string const& original, std::ext::V_opt_string const& filtered);

/**
 * @brief Matches indices from the original vector to the filtered vector.
 * 
 * Similar to the above `match_indices`, but this version accepts non-optional strings
 * in both the original and filtered vectors. If an element in the original vector does not 
 * have a match, the corresponding index is set to -1.
 * 
 * @param original The original vector of strings.
 * @param filtered The filtered vector of strings.
 * @return A vector of matched indices or -1 for non-matches.
*/
std::ext::V_int match_indices(std::ext::V_string const& original, std::ext::V_string const& filtered);

/**
 * @brief Extracts elements from a vector based on specified indices.
 * 
 * @tparam T The type of elements in the vector.
 * @param vec The original vector.
 * @param indices The indices of the elements to extract.
 * @return A vector of optional strings, with nullopt representing NA.
*/
template <typename T>
std::ext::V_opt_string slice(std::vector<T> const& vec, std::ext::V_int const& indices);

/**
 * @brief Checks if any element in the vector is NA (represented as -1).
 * 
 * @param vec The vector to check.
 * @return True if any element is -1, otherwise false.
*/
bool any_isna(std::ext::V_int vec);

/**
 * @brief Checks if elements of vec1 are in vec2.
 * 
 * @param vec1 The first vector of strings.
 * @param vec2 The second vector of strings.
 * @return A boolean vector indicating whether elements of vec1 are found in vec2.
*/
std::ext::V_bool in_op(std::ext::V_string const& vec1, std::ext::V_string const& vec2);

/**
 * @brief Returns a vector of unique elements from the input vector.
 * 
 * @tparam T The type of elements in the vector.
 * @param vec The input vector.
 * @return A vector of unique elements.
*/
template <typename T>
std::vector<T> unique_id(std::vector<T> const& vec);

/**
 * @brief Filters a vector of strings based on a boolean vector.
 * 
 * @param vec1 The vector to filter.
 * @param vec2 The boolean vector indicating which elements to keep.
 * @return A filtered vector of strings.
*/
std::ext::V_string filtered_ids(std::ext::V_string const& vec1, std::ext::V_bool const& vec2);

/**
 * @brief Checks if the input vector of strings is not empty.
 * 
 * @param sample_id The vector of sample IDs to check.
 * @throws std::runtime_error if the vector is empty.
 */
void check_not_empty(std::ext::V_string const& sample_id);


/**
 * @brief Checks if there are any duplicated elements in the vector.
 * 
 * @param vec The vector to check.
 * @return True if there are duplicates, otherwise false.
*/
bool any_duplicated(std::ext::V_string const& vec);

/**
 * @brief Identifies duplicates in a vector of strings.
 * 
 * @param vec The vector to check for duplicates.
 * @return A boolean vector indicating the presence of duplicates.
*/
std::ext::V_bool list_duplicates_bool(std::ext::V_string const& vec);

/**
 * @brief Filters out duplicate rows from a matrix.
 * 
 * This function removes duplicate rows from the input matrix and returns a matrix 
 * with only unique rows, as well as updating the `id_include` vector accordingly.
 * 
 * @param mat The input matrix.
 * @param id_include The vector of IDs to filter.
 * @return A matrix with unique rows.
*/
Mat filter_unique_rows(Mat const& mat, std::ext::V_string& id_include);

/**
 * @brief Applies a function to each column of a matrix and returns a boolean vector.
 * 
 * @param mat The matrix to process.
 * @param func The function to apply to each column.
 * @param threshold The threshold value used in the function.
 * @return A boolean vector indicating which columns meet the criteria.
*/
std::ext::V_bool apply_on_columns(Mat const& mat, const std::function<bool(VectorXd const&)>& func, int threshold);

/**
 * @brief Checks if the number of unique elements in a vector is less than or equal to a threshold.
 * 
 * @param vec The vector to check.
 * @param threshold The threshold value.
 * @return True if the number of unique elements is less than or equal to the threshold, otherwise false.
*/
bool unique_less_equal(DensVec const& vec, int threshold);

/**
 * @brief Generates a strata list mapping observations to their corresponding strata.
 * 
 * @param vec The vector of strata identifiers.
 * @return A map from strata identifiers to lists of observation indices.
 */
std::ext::Map_str_Vint generate_strata_list(std::ext::V_string const& vec);

/**
 * @brief Scales a matrix by centering and/or scaling its columns.
 * 
 * @param mat The matrix to scale.
 * @param center Whether to center the matrix by subtracting the mean of each column.
 * @param scale Whether to scale the matrix by dividing by the standard deviation of each column.
 */
void scale(Mat& mat, bool center = true, bool scale = true);

/**
 * @brief Combines two matrices column-wise.
 * 
 * This function concatenates two matrices `A` and `B` horizontally to create a new matrix.
 * 
 * @param A The first matrix.
 * @param B The second matrix.
 * @return The combined matrix.
 */
Mat cbind(Mat const& A, Mat const& B);

/**
 * @brief Matches included IDs with sample IDs.
 * 
 * This function returns the subset of IDs in `id_include` if they presnt in `sample_id`.
 * 
 * @param id_include A vector of IDs to include.
 * @param sample_id A vector of sample IDs.
 * @return A vector of matched IDs.
 */
std::ext::V_string match_id_include(std::ext::V_string const& id_include, std::ext::V_string const& sample_id);

/**
 * @brief Removes elements equal to -1 from the vector.
 * 
 * @param vec The vector to filter.
 * @return A vector with all elements that are not equal to -1.
*/
std::ext::V_int remove_minus_one(std::ext::V_int vec);

/**
 * @brief Creates strata by concatenating elements in each row of the binary columns of matrix E.
 * 
 * @param Ecat A vector of vector of strings representing categorical data.
 * @return A vector of concatenated strings representing strata.
*/
std::ext::V_string createStrata(std::ext::VV_string const& Ecat);

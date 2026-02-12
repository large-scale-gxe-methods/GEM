#pragma once

#include "SparseInverse.h"
#include "Kinship.h"
#include <optional>
#include <unordered_set>
#include <iterator>
#include <functional>
#include <variant>

const int MAX_N_ITER = 500; 


namespace std
{
    namespace ext
    {
        /**
         * @brief Function signature for the GEM null model fitting routine.
         *
         * This alias defines the callable interface used to invoke the
         * null model fitting implementation (e.g., `fitNullModel2`).
         * It returns results such as fitted values, residuals,
         * regression coefficients, and variance estimates.
         */
        using FitNull_f = std::function<void (int samSize, int numSelCol, int phenoType, double epsilon, 
                    int robust, std::vector<string> covSelHeadersName, std::vector<double> phenodata, 
                    std::vector<double> covdata, std::vector<double>* XinvXTX_ret, vector<double>* miu_ret, 
                    vector<double>* resid_ret, double* sigma2_ret, std::vector<double>& beta_ret,
                    std::vector<double>& Xbeta_ret)>;
        /**
         * @brief Variant type for storing dense or sparse matrices.
         *
         * This type is used to represent matrices that may be stored
         * either in dense format (`DensMat`) or sparse format (`SpaMat`).
         */
        using Matrix_variant = std::variant<DensMat, SpaMat>;
    }
}


/**
 * @brief A structure that contains all required parameters and functions related to the fitting function
 * and return parameters of the method fitglmm_ai in GMMAT class
 * 
 */
struct Fit
{
    DensVec eta;
    DensVec mu;
    DensVec dmu_deta;
    DensMat cov;
    DensVec alpha;
    SpaMat sigma_i; 
    DensVec diag_sigma_i;          
    SpaMat  diag_sigma_i_ZPchol; 
    DensMat sigma_ix;
    std::optional<DensVec> dtau;
    DensVec W;
    /**
     * @brief A function to calculate derivative of mu in respect of eta based on the family type
     * 
     * @param family_t 
     * @param size 
     */
    void calc_dmu_deta(std::string const& family_t, int size);
    /**
     * @brief A function to calculate the square root of W
     * 
     * @return DensVec 
     */
    DensVec calc_sqrtW();
};

/**
 * @brief A structure containing Fit and other variables to pass between GMMAT methods, also return GMMAT object type
 * 
 */
struct Glmmkin
{
    DensVec residuals; // 
    DensVec scaled_residuals;
    std::ext::V_string id_include;
    bool converged;
    Fit fit;
	double sigma2;
    bool run_wb = false;
};


/**
 * @brief A structure to interconnect GEM with GMMAT, especially parameters in Fit structure for the null model.
 * 
 */
struct  GEMFit
{
    std::ext::V_double XinvXTX; 
    std::ext::V_double mu; 
    std::ext::V_double resid; 
    double sigma2; // Keep gf.sigma2 from fitNullModel2
    std::ext::V_double alpha;
    std::ext::V_double eta;
    /**
     * @brief A function to convert GEMFit to Fit structure data type.
     * 
     * @return Fit 
     */
    Fit convert_2_fit();
};

enum ModelType 
{
    LONGITUDINAL_RI,     // Random intercept only
    LONGITUDINAL_RS      // Random slope + intercept
};


/**
 * @brief A class to run association test.
 * 
 */
class GMMAT
{
    public:
        int m_groups = 0;
        int m_n_sel_col;
        std::vector<SparseInverse> m_vkins_sp;
        std::ext::VV_int m_covariance_idx ;
        DensVec m_sqrtW; //sqrtW is the weight matrix
        DensVec m_Y; // working vector
        std::ext::V_int m_group_id;
        std::unordered_map<int, std::ext::V_int> m_group_idx;
        SpaMat m_J;
        SpaMat m_Psi;
        SpaMat m_Z;
        SpaMat m_Zt_Z;
        size_t m_N;
        size_t m_Nobs;
        bool m_dup;
        int m_kins_size;
        DensVec m_tau;
        std::ext::V_int m_fixtau;
        std::ext::V_int m_fixrho;
        std::ext::V_int m_fixrho_idx;
        std::string m_family_t = "binomial";
        std::string m_link = "logit";
        DensVec m_offset;
        DensVec m_y;
        DensMat m_X; 
        DensVec m_rand_slope;
        std::ext::map_str_int m_hdrsMap;
        int m_robust = 0;
        std::ext::V_int m_idxtau;
        std::ext::V_int m_idxtau2;
        ModelType m_modeltype;
        SpaMat m_diag_sigma_im_Z;
        /**
         * @brief Construct the covariance matrix for the random effects.
         *
         * This function builds the random-effect covariance matrix 
         * used in Woodbury longitudinal generalized linear mixed models (GLMMs).
         * ---
         * **1. Random Intercept (RI) model**
         *
         * where:
         * - \f$\Phi_i\f$ are sparse kinship matrices
         * - \f$\theta_i\f$ are variance component parameters stored in `m_tau`
         * - \f$I\f$ is the identity matrix
         *
         * The resulting matrix has dimension \f$N \times N\f$.
         *
         * 
         * **2. Random Intercept + Random Slope (RS) model**
         *
         * where:
         * - \f$\Phi_i\f$ are sparse kinship matrices
         * - \f$\theta_i\f$ are variance component parameters stored in `m_tau`
         * - \f$I\f$ is the identity matrix
         * 
         * The resulting matrix has dimension \f$2N \times 2N\f$.
         * 
         * 
         * **Implementation Notes**
         *
         * - The matrix is assembled efficiently using sparse triplets.
         *
         * @param ng Offset index in `m_tau` where variance components begin.   
         * @param calc_diag_kin To calculate the mean of the diagonal of kinsip   
         * @note The resulting matrix is stored in the member variable `m_Psi`.
         */

        void build_Psi(int const ng, bool calc_diag_kin = false);
        /**
         * @brief Construct the random-effects design matrix.
         *
         * This function builds the sparse design matrix that links
         * random effects to observations in Woodbury longitudinal GLMM models.
         *
         * The structure depends on the model type:
         *
         * ---
         * ## 1. Random Intercept (RI) model
         * where:
         * - \f$J\f$ is the subject-indicator matrix of size \f$N_{obs} \times N\f$
         * - each row assigns an observation to its corresponding subject
         *
         * ---
         * ## 2. Random Intercept + Random Slope (RS) model
         *
         * For `LONGITUDINAL_RS`, the design matrix includes both intercept
         * and slope random effects:
         *
         * where:
         * - the first block represents J
         * - the second block represents random slope effects scaled by J
         *
         * The resulting matrix has dimension:
         *
         * \f$N_{obs} \times 2N\f'
         * ---
         * ## Implementation Notes
         *
         * - The matrix is assembled using sparse triplets for efficiency.
         * - The subject-indicator matrix `m_J` must already be initialized.
         *
         * @note The resulting matrix is stored in the member variable `m_Z`.
         */
        void build_Z();

        /**
         * @brief Fits the null GLMM model using Average Information (AI) algorithm.
         *
         * @param W Working weight vector
         * @return Fit Structure 
         */
        Fit fitglmm_ai(DensVec const& W);
       /**
         * @brief Fits the  GLMM using AI iterations starting from a null model.
         *
         * Iteratively updates:
         *  - fixed-effect coefficients (alpha)
         *  - variance components (tau)
         *  - working response, weights, and residuals
         * until convergence.
         *
         * Supports both Woodbury-based and direct sparse implementations.
         *
         * @param fit_null Initial null model fit.
         * @param maxiter Maximum number of iterations.
         * @param tol Convergence tolerance.
         * @return Glmmkin object.
         */

        Glmmkin glmmkin_ai(Fit fit_null, int maxiter = 500, double tol = 1e-5);
        /**
         * @brief Fits the GLMM using the selected optimization method (currently AI-REML).
         *
         * Prepares group structure and variance components, runs AI-based fitting,
         * and refits automatically if variance estimates hit parameter boundaries.
         *
         * @param fit_null Initial null model fit.
         * @param group_id Group membership for heteroscedastic modeling.
         * @param method Estimation method ("REML" or "ML").
         * @param method_optim Optimization method ("AI").
         * @param maxiter Maximum number of iterations.
         * @param tol Convergence tolerance.
         * @param tau_min Lower bound for variance components.
         * @param tau_max Upper bound for variance components.
         * @param tau_region 
         * @return Fitted GLMM model object.
         */
 
        Glmmkin glmmkin_fit(Fit fit_null, std::ext::V_int group_id, 
                            std::string const method = "REML", 
                            std::string method_optim = "AI", 
                            int maxiter = 500,
                            double tol = 1e-5, double tau_min = 1e-5, 
                            double tau_max = 1e+5, int tau_region = 10);
        /**
         * @brief Initializes and fits a GLMM starting from phenotype and covariate data.
         *
         * @param fit0 Function to compute the initial null model.
         * @param pheno Phenotype data container.
         * @param cov_selected_hdrs Selected covariate column names.
         * @param phenoname Phenotype column name.
         * @param id Sample ID column.
         * @param randomSlopeName Optional random slope variable.
         * @param group Grouping variable for heteroscedastic models.
         * @param method Estimation method ("REML" or "ML").
         * @param method_optim Optimization method ("AI").
         * @param maxiter Maximum number of iterations.
         * @param tol Convergence tolerance.
         * @param tau_min Minimum variance component value.
         * @param tau_max Maximum variance component value.
         * @param tau_region 
         * @return Fitted GLMM object.
         */

        [[nodiscard]] Glmmkin glmmkin_init(std::ext::FitNull_f fit0, Pheno pheno,
                            std::ext::V_string cov_selected_hdrs,
                            std::string phenoname,
                            std::string const& id, 
                            std::ext::V_string int_cov_hdrs_name_new,
                            std::string randomSlopeName,
                            std::string const& groups,
                            int center, int scale,
                            std::string const method = "REML", 
                            std::string method_optim = "AI", 
                            int maxiter = 500,
                            double tol = 1e-5, double tau_min = 1e-5, 
                            double tau_max = 1e+5, int tau_region = 10);
    
    private:
        void fill_mat(const int numRows, std::ext::V_int& indixes_col, int value);
        void fill_J(std::ext::V_string const& sample_ids);
        void set_ai_low_ng(int i, DensVec& score, DensMat& ai, DensVec const& wpy, Fit const& fit, 
                            DensVec const& py, DensVec diagp, DensMat sigma_ixcov);
        void set_ai_high_ng(int i, DensVec& score, DensMat& ai, DensVec const& wpy, Fit const& fit,
                             DensVec const& py, DensMat const& sigma_ixcov, 
                             int ng);
        void set_ai(DensVec& score, DensMat& ai, DensVec const& wpy, Fit const& fit, DensVec const& py, 
                    DensVec diagp, DensMat sigma_ixcov, int ng, int q2);
        bool any_negative();
        bool any_negative(std::ext::V_int vec);
        bool any_nonzero();
        void update_tau(DensVec const& tau0, DensVec const& dtau);
        void update_tau_below_tol(DensVec const& tau0, double tol);
        void update_tau_below_tol2(DensVec const& tau0, double tol);
        void extract_group_idx(std::unordered_set<int> const& group_unique, std::ext::V_int group_id);
        int check_pheno_type();
        /**
         * @brief Calculate the random slope, intercept and their covariance and fill m_covariance_idx
         * 
         * @param kins_size 
         * @param ng 
         */
        void calc_covariance(int &kins_size, int ng);
        void set_dspy(DensMat& dspy, DensVec const& wpy, int ng);
        void set_VZpy(DensMat& VZpy, DensVec const& Zpy, int nk);
        void set_mtau(DensVec V_tr_corr, DensVec tau0, SpaMat& Ztsigma_iZ, 
            DensMat const& Ztsigma_ix, DensMat const& Ztsigma_ixcov, 
            DensVec diagp, int nk, int dimZ, int ng, bool has_random_slope);
        void calc_tr_corr(int i, DensVec &score, DensMat const& Ztsigma_ix, 
            DensMat const& Ztsigma_ixcov, SpaMat& Ztsigma_iZ, 
            int nk, int dimZ, int const& ng, bool has_random_slope);
        void set_score(DensVec &score, DensMat const& Ztsigma_ix, 
            DensMat const& Ztsigma_ixcov, SpaMat& Ztsigma_iZ, 
            DensVec const& diagp, int nk, 
            int dimZ, int const& ng, bool has_random_slope);

        void fill_fixrho_idx(double tol);
        void update_fixtau_fixrho(std::ext::V_int &fixtau_new,std::ext::V_int &fixrho_new, double tol);
        void fill_fixrho_idx0(std::ext::V_int &fixrho_idx0, DensVec const& tau0, double tol);
        void update_mtau_with_m_covariance_idx_tau0(DensVec const& tau0, double tol);
        void update_mtau_with_m_covariance_idx_mfixrho(std::ext::V_int const& idxrho);
        void update_mtau_with_m_covariance_mfixrho_idx_fixrho_idx0(std::ext::V_int const& fixrho_idx0);
        void update_mtau_with_m_covariance_mfixrho_idx();
        void update_mtau_with_m_covariance_idx(double tol);
        std::ext::V_int exclude_idx();
        bool covariate_larger_slope_intercept(double tol);
    };

/**
 * @brief Slice sparse matrix based on rows index
 * 
 * @param spm 
 * @param indices 
 * @return SpaMat 
 */

SpaMat slice_mat(SpaMat const& spm, std::ext::V_int const& indices);

DensMat slice_mat(DensMat const& dm, std::ext::V_int const& indices);

/**
 * @brief Slice dense matrix based on cols and rows index
 * 
 * @param Dense matrix 
 * @param std::vector<int> 
 * @param std::vector<int> 
 * @return Dens matrix 
 */

DensMat slice_mat(DensMat const& dm, std::ext::V_int ind1 , std::ext::V_int ind2);

/**
 * @brief Slice dense matrix based on column and row logical vectors
 * 
 * @param Dense matrix 
 * @param std::vector<bool> 
 * @param std::vector<bool> 
 * @return Dense matrix 
 */
DensMat slice_mat(DensMat const& dm, std::ext::V_bool const& ind1, std::ext::V_bool const& ind2);

/**
 * @brief Slice matrix based on columns and rows variant vectors(either int or bool)
 * 
 * @param Dense matrix 
 * @param std::variant<int, bool> 
 * @param std::variant<int, bool>  
 * @return Dense matrix 
 */
DensMat slice_mat(DensMat const& dm, std::ext::Var_bool_int const& ind1, std::ext::Var_bool_int const& ind2);

/**
 * @brief Slice dense matrix based on rows indexes
 * 
 * @param Dense matrix 
 * @param std::vector<bool> 
 * @return Dense matrix 
 */
DensMat slice_mat(const DensMat& dm, const std::ext::V_bool& indices);
/**
 * @brief Slice dense matrix based on variant row vector(either int or bool)
 * 
 * @param Dense matrix 
 * @param std::variant<int, bool> 
 * @return Dense matrix
 */
DensMat slice_mat(const DensMat& dm, const std::ext::Var_bool_int& indices);

/**
 * @brief Slice dense matrix based on column index vector
 * 
 * @param Dense matrix 
 * @param std::vector<int> 
 * @return Dense matrix 
 */
DensMat slice_mat_cols(DensMat const& dm, std::ext::V_int const& ind2);


/**
 * @brief Slice dense vector based on indexes
 * 
 * @param Dense vector
 * @param std::vector<int> 
 * @return Dense vector 
 */

DensVec slice_vec(DensVec const& dv, std::ext::V_int const& indices);

/**
 * @brief Slice dense vector based on logical vector 
 * 
 * @param Dense vector
 * @param std::vector<bool> 
 * @return Dense vector
 */
DensVec slice_vec(DensVec const& dv, std::ext::V_bool const& indices);

/**
 * @brief Slice dense vector based on variant vector(either int or bool)
 * 
 * @param Dense vector
 * @param std::variant<int, bool> 
 * @return Dense vector
 */
DensVec slice_vec(DensVec const& dv, std::ext::Var_bool_int const& indices);

/**
 * @brief Slice sparse matrix based on cols and rows index
 * 
 * @param Sparse matrix 
 * @param std::vector<int>  
 * @param std::vector<int>  
 * @return sparse matrix
 */
SpaMat slice_mat(SpaMat const& spm, std::ext::V_int ind1, std::ext::V_int ind2, bool check_size = true);

/**
 * @brief Slice sparse matrix based on column and row logical vectors
 * 
 * @param Sparse matrix 
 * @param std::vector<bool>  
 * @param std::vector<bool>  
 * @return sparse matrix
 */
SpaMat slice_mat(SpaMat const& spm, std::ext::V_bool const& ind1, std::ext::V_bool const& ind2, bool check_size = true);

/**
 * @brief Slice sparse matrix based on cols and rows variant(either int or bool)
 * 
 * @param Sparse matrix 
 * @param std::variant<int, bool>  
 * @param std::variant<int, bool>  
 * @return sparse matrix
 */
SpaMat slice_mat(SpaMat const& spm, std::ext::Var_bool_int const& ind1, std::ext::Var_bool_int const& ind2, bool check_size = true);


/**
 * @brief Returns a vector of indices for elements that satisfy the predicate.
 * 
 * @tparam T  type of the vector.
 * @tparam Pred predicate.
 * @param std::vector<T> input vector.
 * @param pred  prdicate to fill  a vector on integer values to return.
 * @return std::ext::V_int, Vector of indices where the predicate is true.
 */
template <typename T, typename Pred>
std::ext::V_int which(std::vector<T> const& vec, Pred pred)
{
    std::vector<int> v_idx;
    for(size_t i{0}; i < vec.size(); ++i)
    {
        if(pred(vec[i]))
        {
            v_idx.push_back(i);
        }
    }
    return v_idx;
}

/**
 * @brief Returns a vector of indices for elements that satisfy the predicate.
 * 
 * @tparam T type of the vector 
 * @tparam Pred predicate
 * @param DensVecInt 
 * @param pred 
 * @return std::ext::V_int 
 */
template<typename Pred>
std::ext::V_int which(DensVecInt const& vec, Pred pred)
{
    std::vector<int> v_idx;
    for(int i{0}; i < vec.size(); ++i)
    {
        if(pred(vec(i)))
        {
            v_idx.push_back(i);
        }
    }
    return v_idx;
}

/**
 * @brief Analogy of crossprod function in R language.
 * it gets two matrix from eigen library of any type (dens or sparse) 
 * and calculates the transpose of the first one and mulitplies it by the second one.
 * 
 * @tparam First  matrix from eigen library
 * @tparam Second  matrix from eigen libray
 * @param first first matrix from eigen library, it can be either dense or sparse.
 * @param second second matrix from eigen library, it can be either dense or sparse.
 * @return auto
 */
template<typename First, typename Second>
auto crossprod(First const& first, Second const& second) {
    return first.transpose() * second;
}

/**
 * @brief Analogy of tcrossprod function in R language.
 * 
 * it gets two matrix from eigen library of any type (dens or sparse) 
 * and mulitplies the first one by the transpose of the second one.
 * 
 * @tparam First matrix from eigen library
 * @tparam Second matrix from eigen libray 
 * @param first  first matrix from eigen library, it can be either dense or sparse. 
 * @param second  second matrix from eigen library, it can be either dense or sparse. 
 * @return auto 
 */
template<typename First, typename Second>
auto tcrossprod(First const& first, Second const& second) {
    return first * second.transpose();
}

/**
 * @brief Removes duplicate values from a vector.
 *
 * Constructs an unordered set from the input vector, keeping only unique elements.
 *
 * @tparam T Element type of the input vector.
 * @param vec Input vector.
 * @return std::unordered_set<T> Set containing the unique elements.
 */
template<typename T>
std::unordered_set<T> unique(std::vector<T> const& vec)
{
    std::unordered_set<T> u_set(vec.begin(), vec.end());
    return u_set;
}

/**
 * @brief Returns a vector of unique elements from the input vector.
 * 
 * @tparam T The type of elements in the vector.
 * @param vec The input vector.
 * @return std::vector<T> A vector of unique elements.
*/
template <typename T>
std::vector<T> unique_id(std::vector<T> const& vec) 
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

inline std::ext::V_int match_indices(std::ext::V_string const& original, std::ext::V_string const& filtered) {
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

#pragma once

#include "SparseInverse.h"
#include "Kinship.h"
// #include "GEM.h"
#include <optional>
#include <unordered_set>
#include <iterator>
#include <functional>
#include <variant>

const int MAX_N_ITER = 500; 

/**
 * @brief To use fitNullModel2 in GEM
 * 
 */
namespace std
{
    namespace ext
    {
        using FitNull_f = std::function<void (int samSize, int numSelCol, int phenoType, double epsilon, 
                    int robust, std::vector<string> covSelHeadersName, std::vector<double> phenodata, 
                    std::vector<double> covdata, std::vector<double>* XinvXTX_ret, vector<double>* miu_ret, 
                    vector<double>* resid_ret, double* sigma2_ret, std::vector<double>& beta_ret,
                    std::vector<double>& Xbeta_ret)>;
        using Matrix_variant = std::variant<Mat, SpaMat>;
    }
}

/**
 * @brief A structure that contains all required parameters and functions related to the fitting function
 * and return parameters of the method fitglmm_ai in GMMAT class
 * 
 */
struct Fit
{
    /**
     * @brief linear predictor, glmm
     * 
     */
    DensVec eta;
    /**
     * @brief conditioanl mean
     * 
     */
    DensVec mu;
    /**
     * @brief derivative of mu in respect of eta
     * 
     */
    DensVec dmu_deta;
    /**
     * @brief 
     * 
     */
    DensMat cov;
    /**
     * @brief  fixed covariate effect.
     * 
     */
    DensVec alpha;
    SpaMat sigma_i;
    DensMat sigma_ix;
    // std::vector<DataFrame> v_models;
    std::optional<DensVec> dtau;
    /**
     * @brief Weight of 
     * 
     */
    DensVec W;
    /**
     * @brief A function for calculating derivative of mu in respect of eta based on the family type
     * 
     * @param family_t : ei
     * @param size 
     */
    void calc_dmu_deta(std::string const& family_t, int size);
    /**
     * @brief A function to cacluate the square root of W
     * 
     * @return DensVec 
     */
    DensVec calc_sqrtW();
    //Fit& operator= (Fit const& fit);
};

/**
 * @brief A structure containg Fit and other variables pass within GMMAT methods and return by GMMAT object
 * 
 */
struct Glmmkin
{
    //int n_pheno, n_group;//might needed
    //SpaMat X;
    DensVec residuals; // 
    DensVec scaled_residuals;
    std::ext::V_string id_include;
    //std::vector<bool> converged;
    bool converged;
    Fit fit;
	double sigma2;

};


/**
 * @brief A structure for interconnect GEM with GMMAT, espcially parameters in Fit sturcture for the null model.
 * 
 */
struct  GEMFit
{
    std::ext::V_double XinvXTX; 
    std::ext::V_double mu; 
    std::ext::V_double resid; 
    double sigma2; // To return the gf.sigma2 from fitnull
    std::ext::V_double alpha;
    std::ext::V_double eta;
    /**
     * @brief A function to convert data and parameters in GEMFit to the ones in Fit structure.
     * 
     * @return Fit 
     */
    Fit convert_2_fit();
};


/**
 * @brief A class to run gene association test.
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
        /**
         * @brief The variance component parameters.
         */
        DensVec m_tau;
        std::ext::V_int m_fixtau;
        std::ext::V_int m_fixrho;
        std::ext::V_int m_fixrho_idx;
        std::string m_family_t = "binomial";
        std::string m_link = "logit";
        DensVec m_offset;
        /**
         * @brief data extracted from pheno files
         * 
         */
        DensVec m_y;
        /**
         * @brief Covariates data extracted from pheno files
         * 
        */
        DensMat m_X;
        /**
         * @brief Data to use for calculating random slope
         * 
        */
        DensVec m_rand_slope;
        std::ext::map_str_int m_hdrsMap;
        int m_robust = 0;
        std::ext::V_int m_idxtau;
        std::ext::V_int m_idxtau2;
        /**
         * @brief 
         * 
         * @param W 
         * @return Fit 
         */
        Fit fitglmm_ai(DensVec const& W);
        /**
         * @brief 
         * 
         * @param fit_null 
         * @param maxiter 
         * @param tol 
         * @return Glmmkin 
         */
        Glmmkin glmmkin_ai(Fit fit_null, int maxiter = 500, double tol = 1e-5);
         /**
          * @brief 
          * 
          * @param fit_null 
          * @param group_id 
          * @param method 
          * @param method_optim 
          * @param maxiter 
          * @param tol 
          * @param tau_min 
          * @param tau_max 
          * @param tau_region 
          * @return Glmmkin 
          */
        Glmmkin glmmkin_fit(Fit fit_null, std::ext::V_int group_id, 
                            std::string const method = "REML", 
                            std::string method_optim = "AI", 
                            int maxiter = 500,
                            double tol = 1e-5, double tau_min = 1e-5, 
                            double tau_max = 1e+5, int tau_region = 10);
        /**
         * @brief 
         * 
         * @param pheno 
         * @param id 
         * @param groups 
         * @param method 
         * @param method_optim 
         * @param maxiter 
         * @param tol 
         * @param tau_min 
         * @param tau_max 
         * @param tau_region 
         * @return Glmmkin 
         */
        [[nodiscard]] Glmmkin glmmkin_final(std::ext::FitNull_f fit0, Pheno pheno,
                            std::ext::V_string covSelectedHeader,
                            std::string phenoname,
                            std::string const& id, 
                            std::string randomSlopeName,
                            std::string const& groups,
                            std::string const method = "REML", 
                            std::string method_optim = "AI", 
                            int maxiter = 500,
                            double tol = 1e-5, double tau_min = 1e-5, 
                            double tau_max = 1e+5, int tau_region = 10);
    
    private:
        void set_ai_low_ng(int i, DensVec& score, DensMat& ai, DensVec const& wpy, Fit const& fit, 
                            DensVec const& py, DensVec diagp, DensMat sigma_ixcov);
        void set_ai_high_ng(int i, DensVec& score, DensMat& ai, DensVec const& wpy, Fit const& fit,
                             DensVec const& py, DensMat sigma_ixcov, int ng);
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
        void calc_rand_effect(int &kins_size, int ng);
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



// //Helper functions declaration
// /**
//  * @brief Convert vector of string to vector of double
//  * 
//  * @param std::vector<string>  
//  * @return std::vector<double> 
//  */
// std::ext::V_double convert_2_vector_of_double(const std::ext::V_string& strings);


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
 * @brief Slice dense matrix based on cols and rows logical vectors
 * 
 * @param Dense matrix 
 * @param std::vector<bool> 
 * @param std::vector<bool> 
 * @return Dens matrix 
 */
DensMat slice_mat(DensMat const& dm, std::ext::V_bool const& ind1, std::ext::V_bool const& ind2);

/**
 * @brief Slice matrix based on cols and rows variant vectors(either int or bool)
 * 
 * @param Dense matrix 
 * @param std::variant<int, bool> > 
 * @param std::variant<int, bool>  
 * @return Dens matrix 
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
 * @param sparse matrix 
 * @param std::vector<int>  
 * @param std::vector<int>  
 * @return sparse matrix
 */
SpaMat slice_mat(SpaMat const& spm, std::ext::V_int ind1, std::ext::V_int ind2, bool check_size = true);

/**
 * @brief Slice sparse matrix based on cols and rows logical vector
 * 
 * @param sparse matrix 
 * @param std::vector<bool>  
 * @param std::vector<boll>  
 * @return sparse matrix
 */
SpaMat slice_mat(SpaMat const& spm, std::ext::V_bool const& ind1, std::ext::V_bool const& ind2, bool check_size = true);

/**
 * @brief Slice sparse matrix based on cols and rows variant(either int or bool)
 * 
 * @param sparse matrix 
 * @param std::variant<int, bool>  
 * @param std::variant<int, bool>  
 * @return sparse matrix
 */
SpaMat slice_mat(SpaMat const& spm, std::ext::Var_bool_int const& ind1, std::ext::Var_bool_int const& ind2, bool check_size = true);


/**
 * @brief A function that returns vector of index for those values that meet the condition in predicate
 * 
 * @tparam T : type of the vector
 * @tparam Pred : predicate
 * @param std::vector<T>  const& : the vector we check the predicate to fill a vector of integer values
 * @param pred : prdicate to fill  a vector on integer values to return 
 * @return std::ext::V_int 
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
 * @brief A function that returns vector of index for those values that meet the condition in predicate
 * 
 * @tparam T : type of the vector 
 * @tparam Pred : predicate
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
 * @tparam First  : matrix from eigen library
 * @tparam Second : matrix from eigen libray
 * @param first  : first matrix from eigen library, it can be either dense or sparse.
 * @param second : second matrix from eigen library, it can be either dense or sparse.
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
 * @tparam First  : matrix from eigen library
 * @tparam Second : matrix from eigen libray 
 * @param first  : first matrix from eigen library, it can be either dense or sparse. 
 * @param second : second matrix from eigen library, it can be either dense or sparse. 
 * @return auto 
 */
template<typename First, typename Second>
auto tcrossprod(First const& first, Second const& second) {
    return first * second.transpose();
}


/**
 * @brief A function to remove the duplication in a given vector
 * 
 * @tparam T : type of the vector in the argument of the function and the return set.
 * @param vec 
 * @return std::unordered_set<T> 
 */
template<typename T>
std::unordered_set<T> unique(std::vector<T> const& vec)
{
    std::unordered_set<T> u_set(vec.begin(), vec.end());
    return u_set;
}

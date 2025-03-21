#include "../include/GMMAT.h"
#include <algorithm>
#include <limits>
#include <cmath>
#include <chrono>
#include <thread>
#include <mutex>

template<typename T>
int sign(T t)
{
    return (t > 0) - (t < 0);
}

std::ext::V_int intersect(std::ext::V_int const& vec1, std::ext::V_int const& vec2)
{
    std::ext::V_int intersect_vec;
    std::ext::V_int sorted_vec1 = vec1;
    std::ext::V_int sorted_vec2 = vec2;

    std::sort(sorted_vec1.begin(), sorted_vec1.end());
    std::sort(sorted_vec2.begin(), sorted_vec2.end());

    std::set_intersection(sorted_vec1.begin(), sorted_vec1.end(), sorted_vec2.begin(),
                         sorted_vec2.end(), std::back_inserter(intersect_vec));
    return intersect_vec;
}

// std::ext::V_int not_intersect(std::ext::V_int const& vec1, std::ext::V_int const& vec2)
// {
//     std::ext::V_int diff_vec;
//     std::ext::V_int sorted_vec1 = vec1;
//     std::ext::V_int sorted_vec2 = vec2;

//     std::sort(sorted_vec1.begin(), sorted_vec1.end());
//     std::sort(sorted_vec2.begin(), sorted_vec2.end());

//     std::set_difference(sorted_vec1.begin(), sorted_vec1.end(), sorted_vec2.begin(),
//                          sorted_vec2.end(), std::back_inserter(diff_vec));
//     return diff_vec;
// }


std::ext::V_int logic_update_fixed_condtion(DensVec dv, double tol)
{
    std::ext::V_int v_ints(dv.size(), 0);

    for(int i {0}; i < dv.size(); ++i)
    {
        if(dv(i) < 1.01 * tol)
        {
            v_ints[i] = 1;
        }
    }
    return v_ints;
}

DensVec conv_vec_vecXd(std::ext::V_double v_std)
{
    Eigen::Map<DensVec> eig_vec(v_std.data(), v_std.size());
    return eig_vec;
}


//Returns a diagonla  matrix
namespace details{
    SpaMat diag(std::ext::V_double vec)
    {
        SpaMat spm(vec.size(), vec.size());
        for(size_t i{0}; i < vec.size(); ++i)
        {
            spm.insert(i, i) = vec[i];
        }
        return spm;
    }
}

//Create matrix X covariate data by adding extra col fill with one and pheno data
DensMat create_covdata(const DataFrame& df) {
    // Get dimensions of the data frame
    int num_rows = df.n_rows();
    int num_cols = df.n_cols();
    
    // Create an Eigen matrix with the same dimensions
   DensMat dmat2ret(num_rows, num_cols + 1);
   for (int i = 0; i < num_rows; ++i) 
   {
        dmat2ret(i, 0) = 1.0;
   }
    // Populate the Eigen matrix from the data frame
    for (int i = 0; i < num_rows; ++i) {
        for (int j = 1; j < num_cols + 1; ++j) {
            
            try
            {
                dmat2ret(i, j) = std::stod(df.m_data.at(df.m_headers[j-1])[i]);  
            }
            catch(std::invalid_argument const&e)
            {
                std::cerr << "value" << df.m_data.at(df.m_headers[j-1])[i] << " is not a valid number \n"; 
            }

              
        }
    } 
    return dmat2ret;
}

//Return map to find headers in matrix for MAGEE
std::ext::map_str_int createHeaderMap(std::ext::V_string const& headers) {
    std::ext::map_str_int headerMap;
    //first column is intercept all filled with one so we start index from 1
    for (int i = 0; i < headers.size(); ++i) 
    {
        headerMap[headers[i]] = i + 1;
    }
    return headerMap;
}


std::ext::V_int conv_stdvs2stdvi(std::ext::V_string const& strings) {
    std::ext::V_int result;
    result.reserve(strings.size());
    for (const std::string& str : strings) {
        result.push_back(std::stoi(str));
    }
    return result;
}


DensVec linkinv(DensVec const& eta, std::string const& family_t, const std::string& link) 
{
    DensVec result(eta.size());
    if (family_t == "gaussian") 
    {
        if (link == "identity") 
        {
            result = eta; // identity link function
        } else if (link == "log") {
            result = eta.array().exp().matrix(); // inverse of the log link function
        } else if (link == "sqrt") 
        {
            result = eta.array().pow(2).matrix(); // inverse of the square root link function
        } else 
        {
            std::cerr << "Unsupported link function for the Gaussian family." << std::endl;
        }
    } 
    else if (family_t == "binomial") 
    {
        if (link == "logit") 
        {
            result = ((eta.array()).exp() / (1.0f + ((eta.array()).exp()))).matrix(); // inverse of the logit link function
        } 
        else 
        {
            std::cerr << "Unsupported link function for the Binomial family." << std::endl;
        }
    } 
    else 
    {
        std::cerr << "Unsupported family type." << std::endl;
    }
    return result;
}


double calc_variance(DensVec const& dv)
{
    double variance = 0.0;
    int size = dv.size();
    if(size < 2)
    {
            fmt::print("Warning: Variance calculation requires at least two elements.\n");
            return variance;
    }

    double mean = dv.mean();

    for(int i{0}; i < size; ++i)
    {
        variance += (dv(i) - mean) * (dv(i) - mean);
    }

    variance /= (size - 1);
    return variance;
}


void Fit::calc_dmu_deta(std::string const& family_t, int size)
{
    dmu_deta.resize(size);
    if(family_t == "binomial")
    {
        for(int i{0}; i < size; ++i)
        {
            dmu_deta(i) = mu(i)*( 1 -mu(i)); 
        }
    }
    else
    {

        dmu_deta = DensVec::Constant(size, 1);
    }
}


DensVec Fit::calc_sqrtW()
{
    return dmu_deta.array().sqrt();
}


std::ext::V_double conv_dm2stdV(DensMat const& dm)
{
    int vec_idx = 0;
    std::ext::V_double v2ret(dm.size());

    for (int row=0; row < dm.rows(); ++row) 
    {
        // v2ret[vec_idx++] = 1.0;
        for (int col=0; col < dm.cols() ; ++col)
        {
            v2ret[vec_idx++] = dm(row,col);
        }
    }
    return v2ret;
}


std::ext::V_double conv_dv2stdVd(DensVec const& dv)
{
    std::ext::V_double v2ret(dv.size());
    for (int i=0; i < dv.size(); ++i) 
    {
            v2ret[i] = dv(i);
    }
    return v2ret;
}

DensVec conv_stdVs2dV(std::ext::V_string const& stdv) {
    DensVec dv2ret(stdv.size());

    for (size_t i = 0; i < stdv.size(); ++i) 
    {
        // Convert each string to a numerical value and assign it to the Eigen vector
        try
        {
            dv2ret(i) = std::stod(stdv[i]);
        }
        catch(std::invalid_argument const& e)
        {
            std::cerr << "value" << stdv[i] << " is not a valid number \n";  
        }
    }

    return dv2ret;
}


DensVec conv_stdV2dv(std::ext::V_double v)
{
    DensVec dv(v.size());
    for(size_t i {0}; i < v.size(); ++i)
    {
        dv(i) = v[i];
    }
    return dv;
}


Fit GEMFit::convert_2_fit()
{
    Fit fit;
    fit.mu = conv_stdV2dv(mu);
    fit.alpha = conv_stdV2dv(alpha);
    fit.eta = conv_stdV2dv(eta);

    return fit;
}


DensVec slice_vec(DensVec const& dv, std::ext::V_int const& indices)
{
    DensVec selected_elements(indices.size());
    for (int j = 0; j < indices.size(); ++j) 
    {
        selected_elements(j) = dv(indices[j]);
    }
    return selected_elements;
}


DensVec slice_vec(DensVec const& dv, std::ext::V_bool const& indices) 
{
    std::vector<int> int_indices;
    for (size_t i = 0; i < indices.size(); ++i) 
    {
        if (indices[i]) 
        {
            int_indices.push_back(i);
        }
    }
    return slice_vec(dv, int_indices);
}


DensVec slice_vec(DensVec const& dv, std::ext::Var_bool_int const& indices) 
{
    return std::visit([&dv](auto&& arg) -> DensVec {
        return slice_vec(dv, arg);
        }, indices);
}


SpaMat slice_mat(SpaMat const& spm, std::ext::V_int indices)
{
    // if(spm.rows() == indices.size()) return spm;
    SpaMat sliced_matrix(indices.size(), spm.cols());
    for(size_t i = 0; i < indices.size(); ++i)
    {
        for(size_t j = 0; j < spm.cols(); ++j)
            sliced_matrix.insert(i,j) = spm.coeff(indices[i],j);
    }
    return sliced_matrix;
}


//Slice dense mat based on rows indixes
DensMat slice_mat(DensMat const& dm, std::ext::V_int const& indices) 
{
    // if (dm.rows() == indices.size()) return dm;
    DensMat sliced_matrix(indices.size(), dm.cols());
    for (size_t i = 0; i < indices.size(); ++i) 
    {
        sliced_matrix.row(i) = dm.row(indices[i]);
    }
    return sliced_matrix;
}

// Function to slice matrix using vector of booleans
DensMat slice_mat(DensMat const& dm, std::ext::V_bool const& indices) 
{
    std::vector<int> int_indices;
    for (size_t i = 0; i < indices.size(); ++i) 
    {
        if (indices[i]) 
        {
            int_indices.push_back(i);
        }
    }
    return slice_mat(dm, int_indices);
}


// Overloaded function using std::variant
DensMat slice_mat(DensMat const& dm, std::ext::Var_bool_int const& indices) 
{
    return std::visit([&dm](auto&& arg) -> DensMat {
        return slice_mat(dm, arg);
    }, indices);
}

// Function to slice SpaMat (sparse matrix) based on row and column indices

SpaMat slice_mat(SpaMat const& spm, std::ext::V_int ind1, std::ext::V_int ind2, bool check_size)
{
    if (check_size && spm.rows() == ind1.size() && spm.cols() == ind2.size())
    {
        return spm;
    }
    
    SpaMat sliced_matrix(ind1.size(), ind2.size());
    
    size_t num_threads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    std::vector<std::ext::VecTuples4spmat> thread_triplets(num_threads);

    // Worker function for slicing
    auto worker = [&](size_t thread_id, size_t start_row, size_t end_row)
    {
        std::ext::VecTuples4spmat local_triplets;
        for (size_t i = start_row; i < end_row; ++i)
        {
            for (size_t j = 0; j < ind2.size(); ++j)
            {
                double value = spm.coeff(ind1[i], ind2[j]);
                if (value != 0.0)
                {
                    local_triplets.emplace_back(i, j, value);
                }
            }
        }
        thread_triplets[thread_id] = std::move(local_triplets);
    };

    // Distribute workload
    size_t chunk_size = ind1.size() / num_threads;
    size_t remainder = ind1.size() % num_threads;
    size_t start_row = 0;

    for (size_t t = 0; t < num_threads; ++t)
    {
        size_t end_row = start_row + chunk_size + (t < remainder ? 1 : 0);
        threads.emplace_back(worker, t, start_row, end_row);
        start_row = end_row;
    }

    // Join threads
    for (auto& thread : threads)
    {
        thread.join();
    }

    // Combine triplets from all threads
    std::ext::VecTuples4spmat all_triplets;
    for (const auto& triplet_vec : thread_triplets)
    {
        all_triplets.insert(all_triplets.end(), triplet_vec.begin(), triplet_vec.end());
    }

    // Build the sparse matrix efficiently
    sliced_matrix.setFromTriplets(all_triplets.begin(), all_triplets.end());

    return sliced_matrix;
}


// Slice mat based on rows and cols logical vector
SpaMat slice_mat(SpaMat const& spm, std::ext::V_bool const& ind1, std::ext::V_bool const& ind2, bool check_size) 
{
    std::ext::V_int int_indices1;
    for (size_t i = 0; i < ind1.size(); ++i) {
        if (ind1[i]) {
            int_indices1.push_back(i);
        }
    }
    std::ext::V_int int_indices2;
    for (size_t i = 0; i < ind2.size(); ++i) {
        if (ind2[i]) {
            int_indices2.push_back(i);
        }
    }
    return slice_mat(spm, int_indices1, int_indices2, check_size);
}

// Overloaded function to handle std::variant<std::vector<int>, std::vector<bool>> for both rows and columns
SpaMat slice_mat(SpaMat const& spm, std::ext::Var_bool_int const& ind1, std::ext::Var_bool_int const& ind2, bool check_size) 
{
    return std::visit([&spm, &ind2, check_size](auto&& arg1) -> SpaMat {
        return std::visit([&spm, &arg1, check_size](auto&& arg2) -> SpaMat {
            return slice_mat(spm, arg1, arg2, check_size);
        }, ind2);
    }, ind1);
}

DensMat slice_mat(DensMat const& dm, std::ext::V_int ind1, std::ext::V_int ind2)
{
    // if(dm.rows() == ind1.size() && dm.cols() == ind2.size()) return dm;
    DensMat sliced_matrix(ind1.size(), ind2.size());
    for(size_t i = 0; i < ind1.size(); ++i)
    {
        for(size_t j = 0; j < ind2.size(); ++j)
        {
            sliced_matrix(i,j) = dm(ind1[i],ind2[j]);
        }
    }
    return sliced_matrix;
}


DensMat slice_mat(DensMat const& dm, std::ext::V_bool const& ind1, std::ext::V_bool const& ind2) 
{
    std::ext::V_int int_indices1;
    for (size_t i = 0; i < ind1.size(); ++i) {
        if (ind1[i]) {
            int_indices1.push_back(i);
        }
    }
    std::ext::V_int int_indices2;
    for (size_t i = 0; i < ind2.size(); ++i) {
        if (ind2[i]) {
            int_indices2.push_back(i);
        }
    }
    return slice_mat(dm, int_indices1, int_indices2);
}

// Overloaded function to handle std::variant<std::vector<int>, std::vector<bool>> for both rows and columns
DensMat slice_mat(DensMat const& dm, std::ext::Var_bool_int const& ind1, std::ext::Var_bool_int const& ind2) {
    return std::visit([&dm, &ind2](auto&& arg1) -> DensMat {
        return std::visit([&dm, &arg1](auto&& arg2) -> DensMat {
            return slice_mat(dm, arg1, arg2);
        }, ind2);
    }, ind1);
}


DensMat slice_mat_cols(DensMat const& dm, std::ext::V_int const& ind2)
{
    if (ind2.empty()) 
    {
        return dm;
    }
    DensMat sliced_matrix(dm.rows(), ind2.size());
    for (size_t j = 0; j < ind2.size(); ++j) 
    {
        sliced_matrix.col(j) = dm.col(ind2[j]);
    }
    return sliced_matrix;
}


bool check_convergence(const DensVec& alpha, const DensVec& alpha0, 
                      const DensVec& tau, const DensVec& tau0, 
                      double tol, size_t& i, int maxiter) 
{
    double max_difference_alpha = ((alpha - alpha0).array().abs() / (alpha.array().abs() + alpha0.array().abs() + tol)).maxCoeff();
    double max_difference_tau   = ((tau - tau0).array().abs() / (tau.  array().abs() + tau0.  array().abs() + tol)).maxCoeff();
    double max_difference = std::max(max_difference_alpha, max_difference_tau);
    if ((2 * max_difference) < tol) 
    {
        return true; // Converged
    }

    if ((tau.array().abs().maxCoeff()) > pow(tol, -2)) {
        std::cerr << "Large variance estimate observed in the iterations, model not converged..." << std::endl;
        i = maxiter;
        return true; // Indicating non-convergence
    }

    return false; // Not converged yet
}




// Find ids for each group, group male=0 and female =1 --> m_group_idx[0]={1,3,5}
void GMMAT::extract_group_idx(std::unordered_set<int> const& group_unique, std::ext::V_int group_id)
{
    for(auto it = group_unique.begin(); it != group_unique.end(); ++it)
    {
        int value = *it;
        m_group_idx[value-1] = which(group_id, [value](int j){return (j==value);});
    }
}

//Check if any negative value exist
bool GMMAT::any_negative()
{
    for(auto i{0}; i < m_tau.size(); ++i)
    {
        if(m_tau(i) < 0) return true;
    }
    return false;
}

bool GMMAT::any_negative(std::ext::V_int vec)
{
    for(auto id : vec)
    {
        if(m_tau(id) < 0) return true;
    }
    return false;
}


//Check if any zero value exist
bool GMMAT::any_nonzero()
{
    return std::any_of(m_fixrho.begin(), m_fixrho.end(), [](int val){ return val != 0 ;});
}


void GMMAT::update_tau(DensVec const& tau0, DensVec const& dtau)
{
    for(size_t i{0}; i< m_idxtau.size(); ++i)
    {
        m_tau(m_idxtau[i]) = tau0(m_idxtau[i]) + dtau(i);
    }
}


void GMMAT::update_tau_below_tol(DensVec const& tau0, double tol)
{
    for (size_t i{0}; i < m_tau.size(); ++i) 
        {
            if(m_tau(i) < tol && tau0(i) < tol) 
            {
                m_tau(i) = 0.0;    
            }
        }
}

void GMMAT::update_tau_below_tol2(DensVec const& tau0, double tol)
{
    for (size_t i{0}; i < m_tau.size(); ++i) 
        {
            if(m_tau(i) < tol) 
            {
                m_tau(i) = 0.0;    
            }
        }
}


void GMMAT::set_ai_low_ng(int i, DensVec& score, DensMat& ai, DensVec const& wpy, Fit const& fit, DensVec const& py, DensVec diagp, DensMat sigma_ixcov)
{
    DensVec results = slice_vec(wpy, m_group_idx[m_idxtau[i]]).cwiseProduct(slice_vec(py, m_group_idx[m_idxtau[i]])) - slice_vec(diagp, m_group_idx[m_idxtau[i]]);
    score(i) = results.sum();
    for(size_t j{0}; j <= i; ++j)
    {
        SpaMat crosspart_fir_a = slice_mat(fit.sigma_i, m_group_idx[m_idxtau[j]], m_group_idx[m_idxtau[i]]);
        DensVec cross_result = crossprod(crosspart_fir_a, slice_vec(wpy, m_group_idx[m_idxtau[j]]));
        double first = (slice_vec(wpy, m_group_idx[m_idxtau[i]]).cwiseProduct(cross_result)).sum();
        DensMat sliced_mat_sigma_ixcov = slice_mat(sigma_ixcov, m_group_idx[m_idxtau[i]]);
        DensVec crosspart_sec_a = crossprod(sliced_mat_sigma_ixcov, slice_vec(wpy, m_group_idx[m_idxtau[i]]));
        DensMat sliced_mat_sigma_ix = slice_mat(fit.sigma_ix, m_group_idx[m_idxtau[j]]);
        DensVec crosspart_sec_b = crossprod(sliced_mat_sigma_ix, slice_vec(wpy, m_group_idx[m_idxtau[j]]));
        double second = (crosspart_sec_a.cwiseProduct(crosspart_sec_b)).sum();
        ai(i,j) = first - second;
        if(j != i) ai(j,i) = ai(i,j);
    }  
}


void GMMAT::set_ai_high_ng(int i, DensVec& score, DensMat& ai, DensVec const& wpy, Fit const& fit, DensVec const& py, DensMat sigma_ixcov, int ng)
{
    SpaMat spmat = m_vkins_sp[m_idxtau[i]-ng].get_spmat();
    DensVec apy = crossprod(spmat, py);
    DensVec first_part = crossprod(fit.sigma_i, apy);   
    DensVec second_part = tcrossprod(fit.sigma_ix, crossprod(sigma_ixcov, apy).transpose());
    DensVec papy = first_part - second_part;
    double first_score = (fit.sigma_i.cwiseProduct(spmat)).sum(); //Sama added .array // let write it without array, sparse matrix does not have array()
    double second_score = (fit.sigma_ix.cwiseProduct(crossprod(spmat, sigma_ixcov))).sum();//Sama added .array
    score(i) = (m_Y.cwiseProduct(papy)).sum() - (first_score - second_score);
    for(size_t j{0}; j <= i; ++j)
    {
        if(m_idxtau[j] < ng)
        {
            ai(i,j) = (slice_vec(wpy, m_group_idx[m_idxtau[j]]).cwiseProduct(slice_vec(papy, m_group_idx[m_idxtau[j]]))).sum();
        }
        else
        {     
            SpaMat spmat_j = m_vkins_sp[m_idxtau[j]-ng].get_spmat();
            ai(i,j) = (py.cwiseProduct(crossprod(spmat_j, papy))).sum();
        }
        if(j != i)
            ai(j,i) = ai(i,j);
    }      
}


void GMMAT::set_ai(DensVec& score, DensMat& ai, DensVec const& wpy, Fit const& fit, DensVec const& py, DensVec diagp, DensMat sigma_ixcov, int ng, int fixtau_0_counts)
{
    for(size_t i{0}; i < fixtau_0_counts; ++i)
    {
        if(m_idxtau[i] < ng)
        {
            set_ai_low_ng(i, score, ai, wpy, fit, py, diagp, sigma_ixcov);
            
        }
        else
        {
            set_ai_high_ng(i, score, ai, wpy, fit, py, sigma_ixcov, ng); 
               
        }       
    }
}


//Calculate the random effect and fill m_covariance_idx
void GMMAT::calc_rand_effect(int &kins_size, int ng)
{
    m_vkins_sp.resize(3*kins_size);
    for(int i{0}; i < kins_size; ++i)
    {
        SparseInverse spi_cov;
        SparseInverse spi_slope;
        SpaMat sp = m_vkins_sp[i].get_spmat();
        SpaMat sp_left = sp * m_rand_slope.asDiagonal();
        SpaMat sp_right = sp_left.transpose();
        SpaMat sp_cov = sp_left + sp_right;
        sp_cov.makeCompressed();
        spi_cov.set_spmat(sp_cov);
        m_vkins_sp[kins_size + i] = spi_cov;
        SpaMat sp_slope = sp_right * m_rand_slope.asDiagonal();
        sp_slope.makeCompressed();
        spi_slope.set_spmat(sp_slope);
        m_vkins_sp[kins_size * 2 + i] = spi_slope;
        m_covariance_idx.push_back({kins_size + i + ng, i + ng, kins_size * 2 + i + ng});
    }
}

//Check constraint of random slope
void GMMAT::fill_fixrho_idx(double tol)
{
    m_fixrho_idx.clear();
    for(int i{0}; i < m_covariance_idx.size(); ++i)
    {
        if(abs(m_tau(m_covariance_idx[i][0])) > (1 - 1.01 * tol) * sqrt(m_tau(m_covariance_idx[i][1]) * m_tau(m_covariance_idx[i][2])))
        {
            m_fixrho_idx.push_back(i);
        }
    }
}


bool GMMAT::covariate_larger_slope_intercept(double tol)
{
    for(int i{0}; i < m_covariance_idx.size(); ++i)
    {        
        if(abs(m_tau(m_covariance_idx[i][0])) > sqrt(m_tau(m_covariance_idx[i][1]) * m_tau(m_covariance_idx[i][2])))
        {
            return true;
        }
    }
    return false;
}

void GMMAT::update_fixtau_fixrho(std::ext::V_int &fixtau_new,std::ext::V_int &fixrho_new, double tol)
{
    if(m_covariance_idx.size() > 0)
    {
        for(int i{0}; i < m_covariance_idx.size(); ++i)
        {
            fixtau_new[m_covariance_idx[i][0]] = 0;
        }
        fixrho_new.clear();
        fixrho_new.resize(m_covariance_idx.size(), 0);
        fill_fixrho_idx(tol);

        for(auto id : m_fixrho_idx)
        {
            fixrho_new[id] = sign(m_tau(m_covariance_idx[id][0]));
        }
    }
}

void GMMAT::fill_fixrho_idx0(std::ext::V_int &fixrho_idx0, DensVec const& tau0, double tol)
{
   for(int i{0}; i < m_covariance_idx.size(); ++i)
    {
        if(abs(tau0(m_covariance_idx[i][0])) > (1 - 1.01 * tol) * sqrt(tau0(m_covariance_idx[i][1]) * tau0(m_covariance_idx[i][2])))
        {
            fixrho_idx0.push_back(i);
        }
    } 
}

void GMMAT::update_mtau_with_m_covariance_idx_tau0(DensVec const& tau0, double tol)
{
    std::ext::V_bool exclude_idx_vec(m_tau.size(), false);

    for(auto elm : m_covariance_idx)
    {
        if(elm[0] < m_tau.size())
        {
            exclude_idx_vec[elm[0]] = true;
        }
    }
    
    for(int i{0}; i < m_tau.size(); ++i)
    {
        if(!exclude_idx_vec[i] && m_tau(i) < tol && tau0(i) < tol)
        {
            m_tau(i) = 0;
        }
    }
    
}

void GMMAT::update_mtau_with_m_covariance_idx(double tol)
{
    std::ext::V_bool exclude_idx_vec(m_tau.size(), false);
    
    for(auto elm : m_covariance_idx)
    {
        if(elm[0] < m_tau.size())
        {
            exclude_idx_vec[elm[0]] = true;
        }
    }
    
    for(int i{0}; i < m_tau.size(); ++i)
    {
        if(!exclude_idx_vec[i] && m_tau(i) < tol)
        {
            m_tau(i) = 0;
        }
    }
}

void GMMAT::update_mtau_with_m_covariance_idx_mfixrho(std::ext::V_int const& idxrho)
{
    
    for(auto idx : idxrho)
    {
        m_tau(m_covariance_idx[idx][0]) = m_fixrho[idx] * sqrt(m_tau(m_covariance_idx[idx][1]) * m_tau(m_covariance_idx[idx][2]));
    }
    
}


void GMMAT::update_mtau_with_m_covariance_mfixrho_idx_fixrho_idx0(std::ext::V_int const& fixrho_idx0)
{
    std::ext::V_int intersect_idx;
    intersect_idx = intersect(m_fixrho_idx, fixrho_idx0);
    for(auto id : intersect_idx)
    {
        m_tau(m_covariance_idx[id][0]) = sign(m_tau(m_covariance_idx[id][0])) * sqrt(m_tau(m_covariance_idx[id][1]) * m_tau(m_covariance_idx[id][2]));
    }
}

void GMMAT::update_mtau_with_m_covariance_mfixrho_idx()
{
    for(auto id : m_fixrho_idx)
    {
        m_tau(m_covariance_idx[id][0]) = sign(m_tau(m_covariance_idx[id][0])) * sqrt(m_tau(m_covariance_idx[id][1]) * m_tau(m_covariance_idx[id][2]));
    }
}

std::ext::V_int GMMAT::exclude_idx()
{
    std::ext::V_bool exclude_idx_bool(m_tau.size(), false);
    std::ext::V_int exclude_idx_vec;

    for(auto elm : m_covariance_idx)
    {
        if(elm[0] < m_tau.size())
        {
            exclude_idx_bool[elm[0]] = true;
        }
    }

    for(int i{0}; i < m_tau.size(); ++i)
    {
        if(!exclude_idx_bool[i])
        {
            exclude_idx_vec.push_back(i);
        }
    }
    return exclude_idx_vec;
}

Fit GMMAT::fitglmm_ai(DensVec const& W)
{
    Fit fit_to_return;
    auto n_rows = m_X.rows();
    auto fixtau_0_counts = std::count(m_fixtau.begin(), m_fixtau.end(), 0); // it is better to find the better name
    auto ng = m_group_idx.size(); 
    std::ext::V_double diag_sigma(n_rows,0); // this is vector of diagonal element of Sigma
    DensVec score; 
    DensVec wpy;
    DensVec py;
    DensVec diagp;
    DensVec apy;
    DensVec papy;

    fit_to_return.W  = W;
    for(int i{0}; i < ng; ++i)
    {
        for (int j : m_group_idx[i]) 
        {
            diag_sigma[j] = m_tau(i) / fit_to_return.W(j);
        }
    }

    // Create diag mat using vector
    SpaMat sigma = details::diag(diag_sigma);  
    int kins_size = m_vkins_sp.size();

    for(int i{0}; i < kins_size; ++i)
    {
        SpaMat curr_kin_spmat = m_vkins_sp[i].get_spmat();
        //if only the upper part of sm1 is stored, convert vec to eigen vec
        sigma = sigma + m_tau(i + ng) * curr_kin_spmat;//for(i in 1:q) Sigma <- Sigma + tau[i+ng]*kins[[i]]
    }
  
    fit_to_return.sigma_i =  SparseInverse::inv_spamat(sigma);
    fit_to_return.sigma_ix = crossprod(fit_to_return.sigma_i, m_X);
    DensMat xsigma_ix = crossprod(m_X, fit_to_return.sigma_ix); 
    fit_to_return.cov = SparseInverse::inv(xsigma_ix);
    DensMat sigma_ixcov = tcrossprod(fit_to_return.sigma_ix, fit_to_return.cov);
    fit_to_return.alpha = crossprod(fit_to_return.cov, crossprod(fit_to_return.sigma_ix, m_Y)); // check in the compile time
    fit_to_return.eta = m_Y - conv_vec_vecXd(diag_sigma).cwiseProduct((crossprod(fit_to_return.sigma_i, m_Y) - tcrossprod(fit_to_return.sigma_ix, fit_to_return.alpha.transpose())));

    if(fixtau_0_counts > 0)
    {
        m_idxtau = which(m_fixtau, [](int i){return (i==0);});
        py = crossprod(fit_to_return.sigma_i, m_Y) - tcrossprod(fit_to_return.sigma_ix, crossprod(sigma_ixcov, m_Y).transpose());
        wpy = py.cwiseQuotient(fit_to_return.W);    
        //int col_n = crossprod(fit_to_return.sigma_ix, sigma_ixcov).cols();
        DensVec diagp = (fit_to_return.sigma_i.diagonal() - ((fit_to_return.sigma_ix.cwiseProduct(sigma_ixcov)).rowwise().sum())).cwiseQuotient(fit_to_return.W);
        DensMat AI = DensMat::Constant(fixtau_0_counts, fixtau_0_counts, std::numeric_limits<double>::quiet_NaN());
        DensVec score = DensVec::Constant(fixtau_0_counts, std::numeric_limits<double>::quiet_NaN());
        set_ai(score, AI, wpy,fit_to_return, py, diagp, sigma_ixcov,ng, fixtau_0_counts);

        // DensVec dv =  AI.fullPivLu().solve(score);
        // fit_to_return.dtau = std::make_optional(dv);

        auto lu_decomp = AI.fullPivLu();
        if (lu_decomp.isInvertible()) 
        {
            DensVec dv = lu_decomp.solve(score);
            fit_to_return.dtau = std::make_optional(dv);
        } 
        else 
        {
            std::cout << "The matrix is not invertible, solve operation failed." << std::endl;
            fit_to_return.dtau = std::nullopt;
            exit(EXIT_FAILURE); 
        }
        return fit_to_return;
    } 
    return fit_to_return;
}

Glmmkin GMMAT::glmmkin_ai(Fit fit_null, int maxiter, double tol)
{
    Glmmkin glmmkin;
    DensVec py;
    DensVec apy;
    DensVec papy; 
    DensVec tau0;
    Fit fit_glmm_ai;
    DensVec alpha0;
    glmmkin.fit.alpha = fit_null.alpha;
    int y_size = m_y.size();
    // std::cout << "Fixed-effect coefficient (alpha):\n" << glmmkin.fit.alpha << '\n';

    if(m_offset.size() < y_size) 
    {
        m_offset = DensVec::Constant(y_size, 0); //we cannot check null in cpp
    }

    m_tau = DensVec::Constant(m_vkins_sp.size() + m_group_idx.size(), 0);
    m_fixtau.resize(m_vkins_sp.size() + m_group_idx.size(), 0);  
    fit_null.calc_dmu_deta(m_family_t, y_size);
    glmmkin.fit = fit_null;
    
    m_Y = fit_null.eta - m_offset + (m_y - glmmkin.fit.mu).cwiseQuotient(glmmkin.fit.dmu_deta); 
    m_sqrtW = glmmkin.fit.calc_sqrtW();
    glmmkin.fit.W = glmmkin.fit.dmu_deta;

    if(m_family_t == "binomial")
    {
        m_tau(0) = 1.0;
        m_fixtau[0] = 1;
    }

    int kins_size = m_vkins_sp.size();
    auto ng = m_group_idx.size(); 
    m_idxtau = which(m_fixtau, [](int i){return (i==0);});
    
    if(m_covariance_idx.size() > 0)
    {
        std::ext::V_int covariance_vec;
        for(auto elm : m_covariance_idx)
        {
            covariance_vec.push_back(elm[0]);
        }
        m_idxtau2 = intersect(covariance_vec, m_idxtau);
    }

    auto fixtau_0_counts = std::count(m_fixtau.begin(), m_fixtau.end(), 0); 

    if(fixtau_0_counts > 0)
    {
        auto tau_value = calc_variance(m_Y) / (kins_size + ng);

        for(int idx : m_idxtau)
        {
             m_tau(idx) = tau_value;//m_fixtau(0) is 1 for binomia so idx!=0 m_tau(0) will be 1 always
        }

        if(m_covariance_idx.size() > 0)
        {
            for(auto id : m_idxtau2)
            {
                m_tau(id) = 0;
            }
        }

        std::ext::V_double diag_sigma(y_size, 0); // vector of diagonal element of Sigma--

        for(int i{0}; i < ng; ++i)
        {
            for (int j : m_group_idx[i]) 
            {
                diag_sigma[j] = m_tau(i) / glmmkin.fit.W(j);
            }
        } 
        
        SpaMat sigma = details::diag(diag_sigma);

        for(int i{0}; i < kins_size; ++i)
        {
            SpaMat curr_kin_spmat = m_vkins_sp[i].get_spmat();
            m_tau(i + ng) = m_tau(i + ng) / (curr_kin_spmat).diagonal().mean();
            sigma = sigma + m_tau(i + ng) * curr_kin_spmat;
        }
        
        glmmkin.fit.sigma_i =  SparseInverse::inv_spamat(sigma);
        glmmkin.fit.sigma_ix = crossprod(glmmkin.fit.sigma_i, m_X);
        DensMat xsigma_ix = crossprod(m_X, glmmkin.fit.sigma_ix);
        glmmkin.fit.cov = SparseInverse::inv(xsigma_ix);
        DensMat sigma_ixcov = tcrossprod(glmmkin.fit.sigma_ix, glmmkin.fit.cov);
        py = crossprod(glmmkin.fit.sigma_i, m_Y) - tcrossprod(glmmkin.fit.sigma_ix, crossprod(sigma_ixcov, m_Y).transpose());
        tau0 = m_tau; 
        
        for(std::size_t i{0}; i < fixtau_0_counts; ++i)
        {
            if(m_idxtau[i] < ng)
            {
                double first_part = tau0(m_idxtau[i]) * tau0(m_idxtau[i]);
                DensVec div = py.cwiseQuotient(m_sqrtW); //cwiseQuotient does elementwise /
                
                double second_part = slice_vec(div, m_group_idx[m_idxtau[i]]).cwiseProduct(slice_vec(div, m_group_idx[m_idxtau[i]])).sum();
                int col_n = (glmmkin.fit.sigma_ix.cwiseProduct(sigma_ixcov)).cols();
                DensVec third_part = glmmkin.fit.sigma_i.diagonal() - (glmmkin.fit.sigma_ix.cwiseProduct(sigma_ixcov)) * DensVec::Ones(col_n);
                DensVec fourth_part = third_part.cwiseQuotient(m_sqrtW);
                double fifth_part = slice_vec(fourth_part, m_group_idx[m_idxtau[i]]).sum();
                double sixth_part = second_part - fifth_part;
                m_tau(m_idxtau[i]) = std::max(0.0, tau0(m_idxtau[i]) + first_part * sixth_part / y_size);
            }
            else
            {   
                apy = crossprod(m_vkins_sp[m_idxtau[i] - ng].get_spmat(), py);
                papy =  crossprod(glmmkin.fit.sigma_i, apy) - tcrossprod(glmmkin.fit.sigma_ix, (crossprod(sigma_ixcov, apy)).transpose());
                double first_part = tau0(m_idxtau[i]) * tau0(m_idxtau[i]); 
                double second_part = (m_Y.cwiseProduct(papy)).sum ();
                double third_part = (glmmkin.fit.sigma_i.cwiseProduct(m_vkins_sp[m_idxtau[i] - ng].get_spmat())).sum();            
                double fourth_part = (glmmkin.fit.sigma_ix.cwiseProduct(crossprod(m_vkins_sp[m_idxtau[i] - ng].get_spmat(), sigma_ixcov))).sum();
                double fifth_part = third_part - fourth_part;            
                double sixth_part =  second_part - fifth_part;
                double seventh_part = tau0(m_idxtau[i]) + (first_part * sixth_part / y_size);            

                if(m_covariance_idx.size() > 0)
                {
                    // std::ext::V_int intersect_vec;
                    // std::ext::V_int not_intersect_vec;
                    // intersect_vec = intersect(m_idxtau, m_idxtau2);
                    // not_intersect_vec = not_intersect(m_idxtau, m_idxtau2);
                    
                    if(std::find(m_idxtau2.begin(), m_idxtau2.end(), m_idxtau[i]) != m_idxtau2.end())
                    {
                        m_tau(m_idxtau[i]) = 0.0;
                    }
                    else
                    {
                        m_tau(m_idxtau[i]) = std::max(0.0, seventh_part);
                    }
                }
                else
                {
                    m_tau(m_idxtau[i]) = std::max(0.0, seventh_part);
                }   
            }
        }  
    }

    size_t i;
    
    for(i = 1; i < maxiter; ++i)
    {
        std::cout << "iteration: " << i << '\n';
        alpha0 = glmmkin.fit.alpha;
        tau0 = m_tau;
        fit_glmm_ai = fitglmm_ai(glmmkin.fit.W); 

        if(fixtau_0_counts > 0)
        {
            if (fit_glmm_ai.dtau.has_value())
            {
                glmmkin.fit.dtau = fit_glmm_ai.dtau;
                update_tau(tau0, glmmkin.fit.dtau.value());

                if(m_covariance_idx.empty())
                {
                    update_tau_below_tol(tau0, tol);
                    while(any_negative())
                    {
                        glmmkin.fit.dtau = glmmkin.fit.dtau.value() / 2;
                        update_tau(tau0, glmmkin.fit.dtau.value());
                        update_tau_below_tol(tau0, tol);
                    }
                    update_tau_below_tol2(tau0, tol);
                }
                else
                {
                    std::ext::V_int fixrho_idx0;
                    std::ext::V_int idxrho;
                    fill_fixrho_idx0(fixrho_idx0, tau0, tol);
                    update_mtau_with_m_covariance_idx_tau0(tau0, tol);

                    if(any_nonzero())
                    {
                        idxrho = which(m_fixrho, [](int i){return i != 0;});
                        update_mtau_with_m_covariance_idx_mfixrho(idxrho);
                    }
                    
                    fill_fixrho_idx(tol);
                    update_mtau_with_m_covariance_mfixrho_idx_fixrho_idx0(fixrho_idx0);
                    std::ext::V_int exclude_idx_vec;
                    exclude_idx_vec = exclude_idx();
                    
                    while(any_negative(exclude_idx_vec) || covariate_larger_slope_intercept(tol))
                    {
                        glmmkin.fit.dtau = glmmkin.fit.dtau.value() / 2;
                        update_tau(tau0, glmmkin.fit.dtau.value());
                        update_mtau_with_m_covariance_idx_tau0(tau0, tol);

                        if(any_nonzero())
                        {
                            idxrho = which(m_fixrho, [](int i){return i != 0;});
                            update_mtau_with_m_covariance_idx_mfixrho(idxrho);
                        }

                        fill_fixrho_idx(tol);
                        update_mtau_with_m_covariance_mfixrho_idx_fixrho_idx0(fixrho_idx0);
                    }

                    update_mtau_with_m_covariance_idx(tol);
                    update_mtau_with_m_covariance_mfixrho_idx();
                }
            }
            else
            {
                std::cerr << "Error: fit_glmm_ai.dtau has no value!" << std::endl;
                // Handle error, possibly return or throw an exception
            }
        } 

        glmmkin.fit.cov = fit_glmm_ai.cov;
        glmmkin.fit.alpha = fit_glmm_ai.alpha;
        glmmkin.fit.eta = fit_glmm_ai.eta + m_offset;
        glmmkin.fit.mu = linkinv(fit_glmm_ai.eta, m_family_t, m_link);
        glmmkin.fit.calc_dmu_deta(m_family_t, y_size);
        m_Y = glmmkin.fit.eta - m_offset + (m_y - glmmkin.fit.mu).cwiseQuotient(glmmkin.fit.dmu_deta);
        m_sqrtW = glmmkin.fit.calc_sqrtW();
        glmmkin.fit.W = glmmkin.fit.dmu_deta;
		glmmkin.fit.sigma_ix = fit_glmm_ai.sigma_ix;
		glmmkin.fit.sigma_i = fit_glmm_ai.sigma_i;
        std::cout << "Variance component estimates (m_tau):\n" << m_tau << '\n';
        std::cout << "Fixed-effect coefficient (alpha):\n" << glmmkin.fit.alpha << '\n';
        if(check_convergence(glmmkin.fit.alpha, alpha0, m_tau, tau0, tol, i, maxiter)) break;
    } 
    
    glmmkin.converged = i < maxiter ? true : false;
    glmmkin.residuals = m_y.cast<double>() - glmmkin.fit.mu;
    
    DensVec res_var = DensVec::Constant(y_size, 1);
    for (int i=0; i< ng; i++)
    {
        auto it = m_group_idx.find(i);
        if (it != m_group_idx.end())
        {
            std::ext::V_int g_idx = it->second;
            for (auto idx : g_idx)
            {
                res_var(idx) = m_tau(i);
            }
        }
    }

    DensVec fit0W = DensVec::Constant(glmmkin.fit.W.size(), 1);
    //fill scaled_residuals
    glmmkin.scaled_residuals =  glmmkin.residuals.array() * fit0W.array() / res_var.array();
    return glmmkin;
}  


Glmmkin GMMAT::glmmkin_fit(Fit fit_null, std::ext::V_int group_id, 
                        std::string const method, 
                        std::string method_optim, 
                        int maxiter,
                        double tol, double tau_min, 
                        double tau_max, int tau_region)
{
    Glmmkin glmmkin;
    if(method_optim == "Brent")
    {
        fmt::println("Error: we do not support Brent");
        exit(EXIT_FAILURE);
    }
    
    if(method_optim == "AI")
    {
        // set of groups_id--> male and female so group_unique size=2{0,1}
        auto group_unique = unique(group_id);
        extract_group_idx(group_unique, group_id);//Update m_group_idx
        int kins_size = m_vkins_sp.size();
        int ng = m_group_idx.size();
       
        std::ext::V_int fixrho_new;
        std::ext::V_int fixrho_old;

        if(m_rand_slope.size() > 0)
        {
            //Calculate random effect and fill 
            calc_rand_effect(kins_size, ng);
            fixrho_old.resize(kins_size, 0);
            kins_size = m_vkins_sp.size();
        }

        std::ext::V_int fixtau_old(kins_size + ng, 0);
        glmmkin = glmmkin_ai(fit_null, maxiter, tol);
        auto fixtau_new = logic_update_fixed_condtion(m_tau, tol);
        //Update fixtau and fixrho
        update_fixtau_fixrho(fixtau_new, fixrho_new, tol);

        while(fixtau_new != fixtau_old || (fixrho_new.size() > 0 && fixrho_new != fixrho_old))
        {
            fmt::print(stderr, "Warning: Variance estimate on the boundary of the parameter space observed, refitting model...\n");
            fixtau_old = fixtau_new;

            if(m_covariance_idx.size() > 0)
            {
                fixrho_old = fixrho_new;
            }
            m_fixtau = fixtau_old;
            m_fixrho = fixrho_old;
            glmmkin = glmmkin_ai(fit_null, maxiter, tol);
            fixtau_new = logic_update_fixed_condtion(m_tau, tol);
            update_fixtau_fixrho(fixtau_new, fixrho_new, tol);
        }

        if(!glmmkin.converged)
        {
            if(ng != 1)
            {
                fmt::print(stderr, "Error: Average Information REML not converged, cannot refit heteroscedastic linear mixed model using Brent or Nelder-Mead methods.\n");
                exit(EXIT_FAILURE);
            }

            if(m_rand_slope.size() > 0)
            {
                fmt::print(stderr, "Error: Average Information REML not converged, cannot refit random slope model for longitudinal data using Brent or Nelder-Mead methods.\n");
                exit(EXIT_FAILURE);
            }

            if(kins_size == 1)
            {
                fmt::print(stderr, "Average Information REML not converged, refitting model using Brent method...\n");
                fmt::print(stderr, "Brent is not available for the time being, stay in touch for updates ;)\n");
                exit(EXIT_FAILURE);
            }
        }
    }
    else
    {
        fmt::print(stderr, "The optimization method is not supported for the time being\n");
        fmt::print(stderr, "Stay in touch for any updates\n");
        exit(EXIT_FAILURE);
    }   
    return glmmkin;
}


Glmmkin GMMAT::glmmkin_final(std::ext::FitNull_f fit0, Pheno pheno,
                            std::ext::V_string cov_selected_hdrs,
                            std::string phenoname,
                            std::string const& id, 
                            std::string randomSlopeName,
                            std::string const& groups,
                            std::string const method, 
                            std::string method_optim, 
                            int maxiter,
                            double tol, double tau_min, 
                            double tau_max, int tau_region)
{
    Glmmkin glmmkin;
    pheno.m_sam_id = id;
    m_y = conv_stdVs2dV(pheno.m_data_frame.get_header(phenoname));
    int y_size = m_y.size();
    std::ext::V_string v_valid_methods {"REML", "ML"};
    std::string rand_slope_hdr = randomSlopeName;

    auto it = std::find(v_valid_methods.begin(), v_valid_methods.end(), method);
    if(it == v_valid_methods.end())
    {
        fmt::print(stderr, "Error: {} is not in GMMAT valid methods (REML, ML)\n", method);
        exit(EXIT_FAILURE);
    }

    if(method ==  "ML" && method_optim == "AI")
    {
        fmt::print(stderr, "Error: {} is not available for {}\n", method, method_optim);
        exit(EXIT_FAILURE);
    }

    if(rand_slope_hdr.size() > 0)
    {
        if(method_optim != "AI")
        {
            fmt::print(stderr, "Error: random slope for longitudinal data is currently only implemented for method.optim \"AI\".");
            exit(EXIT_FAILURE);
        }
        std::ext::V_string slope_temp = pheno.m_data_frame.get_header(rand_slope_hdr);
        m_rand_slope = conv_stdVs2dV(slope_temp);
    }

    Fit fit_null; 
    auto pair = pheno.check_binary(phenoname);
    m_family_t = pair.first;
    m_link = pair.second;
    auto pheno_type = (m_family_t == "binomial") ? 1 : 0;
    // m_X = create_covdata(pheno.m_data_frame.copy_by_hdrs(cov_selected_hdrs));
    GEMFit gf;
    std::ext::V_double pheno_data = conv_dv2stdVd(m_y);
    m_X = create_covdata(pheno.m_data_frame.copy_by_hdrs(cov_selected_hdrs));
    std::ext::V_double cov_data = conv_dm2stdV(m_X); 
    m_n_sel_col = cov_selected_hdrs.size();
    
    fit0(y_size, m_n_sel_col, pheno_type, tol, m_robust, cov_selected_hdrs, pheno_data, cov_data,
                 &gf.XinvXTX, &gf.mu, &gf.resid, &gf.sigma2, gf.alpha, gf.eta); 

    std::cout << std::flush;
    std::cout << "****************************************************************************\n";
    std::cout << "Start association test...\n \n";
    
    fit_null = gf.convert_2_fit(); 
    
    if(pheno.m_data_frame.any_duplicated(pheno.m_sam_id))
    {
        std::cout << "Duplicated id detected...\nAssuming longitudinal data with repeated measures...\n";
        if(!m_vkins_sp[0].kin.m_null_kin) // if there is a kinship file add another matrix
        {
            SparseInverse spi;
            spi.set_spmat(details::diag(std::ext::V_double(y_size, 1)));//Create a diagonal matrix
            m_vkins_sp.emplace_back(spi);
        }

        auto duplicates = pheno.m_data_frame.list_duplicates(pheno.m_sam_id);
        std::ext::IndexMap mapped_indices = m_vkins_sp[0].get_idx_mp();
        SpaMat spi_mat;
        if(!m_vkins_sp[0].kin.m_null_kin)
        {
            spi_mat = m_vkins_sp[1].get_spmat();
        }
        else
        {
            spi_mat = m_vkins_sp[0].get_spmat();
        }
        
        std::ext::VecTuples4spmat triplets;
        for (int k = 0; k < spi_mat.outerSize(); ++k) 
        {
            for (SpaMat::InnerIterator it(spi_mat, k); it; ++it) 
            {
                triplets.emplace_back(it.row(), it.col(), it.value());
            }
        }
        
        for(const auto& dup : duplicates)
        {
            auto v_indices = mapped_indices[dup]; //map id to index
            for(auto& elm1 : v_indices) //loop over duplicated ids to create values-->(id1,id1,1)
            {
                for(auto& elm2 : v_indices)
                {
                    if(spi_mat.coeff(elm1 - 1 , elm2 - 1) == 0)
                    {
                        triplets.emplace_back(elm1 - 1 , elm2 - 1, 1.0);
                    }
                }
            } 
        }

        spi_mat.setFromTriplets(triplets.begin(), triplets.end());
        spi_mat.makeCompressed();
        if(!m_vkins_sp[0].kin.m_null_kin)
        {
            m_vkins_sp[1].set_spmat(spi_mat);
        }
        else
        {
            m_vkins_sp[0].set_spmat(spi_mat);
        }
    }
    else if(m_vkins_sp[0].kin.m_null_kin && rand_slope_hdr.size() > 0)
    {
        fmt::print(stderr, "\"random slope\" ignored for cross-sectional data from unrelated individuals...");
        exit(EXIT_FAILURE);
    }

    std::ext::V_int group_id;
    if(groups.size() == 0)
    {
        group_id = std::ext::V_int (y_size, 1);
    }
    else
    {
        // Convert to vector of int as groups is a vector of string
        group_id = conv_stdvs2stdvi(pheno.m_data_frame.get_header(groups));
        
    }
    glmmkin = glmmkin_fit(fit_null, group_id, method, method_optim, 
                          maxiter, tol, tau_min, tau_max, tau_region);
    
    //To be passed to MAGEE
    glmmkin.id_include = pheno.m_data_frame.get_header(id);
	//To be passed to MAGEE
	glmmkin.sigma2 = gf.sigma2;
	// m_hdrsMap required by MAGEE
    m_hdrsMap = createHeaderMap(cov_selected_hdrs);
    return glmmkin;    
}
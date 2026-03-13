#include "../include/GMMAT.h"
#include <algorithm>
#include <limits>
#include <cmath>
#include <thread>
#include <mutex>

void center_dataframe(
    DataFrame& df,
    const std::vector<std::string>& headers,
    bool center,
    bool scale
)
{
    for (const auto& hdr : headers)
    {
        std::vector<std::string>& col_str = df.m_data[hdr];

        // Convert to double temporarily
        std::vector<double> col(col_str.size());
        for (size_t i = 0; i < col.size(); ++i)
            col[i] = std::stod(col_str[i]);

        double mean = 0.0;
        double sd = 1.0;

        if (center || scale)
        {
            for (double v : col)
                mean += v;
            mean /= static_cast<double>(col.size());
        }

        if (scale)
        {
            double var = 0.0;
            for (double v : col)
                var += (v - mean) * (v - mean);

            sd = std::sqrt(var / (col.size() - 1.0));
            if (sd == 0.0)
                sd = 1.0;
        }

        // Apply transform
        for (size_t i = 0; i < col.size(); ++i)
        {
            if (center)
                col[i] -= mean;

            if (scale)
                col[i] /= sd;

            col_str[i] = std::to_string(col[i]);
        }
    }
}


DensMat apply_wb(Eigen::Ref<const Eigen::MatrixXd> const& mat, DensVec const& diag_sigma_i, SpaMat const& diag_sigma_i_ZPchol)
{
    DensMat sigma_imat = (mat.array().colwise() * diag_sigma_i.array()).matrix() - (diag_sigma_i_ZPchol * crossprod(diag_sigma_i_ZPchol,mat));
    return sigma_imat;
}

void GMMAT::fill_mat(const int numRows, std::ext::V_int& indixes_col, int value)
{
    std::ext::VecTriple_i tripletList;

    for (int i = 0; i < numRows; ++i) 
    {
        if (indixes_col[i] != -1)
        {
            tripletList.emplace_back(i, indixes_col[i], value); 
        }
    } 

    m_J.setFromTriplets(tripletList.begin(), tripletList.end());   
}


void GMMAT::fill_J(std::ext:: V_string const& sample_ids)
{
    std::ext::V_string unique_sample_ids = unique_id(sample_ids);
    std::ext::V_int indixes_col = match_indices(sample_ids, unique_sample_ids); 
    m_J.resize(sample_ids.size(), unique_sample_ids.size());
    fill_mat(sample_ids.size(), indixes_col, 1);
}

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

//Returns a diagonal matrix
namespace details{
    SpaMat diag(DensVec vec)
    {
        SpaMat spm(vec.size(), vec.size());
        for(size_t i{0}; i < vec.size(); ++i)
        {
            spm.insert(i, i) = vec[i];
        }
        return spm;
    }

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

//Return a map to find headers for MAGEE
std::ext::map_str_int createHeaderMap(std::ext::V_string const& headers) 
{
    std::ext::map_str_int headerMap;
    // First column is intercept all filled with one
    for (int i = 0; i < headers.size(); ++i) 
    {
        headerMap[headers[i]] = i + 1;
    }
    return headerMap;
}

std::ext::V_int conv_stdvs2stdvi(std::ext::V_string const& strings) 
{
    std::ext::V_int result;
    result.reserve(strings.size());
    for (const std::string& str : strings) 
    {
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
            result = eta; 
        } 
        else if (link == "log") 
        {
            result = eta.array().exp().matrix(); 
        } 
        else if (link == "sqrt") 
        {
            result = eta.array().pow(2).matrix(); 
        } 
        else 
        {
            std::cerr << "Unsupported link function for the Gaussian family." << "\n";
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
        std::cerr << "Unsupported family type." << "\n";
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


SpaMat slice_mat(SpaMat const& spm, std::ext::V_int const& indices)
{
    SpaMatRow A = spm;  // Row-major copy for fast row iteration

    std::ext::V_int new_row(A.rows(), -1);
    for (int i = 0; i < indices.size(); ++i) new_row[indices[i]] = i;

    std::ext::VecTuples4spmat trip;
    trip.reserve(A.nonZeros()); // 

    for (int r_old : indices)
    {
        int r_new = new_row[r_old];
        for (SpaMatRow::InnerIterator it(A, r_old); it; ++it)
        {
            trip.emplace_back(r_new, it.col(), it.value());
        }
    }

    SpaMat out(indices.size(), spm.cols());
    out.setFromTriplets(trip.begin(), trip.end());
    out.makeCompressed();
    return out;
}

//Slice dense mat based on rows indixes
DensMat slice_mat(DensMat const& dm, std::ext::V_int const& indices) 
{
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
SpaMat slice_mat(SpaMat const& spm, std::ext::V_int rows, std::ext::V_int cols, bool check_size)
{
    if (check_size &&
        rows.size() == static_cast<size_t>(spm.rows()) &&
        cols.size() == static_cast<size_t>(spm.cols())) {
        return spm; // identity selection
    }

    std::ext::V_int row_map(spm.rows(), -1);
    for (size_t i = 0; i < rows.size(); ++i)
        row_map[rows[i]] = static_cast<int>(i);

    SpaMat out(rows.size(), cols.size());

    size_t num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 1;
    num_threads = std::min(num_threads, std::max<size_t>(cols.size(), 1));

    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    std::vector<std::ext::VecTuples4spmat> per_thread_triples(num_threads);

    auto worker = [&](size_t t, size_t c0, size_t c1)
    {
        auto& tripel = per_thread_triples[t];

        size_t nz          = static_cast<size_t>(spm.nonZeros());
        size_t total_cols  = static_cast<size_t>(spm.cols());
        size_t denom       = std::max<size_t>(1, total_cols);
        size_t span_cols   = (c1 - c0);
        size_t avg_per_col = nz / denom;
        tripel.reserve(span_cols * (avg_per_col + 1));

        for (size_t jj = c0; jj < c1; ++jj) {
            int old_j = cols[jj];
            for (SpaMat::InnerIterator it(spm, old_j); it; ++it) {
                int i_old = it.row();
                int i_new = row_map[i_old];
                if (i_new != -1) {
                    tripel.emplace_back(i_new, static_cast<int>(jj), it.value());
                }
            }
        }
    };

    // Split columns across threads
    size_t chunk = cols.size() / num_threads;
    size_t rem   = cols.size() % num_threads;
    size_t start = 0;
    for (size_t t = 0; t < num_threads; ++t) 
    {
        size_t end = start + chunk + (t < rem ? 1 : 0);
        threads.emplace_back(worker, t, start, end);
        start = end;
    }
    for (auto& th : threads) th.join();

    // Merge triplets and build
    std::ext::VecTuples4spmat triples;
    size_t total = 0; for (auto& v : per_thread_triples) total += v.size();
    triples.reserve(total);
    for (auto& v : per_thread_triples) {
        triples.insert(triples.end(),
                     std::make_move_iterator(v.begin()),
                     std::make_move_iterator(v.end()));
    }

    out.setFromTriplets(triples.begin(), triples.end());
    return out;
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
    for (size_t i = 0; i < ind1.size(); ++i) 
    {
        if (ind1[i]) {
            int_indices1.push_back(i);
        }
    }
    std::ext::V_int int_indices2;
    for (size_t i = 0; i < ind2.size(); ++i) 
    {
        if (ind2[i]) 
        {
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
    double max_difference_tau   = ((tau - tau0).array().abs() / (tau.array().abs() + tau0.  array().abs() + tol)).maxCoeff();
    double max_difference = std::max(max_difference_alpha, max_difference_tau);
    if ((2 * max_difference) < tol) 
    {
        return true; // Converged
    }
    if ((tau.array().abs().maxCoeff()) > pow(tol, -2)) {
        std::cerr << "Large variance estimate observed in the iterations, model not converged..." << "\n";
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

void GMMAT::set_dspy(DensMat& dspy, DensVec const& wpy, int ng)//Mask for groups
{
    dspy.setZero();

    for (int i = 0; i < ng; ++i)
    {
        for (int j : m_group_idx[i])   // Each group_idx[i] is a vector<int> of row indices
        {
            dspy(j, i) = wpy(j);
        }
    }
}


void GMMAT::set_VZpy(DensMat& VZpy, DensVec const& Zpy, int nk)
{
    bool has_random_slope = (m_rand_slope.size() > 0);
    int dimZ = m_Z.cols();
    int n = dimZ / (has_random_slope ? 2 : 1);
    DensVec mask10, mask01;

    for (int i = 0; i < nk; ++i)
    {
        bool is_identity = (i == nk - 1);

        //  No random slope   → (dimZ = n)
        if (!has_random_slope)
        {
            if (is_identity) 
            {
                // Identity block: VZpy[,i] = Zpy (n×1)
                VZpy.col(i) = Zpy;
            } 
            else 
            {
                // kinship block: VZpy[,i] = K * Zpy
                SpaMat &K = m_vkins_sp[i].get_spmat();
                VZpy.col(i) = K * Zpy;   // n×1
            }
            continue;   
        }
        // With random slope (dimZ = 2n)
        if (is_identity)
        {
            // Identity matrix case
            // Intercept block
            VZpy.col(i).head(n) = Zpy.head(n);       // top half
            VZpy.col(i).tail(n).setZero();           // bottom half
            // covariance block
            VZpy.block(0, i + nk, n, 1) = Zpy.segment(n, n);
            VZpy.block(n, i + nk, n, 1) = Zpy.segment(0, n);
            // variance block
            VZpy.col(i + 2*nk).head(n).setZero();    // top half
            VZpy.col(i + 2*nk).tail(n) = Zpy.tail(n);// bottom half
        }
        else
        {
            // Kinship exist
            SpaMat &K = m_vkins_sp[i].get_spmat();
            // Split Zpy into upper/lower halves
            DensVec Zpy_top    = Zpy.segment(0, n);
            DensVec Zpy_bottom = Zpy.segment(n, n);
            // Compute kinsZpy 
            DensVec kinsZpy(dimZ);
            kinsZpy.segment(0, n) = K * Zpy_top;
            kinsZpy.segment(n, n) = K * Zpy_bottom;
            // Intercept block
            VZpy.col(i).head(n) = kinsZpy.head(n);       // top half
            VZpy.col(i).tail(n).setZero();           // bottom half
            // Covariance block
            VZpy.block(0, i + nk, n, 1) = kinsZpy.segment(n, n);  // bottom → top
            VZpy.block(n, i + nk, n, 1) = kinsZpy.segment(0, n);  // top → bottom

            // variance block
            VZpy.col(i + 2*nk).head(n).setZero();      // top half
            VZpy.col(i + 2*nk).tail(n) = kinsZpy.tail(n);           // bottom half
        }
    }
}

static void ensure_compressed(SpaMat& M) //static for ODR
{
    if (!M.isCompressed()) M.makeCompressed();
}


double hadamard_sum_block(SpaMat& Z, const SpaMat& K, int n, int rowOff, int colOff)
{
    ensure_compressed(Z);
    double s = 0.0;

    for (int c = 0; c < n; ++c)
    {
        const int zcol = colOff + c;

        SpaMat::InnerIterator itZ(Z, zcol);
        SpaMat::InnerIterator itK(K, c);

        while (itZ && itK)
        {
            const int rZ = itZ.row();
            if (rZ < rowOff) { ++itZ; continue; }
            if (rZ >= rowOff + n) break; 

            const int rK_abs = rowOff + itK.row(); // K row mapped into Z's block rows

            if (rZ == rK_abs) 
            {
                s += itZ.value() * itK.value();
                ++itZ; ++itK;
            } 
            else if (rZ < rK_abs) 
            {
                ++itZ;
            } 
            else 
            {
                ++itK;
            }
        }
    }
    return s;
}


void GMMAT::calc_tr_corr(int i, DensVec &score, DensMat const& Ztsigma_ix, 
    DensMat const& Ztsigma_ixcov, SpaMat& Ztsigma_iZ, 
    int nk, int dimZ, int const& ng, bool has_random_slope)
    {
        int kin_index = i - ng; 
    // ----- identity matrix case -----
        if (kin_index == nk - 1)
        {
            if (!has_random_slope)
            {
                double tr = Ztsigma_iZ.diagonal().sum();   

                double corr = Ztsigma_ix.cwiseProduct(Ztsigma_ixcov).sum();
                score(i) -= (tr - corr);
            }
            if (has_random_slope)
            {
                int nind = dimZ / 2; 
                DensVec diagZ = Ztsigma_iZ.diagonal();     // length 2n (DensVec)
                double tr11 = diagZ.head(nind).sum();
                
                double corr11 = Ztsigma_ix.topRows(nind)
                        .cwiseProduct(Ztsigma_ixcov.topRows(nind)).sum();

                score(i) -= (tr11 - corr11);

                double tr12 = 0.0;
                for (int k = 0; k < nind; ++k)
                {
                    tr12 += Ztsigma_iZ.coeff(k, nind + k);
                }

                double corr12 =
                    Ztsigma_ix.topRows(nind)
                        .cwiseProduct(Ztsigma_ixcov.bottomRows(nind)).sum();
                score(i + nk) -= 2.0 * (tr12 - corr12);
                
                double tr22 = diagZ.tail(nind).sum();
                double corr22 =
                    Ztsigma_ix.bottomRows(nind)
                        .cwiseProduct(Ztsigma_ixcov.bottomRows(nind)).sum();

                score(i + 2*nk) -= (tr22 - corr22);
            }
        }

        // ----- kinship matrix case -----
        else if (kin_index >= 0 && kin_index < nk - 1)
        {
            const SpaMat& K =  m_vkins_sp[kin_index].get_spmat();  // n × n
        if (!has_random_slope)
        {
            DensMat kinsZtsigma_ixcov = K * Ztsigma_ixcov;   // n×p
            double tr = hadamard_sum_block(Ztsigma_iZ, K, dimZ, 0, 0);
            double corr = Ztsigma_ix.cwiseProduct(kinsZtsigma_ixcov).sum();
            score(i) -= (tr - corr);
        }
        if (has_random_slope)
        {
            int nind = dimZ / 2;              
            int p = Ztsigma_ix.cols();

            DensMat top = K * Ztsigma_ixcov.topRows(nind);      // n×p
            DensMat bot = K * Ztsigma_ixcov.bottomRows(nind);   // n×p
            DensMat kinsZtsigma_ixcov(dimZ, p);
            kinsZtsigma_ixcov.topRows(nind)    = top;
            kinsZtsigma_ixcov.bottomRows(nind) = bot;

            // score(i): top-left block
            double tr11   = hadamard_sum_block(Ztsigma_iZ, K, nind, 0, 0);
            double corr11 = Ztsigma_ix.topRows(nind)
                                .cwiseProduct(kinsZtsigma_ixcov.topRows(nind))
                                .sum();
            score(i) -= (tr11 - corr11);

            // score(i+nk): (Z12 + Z21) block
            double tr12 = hadamard_sum_block(Ztsigma_iZ, K, nind, 0, nind)
                        + hadamard_sum_block(Ztsigma_iZ, K, nind, nind, 0);

            double corr12 = Ztsigma_ix.topRows(nind)
                            .cwiseProduct(kinsZtsigma_ixcov.bottomRows(nind))
                            .sum();
            double corr21 = Ztsigma_ix.bottomRows(nind)
                            .cwiseProduct(kinsZtsigma_ixcov.topRows(nind)).sum();
                            
            score(i + nk) -= (tr12 - corr12 - corr21);

            // score(i+2nk): bottom-right block
            double tr22   = hadamard_sum_block(Ztsigma_iZ, K, nind, nind, nind);
            double corr22 = Ztsigma_ix.bottomRows(nind)
                            .cwiseProduct(kinsZtsigma_ixcov.bottomRows(nind)).sum();
            score(i + 2*nk) -= (tr22 - corr22);
        }
    }
}


void GMMAT::set_mtau(DensVec V_tr_corr, DensVec tau0, SpaMat& Ztsigma_iZ,
    DensMat const& Ztsigma_ix, DensMat const& Ztsigma_ixcov, 
    DensVec diagp, int nk, int dimZ, int ng, bool has_random_slope)
{
    int size_m_tau = nk;
    if(has_random_slope)
    {
        size_m_tau = (nk * 3) + ng;
    }
    else 
    {
        size_m_tau = nk + ng;
    }

    for(std::size_t i{0}; i < size_m_tau; ++i)
    {         
        double tau0_sq = tau0(i) * tau0(i);  
        if(i < ng)
        {
            DensVec diagp_div_w = diagp.cwiseQuotient(m_sqrtW);
            double diagp_sl = slice_vec(diagp_div_w, m_group_idx[i]).sum();
            double wpy_diagp = V_tr_corr(i) - diagp_sl; //V_tr_corr -> wpy pow 2
            m_tau(i) = std::max(0.0, tau0(i) + tau0_sq * wpy_diagp / m_Nobs);
        }
        else
        {
            calc_tr_corr(i, V_tr_corr, Ztsigma_ix, Ztsigma_ixcov, Ztsigma_iZ, 
                    nk, dimZ, ng, has_random_slope);
            double score_tau = tau0(i) + (tau0_sq * V_tr_corr(i) / m_Nobs);  
            if(m_covariance_idx.size() > 0)
            {
                if(std::find(m_idxtau2.begin(), m_idxtau2.end(), i) != m_idxtau2.end())
                {
                    m_tau(i) = 0.0;
                }
                else
                {
                    m_tau(i) = std::max(0.0, score_tau);
                }
            }
            else
            {
                m_tau(i) = std::max(0.0, score_tau);
            }   
        }
    }
}


void GMMAT::set_score(DensVec &score, DensMat const& Ztsigma_ix, 
    DensMat const& Ztsigma_ixcov, SpaMat& Ztsigma_iZ, 
    DensVec const& diagp, int nk,
    int dimZ, int const& ng, bool has_random_slope)
{
    int p = Ztsigma_ix.cols();
    //residuals for ngs
    for (int i = 0; i < ng + nk; ++i)
    {
        if (i < ng)
        {
            double sum_diag = 0.0;
            for (int idx : m_group_idx[i])
                sum_diag += diagp(idx);

            score(i) -= sum_diag;
        }
        else
        {
            calc_tr_corr(i, score, Ztsigma_ix, Ztsigma_ixcov, Ztsigma_iZ, 
                nk, dimZ, ng, has_random_slope);
        }
    }
}

void GMMAT::set_ai_low_ng(int i, DensVec& score, DensMat& ai, DensVec const& wpy, Fit const& fit, DensVec const& py, DensVec diagp, DensMat sigma_ixcov)
{
    DensVec results = slice_vec(wpy, m_group_idx[m_idxtau[i]]).cwiseProduct(slice_vec(py, m_group_idx[m_idxtau[i]])) - slice_vec(diagp, m_group_idx[m_idxtau[i]]);
    score(i) = results.sum();
    for(size_t j{0}; j <= i; ++j)
    {
        SpaMat sigma_i_ij = slice_mat(fit.sigma_i, m_group_idx[m_idxtau[j]], m_group_idx[m_idxtau[i]]);
        DensVec sigma_i_ij_wPy_j = crossprod(sigma_i_ij, slice_vec(wpy, m_group_idx[m_idxtau[j]]));
        double main_term = (slice_vec(wpy, m_group_idx[m_idxtau[i]]).cwiseProduct(sigma_i_ij_wPy_j)).sum();
        DensMat sigma_ixcovi = slice_mat(sigma_ixcov, m_group_idx[m_idxtau[i]]);
        DensVec sigma_ixcoviwpyi = crossprod(sigma_ixcovi, slice_vec(wpy, m_group_idx[m_idxtau[i]]));
        DensMat sigma_ixj = slice_mat(fit.sigma_ix, m_group_idx[m_idxtau[j]]);
        DensVec sigma_ixjwpyj = crossprod(sigma_ixj, slice_vec(wpy, m_group_idx[m_idxtau[j]]));
        double correction_term = (sigma_ixcoviwpyi.cwiseProduct(sigma_ixjwpyj)).sum();
        ai(i,j) = main_term - correction_term;
        if(j != i) 
        {
            ai(j,i) = ai(i,j);
        }       
    }  
}


void GMMAT::set_ai_high_ng(int i, DensVec& score, DensMat& ai, DensVec const& wpy,
     Fit const& fit, DensVec const& py, DensMat const& sigma_ixcov, int ng)
{
    SpaMat A = m_vkins_sp[m_idxtau[i]-ng].get_spmat();
    DensVec apy = crossprod(A, py);
    DensVec papy_main = crossprod(fit.sigma_i, apy);   
    DensVec papy_correction = tcrossprod(fit.sigma_ix, crossprod(sigma_ixcov, apy).transpose());
    DensVec papy = papy_main - papy_correction;
    double sigma_ia = (fit.sigma_i.cwiseProduct(A)).sum(); 
    double sigma_iax = (fit.sigma_ix.cwiseProduct(crossprod(A, sigma_ixcov))).sum();
    score(i) = (m_Y.cwiseProduct(papy)).sum() - (sigma_ia - sigma_iax);
    for(size_t j{0}; j <= i; ++j)
    {
        if(m_idxtau[j] < ng)
        {
            ai(i,j) = (slice_vec(wpy, m_group_idx[m_idxtau[j]]).cwiseProduct(slice_vec(papy, m_group_idx[m_idxtau[j]]))).sum();
        }
        else
        {     
            SpaMat B = m_vkins_sp[m_idxtau[j]-ng].get_spmat();
            ai(i,j) = (py.cwiseProduct(crossprod(B, papy))).sum();
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
void GMMAT::calc_covariance(int &kins_size, int ng)
{
    if(m_vkins_sp[0].m_ratio_Nob2N < m_vkins_sp[0].m_thr)
    {
        m_vkins_sp.resize(3*kins_size);
    }
    for(int i{0}; i < kins_size; ++i)
    {
        if(m_vkins_sp[0].m_ratio_Nob2N < m_vkins_sp[0].m_thr)
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
        }
        m_covariance_idx.push_back({kins_size + i + ng, i + ng, kins_size * 2 + i + ng});
    }
}


void GMMAT::build_Psi(int const ng, bool calc_diag_kin)
{
    // ----------longitudinal (RI only) ----------
    if (m_modeltype == LONGITUDINAL_RI)
    {
        m_Psi.resize(m_N, m_N);// N in N mat
        const double tauI_eps = 1e-6;
        
        std::ext::VecTuples4spmat triples;
        
        // Reserve: sum of all Φ_i nnz + identity diagonal
        size_t total_nnz = m_N;

        int K = 0;
        for (auto& ks : m_vkins_sp) 
        {
            if (!ks.kin.m_null_kin) 
            {   
                ++K;
                total_nnz += ks.get_spmat().nonZeros();
            }
        }

        triples.reserve(total_nnz);

        // ---- Loop over kinship matrices ----
        for (int i = 0; i < K; i++)
        {
            SpaMat phi = m_vkins_sp[i].get_spmat();

            double th = std::max(m_tau[ng + i], tauI_eps);   // θ2, θ3, ..., θ(K+1)

            for (int col = 0; col < phi.outerSize(); ++col)
            {
                for (SpaMat::InnerIterator it(phi, col); it; ++it)
                {
                    triples.emplace_back(it.row(), it.col(), th * it.value());
                }
            }
        }

        // ---- Add θ(identity) * I ----

        double thI = std::max(m_tau[ng + K], tauI_eps);   // last one → θ(ng+K)
        for (int j = 0; j < m_N; j++)
            triples.emplace_back(j, j, thI);

        m_Psi.setFromTriplets(triples.begin(), triples.end());
        return;
    }

    // ---------- longitudinal + slope ----------
    // Reserve triplets: each B_i contributes ~ (Phi_i.nnz * 4)
    DensVec diag_block1;
    DensVec diag_block4;
    DensVec diag_block2;
    DensVec diag_kin;
    double mean_diag_block1 = 1.0; 
    double mean_diag_block4 = 1.0;
    double mean_diag_block2 = 1.0;
    
    if(calc_diag_kin)
    {
        diag_block1.resize(m_N);
        diag_block4.resize(m_N);
        diag_block2.resize(m_N);

        for (int i = 0; i < m_N; ++i)
        {
            // Block (1,1) diagonal
            diag_block1(i) = m_Zt_Z.coeff(i, i);

            // Block (2,2) diagonal
            diag_block4(i) = m_Zt_Z.coeff(m_N + i, m_N + i);

            // Block (1,2) diagonal
            diag_block2(i) = m_Zt_Z.coeff(i, m_N + i) * 2;
        }
    }
    
    
    size_t total_nnz = 0;
    int K = 0;
    for (auto& ks : m_vkins_sp) 
    {
        if (!ks.kin.m_null_kin) 
        {   
            ++K;
            total_nnz += ks.get_spmat().nonZeros();
        }
    }
    int kins_size = K + 1;          // plus identity block
    total_nnz += m_N;               // identity block
    
    m_Psi.resize(2*m_N, 2*m_N);
    std::ext::VecTuples4spmat triples;
    triples.reserve(4 * total_nnz);

    // ---------- loop over kinship matrixes and I matrix ----------
    for (int i = 0; i <= K; i++)
    {
        const SpaMat *phi_ptr;

        if (i < K)
        {
            phi_ptr = &m_vkins_sp[i].get_spmat(); 
            if(calc_diag_kin)
            {
                diag_kin = phi_ptr->diagonal();
                mean_diag_block1 = diag_block1.dot(diag_kin) / m_Nobs;
                mean_diag_block2 =  diag_block2.dot(diag_kin) / m_Nobs;
                mean_diag_block4 =  diag_block4.dot(diag_kin) / m_Nobs;
            }
        }
        else
        {
            phi_ptr = nullptr;    // identity I
            if(calc_diag_kin)
            {
                mean_diag_block1 = diag_block1.sum() / m_Nobs;
                mean_diag_block2 =  diag_block2.sum() / m_Nobs;
                mean_diag_block4 =  diag_block4.sum() / m_Nobs;
            }
        }
        // θ-indexing for block i
        m_tau[ng + i] = m_tau[ng + i] / mean_diag_block1;
        m_tau[ng + i + kins_size] = m_tau[ng + i + kins_size] / mean_diag_block2;
        m_tau[ng + i + 2 * kins_size] = m_tau[ng + i + 2 * kins_size] / mean_diag_block4;

        double th11 = m_tau[ng + i] ;
        double th12 = m_tau[ng + i + kins_size] ;
        double th22 = m_tau[ng + i + 2 * kins_size] ;

        // small det near zero  -> singular mat 
        const double eps = 1e-10;     // min variance
        const double rho = 0.999;      // <-- margin 

        th11 = std::max(th11, eps);
        th22 = std::max(th22, eps);
        
        double max_abs_th12 = rho * std::sqrt(th11 * th22);
        th12 = std::clamp(th12, -max_abs_th12, max_abs_th12);

        // block (1,1): th11 * Φ_i or th11*I
        if (phi_ptr)
        {
            const SpaMat &phi = *phi_ptr;
            for (int col = 0; col < phi.outerSize(); ++col)
                for (SpaMat::InnerIterator it(phi, col); it; ++it)
                    triples.emplace_back(it.row(), it.col(), th11 * it.value());
        }
        else
        {
            for (int j = 0; j < m_N; j++)
                triples.emplace_back(j, j, th11);
        }
        
        // block (1,2): th12 * Φ_i or th12*I
        if (phi_ptr)
        {
            const SpaMat &phi = *phi_ptr;
            for (int col = 0; col < phi.outerSize(); ++col)
                for (SpaMat::InnerIterator it(phi, col); it; ++it)
                    triples.emplace_back(it.row(), m_N + it.col(), th12 * it.value());
        }
        else
        {
            for (int j = 0; j < m_N; j++)
                triples.emplace_back(j, m_N + j, th12);
        }

        // block (2,1): same as (1,2)
        if (phi_ptr)
        {
            const SpaMat &phi = *phi_ptr;
            for (int col = 0; col < phi.outerSize(); ++col)
                for (SpaMat::InnerIterator it(phi, col); it; ++it)
                    triples.emplace_back(m_N + it.row(), it.col(), th12 * it.value());
        }
        else
        {
            for (int j = 0; j < m_N; j++)
                triples.emplace_back(m_N + j, j, th12);
        }
        
        // block (2,2): th22 * Φ_i or th22*I
        if (phi_ptr)
        {          
            const SpaMat &phi = *phi_ptr;
            
            for (int col = 0; col < phi.outerSize(); ++col)
                for (SpaMat::InnerIterator it(phi, col); it; ++it)
                    triples.emplace_back(m_N + it.row(), m_N + it.col(), th22 * it.value());
            
        }
        else
        {
            for (int j = 0; j < m_N; j++)
                triples.emplace_back(m_N + j, m_N + j, th22);
        }
    }
    
    m_Psi.setFromTriplets(triples.begin(), triples.end());   
}

void GMMAT::build_Z()
{
    m_Z.resize(m_Nobs, 0);

    std::ext::VecTuples4spmat triples;

    switch (m_modeltype)
    {
        case LONGITUDINAL_RI:
        {
            m_Z = m_J;    // Assumes m_J is already Nobs x N
            return;
        }

        case LONGITUDINAL_RS:
        {
            // Z = [ J | diag(E)J ]  (Nobs x 2N)
            m_Z.resize(m_Nobs, 2 * m_N);
            triples.reserve(2 * m_J.nonZeros());

            // Block 1: J
            for (int col = 0; col < m_J.outerSize(); ++col)
                for (SpaMat::InnerIterator it(m_J, col); it; ++it)
                    triples.emplace_back(it.row(), it.col(), it.value());

            // Block 2: diag(E) * J
            for (int col = 0; col < m_J.outerSize(); ++col)
                for (SpaMat::InnerIterator it(m_J, col); it; ++it)
                {
                    const auto i = it.row();
                    triples.emplace_back(i, m_N + col, m_rand_slope(i) * it.value());
                }
            m_Z.setFromTriplets(triples.begin(), triples.end());
            return;
        }

        default:
        {
            // keep m_Z as (Nobs x 0)
            return;
        }
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
    int fixtau_0_counts = std::count(m_fixtau.begin(), m_fixtau.end(), 0); 
    int ng = m_group_idx.size(); 
    DensVec diag_sigma = DensVec::Zero(m_Nobs); 
    DensVec score; 
    DensVec wpy;
    DensVec py;
    DensVec diagp;
    DensVec apy;
    DensVec papy;
    fit_to_return.W  = W;
    if(m_family_t == "binomial")
    {
        m_tau(0) = 1.0;
    }
    for(int i{0}; i < ng; ++i)
    {
        for (int j : m_group_idx[i]) 
        {
            diag_sigma[j] = m_tau(i) / fit_to_return.W(j);
        }
    }

    if(m_vkins_sp[0].m_ratio_Nob2N >= m_vkins_sp[0].m_thr)
    {
        fit_to_return.sigma_i.resize(0,0);
        int nk = 0;
        if(m_rand_slope.size() > 0)
        {
            nk =  m_kins_size / 3; // kinsips + I matrixes
        }
        else
        {
             nk =  m_kins_size; // kinsips + I matrixes
        }

        fit_to_return.diag_sigma_i = diag_sigma.cwiseInverse();

        m_diag_sigma_im_Z = details::diag(fit_to_return.diag_sigma_i) * m_Z; 
        SpaMat Ztdiag_sigma_iZ = crossprod(m_Z, m_diag_sigma_im_Z); //2N in 2N or N in N
        build_Psi(ng);
        SpaMat Pchol = SparseInverse::inv_spamat_chol(SparseInverse::inv_spamat_chol(m_Psi, true) + Ztdiag_sigma_iZ); //2N in 2N or N in N
        fit_to_return.diag_sigma_i_ZPchol = m_diag_sigma_im_Z * Pchol; //Nobs * 2N sparse or in N
        DensMat YX(m_Y.rows(), m_Y.cols() + m_X.cols());
        YX << m_Y, m_X;      // [Y | X]
        DensMat sigma_im_YX = apply_wb(YX, fit_to_return.diag_sigma_i, fit_to_return.diag_sigma_i_ZPchol); //Nobs in (1+P)
        DensVec sigma_im_Y = sigma_im_YX.col(0);

        // All remaining columns → Sigma_iX
        fit_to_return.sigma_ix = sigma_im_YX.rightCols(sigma_im_YX.cols() - 1);
        DensMat xsigma_ix = crossprod(m_X, fit_to_return.sigma_ix); 
        fit_to_return.cov = SparseInverse::inv(xsigma_ix); //P in P
        DensMat sigma_ixcov = tcrossprod(fit_to_return.sigma_ix, fit_to_return.cov);//Nobs in P
        fit_to_return.alpha = crossprod(fit_to_return.cov, crossprod(fit_to_return.sigma_ix, m_Y)); //P * Nobs
        py = sigma_im_Y - fit_to_return.sigma_ix * fit_to_return.alpha; //Nobs * 1
        fit_to_return.eta = m_Y - diag_sigma.cwiseProduct(py);
        DensMat Ztsigma_ix = crossprod(m_Z, fit_to_return.sigma_ix); //2N in P or N in P
        DensMat  Ztsigma_ixcov = tcrossprod(Ztsigma_ix, fit_to_return.cov); //2N in P or N in P
        if(fixtau_0_counts > 0)
        {
            m_idxtau = which(m_fixtau, [](int i){return (i==0);});
            wpy = py.cwiseQuotient(fit_to_return.W); //Nobs * 1
            DensVec Zpy = crossprod(m_Z, py); //2N vec or N if no random slope
            SpaMat Ztdiag_sigma_iZ_Pchol = crossprod(m_Z, fit_to_return.diag_sigma_i_ZPchol); //2n x 2n sparse
            SpaMat Ztsigma_iZ = Ztdiag_sigma_iZ -tcrossprod(Ztdiag_sigma_iZ_Pchol, Ztdiag_sigma_iZ_Pchol);
            DensVec diagp;
            {
                DensVec rowsum = DensVec::Zero(fit_to_return.diag_sigma_i_ZPchol.rows());
    
                for (int k=0; k < fit_to_return.diag_sigma_i_ZPchol.outerSize(); ++k)
                {
                    for (SpaMat::InnerIterator it(fit_to_return.diag_sigma_i_ZPchol, k); it; ++it)
                    {
                        rowsum(it.row()) += it.value() * it.value();
                    }
                }
                DensVec sigma_i_diag = fit_to_return.diag_sigma_i - rowsum;// To cal diag(sigma_i)
                diagp = (sigma_i_diag - ((fit_to_return.sigma_ix.cwiseProduct(sigma_ixcov)).rowwise().sum())).cwiseQuotient(fit_to_return.W); //Nobs in 1;
            }
            
            DensVec score;
            DensMat dspy(m_Nobs, m_kins_size + ng);
            set_dspy(dspy, wpy, ng);
            DensMat VZpy(m_Z.cols(), m_kins_size); //2N in nk*3 or nk
            VZpy.setConstant(std::numeric_limits<double>::quiet_NaN());
            set_VZpy(VZpy, Zpy, nk);
            dspy.rightCols(dspy.cols() - ng) = m_Z * VZpy; //Nobs nk*3
            DensMat pdspy = dspy.array().colwise() / diag_sigma.array();
            pdspy -= fit_to_return.diag_sigma_i_ZPchol * crossprod(fit_to_return.diag_sigma_i_ZPchol, dspy); //Nobs in nk*3 or nk
            pdspy -= fit_to_return.sigma_ix * crossprod(sigma_ixcov, dspy); //Nobs in nk*3 or nk
            DensMat AI = crossprod(dspy, pdspy); //nk in nk or nk *3
            score = crossprod(dspy, py); //nk or nk *3
            bool has_random_slope = (m_Z.cols() == 2*m_N);   // or (m_rand_slope.size() > 0)
            int dimZ = m_Z.cols();    // N or 2N

            set_score(score, Ztsigma_ix, Ztsigma_ixcov, Ztsigma_iZ, diagp,
                nk, dimZ, ng, has_random_slope);

            const int k = (m_idxtau.size());
            DensMat A_sub(k, k);
            DensVec b_sub(k);
 
            for (int a = 0; a < k; ++a)
            {
                const int ia = m_idxtau[a];
                b_sub(a) = score(ia);

                for (int b = 0; b < k; ++b)
                {
                    const int ib = m_idxtau[b];
                    A_sub(a, b) = AI(ia, ib);
                }
            }

            auto lu_decomp = A_sub.fullPivLu(); 

            if (lu_decomp.isInvertible()) 
            {
                DensVec dv = lu_decomp.solve(b_sub);
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
    }
    else
    {
        fit_to_return.diag_sigma_i.resize(0);
        fit_to_return.diag_sigma_i_ZPchol.resize(0,0);
        // Create diagonal mat using vector
        SpaMat sigma = details::diag(diag_sigma);  
        int kins_size = m_vkins_sp.size();

        for(int i{0}; i < kins_size; ++i)
        {
            SpaMat curr_kin_spmat = m_vkins_sp[i].get_spmat();
            sigma = sigma + m_tau(i + ng) * curr_kin_spmat;
        }
    
        fit_to_return.sigma_i =  SparseInverse::inv_spamat_chol(sigma, true);
        fit_to_return.sigma_ix = crossprod(fit_to_return.sigma_i, m_X);
        DensMat xsigma_ix = crossprod(m_X, fit_to_return.sigma_ix); 
        fit_to_return.cov = SparseInverse::inv(xsigma_ix);
        DensMat sigma_ixcov = tcrossprod(fit_to_return.sigma_ix, fit_to_return.cov);
        fit_to_return.alpha = crossprod(fit_to_return.cov, crossprod(fit_to_return.sigma_ix, m_Y)); 
        fit_to_return.eta = m_Y - diag_sigma.cwiseProduct((crossprod(fit_to_return.sigma_i, m_Y) - tcrossprod(fit_to_return.sigma_ix, fit_to_return.alpha.transpose())));
        
        if(fixtau_0_counts > 0)
        {
            m_idxtau = which(m_fixtau, [](int i){return (i==0);});
            py = crossprod(fit_to_return.sigma_i, m_Y) - tcrossprod(fit_to_return.sigma_ix, crossprod(sigma_ixcov, m_Y).transpose()); //Nobs in 1
            wpy = py.cwiseQuotient(fit_to_return.W); 
            DensVec diagp = (fit_to_return.sigma_i.diagonal() - ((fit_to_return.sigma_ix.cwiseProduct(sigma_ixcov)).rowwise().sum())).cwiseQuotient(fit_to_return.W);
            DensMat AI = DensMat::Constant(fixtau_0_counts, fixtau_0_counts, std::numeric_limits<double>::quiet_NaN());
            DensVec score = DensVec::Constant(fixtau_0_counts, std::numeric_limits<double>::quiet_NaN());
            set_ai(score, AI, wpy,fit_to_return, py, diagp, sigma_ixcov,ng, fixtau_0_counts);
            auto lu_decomp = AI.fullPivLu();

            if (lu_decomp.isInvertible()) 
            {
                DensVec dv = lu_decomp.solve(score);
                fit_to_return.dtau = std::make_optional(dv);
            } 
            else 
            {
                std::cout << "The matrix is not invertible, solve operation failed." << "\n";
                fit_to_return.dtau = std::nullopt;
                exit(EXIT_FAILURE); 
            }
            return fit_to_return;
        }
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
    DensVec wpy;

    int kins_size = 0;
    glmmkin.fit.alpha = fit_null.alpha;
    
    if(m_offset.size() < m_Nobs) 
    {
        m_offset = DensVec::Constant(m_Nobs, 0); 
    }
    m_tau = DensVec::Constant(m_kins_size + m_group_idx.size(), 0);
  
    m_fixtau.resize(m_kins_size + m_group_idx.size(), 0);  
    fit_null.calc_dmu_deta(m_family_t, m_Nobs);
    glmmkin.fit = fit_null;
    m_Y = fit_null.eta - m_offset + (m_y - glmmkin.fit.mu).cwiseQuotient(glmmkin.fit.dmu_deta); 
    m_sqrtW = glmmkin.fit.calc_sqrtW();
    glmmkin.fit.W = glmmkin.fit.dmu_deta;
    if(m_family_t == "binomial")
    {
        m_tau(0) = 1.0;
        m_fixtau[0] = 1;
    }

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
        auto tau_value = calc_variance(m_Y) / (m_kins_size + ng);
        
        for(int idx : m_idxtau)
        {
            m_tau(idx) = tau_value;//m_fixtau(0) is 1 for binomia  -> m_idxtau(0)=1, idx!=0  -> m_tau(0) will be 1 always
        }

        if(m_covariance_idx.size() > 0)
        {
            for(auto id : m_idxtau2)
            {
                m_tau(id) = 0;
            }
        }

        DensVec diag_sigma = DensVec::Zero(m_Nobs); 
        for(int i{0}; i < ng; ++i)
        {
            for (int j : m_group_idx[i]) 
            {
                diag_sigma[j] = m_tau(i) / glmmkin.fit.W(j);
            }
        } 

        if(m_vkins_sp[0].m_ratio_Nob2N >= m_vkins_sp[0].m_thr)
        {
            glmmkin.fit.sigma_i.resize(0,0);
            int nk = 0;
            if(m_rand_slope.size() > 0)
            {
                nk =  m_kins_size / 3; // kinsips + I matrixes
            }
            else
            {
                nk =  m_kins_size; // kinsips + I matrixes
            }
            glmmkin.fit.diag_sigma_i = diag_sigma.cwiseInverse();
            build_Z();
            m_diag_sigma_im_Z = details::diag(glmmkin.fit.diag_sigma_i) * m_Z; 
            SpaMat Ztdiag_sigma_iZ = crossprod(m_Z, m_diag_sigma_im_Z); //2N in 2N
            
            // Calculate Zt_Z to do -->  m_tau(i + ng) / (curr_kin_spmat).diagonal().mean();
            m_Zt_Z.resize(m_Z.cols(), m_Z.cols());
            m_Zt_Z = m_Z.transpose() * m_Z;
            bool calc_diag_kin = true;
            build_Psi(ng, calc_diag_kin);

            SpaMat Pchol = SparseInverse::inv_spamat_chol(SparseInverse::inv_spamat_chol(m_Psi, true) + Ztdiag_sigma_iZ); //2N in 2N or N in N
            glmmkin.fit.diag_sigma_i_ZPchol = m_diag_sigma_im_Z * Pchol; //Nobs * 2N sparse or in N
            DensMat YX(m_Y.rows(), m_Y.cols() + m_X.cols());
            YX << m_Y, m_X;      // [Y | X]
            DensMat sigma_im_YX = apply_wb(YX, glmmkin.fit.diag_sigma_i, glmmkin.fit.diag_sigma_i_ZPchol); //Nobs in (1+P)
            DensVec sigma_im_Y = sigma_im_YX.col(0);

            glmmkin.fit.sigma_ix = sigma_im_YX.rightCols(sigma_im_YX.cols() - 1);
            DensMat xsigma_ix = crossprod(m_X, glmmkin.fit.sigma_ix); 
            glmmkin.fit.cov = SparseInverse::inv(xsigma_ix); //P in P
            DensMat sigma_ixcov = tcrossprod(glmmkin.fit.sigma_ix, glmmkin.fit.cov); //Nobs in P
            py = sigma_im_Y - glmmkin.fit.sigma_ix * glmmkin.fit.alpha; //Nobs * 1
            wpy = py.cwiseQuotient(m_sqrtW); 
            tau0 = m_tau;

            DensMat dspy(m_Nobs, m_kins_size + ng);
            set_dspy(dspy, wpy, ng);
            DensMat VZpy(m_Z.cols(), m_kins_size); //2N in nk*3
            DensVec Zpy = crossprod(m_Z, py);

            VZpy.setConstant(std::numeric_limits<double>::quiet_NaN());
            set_VZpy(VZpy, Zpy, nk);
            dspy.rightCols(dspy.cols() - ng) = m_Z * VZpy; 
            DensVec dspypy(dspy.cols());
            dspypy.setZero();
            dspypy.tail(dspy.cols() - ng) = crossprod(dspy.rightCols(dspy.cols() - ng), py);
            dspypy.head(ng) = dspy.leftCols(ng).array().square().colwise().sum().transpose().matrix();
            DensVec diagp;
            {
                DensVec rowsum = DensVec::Zero(glmmkin.fit.diag_sigma_i_ZPchol.rows());
    
                for (int k=0; k < glmmkin.fit.diag_sigma_i_ZPchol.outerSize(); ++k)
                {
                    for (SpaMat::InnerIterator it(glmmkin.fit.diag_sigma_i_ZPchol, k); it; ++it)
                    {
                        rowsum(it.row()) += it.value() * it.value();
                    }
                }
                DensVec sigma_i_diag = glmmkin.fit.diag_sigma_i - rowsum;
                int col_n = sigma_ixcov.cols();
                diagp = sigma_i_diag - (glmmkin.fit.sigma_ix.cwiseProduct(sigma_ixcov)) * DensVec::Ones(col_n);; //Nobs in 1;
            }

            SpaMat Ztdiag_sigma_iZ_Pchol = crossprod(m_Z, glmmkin.fit.diag_sigma_i_ZPchol); //2N x 2N sparse
            SpaMat Ztsigma_iZ = Ztdiag_sigma_iZ -tcrossprod(Ztdiag_sigma_iZ_Pchol, Ztdiag_sigma_iZ_Pchol);
            DensMat Ztsigma_ix = crossprod(m_Z, glmmkin.fit.sigma_ix); //2N in P or N in P
            DensMat  Ztsigma_ixcov = tcrossprod(Ztsigma_ix, glmmkin.fit.cov); //2N in P or N in P
            bool has_random_slope = (m_Z.cols() == 2*m_N);   // or (m_rand_slope.size() > 0)
            int dimZ = m_Z.cols(); 
            set_mtau(dspypy, tau0, Ztsigma_iZ, Ztsigma_ix, Ztsigma_ixcov, 
                    diagp, nk, dimZ, ng, has_random_slope); 
        }
        else
        {    
            glmmkin.fit.diag_sigma_i.resize(0);
            glmmkin.fit.diag_sigma_i_ZPchol.resize(0,0);
            SpaMat sigma = details::diag(diag_sigma); 
            for(int i{0}; i < m_kins_size; ++i)
            {
                SpaMat curr_kin_spmat = m_vkins_sp[i].get_spmat();
                m_tau(i + ng) = m_tau(i + ng) / (curr_kin_spmat).diagonal().mean();
                sigma = sigma + m_tau(i + ng) * curr_kin_spmat;
            }

            glmmkin.fit.sigma_i =  SparseInverse::inv_spamat_chol(sigma, true);
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
                    DensVec wpy = py.cwiseQuotient(m_sqrtW); // elementwise 
                    double tau0_sq = tau0(m_idxtau[i]) * tau0(m_idxtau[i]);                    
                    double wpyi = slice_vec(wpy, m_group_idx[m_idxtau[i]]).cwiseProduct(slice_vec(wpy, m_group_idx[m_idxtau[i]])).sum();
                    int col_n = sigma_ixcov.cols();
                    DensVec diagp = glmmkin.fit.sigma_i.diagonal() - (glmmkin.fit.sigma_ix.cwiseProduct(sigma_ixcov)) * DensVec::Ones(col_n);
                    DensVec diagpw = diagp.cwiseQuotient(m_sqrtW);
                    double diagpwi = slice_vec(diagpw, m_group_idx[m_idxtau[i]]).sum();
                    double score_component = wpyi - diagpwi;
                    m_tau(m_idxtau[i]) = std::max(0.0, tau0(m_idxtau[i]) + tau0_sq * score_component / m_Nobs);
                }
                else
                {   
                    apy = crossprod(m_vkins_sp[m_idxtau[i] - ng].get_spmat(), py);
                    papy =  crossprod(glmmkin.fit.sigma_i, apy) - tcrossprod(glmmkin.fit.sigma_ix, (crossprod(sigma_ixcov, apy)).transpose());
                    double tau0_sq = tau0(m_idxtau[i]) * tau0(m_idxtau[i]); 
                    double ypapy = (m_Y.cwiseProduct(papy)).sum ();
                    double sigma_ia = (glmmkin.fit.sigma_i.cwiseProduct(m_vkins_sp[m_idxtau[i] - ng].get_spmat())).sum();            
                    double sigma_ixaxcov = (glmmkin.fit.sigma_ix.cwiseProduct(crossprod(m_vkins_sp[m_idxtau[i] - ng].get_spmat(), sigma_ixcov))).sum();
                    double sigma_i_corrected = sigma_ia - sigma_ixaxcov;            
                    double score_component =  ypapy - sigma_i_corrected;
                    double score_tau_i = tau0(m_idxtau[i]) + (tau0_sq * score_component / m_Nobs);            

                    if(m_covariance_idx.size() > 0)
                    {
                        if(std::find(m_idxtau2.begin(), m_idxtau2.end(), m_idxtau[i]) != m_idxtau2.end())
                        {
                            m_tau(m_idxtau[i]) = 0.0;
                        }
                        else
                        {
                            m_tau(m_idxtau[i]) = std::max(0.0, score_tau_i);
                        }
                    }
                    else
                    {
                        m_tau(m_idxtau[i]) = std::max(0.0, score_tau_i);
                    }   
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
                exit(EXIT_FAILURE);
            }
        } 

        glmmkin.fit.cov = fit_glmm_ai.cov;
        glmmkin.fit.alpha = fit_glmm_ai.alpha;
        glmmkin.fit.eta = fit_glmm_ai.eta + m_offset;
        glmmkin.fit.mu = linkinv(fit_glmm_ai.eta, m_family_t, m_link);
        glmmkin.fit.calc_dmu_deta(m_family_t, m_Nobs);
        m_Y = glmmkin.fit.eta - m_offset + (m_y - glmmkin.fit.mu).cwiseQuotient(glmmkin.fit.dmu_deta);
        m_sqrtW = glmmkin.fit.calc_sqrtW();
        glmmkin.fit.W = glmmkin.fit.dmu_deta;
		glmmkin.fit.sigma_ix = fit_glmm_ai.sigma_ix;
		glmmkin.fit.sigma_i = fit_glmm_ai.sigma_i;
        glmmkin.fit.diag_sigma_i = fit_glmm_ai.diag_sigma_i;
        glmmkin.fit.diag_sigma_i_ZPchol = fit_glmm_ai.diag_sigma_i_ZPchol;
        std::cout << "Variance component estimates (m_tau):\n" << m_tau << '\n';
        std::cout << "Fixed-effect coefficient (alpha):\n" << glmmkin.fit.alpha << '\n';
        if(check_convergence(glmmkin.fit.alpha, alpha0, m_tau, tau0, tol, i, maxiter)) break;
    } 
    
    glmmkin.converged = i < maxiter ? true : false;
    glmmkin.residuals = m_y.cast<double>() - glmmkin.fit.mu;
    
    DensVec res_var = DensVec::Constant(m_Nobs, 1);
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
    //Fill scaled_residuals
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
    
    bool run_wb = m_vkins_sp[0].m_ratio_Nob2N >= m_vkins_sp[0].m_thr;
    if(method_optim == "AI")
    {
        // Set of groups_id--> male and female so group_unique size=2{0,1}
        auto group_unique = unique(group_id);
        extract_group_idx(group_unique, group_id); //Update m_group_idx
        int kins_size = 0;
        int ng = m_group_idx.size();
       
        std::ext::V_int fixrho_new;
        std::ext::V_int fixrho_old;
        if(run_wb)
        {
            for (auto const& ks : m_vkins_sp) 
            {
                if (!ks.kin.m_null_kin) 
                {   
                    ++kins_size;
                }
            }
            ++kins_size; //I matrix
        }
        else
        {
            kins_size = m_vkins_sp.size();
        }

        m_kins_size = kins_size;

        if(m_rand_slope.size() > 0)
        {
            //Fill m_covariance_idx
            calc_covariance(kins_size, ng);
            fixrho_old.resize(kins_size, 0);
            if(m_vkins_sp[0].m_ratio_Nob2N >= m_vkins_sp[0].m_thr)
            {
                m_kins_size = kins_size * 3; // To adjust tau size
            }
            else
            {
                m_kins_size = m_vkins_sp.size();
            }
        }

        std::ext::V_int fixtau_old(m_kins_size + ng, 0);
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
                fmt::print(stderr, "Brent is not available for the time being, stay in touch for updates.\n");
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
    glmmkin.run_wb = run_wb;
    return glmmkin;
}


Glmmkin GMMAT::glmmkin_init(std::ext::FitNull_f fit0, Pheno pheno,
                            std::ext::V_string cov_selected_hdrs_new,
                            std::string phenoname,
                            std::string const& id, 
                            std::ext::V_string int_cov_hdrs_name_new,
                            std::string randomSlopeName,
                            std::string const& group,
                            int center, int scale,
                            std::string const method, 
                            std::string method_optim, 
                            int maxiter,
                            double tol, double tau_min, 
                            double tau_max, int tau_region)
{
    Glmmkin glmmkin;
    pheno.m_sam_id = id;
    m_y = conv_stdVs2dV(pheno.m_data_frame.get_header(phenoname));
    m_Nobs = m_y.size();
    DensMat m_X_org = create_covdata(pheno.m_data_frame.copy_by_hdrs(cov_selected_hdrs_new));
    // Centering 
    if (center == 2)
    {
        center_dataframe(pheno.m_data_frame, int_cov_hdrs_name_new, true, scale);

        std::cout << "GEM centered only the interaction covariate(s): ";
        for (size_t i = 0; i < int_cov_hdrs_name_new.size(); ++i)
        {
            std::cout << int_cov_hdrs_name_new[i];
            if (i != int_cov_hdrs_name_new.size() - 1)
                std::cout << ",";
        }
        std::cout << "." << "\n";
        std::cout << "*********************************************************\n";
    }
    if (center == 1 || scale == 1)
    {
        center_dataframe(pheno.m_data_frame, cov_selected_hdrs_new, center, scale);

        std::cout << "Warning:" << std::endl;
        std::cout << "All interaction covariates, exposure and covariates were centered. "
            << "Meta-analysis is not recommended using centered results." << std::endl;
        std::cout << "*********************************************************\n";
    }
    if (center == 0)
    {
        std::cout<<"None of the interaction covariates, exposure and covariates were centered."<< std::endl;
        if (int_cov_hdrs_name_new.size() > 0)
        {
            std::cout<< "Warning:"<< std::endl;
            std::cout<< "It is strongly recommended to center all interaction covariates (program default) for better interpretation of the joint test for genetic main effects and gene-exposure interactions"<< "\n";
            std::cout << "*********************************************************\n";
        }
    }
    m_X = create_covdata(pheno.m_data_frame.copy_by_hdrs(cov_selected_hdrs_new));
    
    std::ext::V_string v_valid_methods {"REML", "ML"};
    std::string rand_slope_hdr = randomSlopeName;

    auto it = std::find(v_valid_methods.begin(), v_valid_methods.end(), method);
    if(it == v_valid_methods.end())
    {
        fmt::print(stderr, "Error: {} is not in GMMAT valid methods (REML, ML)\n", method);
        std::exit(EXIT_FAILURE);
    }

    if(method ==  "ML" && method_optim == "AI")
    {
        fmt::print(stderr, "Error: {} is not available for {}\n", method, method_optim);
        std::exit(EXIT_FAILURE);
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

    if(pheno_type == 1 && group.size() > 0)
    {
        std::cerr << "Error: 'group' option cannot be used with a binary phenotype ";
        exit(EXIT_FAILURE); 
    }
    
    GEMFit gf;
    std::ext::V_double cov_data = conv_dm2stdV(m_X); 
    std::ext::V_double pheno_data = conv_dv2stdVd(m_y);
    m_n_sel_col = cov_selected_hdrs_new.size();
    fit0(m_Nobs, m_n_sel_col, pheno_type, tol, m_robust, cov_selected_hdrs_new, pheno_data, cov_data,
                 &gf.XinvXTX, &gf.mu, &gf.resid, &gf.sigma2, gf.alpha, gf.eta); 


    std::cout << std::flush;
    std::cout << "****************************************************************************\n";
    std::cout << "Start Fitting the null model...\n \n";
    
    fit_null = gf.convert_2_fit(); 
    std::ext::V_string  sample_ids = pheno.m_data_frame.get_header(id);
    m_dup = pheno.m_data_frame.any_duplicated(pheno.m_sam_id);
    if(m_dup && m_vkins_sp[0].m_ratio_Nob2N >= m_vkins_sp[0].m_thr)
    {
        std::cout << "Duplicate IDs detected...\nAssuming longitudinal data with repeated measures...\n";
        fill_J(sample_ids);
        m_N = m_J.cols(); 

        if(rand_slope_hdr.size() > 0)
        {
            m_modeltype = LONGITUDINAL_RS;
        }
        else
        {
             m_modeltype = LONGITUDINAL_RI;
        }
    }
    else if(m_dup && m_vkins_sp[0].m_ratio_Nob2N < m_vkins_sp[0].m_thr)
    {
        std::cout << "Duplicate IDs detected...\nAssuming longitudinal data with repeated measures...\n";
        if(!m_vkins_sp[0].kin.m_null_kin) // If there is a kinship file add another matrix
        {
            SparseInverse spi;
            spi.set_spmat(details::diag(std::ext::V_double(m_Nobs, 1)));//Create a diagonal matrix
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
            auto v_indices = mapped_indices[dup]; //Map id to index
            for(auto& elm1 : v_indices) //Loop over duplicated ids to create values-->(id1,id1,1)
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

    std::ext::V_int group_id;
    if(group.size() == 0)
    {
        group_id = std::ext::V_int (m_Nobs, 1);
    }
    else
    {
        group_id = conv_stdvs2stdvi(pheno.m_data_frame.get_header(group));
    }
    
    glmmkin = glmmkin_fit(fit_null, group_id, method, method_optim, 
                          maxiter, tol, tau_min, tau_max, tau_region);
    glmmkin.id_include = pheno.m_data_frame.get_header(id);
	glmmkin.sigma2 = gf.sigma2;
    m_hdrsMap = createHeaderMap(cov_selected_hdrs_new);
    m_X = m_X_org;
    return glmmkin;    
}
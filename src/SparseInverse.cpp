#include "../include/SparseInverse.h"
#include <fstream> 
#include <algorithm>
#include <iterator>
#include <unordered_set>
#include <chrono>


SparseInverse::SparseInverse(std::string kin_add, std::string pheno_add, char kin_delim, double kin_diag,
                             char pheno_delim, std::string &m_sam_id, std::ext::V_string &m_v_hdrs, 
                             std::ext::V_string &bgen_sample_id, std::string missing_key) 
                             : kin_delim(kin_delim), pheno_delim(pheno_delim)
{
    pheno.m_sam_id = m_sam_id;
    pheno.m_v_hdrs = m_v_hdrs;
    pheno.m_data_frame.m_geno_ids = bgen_sample_id;
    pheno.m_data_frame.m_missing_key = missing_key;
    kin.m_data_frame.m_missing_key = missing_key;
    kin.m_diag = kin_diag;
    if (kin_add.size() > 0)
    {
        kin.set_path(kin_add);
    }
    else
    {
        kin.m_null_kin = true;
    }
    pheno.set_path(pheno_add);
    set_spmat();
}


void SparseInverse::set_idx_mp(std::ext::V_string v_strs)
{

    for(unsigned int i{0}; i < v_strs.size(); ++i)
    {
        //auto vst = "\"" + v_strs[i] + "\"";
        m_idx_mp[v_strs[i]].push_back(i+1);
    }
  
}


std::ext::IndexMap SparseInverse::get_idx_mp()
{
    return m_idx_mp;
}


SpaMat& SparseInverse::get_spmat()
{
    return m_spmat;
}


//Helper function to generate unique value
uint64_t combine_indices(long long a, long long b) {
    if (a > b) std::swap(a, b);  // Ensure a <= b for order independence
    return ((a + b) * (a + b + 1)) / 2 + b;
}


std::ext::VecTuples4spmat SparseInverse::create_tuple4spmat()
{
    std::ext::VecTuples4spmat vt4spmat;
    auto path = pheno.get_path();
    pheno.read_file(path, pheno_delim);
    //Match genofile sample IDs
    size_t unique_phenoIDs_beforematch = pheno.m_data_frame.size_wo_duplicates(pheno.m_sam_id);
    fmt::println("The number of unique IDs in Phenotype file before matching IDs is: {}", unique_phenoIDs_beforematch);
    pheno.m_data_frame.match_genoids(pheno.m_sam_id, pheno.m_v_hdrs);
    size_t unique_phenoIDs = pheno.m_data_frame.size_wo_duplicates(pheno.m_sam_id);
    fmt::println("The number of unique IDs in Phenotype file after matching IDs is: {}", unique_phenoIDs);
    fmt::println("****************************************************************************");
    //Map pheno sample ids to int to be used as matrix indices
    set_idx_mp(pheno.m_data_frame.m_data[pheno.m_sam_id]);
    if(kin.m_null_kin)
    {
        for(unsigned int i {0}; i < pheno.size(); ++i)
        {
            vt4spmat.emplace_back(std::ext::Triplet_d(i, i, 1.0)); 
        }
    }

    if(!kin.m_null_kin)
    {
        kin.read_file(kin.m_path, kin_delim);
        std::unordered_set<uint64_t> added_pairsdiag;
        std::unordered_set<uint64_t> added_pairs;

        if (pheno.m_data_frame.any_duplicated(pheno.m_sam_id))
        {
            // std::cout << "Duplicated id detected... \nAssuming longitudinal data with repeated measures...\n";
            auto duplicates = pheno.m_data_frame.list_duplicates(pheno.m_sam_id);
            for( auto dup : duplicates)
            {
                auto v_indices = m_idx_mp[dup]; //map id to index
            
                for(auto& id1 : v_indices) //loop over duplicated ids to create values-->(id1,id1,1)
                {
                    for(auto& id2 : v_indices)
                    {
                        uint64_t combined_key = combine_indices(id1, id2);
                        if(added_pairsdiag.insert(combined_key).second)
                        {
                            vt4spmat.emplace_back(std::ext::Triplet_d(id1 - 1 , id2 - 1, kin.m_diag));
                            if (id1 != id2) 
                            {
                                vt4spmat.emplace_back(std::ext::Triplet_d(id2 - 1, id1 - 1, kin.m_diag));
                            }
                        }
                    }
                } 
            }
        }
        else
        {
            // auto id0 = kin.m_data_frame.m_data[kin.m_data_frame.m_headers[0]];
            // auto id1 = kin.m_data_frame.m_data[kin.m_data_frame.m_headers[1]];
            // std::set<std::string> set_id0(id0.begin(), id0.end());
            // std::set<std::string> set_id1(id1.begin(), id1.end());
            // std::set<std::string> combined_ids(id0.begin(), id0.end());
            // combined_ids.insert(id1.begin(), id1.end());
            // std::unordered_set<int> mapped_ids;         
            // for(auto id : combined_ids)
            // {
            //     id.erase(std::remove(id.begin(), id.end(), '\"'), id.end());
            //     if (m_idx_mp.find(id) != m_idx_mp.end()) 
            //     {
            //         auto ids = m_idx_mp[id];
            //         for(auto& diag : ids)
            //         {
            //             mapped_ids.insert(diag);
            //         }
            //     }
            // }
            for(unsigned int i {0}; i < pheno.size(); ++i)
            {
                vt4spmat.emplace_back(std::ext::Triplet_d(i, i, kin.m_diag)); 
            }
            // if(mapped_ids.size() == 0)
            // {
            //     std::cerr << "Warning: There is no matching ID in kinship file and Phenotype file\n";
            //     std::exit(EXIT_FAILURE);
            // }
        }

        for(unsigned int i{0}; i < kin.size(); ++i)
        {
            auto val0 = kin.m_data_frame.m_data[kin.m_data_frame.m_headers[0]][i];//Retutn first ID in kin
            auto val1 = kin.m_data_frame.m_data[kin.m_data_frame.m_headers[1]][i];//Retutn second ID in kin

            val0.erase(std::remove(val0.begin(), val0.end(), '\"'), val0.end());
            val1.erase(std::remove(val1.begin(), val1.end(), '\"'), val1.end());
            // Map indxes to be meaningful in matrix if the sample exist in pheno it return corresponding index
            std::ext::V_int v_v0_idx;
            std::ext::V_int v_v1_idx;
            
            if (m_idx_mp.find(val0) != m_idx_mp.end() && m_idx_mp.find(val1) != m_idx_mp.end()) 
            {
                v_v0_idx = m_idx_mp[val0];
                v_v1_idx = m_idx_mp[val1];
            }

            if(v_v0_idx.size() > 0 && v_v1_idx.size() > 0)
            {
                for (int v0_idx : v_v0_idx) 
                {
                    for (int v1_idx : v_v1_idx) 
                    {
                        uint64_t combined_key = combine_indices(v0_idx, v1_idx);
                        //Check if pairs are not repeated and kinship's ids are in the phenodata
                        if (!is_missing(v0_idx, v1_idx) && added_pairs.insert(combined_key).second) 
                        {
                            double kinship_val = std::stod(kin.m_data_frame.m_data[kin.m_data_frame.m_headers[2]][i]);
                            vt4spmat.emplace_back(std::ext::Triplet_d(v0_idx - 1, v1_idx - 1, kinship_val));
                            
                            if (v0_idx != v1_idx) 
                            {
                                vt4spmat.emplace_back(std::ext::Triplet_d(v1_idx - 1, v0_idx - 1, kinship_val));
                            }
                        }
                    }
                }
            }
        }
    }
    return vt4spmat;
}

bool SparseInverse::is_missing(int id1, int id2)
{
    return !(id1 != 0 && id2 != 0);
}

void SparseInverse::set_spmat()
{
    auto vt = create_tuple4spmat();
    auto mat_size = pheno.size();
    m_spmat.resize(mat_size, mat_size);
    m_spmat.setFromTriplets(vt.begin(), vt.end());
    m_spmat.finalize();
    m_spmat.makeCompressed();
}

void SparseInverse::set_spmat(SpaMat sm)
{
    m_spmat = sm;
}

SpaMat SparseInverse::inv_spamat()
{
    Eigen::SimplicialLDLT<SpaMat, Eigen::Upper> solver;
    solver.compute(m_spmat);
    if (solver.info() != Eigen::Success) 
    { 
        std::cerr << "Decomposition failed!" << std::endl; 
        exit(EXIT_FAILURE);
    }
    std::cout << "Solver computed finished" << std::endl;
    SpaMat I(m_spmat.rows(), m_spmat.rows()); 
    I.setIdentity();
    return solver.solve(I);
}

// Function to convert an Eigen SparseMatrix to SuiteSparse cholmod_sparse
cholmod_sparse* convertEigenToSuiteSparse(SpaMat const& sm, cholmod_common *cm) 
{
    // Retrieve matrix properties
    int rows = sm.rows();
    int cols = sm.cols();
    int nnz = sm.nonZeros();  

    // Allocate memory for cholmod_sparse matrix
    cholmod_sparse* ssm = cholmod_allocate_sparse(
        rows,                  
        cols,                  
        nnz,                   
        1,                    
        1,                     
        1,                     
        CHOLMOD_REAL,          
        cm                     
    );

    // Check if allocation was successful
    if (ssm == nullptr) 
    {
        std::cerr << "Failed to allocate cholmod_sparse matrix." << std::endl;
        return nullptr;
    }

    // Copy data from Eigen matrix to cholmod_sparse
    std::memcpy(ssm->p, sm.outerIndexPtr(), (cols + 1) * sizeof(int));
    std::memcpy(ssm->i, sm.innerIndexPtr(), nnz * sizeof(int));
    std::memcpy(ssm->x, sm.valuePtr(), nnz * sizeof(double));

    // Set other fields 
    ssm->nzmax = nnz;            
    ssm->nrow = rows;            
    ssm->ncol = cols;           
    // ssm->xtype = CHOLMOD_REAL;   
    return ssm;
}

// Convert cholmod to Csparse
cs_di* cholmodSparseToCsparse(cholmod_sparse* L) 
{
    int n_row = L->nrow;  // Number of rows// square matrix
    int nzmax = L->nzmax;  // Maximum number of non-zero elements
    cs_di* csMatrix = cs_di_spalloc(n_row, n_row, nzmax, 1, 0);  // Allocate CSparse matrix
    if (!csMatrix) 
    {
        std::cerr << "Memory allocation for CSparse matrix failed.\n";
        return nullptr;
    }

    int* Lp = static_cast<int*>(L->p);
    int* Li = static_cast<int*>(L->i);
    double* Lx = static_cast<double*>(L->x);

    // Copy column pointers
    // p array has size n_col+1 or n_row+1, square matrix
    for (int i = 0; i <= n_row; ++i) 
    {
        csMatrix->p[i] = Lp[i];
        if (Lp[i] < 0 || Lp[i] > nzmax)
        {
            std::cerr << "Invalid column pointer in CSparse matrix: Lp[" << i << "] = " << Lp[i] << "\n";
        }
    }

    // Copy row indices and values
    for (int i = 0; i < nzmax; ++i) 
    {
        csMatrix->i[i] = Li[i];
        csMatrix->x[i] = Lx[i];
        if (Li[i] < 0 || Li[i] >= n_row) 
        {
            std::cerr << "Invalid row index in CSparse matrix: Li[" << i << "] = " << Li[i] << "\n";
        }
    }

    return csMatrix;
}


cs_di* createSparseIdentity(int m) 
{
    // Allocate space for an identity matrix in sparse format
    cs_di* csMatrix = cs_di_spalloc(m, m, m, 1, 0); // m non-zero entries, since it's an identity matrix
    if (!csMatrix) 
    {
        std::cerr << "Memory allocation for sparse identity matrix failed.\n";
        return nullptr;
    }

    // Set up the identity matrix in compressed column storage (CCS) format
    for (int i = 0; i < m; ++i) {
        csMatrix->p[i] = i;      // Column pointers: each column starts at index i
        csMatrix->i[i] = i;      // Row indices: diagonal element at row i
        csMatrix->x[i] = 1.0;    // Value: all diagonal elements are 1
    }
    csMatrix->p[m] = m;       // p array has m + 1, the last column pointer equals the total number of non-zero elements

    return csMatrix;
}


SpaMat SparseInverse::inv_spamat(SpaMat const& sm) 
{
    cholmod_common cm;
    cholmod_start(&cm);
    cm.supernodal = CHOLMOD_SIMPLICIAL;  // Use simplicial factorization
    // cm.supernodal = CHOLMOD_SUPERNODAL;
    cm.method[0].ordering = CHOLMOD_NATURAL; // Use natural ordering (no permutation)
    cm.postorder = 0;
    cm.nmethods = 1;                 
    cm.final_ll = 1;     //if true, simplicial factors are converted to LL                
    cm.final_pack = 1;
    // cm.final_asis = 0;
    // cm.final_monotonic = 1;
    // Convert Eigen sparse matrix to CHOLMOD format
    cholmod_sparse* spMatrix = convertEigenToSuiteSparse(sm, &cm); 
    if (spMatrix == nullptr) {
        std::cerr << "Conversion to SuiteSparse format failed.";
        cholmod_finish(&cm);
        throw std::runtime_error("Conversion to SuiteSparse format failed.\n");
    }

    cholmod_factor* factor = cholmod_analyze(spMatrix, &cm);
    if (factor == nullptr) 
    {
        std::cerr << "CHOLMOD analysis failed.\n";
        cholmod_free_factor(&factor, &cm);// FIXEME do we need it?
        cholmod_free_sparse(&spMatrix, &cm);
        cholmod_finish(&cm);
        throw std::runtime_error("CHOLMOD analysis failed.\n");
    }

    cm.final_asis = 0;
	// cm.final_super = 0;
	cm.final_monotonic = 1;
    int status = cholmod_factorize(spMatrix, factor, &cm);
	
   if (!status) 
   {
        std::cerr << "CHOLMOD factorization failed.\n";
        cholmod_free_factor(&factor, &cm);
        cholmod_free_sparse(&spMatrix, &cm);
        cholmod_finish(&cm);
        throw std::runtime_error("CHOLMOD factorization failed.\n");
    }

    // Convert CHOLMOD factor to CSparse format
    cholmod_sparse* L = cholmod_factor_to_sparse(factor, &cm);//Suitesparse func
    cs_di* csFactor = cholmodSparseToCsparse(L);
    cholmod_free_sparse(&L, &cm);
 
    if (csFactor == nullptr) {
        std::cerr << "Conversion to CSparse format failed.";
        cholmod_free_factor(&factor, &cm);
        cholmod_free_sparse(&spMatrix, &cm);
        cholmod_finish(&cm);
        throw std::runtime_error("Conversion to CSparse format failed.\n");
    }

    int n_cols = sm.cols();
    SpaMat eigenInverse(n_cols, n_cols);
    std::vector<Eigen::Triplet<double>> tripletList;
    // Workspace for cs_di_spsolve
    int* xi = (int*)malloc(2 * n_cols * sizeof(int));
    double* x = (double*)malloc(n_cols * sizeof(double));
    if (!xi || !x) {
        std::cerr << "Memory allocation failed.";
        cs_di_spfree(csFactor);
        cholmod_free_factor(&factor, &cm);
        cholmod_free_sparse(&spMatrix, &cm);
        free(xi);
        free(x);
        throw std::runtime_error("Memory allocation failed.\n");
    }
    cs_di* csb = createSparseIdentity(n_cols);//or nnz
    // Solve the system for each column 
    auto start1 = std::chrono::high_resolution_clock::now();
    for (int col = 0; col < n_cols; ++col) 
    {
        int top = cs_di_spsolve(csFactor, csb, col, xi, x, nullptr, 1); // 0 indicates upper triangular csFactor should be llt or supernudal
        if (top == -1) 
        {
            std::cerr << "CSparse solve failed for column " << col << ".";
            cs_di_spfree(csFactor); 
            cs_di_spfree(csb);
            cholmod_free_factor(&factor, &cm);
            cholmod_free_sparse(&spMatrix, &cm);
            free(xi);
            free(x);
            cholmod_finish(&cm);
            throw std::runtime_error("CSparse solve failed.\n");
        }
        // Store the result in the triplet list
        for (int i = top; i < n_cols; ++i)
        {
            tripletList.emplace_back(xi[i], col, x[xi[i]]); 
        }
    } 
        
    // Populate the Eigen sparse inverse matrix
    eigenInverse.setFromTriplets(tripletList.begin(), tripletList.end());
    SpaMat eigenInverseFinal = eigenInverse.transpose() * eigenInverse; //To calculate upper and lower
    cs_di_spfree(csb);
    cs_di_spfree(csFactor);
    cholmod_free_factor(&factor, &cm);
    cholmod_free_sparse(&spMatrix, &cm);
    free(xi);
    free(x);
    // Finish CHOLMOD
    cholmod_finish(&cm);
    return eigenInverseFinal;
}


DensMat SparseInverse::inv(DensMat const& dm)
{    
    Eigen:: LDLT<DensMat> ldlt;
    
    ldlt.compute(dm);
    
    if (ldlt.info() != Eigen::Success)
    {
        std::cerr << "Error on LDLT!" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    DensMat I(dm.rows(), dm.rows()); 
    
    I.setIdentity();
    
    return ldlt.solve(I);
}
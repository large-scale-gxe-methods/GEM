#pragma once
// #define ARMA_USE_LAPACK
// #define ARMA_USE_BLAS
// #define DARMA_DONT_USE_WRAPPER
#include <armadillo>
// #define EIGEN_USE_MKL_ALL
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <cholmod.h>
#include <cs.h> 
#include "Pheno.h"
#include "Kinship.h"
#include <tuple>
#include <iostream>


namespace std
{
    namespace ext
    {
        using IndexMap = unordered_map<string, V_int>;
        using IndexMapRev = unordered_map<int, string>;
        using Triplet_d = Eigen::Triplet<double>;
        using Triplet_i = Eigen::Triplet<int>;
        using VecTuples4spmat = vector<Triplet_d>;
        using VecTriple_i = vector<Triplet_i>;


        // a template helper function and functions for printing tuples
        template<typename Tuple, size_t... Indices>
        void tuplePrint(const Tuple& t, std::index_sequence<Indices...>)
        {
            ((cout << std::get<Indices>(t) << " "), ...);
        }

        template<typename... Args>
        void tuplePrint(const std::tuple<Args...>& t)
        {
            tuplePrint(t, std::index_sequence_for<Args...>());
            cout << "\n";
        }
        
    }
}

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (std::pair<T1,T2> const& pair) const 
    {
        auto h1 = std::hash<T1>{}(pair.first);
        auto h2 = std::hash<T2>{}(pair.second);
        return h1 ^ h2; 
    }
};

using SpaMat = Eigen::SparseMatrix<double, Eigen::ColMajor>;
// using SpaMat = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;
using Mat = Eigen::MatrixXd;
using DensMat = Eigen::MatrixXd;
using DensVec = Eigen::VectorXd;
using DensVecInt = Eigen::VectorXi;
using DenseMatInt = Eigen::MatrixXi;


class SparseInverse
{
    public:
        Pheno pheno;
        Kinship kin;
        char pheno_delim;
        char kin_delim; 

        SparseInverse() =  default;
        SparseInverse(std::string kin_add, std::string pheno_add, char kin_delim,
                      double kin_diag, char pheno_delim, std::string &m_sam_id,
                      std::ext::V_string &m_v_hdrs, std::ext::V_string &bgen_sample_id,
                      std::string phenoMissingKey) ;
        void set_idx_mp(std::ext::V_string v_strs);
        std::ext::IndexMap get_idx_mp();        
        void set_spmat();
        void set_spmat(SpaMat sm);
        SpaMat& get_spmat();
        SpaMat inv_spamat();
        //define static as we might don't create an object
        static SpaMat inv_spamat(SpaMat const& sm);
        static DensMat inv(DensMat const &dm);

    private:
       std::ext::IndexMap m_idx_mp;
       SpaMat m_spmat;
    //    std::vector<SpaMat> m_vspmat;
       //std::ext::V_string get_pheno_samid();
       [[nodiscard]] std::ext::VecTuples4spmat create_tuple4spmat();
       [[nodiscard]] bool is_missing(int id1, int id2);
};
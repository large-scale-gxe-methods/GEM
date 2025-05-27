#include "MAGEE.h"
#include "declars.h"
#include "ReadBGEN.h"
#include "../thirdparty/zstd-1.5.5/lib/zstd.h"
#include "../thirdparty/libdeflate-1.18/libdeflate.h"
#include "../thirdparty/plink-2.0/plink2_bits.h"
#include "../thirdparty/plink-2.0/plink2_base.h"
#include "../thirdparty/plink-2.0/pgenlib_misc.h"
#include "../thirdparty/plink-2.0/pgenlib_read.h"
#include "../thirdparty/plink-2.0/pgenlib_ffi_support.h"
#include <stdexcept>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <zlib.h>
#include <memory>
#include <stdint.h>
#include <limits.h>

#ifndef DBL_EPSILON
#define DBL_EPSILON 2.2204460492503131e-16
#endif

using uint = unsigned int;
using uchar = unsigned char;
using ushort = unsigned short;
using llui = long long unsigned int;


template <typename Container, typename T>
std::vector<int> find_indx(const Container &container, T value)
{
    std::vector<int> ret;
    for (auto i = 0; i < container.size(); ++i)
    {
        if (container[i] == value)
        {
            ret.push_back(i);
        }
    }
    return ret;
}

template <typename Pred>
std::vector<int> filter_elements(const arma::vec &container, const arma::Col<int> &condition, Pred pred)
{
    std::vector<int> ret;

    for (int i = 0; i < condition.size() && i < container.size(); ++i)
    {
        if (pred(condition(i)))
        {
            ret.push_back(container(i));
        }
    }
    return ret;
}

double chi_square_CDF(double x, double df, bool lower_tail, bool log_p)
{
    // Check if x is NaN or less than or equal to zero
    if (std::isnan(x) || x <= 0.0)
    {
        // Handle the invalid input as needed:
        return std::numeric_limits<double>::quiet_NaN();
    }
    // Proceed with the chi-square CDF calculation
    boost::math::chi_squared dist(df);
    double cdf_value;

    try
    {
        if (lower_tail)
        {
            cdf_value = boost::math::cdf(dist, x);
        }
        else
        {
            cdf_value = boost::math::cdf(boost::math::complement(dist, x));
        }

        if (log_p)
        {
            return std::log(cdf_value);
        }
        else
        {
            return cdf_value;
        }
    }
    catch (const std::domain_error& e)
    {
        // Handle any domain errors thrown by Boost
        std::cerr << "Domain error in chi square CDF: " << e.what() << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
}

void glmm_gei(std::string snpID, arma::mat &G, arma::uvec &snp_skip, size_t &npbidx, size_t npb, size_t n, 
                        int ei, int qi, const Magee_Arma &null_obj, std::ofstream &writefile, 
                        bool meta_output, std::ext::V_string &tmpout, uint m, uint end) 
{
	if ((m == end) || (npbidx == npb))
		{
			if (npbidx != npb)
			{
				G.reshape(n, npbidx);
				snp_skip = snp_skip.rows(0,npbidx-1);
			}
			arma::uvec snp_idx = find(snp_skip == 0);
			G = G.cols(snp_idx);
			int ng = G.n_cols;
			arma::mat IV_U;
			arma::mat IV_U1;
			arma::mat IV_V_i;
			arma::mat IV_V_i1;
			arma::mat IV_E_i;
			arma::mat IV_GE_i;
			arma::mat STAT_JOINT_tmp;
			arma::vec BETA_MAIN;
			arma::vec STAT_INT;
			arma::vec STAT_JOINT;
			arma::vec SE_MAIN;
			arma::mat BETA_INT;
			arma::mat BETA_INT1;
			arma::vec PVAL_MAIN(ng);
			PVAL_MAIN.fill(arma::datum::nan);
			arma::vec PVAL_INT(ng);
			PVAL_INT.fill(arma::datum::nan);
			arma::vec PVAL_JOINT(ng);
			PVAL_JOINT.fill(arma::datum::nan);
			size_t ngei1 = ng * (ei + 1);
			if (G.n_cols != 0)
			{
				arma::sp_mat Gsp(G);
				arma::sp_mat PG;
				arma::vec U;
				arma::mat GPG;
				
				U = G.t() * null_obj.Jres;// 1*1 if G be a vec
				arma::sp_mat Gsigma_ixJ =  Gsp.t() * null_obj.sigma_ixJ;
				PG = (null_obj.sigma_iJJ.t() * Gsp) - (null_obj.sigma_ixJ * (Gsigma_ixJ * null_obj.cov.t()).t());//P projection matrix
				GPG = (G.t() * PG) % kron(arma::ones<arma::mat>(1, 1), arma::mat(ng, ng, arma::fill::eye));         
				arma::mat GPG_i;
				bool is_non_singular = inv(GPG_i, GPG);
				if (!is_non_singular) 
				{
					GPG_i = pinv(GPG);
				}
				arma::mat V_i;
				V_i = diagvec(GPG_i);           
				arma::vec V_MAIN_adj = diagvec(GPG_i);
				V_MAIN_adj = V_MAIN_adj.rows(0, ng-1);            
				arma::vec BETA_MAIN_adj = GPG_i.t() * U;
				BETA_MAIN_adj= BETA_MAIN_adj.rows(0, ng-1);
				arma::vec STAT_MAIN_adj(ng);
				STAT_MAIN_adj.fill(arma::datum::nan);

				for (size_t s = 0; s < V_MAIN_adj.size(); s++) 
				{
					if (V_MAIN_adj[s] > 0) 
					{
						STAT_MAIN_adj[s] = (BETA_MAIN_adj[s] * BETA_MAIN_adj[s]) / V_MAIN_adj[s];
					}
				} 

				BETA_MAIN = V_i % U.rows(0,ng-1);
				SE_MAIN = sqrt(V_i);
				arma::vec STAT_MAIN = BETA_MAIN % U.rows(0,ng-1);

				for (size_t s = 0; s < STAT_MAIN.size(); s++) 
				{
					if (STAT_MAIN[s] > 0) 
					{
						PVAL_MAIN[s] = chi_square_CDF(STAT_MAIN[s], 1, 0, 0);
					}
				}
				
				arma::mat Hv(ei+qi+1, ei+qi+1, arma::fill::zeros);
				arma::mat Gtblock = kron(arma::mat(ei+qi+1, ei+qi+1, arma::fill::eye), G.t());
				// arma::sp_mat Gtblock = kron(arma::speye<arma::sp_mat>(ei+qi+1, ei+qi+1), arma::sp_mat(G).t());
				arma::mat GblockXi = Gtblock * null_obj.Xi;
				// arma::sp_mat GblockXi = Gtblock * null_obj.Xi;
				Hv = Gtblock * null_obj.Psi * Gtblock.t() - (GblockXi * null_obj.cov.t() * GblockXi.t());
				Hv = Hv % arma::kron(arma::ones(ei+qi+1, ei+qi+1), arma::mat(ng, ng, arma::fill::eye));
				bool is_non_singular_Hv = inv(IV_V_i1, Hv);
				if (!is_non_singular_Hv) 
				{
					IV_V_i1 = arma::pinv(Hv);
				}

				arma::mat cross_Eres_G;
				cross_Eres_G = kron(arma::mat(ei+qi+1, ei+qi+1, arma::fill::eye), G.t()) * null_obj.JEresblock;//G.t() * cross_K1_J;
				arma::mat IV_U_kron1 = arma::kron(arma::ones<arma::mat>(ei + qi + 1, 1), arma::mat(ng, ng, arma::fill::eye));
				IV_U1 = IV_U_kron1.each_col() % cross_Eres_G;// V_U_kron1.each_col() % cross_K1_J_G //cross_K_res;
				BETA_INT1 = (IV_V_i1 * IV_U1);
				bool is_non_singular_IV_V1 = arma::inv(IV_E_i,IV_V_i1(arma::span(ng,ngei1-1), arma::span(ng,ngei1-1)));
				if (!is_non_singular_IV_V1) 
				{
					IV_E_i = arma::pinv(IV_V_i1(arma::span(ng,ngei1-1), arma::span(ng,ngei1-1))); 
				}

				IV_U = IV_E_i* BETA_INT1.rows(ng,ngei1-1);
				STAT_INT = diagvec(IV_U.t()*BETA_INT1.rows(ng,ngei1-1));

				try 
				{
					IV_GE_i = arma::inv(IV_V_i1(arma::span(0,ngei1-1), arma::span(0,ngei1-1))); 
					STAT_JOINT_tmp=  IV_GE_i*BETA_INT1.rows(0,ngei1-1);
					STAT_JOINT = diagvec(STAT_JOINT_tmp.t()*BETA_INT1.rows(0,ngei1-1));
				
					for (size_t s = 0; s < STAT_INT.size(); s++) 
					{
						// PVAL_INT[s] = chi_square_CDF(STAT_INT[s], ei, 0, 0);
						if (arma::is_finite(PVAL_MAIN[s])) 
						{
							PVAL_JOINT[s] = chi_square_CDF(STAT_JOINT[s], 1+ei, 0, 0);
						}
					} 
				} 

				catch (std::runtime_error const& error) 
				{
					std::cout << "Warning: A singular matrix was observed for snpID: "<< snpID << "\n";
					for (size_t s = 0; s < STAT_INT.size(); s++) 
					{
						PVAL_JOINT[s] = DBL_EPSILON;
					}
				}

				for (size_t s = 0; s < STAT_INT.size(); s++)
				{
					PVAL_INT[s] = chi_square_CDF(STAT_INT[s], ei, 0, 0);
				}
			}

			arma::uvec b_idx1 = arma::regspace<arma::uvec>(0, ng, (ei+qi) * ng);
			int ng_j = 0;
			//Write in the file

			for (size_t j = 0; j < npbidx; ++j)
			{
				if (snp_skip[j] == 1)
				{ 
					continue;
				}
				else
				{
					writefile << tmpout[j] <<  BETA_MAIN[ng_j] << "\t" << SE_MAIN[ng_j] << "\t";
					
					if (meta_output)
					{
						// Beta Int the diaganol of BETA_INT1
						for (int b = 0; b < ei + qi + 1; b++)
						{
							int row = b_idx1[b] + ng_j;
							writefile << BETA_INT1(row, ng_j) << "\t";
						}
						// Var (the diaganol of IV_V_i1)
						for (int b = 0; b < ei + qi + 1; b++)
						{
							int col = b_idx1[b] + ng_j;
							for (int d = 0; d < ei + qi + 1; d++)
							{
								if (b == d)
								{
									int row = b_idx1[d] + ng_j;
									writefile << std::sqrt(IV_V_i1(row, col)) << "\t";
								}
							}
						}

						// Cov (the lower triangle elements of IV_V_i1)
						for (int b = 0; b < ei + qi + 1; b++)
						{
							int col = b_idx1[b] + ng_j;
							for (int d = 0; d < ei + qi + 1; d++)
							{
								if (d > b)
								{
									int row = b_idx1[d] + ng_j;
									writefile << IV_V_i1(row, col) << "\t";
								}
							}
						}
					}
					else
					{
						int ncolE = ei + qi + 1;
						arma::mat split_mat(ncolE,ncolE);

						for (int i=0; i<ncolE; i++) 
						{
							for (int j=0; j<ncolE; j++)
							split_mat(i,j) = i+1-ncolE+ncolE*(j+1);
						}
						
						split_mat = split_mat(arma::span(1,ei), arma::span(1,ei));

						if (split_mat.size() == 1)
						{
							for (int b = 0; b < ei + 1; b++)
							{
								int row = b_idx1[b] + ng_j;
								// NOT the first ng row

								if (row > ng - 1)
								{

									writefile << BETA_INT1(row, ng_j) << "\t";
								}
							}
							// Var (the diaganol elements)
							for (int b = 0; b < ei + 1; b++)
							{
								int col = b_idx1[b] + ng_j;
								for (int d = 0; d < ei + 1; d++)
								{
									if (b == d)
									{
										int row = b_idx1[d] + ng_j;
										// NOT the first ng row or first ng col
										if (row > ng - 1 && col > ng - 1)
										{
											writefile << std::sqrt(IV_V_i1(row, col)) << "\t";
										}
									}
								}
							}
						}

						else
						{
							for (int b = 0; b < ei + 1; b++)
							{
								int row = b_idx1[b] + ng_j;
								// NOT the first ng row
								if (row > ng - 1)
								{
									writefile << BETA_INT1(row, ng_j) << "\t";
								}
							}

							for (int b = 0; b < ei + 1; b++)
							{
								int col = b_idx1[b] + ng_j;
								for (int d = 0; d < ei + 1; d++)
								{
									if (b == d)
									{
										int row = b_idx1[d] + ng_j;
										// NOT the first ng row or first ng col
										if (row > ng - 1 && col > ng - 1)
										{
											writefile << std::sqrt(IV_V_i1(row, col)) << "\t";
										}
									}
								}
							}

							for (int b = 0; b < ei + 1; b++)
							{
								int col = b_idx1[b] + ng_j;
								for (int d = 0; d < ei + 1; d++)
								{
									if (d > b)
									{
										int row = b_idx1[d] + ng_j;
										if (row > ng - 1 && col > ng - 1)
										{
											writefile << IV_V_i1(row, col) << "\t";
										}
									}
								}
							}
						}
					}
					writefile << PVAL_MAIN[ng_j] << "\t" << PVAL_INT[ng_j] << "\t" << PVAL_JOINT[ng_j] << "\n"; 
					ng_j++;
				}
			}

			npbidx = 0;
			snp_skip.zeros();
			G.reshape(n, npb);
		}

		if ((m) % 10000 == 0)
		{
			writefile << std::flush;
		} 

}

void glmm_gei_bgen13(Magee_Arma const& null_obj, string const &bgenfile,
					string const &outfile, double minmaf, double missrate,
					size_t npb, int ei, int qi, std::ext::Map_str_Vint const &strata_list, 
					uint begin, uint end, long long unsigned int byte, uint Nbgen,
					uint compression, bool meta_output)
{
	// bool isDupeID = null_obj.dupflag;
    std::ext::V_int select = null_obj.select;
    int strataList_size = strata_list.size();
    bool skip_strata = strata_list.empty();
    double maxmaf = 1 - minmaf;

	std::string output = outfile + "_bin_" + std::to_string(begin) + ".tmp";
	std::ofstream writefile(output, std::ios::binary);
    if (!writefile) {
        std::cerr << "Unable to open file for appending: " << outfile << std::endl;
        std::exit(EXIT_FAILURE);
    }

	size_t n = null_obj.n;
	size_t n_obs = null_obj.n_obs;

	arma::mat G(n, npb);
	arma::vec g(n);
    arma::uvec gmiss(n);
	arma::uvec snp_skip = arma::zeros<arma::uvec>(npb);
	// arma::mat G(n_obs, npb);
	std::ext::V_string tmpout(npb);

	double gmean, gsqmean, geno, gmax, gmin;
	size_t ncount, nmiss, npbidx = 0;
	
	struct libdeflate_decompressor *decompressor = libdeflate_alloc_decompressor();

	uint maxLA = 65536;
	std::vector<uchar> zBuf12;
	std::vector<uchar> shortBuf12;

	std::vector<char> snpID(maxLA + 1);
	std::vector<char> rsID(maxLA + 1);
	std::vector<char> chrStr(maxLA + 1);
	std::vector<char> allele1(maxLA + 1);
	std::vector<char> allele0(maxLA + 1);

	arma::vec g2(Nbgen);
	g2.fill(arma::datum::nan);
	
	FILE *fp = fopen(bgenfile.c_str(), "rb"); 

	fseek(fp, byte, SEEK_SET); 
	int ret;

	arma::mat strata_AF(end +1, strataList_size);
	arma::mat strata_Var(end +1, strataList_size);
	arma::mat strata_N(end +1, strataList_size);

	for (uint m = begin; m <= end; ++m)
	{
		std::ostringstream writeout;
		ushort LS;
		ret = fread(&LS, 2, 1, fp);
		snpID.resize(LS);
		ret = fread(&snpID[0], 1, LS, fp);
		std::string str_snpID = "NA";
		if (LS != 0)
		{
			str_snpID = std::string(&snpID[0]);
		}
		ushort LR;
		ret = fread(&LR, 2, 1, fp);
		rsID.resize(LR);
		ret = fread(&rsID[0], 1, LR, fp);

		ushort LC;
		ret = fread(&LC, 2, 1, fp);
		chrStr.resize(LC);
		ret = fread(&chrStr[0], 1, LC, fp);

		uint physpos;
		ret = fread(&physpos, 4, 1, fp);
		std::string physpos_tmp = std::to_string(physpos);

		ushort LKnum;
		ret = fread(&LKnum, 2, 1, fp);
		
		if (LKnum != 2)
		{
			std::cerr << "\nERROR: " << str_snpID << " is a non-bi-allelic variant with " << LKnum << " alleles. Please filter these variants for now.\n\n";
			std::exit(EXIT_FAILURE);
		}

		uint LA;
		ret = fread(&LA, 4, 1, fp);
		allele1.resize(LA);
		ret = fread(&allele1[0], 1, LA, fp);

		uint LB;
		ret = fread(&LB, 4, 1, fp);
		allele0.resize(LB);
		ret = fread(&allele0[0], 1, LB, fp);

		uint cLen;
		ret = fread(&cLen, 4, 1, fp);
		uchar *bufAt;
		
		if (compression == 1)
		{
			zBuf12.resize(cLen - 4);
			uint dLen;
			ret = fread(&dLen, 4, 1, fp);
			ret = fread(&zBuf12[0], 1, cLen - 4, fp);
			shortBuf12.resize(dLen);
			uLongf destLen = dLen;
			// The function returns an enum value of type libdeflate_result, which indicates the result of the decompression operation.
			if (libdeflate_zlib_decompress(decompressor, &zBuf12[0], cLen - 4, &shortBuf12[0], destLen, NULL) != LIBDEFLATE_SUCCESS)
			{
				std::cerr << "\nERROR: Decompressing " << str_snpID << " block failed with libdeflate.\n\n";
				//throw std::runtime_error("Decompressing " + std::string(rsID.get()) + " genotype block failed with libdeflate.");
				std::exit(EXIT_FAILURE);
			}
			bufAt = &shortBuf12[0];
		}
		else if (compression == 2)
		{
			zBuf12.resize(cLen - 4);
			uint dLen;
			ret = fread(&dLen, 4, 1, fp);
			ret = fread(&zBuf12[0], 1, cLen - 4, fp);
			shortBuf12.resize(dLen);

			uLongf destLen = dLen;
			size_t ret = ZSTD_decompress(&shortBuf12[0], destLen, &zBuf12[0], cLen - 4);
			if (ret > destLen)
			{
				if (ZSTD_isError(ret))
				{
					std::cerr << "Error reading bgen file: Decompressing genotype block failed. \n" << "ZSTD ERROR: " << ZSTD_getErrorName(ret);
					std::exit(EXIT_FAILURE);
				}
			}
			bufAt = &shortBuf12[0];
		}
		else
		{
			zBuf12.resize(cLen);
			ret = fread(&zBuf12[0], 1, cLen, fp);
			bufAt = &zBuf12[0];
		}
		
		uint N;
		std::memcpy(&N, bufAt, sizeof(int32_t));

		if (N != Nbgen)
		{
			std::cerr << "\nERROR: " << str_snpID << " number of samples (" << N << ") with genotype probabilties does not match number of samples specified in BGEN file (" << Nbgen << ").\n\n";
			std::exit(EXIT_FAILURE);
		}
		ushort K;
		std::memcpy(&K, &(bufAt[4]), sizeof(int16_t));
		if (K != 2)
		{
			std::cerr << "\nERROR: There are SNP(s) with more than 2 alleles (non-bi-allelic). Currently unsupported. \n\n";
			std::exit(EXIT_FAILURE);
		}
		
		const uint min_ploidy = bufAt[6];
		if (min_ploidy != 2)
		{
			std::cerr << "\nERROR: " << str_snpID << " has minimum ploidy " << min_ploidy << ". Currently unsupported. \n\n";
			std::exit(EXIT_FAILURE);
		}
		const uint max_ploidy = bufAt[7];
		if (max_ploidy != 2)
		{
			std::cerr << "\nERROR: " << str_snpID << " has minimum ploidy " << max_ploidy << ". Currently unsupported. \n\n";
			std::exit(EXIT_FAILURE);
		}
		const unsigned char *missing_and_ploidy_info = &(bufAt[8]);
		const unsigned char *probs_start = &(bufAt[10 + N]);
		const uint is_phased = probs_start[-2]; 
		if (is_phased != 1 && is_phased != 0)	
		{
			std::cerr << "\nERROR: " << str_snpID << " has phased value of " << is_phased << ". This must be 0 or 1. \n\n";
			std::exit(EXIT_FAILURE);
		}

		const uint B = probs_start[-1];

		if (B != 8 && B != 16 && B != 24 && B != 32)
		{
			std::cerr << "Error reading bgen file: Bits to store probabilities must be 8, 16, 24, or 32. \n";
			std::exit(EXIT_FAILURE);
		}

		const uintptr_t numer_mask = (1U << B) - 1;
		const uintptr_t probs_offset = B / 8;

		gmean = 0.0;
		//double mac = 0.0;
		gsqmean = 0.0;
		gmax = -100.0;
		gmin = 100.0;
		//double rsq = 0.0;
		nmiss = 0;
		ncount = 0;

		if (!is_phased)
		{
			for (size_t i = 0; i < N; i++)
			{
				const uint missing_and_ploidy = missing_and_ploidy_info[i];
				uintptr_t numer_aa;
				uintptr_t numer_ab;

				if (missing_and_ploidy == 2)
				{
					Bgen13GetTwoVals(probs_start, B, probs_offset, &numer_aa, &numer_ab);
					probs_start += (probs_offset * 2);
				}
				else if (missing_and_ploidy == 130)
				{
					probs_start += (probs_offset * 2);
					if (select[ncount] >= 0)
					{
						gmiss(select[ncount]) = 1;
						nmiss++;
					}
					ncount++;
					continue;
				}
				else
				{
					std::cerr << "Error reading bgen file: Ploidy value " << missing_and_ploidy << " is unsupported. Must be 2 or 130. \n";
					std::exit(EXIT_FAILURE);
				}
				
				if (select[ncount] >= 0)
				{
					double p11 = numer_aa / double(1.0 * numer_mask);
					double p10 = numer_ab / double(1.0 * numer_mask);
					// geno = 2 * p11 + p10;
					geno = 2 * (1 - p11 - p10) + p10;
					gmiss(select[ncount]) = 0;
					g[select[ncount]] = geno;
					g2[select[ncount]] = geno * geno;
					
					gmean += geno;
					gsqmean += geno * geno;

					if (geno > gmax)
					{
						gmax = geno;
					}
					if (geno < gmin)
					{
						gmin = geno;
					}
				}
				ncount++;
			} 
		} //end of if is_phased
		else
		{
			for (size_t i = 0; i < N; i++)
			{
				const uint missing_and_ploidy = missing_and_ploidy_info[i];
				uintptr_t numer_aa;
				uintptr_t numer_ab;

				if (missing_and_ploidy == 2)
				{
					Bgen13GetTwoVals(probs_start, B, probs_offset, &numer_aa, &numer_ab);
					probs_start += (probs_offset * 2);
				}
				else if (missing_and_ploidy == 130)
				{
					probs_start += (probs_offset * 2);
					if (select[ncount] >= 0)
					{
						gmiss(select[ncount]) = 1;
						nmiss++;
					}
					ncount++;
					continue;
				}
				else
				{
					std::cerr << "Error reading bgen file: Ploidy value " << missing_and_ploidy << " is unsupported. Must be 2 or 130. \n";
					// return R_NilValue;
					exit(EXIT_FAILURE);
				}

				if (select[ncount] >= 0)
				{
					double p11 = numer_aa / double(1.0 * numer_mask);
					double p10 = numer_ab / double(1.0 * numer_mask);
					geno = double(1.0 * missing_and_ploidy) - (p11 + p10); // help to calculate mac
					gmiss(select[ncount]) = 0;
					g[select[ncount]] = geno;
					g2[select[ncount]] = geno * geno;
					gmean += geno;
					gsqmean += geno * geno;

					if (geno > gmax)
					{
						gmax = geno;
					}
					if (geno < gmin)
					{
						gmin = geno;
					}
				}
				ncount++;
			}
		} //  end of else isphased
		// mac = gmean;
		// if (static_cast<double>(n - nmiss) < mac)
		// {
		// 	mac = static_cast<double>(n - nmiss) * 2.0 - mac;
		// }
		gmean /= static_cast<double>(n - nmiss);
		gsqmean /= static_cast<double>(n - nmiss);
		//rsq = (gsqmean - gmean * gmean) * static_cast<double>(n - nmiss) / static_cast<double>(n - nmiss - 1) / (gmean * (1.0 - gmean / 2.0));
		double var = (gsqmean - gmean * gmean) * static_cast<double>(n - nmiss) / static_cast<double>(n - nmiss - 1);
		
		if (skip_strata)
		{
			writeout << std::string (snpID.begin(), snpID.end()) << "\t" << std::string (rsID.begin(), rsID.end()) << "\t" << std::string (chrStr.begin(), chrStr.end()) << "\t" << physpos_tmp << "\t" << std::string (allele1.begin(), allele1.end()) << "\t" << std::string (allele0.begin(), allele0.end()) << "\t" << (n - nmiss) << "\t" << gmean / 2.0 << "\t" << var << "\t";
		}
		else
		{
			writeout << std::string (snpID.begin(), snpID.end()) << "\t" << std::string (rsID.begin(), rsID.end()) << "\t" << std::string (chrStr.begin(), chrStr.end()) << "\t" << physpos_tmp << "\t" << std::string (allele1.begin(), allele1.end()) << "\t" << std::string (allele0.begin(), allele0.end()) << "\t" << (n - nmiss) << "\t" << gmean / 2.0 << "\t" << var << "\t";
			std::vector<double> strata_range(strataList_size);
			int strata_cnt = 0;
			arma::uvec strata_gmiss;
			arma::vec strata_g;
			arma::vec strata_g2;

			for (const auto &strata : strata_list)
			{
				std::ext::V_int vec = strata.second;
				arma::uvec strata_tmp = arma::conv_to<arma::uvec>::from(vec);
				strata_gmiss = gmiss.elem(strata_tmp);
				strata_g = g.elem(strata_tmp);
				strata_g2 = g2.elem(strata_tmp);//To calc var
				strata_AF(m, strata_cnt) = mean(strata_g.elem(find(strata_gmiss == 0))) / 2.0;
				arma::vec tmp = strata_g.elem(find(strata_gmiss == 0));
				arma::vec tmp2 = strata_g2.elem(find(strata_gmiss == 0)); //To calc var
				strata_Var(m, strata_cnt) = (mean(tmp2) - (mean(tmp) * mean(tmp))) * static_cast<double>(tmp.n_elem) / static_cast<double>(tmp.n_elem -1);//(gsqmean - gmean * gmean) * static_cast<double>(n - nmiss) / static_cast<double>(n - nmiss - 1);
				strata_N(m, strata_cnt) = tmp.n_elem;
		
				strata_cnt++;
			}

			for (int strata_idx = 0; strata_idx < strataList_size; strata_idx++)
			{
				writeout << strata_N(m, strata_idx) << "\t";
				writeout << strata_AF(m, strata_idx) << "\t";
				writeout << strata_Var(m, strata_idx) << "\t";
			}
		}

		for (size_t j = 0; j < n; ++j)
		{
			if (gmiss(j) == 1)
			{
				g[j] = gmean;
			}
		}

		double AF = gmean / 2.0; // convert mean to allele freq
		
		if ((static_cast<double>(nmiss) / n > missrate) || ((AF < minmaf) || (AF > maxmaf)))
		{ 
			snp_skip[npbidx] = 1;
		}
		else
		{
			G.col(npbidx) = g; 
		}

		tmpout[npbidx] = writeout.str();
		writeout.clear();
		npbidx++;

		glmm_gei(std::string (snpID.begin(), snpID.end()), G, snp_skip, npbidx, npb, n, ei, qi, null_obj, 
				writefile, meta_output, tmpout, m, end);


		if ((m) % 10000 == 0)
		{
			writefile << std::flush;
		} 
	}

	if (end % 10000 != 0)
	{
		writefile << std::flush;
	}

	writefile.close();
	writefile.clear();
	libdeflate_free_decompressor(decompressor);
	fclose(fp);
}


void glmm_gei_pgen13(Magee_Arma const& null_obj, string const &pgenfile, 
					std::string pvarFile, string const &outfile, double minmaf,
					double missrate, size_t npb, int ei, int qi, 
					std::ext::Map_str_Vint const &strata_list, 
					uint begin, uint end, std::ext::V_lluint pgenPos, 
					bool filterVariants, int pvarLength, int pvarLast,
					std::ext::V_int pvarIndex, bool meta_output)
{
    std::ifstream fIDMat;
    fIDMat.open(pvarFile);
    std::string IDline;
    std::string tmpvalue;
    std::ext::V_string tmpvalues;
    int prev = fIDMat.tellg();
	std::ext::V_string geno_snpid(npb);
	// bool isDupeID = null_obj.dupflag;
    std::ext::V_int select = null_obj.select;
    int strataList_size = strata_list.size();
    bool skip_strata = strata_list.empty();
    double maxmaf = 1 - minmaf;
	size_t n = null_obj.n;
	size_t n_obs = null_obj.n_obs;
	double gsqmean;
	arma::mat G(n, npb);
	arma::vec g(n);
	arma::vec g2(n);
    arma::uvec gmiss(n);
	arma::uvec snp_skip = arma::zeros<arma::uvec>(npb);
	std::ext::V_string tmpout(npb);
	double gmean, geno, gmax, gmin;
	size_t ncount, nmiss, npbidx = 0;
	
	arma::mat strata_AF(end +1, strataList_size);
	arma::mat strata_N(end +1, strataList_size);
	arma::mat strata_Var(end +1, strataList_size);

	std::string output = outfile + "_bin_" + std::to_string(begin) + ".tmp";
	std::ofstream writefile(output, std::ios::binary);

    if (!writefile) {
        std::cerr << "Unable to open file for appending: " << outfile << std::endl;
        std::exit(EXIT_FAILURE);
    }


	while (getline(fIDMat, IDline)) 
	{
        std::istringstream iss(IDline);
        while (getline(iss, tmpvalue, '\t')) 
		{
            tmpvalue.erase(std::remove(tmpvalue.begin(), tmpvalue.end(), '\r'), tmpvalue.end());
            tmpvalues.push_back(tmpvalue);
        }
        if (tmpvalues[0].rfind("##", 0) != 0) 
		{
            break;
        }
        prev = fIDMat.tellg();
        tmpvalues.clear();
    }

    if (tmpvalues[0].compare("#CHROM") != 0) 
	{
        fIDMat.seekg(prev);
    }

    uint32_t skipIndex = 0;
    if (!filterVariants) 
	{
        while (skipIndex != begin) 
		{
            getline(fIDMat, IDline);
            skipIndex++;
        }
    }
    else {
        while (skipIndex != pgenPos[begin]) 
		{
            getline(fIDMat, IDline);
            skipIndex++;
        }
    }

    const char* geno_filename = pgenfile.c_str();
    plink2::PgenFileInfo _info_ptr;
    plink2::PreinitPgfi(&_info_ptr);
    plink2::PgenHeaderCtrl header_ctrl;

    uintptr_t pgfi_alloc_cacheline_ct;
    char errstr_buf[plink2::kPglErrstrBufBlen];
    if (PgfiInitPhase1(geno_filename, geno_filename, UINT32_MAX, UINT32_MAX, &header_ctrl, &_info_ptr, &pgfi_alloc_cacheline_ct, errstr_buf) != plink2::kPglRetSuccess) 
	{
        throw std::runtime_error(errstr_buf);
    }

    const uint32_t raw_variant_ct = _info_ptr.raw_variant_ct;
    const uint32_t file_sample_ct = _info_ptr.raw_sample_ct;

    unsigned char* pgfi_alloc = nullptr;
    if (plink2::cachealigned_malloc(pgfi_alloc_cacheline_ct * plink2::kCacheline, &pgfi_alloc)) 
	{
        cerr << "Out of memory" << endl;
    }

    uint32_t max_vrec_width;
    uintptr_t pgr_alloc_cacheline_ct;
    if (PgfiInitPhase2(header_ctrl, 1, 0, 0, 0, raw_variant_ct, &max_vrec_width, &_info_ptr, pgfi_alloc, &pgr_alloc_cacheline_ct, errstr_buf)) 
	{
        if (pgfi_alloc && (!_info_ptr.vrtypes)) 
		{
            plink2::aligned_free(pgfi_alloc);
        }
        throw std::runtime_error(errstr_buf);
    }
    
    plink2::PgenVariant _pgv;
    plink2::PgenReader _state_ptr;
    plink2::PreinitPgr(&_state_ptr);
    plink2::PgrSetFreadBuf(nullptr, &_state_ptr);
    const uintptr_t pgr_alloc_main_byte_ct = pgr_alloc_cacheline_ct * plink2::kCacheline;
    const uintptr_t sample_subset_byte_ct = plink2::DivUp(file_sample_ct, plink2::kBitsPerVec) * plink2::kBytesPerVec;
    const uintptr_t cumulative_popcounts_byte_ct = plink2::DivUp(file_sample_ct, plink2::kBitsPerWord * plink2::kInt32PerVec) * plink2::kBytesPerVec;
    const uintptr_t genovec_byte_ct = plink2::DivUp(file_sample_ct, plink2::kNypsPerVec) * plink2::kBytesPerVec;
    const uintptr_t dosage_main_byte_ct = plink2::DivUp(file_sample_ct, (2 * plink2::kInt32PerVec)) * plink2::kBytesPerVec;
    uintptr_t multiallelic_hc_byte_ct = 0;

    unsigned char* pgr_alloc;
    if (plink2::cachealigned_malloc(pgr_alloc_main_byte_ct + (2 * plink2::kPglNypTransposeBatch + 5) * sample_subset_byte_ct + cumulative_popcounts_byte_ct + (1 + plink2::kPglNypTransposeBatch) * genovec_byte_ct + multiallelic_hc_byte_ct + dosage_main_byte_ct + plink2::kPglBitTransposeBufbytes + 4 * (plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 8), &pgr_alloc)) 
	{
        cerr << "Out of memory" << endl;
    }

    plink2::PglErr reterr = PgrInit(geno_filename, max_vrec_width, &_info_ptr, &_state_ptr, pgr_alloc);
    if (reterr) 
	{
        throw std::runtime_error("Out of memory.");
    }

    unsigned char* pgr_alloc_iter = &(pgr_alloc[pgr_alloc_main_byte_ct]);
    uintptr_t* _subset_include_vec = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    uintptr_t* _subset_include_interleaved_vec = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    _subset_include_interleaved_vec[-1] = 0;

    pgr_alloc_iter = &(pgr_alloc_iter[cumulative_popcounts_byte_ct]);
    _pgv.genovec = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[genovec_byte_ct]);

    _pgv.phasepresent = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    _pgv.phaseinfo = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    _pgv.dosage_present = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    _pgv.dosage_main = reinterpret_cast<uint16_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[dosage_main_byte_ct]);

    uint32_t _subset_size = file_sample_ct;
    plink2::PgrSampleSubsetIndex _subset_index;
    pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglBitTransposeBufbytes]);

    // int variant_index = 0;
    // int keepIndex = 0;
    vector<double> buf(file_sample_ct);
    
	for (uint m = begin; m <= end; ++m)
	{
		std::ostringstream writeout;
		
		if (m >= _info_ptr.raw_variant_ct) 
		{
			char errstr_buf[256];
			sprintf(errstr_buf, "variant_num out of range (%d; must be 1..%u)", m + 1, _info_ptr.raw_variant_ct);
			cerr << errstr_buf << "\n";
		}

		uint32_t dosage_ct;
		string value;
		std::ext::V_string values;
		if (!filterVariants) 
		{
			reterr = plink2::PgrGet1D(_subset_include_vec, _subset_index, _subset_size, m, 1, &_state_ptr, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, &dosage_ct);
			getline(fIDMat, IDline);
			std::istringstream iss(IDline);
			while (getline(iss, value, '\t')) 
			{
				values.push_back(value);					
			}
		}
		else 
		{
			while (skipIndex != pgenPos[m]) 
			{
				getline(fIDMat, IDline);
				skipIndex++;
			}
			getline(fIDMat, IDline);
			skipIndex++;
			std::istringstream iss(IDline);
			while (getline(iss, value, '\t')) 
			{
				values.push_back(value);
			}
			reterr = plink2::PgrGet1D(_subset_include_vec, _subset_index, _subset_size, pgenPos[m], 1, &_state_ptr, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, &dosage_ct);
		}

		plink2::Dosage16ToDoubles(plink2::kGenoDoublePairs, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, _subset_size, dosage_ct, &buf[0]);
		
		gmean = 0.0;
		gsqmean = 0.0;
		gmax = -100.0;
		gmin = 100.0;
		nmiss = 0;
		ncount = 0;
		int idx_k = 0;

		for (uint32_t ct = 0; ct < file_sample_ct; ct++) 
		{
			if (buf[ct] == -9.0) 
			{
				if (select[idx_k] >= 0)
				{
					// missingIndex.push_back(idx_k);
					gmiss(select[idx_k]) = 1;
					nmiss++;
					// idx_k++;
				}
				idx_k++;
				continue;
			}

			if (select[idx_k] >= 0) 
			{
				geno = buf[ct];
				gmiss(select[idx_k]) = 0;
				g[select[idx_k]] = geno;
				g2[select[idx_k]] = geno * geno;
				gmean += geno;
				gsqmean += geno * geno;

				if (geno > gmax)
				{
					gmax = geno;
				}
				if (geno < gmin)
				{
					gmin = geno;
				}
			}
			idx_k++;
		}

		gmean /= static_cast<double>(n - nmiss);
		gsqmean /= static_cast<double>(n - nmiss);
		//rsq = (gsqmean - gmean * gmean) * static_cast<double>(n - nmiss) / static_cast<double>(n - nmiss - 1) / (gmean * (1.0 - gmean / 2.0));
		double var = (gsqmean - gmean * gmean) * static_cast<double>(n - nmiss) / static_cast<double>(n - nmiss - 1);
		double AF = gmean / 2.0;
		double percMissing = nmiss / (n * 1.0);

		for (size_t j = 0; j < n; ++j)
		{
			if (gmiss(j) == 1)
			{
				g[j] = gmean;
			}
		}

		if ((static_cast<double>(nmiss) / n > missrate) || ((AF < minmaf) || (AF > maxmaf)))
		{ 
			snp_skip[npbidx] = 1;
		}
		else
		{
			G.col(npbidx) = g; 
		}

		std::string tmpString = "";
		std::string snpID;
		values[pvarLast].erase(std::remove(values[pvarLast].begin(), values[pvarLast].end(), '\r'), values[pvarLast].end());

		for (int p = 0; p < pvarLength; p++) 
		{
			tmpString = tmpString + values[pvarIndex[p]] + "\t";
			if(p == 0)
			{
				snpID = values[pvarIndex[p]];
			}
		}

		geno_snpid[npbidx] = tmpString + std::to_string(n - nmiss);
		
		
		if (skip_strata)
		{
			writeout << geno_snpid[npbidx] << "\t" << gmean / 2.0 << "\t" << var << "\t";
		}
		else
		{
			writeout << geno_snpid[npbidx] << "\t" << gmean / 2.0 << "\t" << var << "\t";
			std::vector<double> strata_range(strataList_size);
			int strata_cnt = 0;
			arma::uvec strata_gmiss;
			arma::vec strata_g;
			arma::vec strata_g2;

			for (const auto &strata : strata_list)
			{
				std::ext::V_int vec = strata.second;
				arma::uvec strata_tmp = arma::conv_to<arma::uvec>::from(vec);
				strata_gmiss = gmiss.elem(strata_tmp);
				strata_g = g.elem(strata_tmp);
				strata_g2 = g2.elem(strata_tmp);
				strata_AF(m, strata_cnt) = mean(strata_g.elem(find(strata_gmiss == 0))) / 2.0;
				arma::vec tmp = strata_g.elem(find(strata_gmiss == 0));
				arma::vec tmp2 = strata_g2.elem(find(strata_gmiss == 0)); //To calc var
				strata_Var(m, strata_cnt) = (mean(tmp2) - (mean(tmp) * mean(tmp))) * tmp.n_elem/ static_cast<double>(tmp.n_elem -1);//(gsqmean - gmean * gmean) * static_cast<double>(n - nmiss) / static_cast<double>(n - nmiss - 1);
				strata_N(m, strata_cnt) = tmp.n_elem;
				strata_cnt++;
			}
		
			for (int strata_idx = 0; strata_idx < strataList_size; strata_idx++)
			{
				writeout << strata_N(m, strata_idx) << "\t";
				writeout << strata_AF(m, strata_idx) << "\t";
				writeout << strata_Var(m, strata_idx) << "\t";
			}
		}

		
		tmpout[npbidx] = writeout.str();
		writeout.clear();
		npbidx++;
		
		glmm_gei(snpID, G, snp_skip, npbidx, npb, n, ei, qi, null_obj, 
				writefile, meta_output, tmpout, m, end);

		if ((m) % 10000 == 0)
		{
			writefile << std::flush;
		} 
	}

	if (end % 10000 != 0)
	{
		writefile << std::flush;
	}

	writefile.close();
	writefile.clear();
	fIDMat.close();
}


void glmm_gei_bed13(Magee_Arma const& null_obj, string const &bedfile, 
					std::string bimFile, string const &outfile, double minmaf,
					double missrate, size_t npb, int ei, int qi, 
					std::ext::Map_str_Vint const &strata_list, 
					uint begin, uint end, std::ext::V_lluint bedPos, 
					bool filterVariants, char bimDelim, int bimLast, 
					uint32_t n_samples, bool meta_output)
{
	std::ifstream fIDMat;
    fIDMat.open(bimFile);
    std::string IDline;
	std::string geno_snpid;
    std::ext::V_int select = null_obj.select;
    int strataList_size = strata_list.size();
    bool skip_strata = strata_list.empty();
    double maxmaf = 1 - minmaf;
	double gsqmean;
	std::string output = outfile + "_bin_" + std::to_string(begin) + ".tmp";
	std::ofstream writefile(output, std::ios::binary);

    if (!writefile) 
	{
        std::cerr << "Unable to open file for appending: " << outfile << std::endl;
        std::exit(EXIT_FAILURE);
    }

	size_t n = null_obj.n;
	size_t n_obs = null_obj.n_obs;

	arma::mat G(n, npb);
	arma::vec g(n);
	arma::vec g2(n);
    arma::uvec gmiss(n);
	arma::uvec snp_skip = arma::zeros<arma::uvec>(npb);
	std::ext::V_string tmpout(npb);
	double gmean, geno, gmax, gmin;
	size_t ncount, nmiss, npbidx = 0;
	arma::mat strata_AF(end +1, strataList_size);
	arma::mat strata_N(end +1, strataList_size);
	arma::mat strata_Var(end +1, strataList_size);
	uint32_t skipIndex = 0;
	
    
    if (!filterVariants) 
	{
        while (skipIndex != begin) 
		{
            getline(fIDMat, IDline);
            skipIndex++;
        }
    }
    else 
	{
        while (skipIndex != bedPos[begin]) 
		{
            getline(fIDMat, IDline);
            skipIndex++;
        }       
    }
	   
    std::ifstream readbedfile(bedfile.c_str(), std::ios::binary);
    uint nblocks = (n_samples + 3) / 4, pos;
    unsigned char temp[2];
    unsigned char* buffer = new unsigned char[nblocks];

	for (uint m = begin; m <= end; ++m)
	{
		std::ostringstream writeout;
		std::string value;
		std::vector <string> values;
		if (!filterVariants) 
		{
			readbedfile.seekg((std::streamoff)m * nblocks + 3, readbedfile.beg);
			getline(fIDMat, IDline);             
			std::istringstream iss(IDline);
			while (getline(iss, value, bimDelim)) 
			{
				values.push_back(value);  
			}			
		}
		else 
		{
			while (skipIndex != bedPos[m]) 
			{
				getline(fIDMat, IDline);
				skipIndex++;
			}
			readbedfile.seekg((std::streamoff)bedPos[m] * nblocks + 3, readbedfile.beg);                
			getline(fIDMat, IDline);
			skipIndex++;
			std::istringstream iss(IDline);
			while (getline(iss, value, bimDelim)) 
			{
				values.push_back(value);
			}
		}

		readbedfile.read((char*)buffer, nblocks);
		gmean = 0.0;
		gsqmean = 0.0;
		gmax = -100.0;
		gmin = 100.0;
		nmiss = 0;
		ncount = 0;
		int idx_k = 0;

		for (size_t block = 0; block < nblocks; block++) 
		{
            pos = 0;

            for (int i = 0; i < 4; i++) 
			{
				if ((ncount == n_samples) && (block == nblocks - 1)) 
				{
					break;
				}
				for (size_t l = 0; l < 2; ++l) 
				{
					temp[l] = (buffer[block] >> pos) & 1;
					pos++;
				}
				
				if (select[idx_k] == -1) 
				{
					ncount++;
					idx_k++;
					continue;
				}
			
				if (temp[0] == 0 && temp[1] == 0) 
				{
					geno = 2.0;
				}
				else if (temp[0] == 1 && temp[1] == 1) 
				{
					geno = 0.0;
				}
				else if (temp[0] == 0 && temp[1] == 1) 
				{
					geno = 1.0;
				}
				else 
				{
					// missingIndex.push_back(idx_k);
					gmiss(select[idx_k]) = 1;
					nmiss++;
					idx_k++;
					ncount++;
					continue;
				}

				gmean += geno;
				gsqmean += geno * geno;

				if (geno > gmax)
				{
					gmax = geno;
				}
				if (geno < gmin)
				{
					gmin = geno;
				}

				g(select[idx_k]) = geno;
				g2(select[idx_k]) = geno * geno;
				gmiss(select[idx_k]) = 0;
				idx_k++;
				ncount++;
			}
        }
		
		gmean /=  static_cast<double>(n - nmiss);
		gsqmean /= static_cast<double>(n - nmiss);
		//rsq = (gsqmean - gmean * gmean) * static_cast<double>(n - nmiss) / static_cast<double>(n - nmiss - 1) / (gmean * (1.0 - gmean / 2.0));
		double var = (gsqmean - gmean * gmean) * static_cast<double>(n - nmiss) / static_cast<double>(n - nmiss - 1);
		double AF = gmean / 2.0;
		double percMissing = nmiss / (n * 1.0);
		
		for (size_t j = 0; j < n; ++j)
		{
			if (gmiss(j) == 1)
			{
				g[j] = gmean;
			}
		}

		if ((static_cast<double>(nmiss) / n > missrate) || ((AF < minmaf) || (AF > maxmaf)))
		{
            snp_skip[npbidx] = 1;
			continue;                           
        }
		else
		{
			G.col(npbidx) = g; 
		}

		values[bimLast].erase(std::remove(values[bimLast].begin(), values[bimLast].end(), '\r'), values[bimLast].end());
        geno_snpid = values[1] + "\t" + values[0] + "\t" + values[bimLast - 2] + "\t" + values[bimLast] + "\t" + values[bimLast - 1] + "\t" + std::to_string(n - nmiss);

		if (skip_strata)
		{
			writeout << geno_snpid << "\t" << gmean / 2.0 << "\t" << var << "\t";
		}
		else
		{
			writeout << geno_snpid << "\t" << gmean / 2.0 << "\t" << var << "\t";
			std::vector<double> strata_range(strataList_size);
			int strata_cnt = 0;
			arma::uvec strata_gmiss;
			arma::vec strata_g;
			arma::vec strata_g2;

			for (const auto &strata : strata_list)
			{
				std::ext::V_int vec = strata.second;
				arma::uvec strata_tmp = arma::conv_to<arma::uvec>::from(vec);
				strata_gmiss = gmiss.elem(strata_tmp);
				strata_g = g.elem(strata_tmp);
				strata_g2 = g2.elem(strata_tmp);
				strata_AF(m,strata_cnt) = mean(strata_g.elem(find(strata_gmiss == 0))) / 2.0;
				arma::vec tmp = strata_g.elem(find(strata_gmiss == 0));
				arma::vec tmp2 = strata_g2.elem(find(strata_gmiss == 0));
				strata_Var(m, strata_cnt) = (mean(tmp2) - (mean(tmp) * mean(tmp))) * tmp.n_elem/ static_cast<double>(tmp.n_elem -1);//(gsqmean - gmean * gmean) * static_cast<double>(n - nmiss) / static_cast<double>(n - nmiss - 1);
				strata_N(m,strata_cnt) = tmp.n_elem;
				strata_cnt++;
			}
		
			for (int strata_idx = 0; strata_idx < strataList_size; strata_idx++)
			{
				writeout << strata_N(m, strata_idx) << "\t";
				writeout << strata_AF(m, strata_idx) << "\t";
				writeout << strata_Var(m, strata_idx) << "\t";
			}
		}

		tmpout[npbidx] = writeout.str();
				
		writeout.clear();
		npbidx++;

		glmm_gei(values[1], G, snp_skip, npbidx, npb, n, ei, qi, null_obj, 
				writefile, meta_output, tmpout, m, end);

		if ((m) % 10000 == 0)
		{
			writefile << std::flush;
		} 
	}

	delete[] buffer;
	buffer = nullptr;
	
	if (end % 10000 != 0)
	{
		writefile << std::flush;
	}

	writefile.close();
	writefile.clear();
	fIDMat.close();
}


#pragma once
#include "GMMAT.h"
#include "MAGEE.h"
#include "declars.h"
#include "Logger.h"


/**
 * @brief Extract the output/log file name from command-line arguments.
 *
 * This function scans the command-line inputs and looks for the option
 * `--out`. If found, the following argument is returned as the log/output
 * file name.
 *
 * If the option is not provided, an empty string (fallback) is returned.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line argument strings.
 *
 * @return std::string The file name specified after `--out`, or an empty string.
 */
std::string get_log_name(int argc, char* argv[]);

/**
 * @brief Determine whether the phenotype is binary or continuous.
 *
 * This function inspects the phenotype values across all provided sample IDs.
 * It counts the number of unique phenotype categories:
 *
 * - If only 2 unique values exist → phenotype is treated as binary.
 * - If more than 2 unique values exist → phenotype is treated as continuous.
 *
 * If all phenotype values are identical, an error is printed.
 *
 * When binary, the logistic regression convergence threshold is also reported.
 *
 * @param phenoMap Map from sample ID to phenotype vector.
 * @param sampleID Vector of sample IDs in analysis order.
 * @param epsilon Convergence tolerance used in logistic regression.
 *
 * @return int Returns 1 if phenotype is binary, 0 if continuous.
 */
int  checkBinary(unordered_map<string, vector<vector<string>>> phenoMap, vector<string> sampleID, double epsilon);


void center(int center, int scale, int samSize, int numSelCol, vector<double> covdata, vector<double>* covdata_ret);

/**
 * @brief Fit null regression model and return additional intermediate outputs.
 *
 * This is an extended version of fitNullModel().
 *
 * In addition to standard null model outputs, it also returns:
 *
 * - beta estimates (alpha)
 * - linear predictor X*beta (eta)
 *
 * These intermediate quantities are required for downstream
 * GMMAT-style score statistic calculations.
 *
 * Supports both:
 * - Linear regression 
 * - Logistic regression 
 *
 * @param samSize 
 * @param numSelCol 
 * @param phenoType 
 * @param epsilon 
 * @param robust 
 * @param covSelHeadersName 
 * @param phenodata 
 * @param covdata 
 * @param XinvXTX_ret 
 * @param miu_ret 
 * @param resid_ret 
 * @param sigma2_ret 
 * @param beta_ret 
 * @param Xbeta_ret 
 */
void fitNullModel2(int samSize, int numSelCol, int phenoType, double epsilon, 
                    int robust, std::vector<string> covSelHeadersName, std::vector<double> phenodata, 
                    std::vector<double> covdata, std::vector<double>* XinvXTX_ret, vector<double>* miu_ret, 
                    vector<double>* resid_ret, double* sigma2_ret, std::vector<double>& beta_ret,
                    std::vector<double>& Xbeta_ret);



/**
 * @brief Fit the null regression model (linear or logistic).
 *
 * This function fits a regression model without genotype effects:
 *
 * - Linear regression for continuous traits
 * - Logistic regression for binary traits
 *
 * Outputs include:
 * - Residual vector
 * - Dispersion estimate (sigma²)
 * - Mean response (mu)
 * - Precomputed matrix X 
 * If logistic regression fails to converge after MAX_ITER,
 * the program terminates with an error.
 *
 * @param samSize 
 * @param numSelCol 
 * @param phenoType 
 * @param epsilon 
 * @param robust 
 * @param covSelHeadersName 
 * @param phenodata 
 * @param covdata 
 * @param XinvXTX_ret 
 * @param miu_ret 
 * @param resid_ret 
 * @param sigma2_ret 
 */
void fitNullModel(int samSize, int numSelCol, int phenoType, double epsilon, 
                    int robust, std::vector<string> covSelHeadersName, std::vector<double> phenodata, 
                    std::vector<double> covdata, std::vector<double>* XinvXTX_ret, vector<double>* miu_ret, 
                    vector<double>* resid_ret, double* sigma2_ret);

/**
 * @brief Print regression coefficients and covariance matrix of estimates.
 *
 *  Function outputs:
 *
 * - Estimated regression coefficients (beta)
 * - Standard errors
 * - Z-values
 * - Chi-square p-values
 *
 * It also prints the full variance-covariance matrix.
 *
 * Used after fitting the null model in linear or logistic regression.
 *
 * @param numCovs Number of covariates.
 * @param covNames Names of covariates.
 * @param covVarMat Variance-covariance matrix.
 * @param beta Estimated regression coefficients.
 * @param phenoType Phenotype type (1=binary, 0=continuous).
 * @param samSize Number of samples.
 */
void printCovVarMat(int numCovs, vector<string> covNames, double* covVarMat, double* beta, int phenoType, int samSize);


/**
 * @brief Write the header line for the GWAS output results file.
 *
 * This function generates the correct column header depending on:
 *
 * The header includes:
 * - Variant information (CHR, POS, alleles, AF)
 * - Marginal effects
 * - Interaction effects
 * - Standard errors
 * - Covariance terms
 * - P-values 
 *
 * @param useBgen 
 * @param numExpSelCol_new 
 * @param Sq1 
 * @param covNames 
 * @param output 
 * @param outStyle 
 * @param robust 
 * @param sigma2 
 * @param binE 
 */
void printOutputHeader(bool useBgen, int numExpSelCol_new, int Sq1, vector<string> covNames, string output, string outStyle, 
                       int robust, double sigma2, BinE binE);

            


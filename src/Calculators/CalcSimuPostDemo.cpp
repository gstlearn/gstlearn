/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "Calculators/CalcSimuPostDemo.hpp"

CalcSimuPostDemo::CalcSimuPostDemo()
    : CalcSimuPost()
{
}

CalcSimuPostDemo::~CalcSimuPostDemo()
{
}

/**
 * Returns the number of variables after transformation
 * In the current version, this number is set to 2
 * @return
 */
int CalcSimuPostDemo::getTransfoNvar() const
{
  return 2;
}

/**
 * Perform the Transformation to convert the multivariate input vector 'tabin'
 * into a multivariate output vector.
 * @param Z_n_k_s  Input information (Dimension: 's' x 'n')
 *
 * @return Output information (Dimension: 's' x 'p')
 *
 * @remark In the current version, the Transformation consists in calculating:
 * - the sum of all values
 * - the standard deviation of all values
 */
VectorVectorDouble CalcSimuPostDemo::transformFunction(const VectorVectorDouble& Z_n_k_s) const
{
  int nvarout = getTransfoNvar();
  int nsample = (int) Z_n_k_s.size();
  VectorVectorDouble Y_p_k_s(nsample);

  for (int is = 0; is < nsample; is++)
    Y_p_k_s[is].resize(nvarout);

  for (int is = 0; is < nsample; is++)
    Y_p_k_s[is][0] = VH::cumul(Z_n_k_s[is]);

  for (int is = 0; is < nsample; is++)
    Y_p_k_s[is][1] = VH::stdv(Z_n_k_s[is]);

  return Y_p_k_s;
}

/**
 *
 * @param dbin Input data base
 * @param dbout Output data base (must be a Grid)
 * @param names Vector of simulation names
 * @param flag_match True if the ranks of simulations must match; otherwise: product
 * @param upscale Option within EPostUpscale
 * @param stats Vector of options within EPostStat
 * @param verbose Verbose flag
 * @param rank_check Rank (1-based) of the target element to be checked (0: None; -1: All)
 * @param namconv Naming convention
 * @return Error code
 *
 * @remark
 * - N : number of variables in 'dbin' defined by 'names' (index 'n')
 * - k_1, ..., k_N : number of outcomes for each variable
 * - K : number of multivariate simulations according to 'flag_match' (index 'k')
 * - Transformation function 'F' from Z in R_N to Y in R_P (index 'p')
 *
 * Description of the Flow Chart:
 *
 * -# Loop over cells of 'dbout' (index 'C')
 * -# Loop over the simulations (index 'k')
 * -# Find the active samples of 'dbin' in 'C' (index 's') and build table of Z_n^{k}(s)
 * -# Apply the Transform function to table: Y_p^{k}(s)
 * -# Upscale to the target cell according to upscaling rule '_upscale': up_Y_p^{k}(C)
 * -# Compute statistics according to stat rule '_stats'
 */
int simuPostDemo(Db *dbin,
                 DbGrid *dbout,
                 const VectorString &names,
                 bool flag_match,
                 const EPostUpscale &upscale,
                 const std::vector<EPostStat> &stats,
                 bool verbose,
                 int rank_check,
                 const NamingConvention &namconv)
{
  CalcSimuPostDemo calcul;
  calcul.setDbin(dbin);
  calcul.setDbout(dbout);
  calcul.setNames(names);
  calcul.setUpscale(upscale);
  calcul.setStats(stats);
  calcul.setFlagMatch(flag_match);
  calcul.setVerbose(verbose);
  calcul.setRankCheck(rank_check);
  calcul.setNamingConvention(namconv);

  int error = (calcul.run()) ? 0 : 1;
  return error;
}

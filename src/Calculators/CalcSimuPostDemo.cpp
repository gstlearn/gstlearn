/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
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
int CalcSimuPostDemo::_getTransfoNvar() const
{
  return 2;
}

/**
 * Perform the Transformation to convert one multivariate input vector 'tabin'
 * into one multivariate output vector.
 * @param Z_n_k_s  Input information (Dimension: 'n')
 * @param Y_p_k_s  Output information (Dimension: 'p')
 *
 * @remark In the current version, the Transformation consists in calculating:
 * - the mean of all values
 * - the standard deviation of all values
 */
void CalcSimuPostDemo::_transformFunction(const VectorDouble& Z_n_k_s, VectorDouble& Y_p_k_s) const
{
  Y_p_k_s[0] = VH::mean(Z_n_k_s);
  Y_p_k_s[1] = VH::stdv(Z_n_k_s);
}

/**
 * This is a particular use of Simulation Post-Processing functions.
 *
 * Its specificity is an embedded transformation function which transforms the multivariate simulation vector
 * (combining the one outcome of each variable, and for a given sample)
 * This transformation (which is provided as an example here):
 *
 * - uses the multivariate simulation vector in input (N elements)
 * - produces an output vector (Dimension 2) with the *mean* and the *standard deviation* of the values of the input vector
 *
 * For a detailed list of arguments, see \link CalcSimuPost.cpp simuPost \endlink
 *
 */
int simuPostDemo(Db *dbin,
                 DbGrid *dbout,
                 const VectorString &names,
                 bool flag_match,
                 const EPostUpscale &upscale,
                 const std::vector<EPostStat> &stats,
                 bool verbose,
                 const VectorInt& check_targets,
                 int check_level,
                 const NamingConvention &namconv)
{
  CalcSimuPostDemo calcul;
  calcul.setDbin(dbin);
  if (dbout != nullptr)
  {
    calcul.setFlagUpscale(true);
    calcul.setDbout(dbout);
  }
  calcul.setNames(names);
  calcul.setUpscale(upscale);
  calcul.setStats(stats);
  calcul.setFlagMatch(flag_match);
  calcul.setVerbose(verbose);
  calcul.setCheckTargets(check_targets);
  calcul.setCheckLevel(check_level);
  calcul.setNamingConvention(namconv);

  int error = (calcul.run()) ? 0 : 1;
  return error;
}

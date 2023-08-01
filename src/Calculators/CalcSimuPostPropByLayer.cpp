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
#include <Calculators/CalcSimuPostPropByLayer.hpp>

CalcSimuPostPropByLayer::CalcSimuPostPropByLayer()
    : CalcSimuPost(),
      _dbgrid(nullptr)
{
}

CalcSimuPostPropByLayer::~CalcSimuPostPropByLayer()
{
}

bool CalcSimuPostPropByLayer::_check()
{
  // Be compliant with the general constraints of SimuPost

  if (! CalcSimuPost::_check()) return false;

  // The upscaling is compulsory here
  if (! _getFlagUpscale())
  {
    messerr("The output 'Db' (organized as a Grid is compulsory");
    return false;
  }

  // This version has some restrictions:
  // - It has been coded for a space dimension of output Grid limited to 3

  int ndim_out = getDbout()->getNDim();
  if (ndim_out > 3)
  {
    messerr("The current version has been coded for a Space Dimension of 'dbout' (%d) limited to 3",
            ndim_out);
    return false;
  }
  _dbgrid = dynamic_cast<const DbGrid*>(getDbout());
  if (_dbgrid->isGridRotated())
  {
    messerr("The current version has been coded for a non-rotated output grid");
    return false;
  }
  return true;
}

/**
 * Returns the number of variables after transformation
 * This number is one more than the number of input simulated variables
 * @return Number of output variables
 */
int CalcSimuPostPropByLayer::_getTransfoNvar() const
{
  return _getNVar() + 1;
}

/**
 * Perform the Transformation to convert one multivariate input vector 'tabin'
 * into one multivariate output vector.
 *
 * @param Z_n_k_s  Input information (Dimension: 'n')
 * @param Y_p_k_s  Output information (Dimension: 'p')
 *
 */
void CalcSimuPostPropByLayer::_transformFunction(const VectorDouble& Z_n_k_s, VectorDouble& Y_p_k_s) const
{
  int nlayer = (int) Z_n_k_s.size();
  int ndim_out = getDbout()->getNDim();
  int iechout = _getIechout();

  double z_ref  = _dbgrid->getCoordinate(iechout, ndim_out-1);
  double h_max  = _dbgrid->getDX(ndim_out - 1);
  double z_base = z_ref - h_max / 2.;

  /* previous implementation
  double previous = 0.;
  double cote = 0.;
  for (int ilayer = 0; ilayer < nlayer; ilayer++)
  {
    cote += Z_n_k_s[ilayer];
    Y_p_k_s[ilayer] = MIN(MAX(cote - z_base, 0.), h_max) - previous;
    previous = Y_p_k_s[ilayer];
  }
  Y_p_k_s[nlayer] = h_max - Y_p_k_s[nlayer-1];

   */

  // from bottom to top
  double cote = 0.;
  for (int ilayer = 0; ilayer < nlayer; ilayer++)
  {
    cote += Z_n_k_s[ilayer];
    Y_p_k_s[ilayer] = MIN(MAX(cote - z_base, 0.), h_max);
  }
  Y_p_k_s[nlayer] = h_max;

  // from top to bottom
  for (int ilayer = nlayer; ilayer > 0; ilayer--)
    {
      Y_p_k_s[ilayer] -= Y_p_k_s[ilayer-1];
    }

  // Normalize by the extension of the cell
  for (int ilayer = 0; ilayer <= nlayer; ilayer++)
    Y_p_k_s[ilayer] /= h_max;
}

/**
 * This is a particular use of Simulation Post-Processing functions. Its specificity comes from its transformation function.
 *
 * It is assumed that each input variable corresponds to the thickness of ordered layers.
 * This function receives a vector of multivariate information (for each combination of simulation outcome
 * and for each sample of the input 'db').
 *
 * If N designates the number of elements of the input vector, the transformation returns a vector of N+1 elements.
 * This vector corresponds to the percentage of each layer within the target block
 *
 * @remarks: Coding Rule
 * @remarks: - code = 0 -Inf            < z <= H1
 * @remarks: - code = 1   H1            < z <= H1 + H2
 * @remarks: - code = 2   H1 + H2       < z <= H1 + H2 + H3
 * @remarks: - code = 3   H1 + H2 + ... < z <= + Inf
 *
 * For a detailed list of arguments, see \link CalcSimuPost.cpp simuPost \endlink
 */

int simuPostPropByLayer(Db *dbin,
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
  CalcSimuPostPropByLayer calcul;
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

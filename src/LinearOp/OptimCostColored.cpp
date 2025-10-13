/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "LinearOp/OptimCostColored.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "LinearOp/OptimCostBinary.hpp"

namespace gstlrn
{
OptimCostColored::OptimCostColored()
  : OptimCostBinary()
  , _nprop(0)
  , _splits()
  , _meanProps()
{
}

OptimCostColored::OptimCostColored(Id nprop,
                                   PrecisionOp* pmat,
                                   const ProjMatrix* projdata,
                                   const ProjMatrix* projseis,
                                   const VectorDouble& propseis,
                                   const VectorDouble& varseis)
  : OptimCostBinary()
  , _nprop(0)
  , _splits()
  , _meanProps()
{
  reset(nprop, pmat, projdata, projseis, propseis, varseis);
}

OptimCostColored::OptimCostColored(const OptimCostColored& m)
  : OptimCostBinary(m)
  , _nprop(m._nprop)
  , _splits(m._splits)
  , _meanProps(m._meanProps)
{
}

OptimCostColored& OptimCostColored::operator=(const OptimCostColored& m)
{
  if (this != &m)
  {
    _nprop     = m._nprop;
    _splits    = m._splits;
    _meanProps = m._meanProps;
  }
  return *this;
}

OptimCostColored::~OptimCostColored()
{
}

/*****************************************************************************/
/*!
**  Initialize the Hessian Operator
**
** \param[in]  nprop    Number of different proportions (or facies)
** \param[in]  pmat     The precision matrix to be optimized
** \param[in]  projdata The Projection operator between Data and Meshing
** \param[in]  projseis The Projection operator between Seismic and Meshing
** \param[in]  propseis Array of facies proportions
** \param[in]  varseis  Array of variance attached to the seismic
**
*****************************************************************************/
void OptimCostColored::reset(Id nprop,
                             PrecisionOp* pmat,
                             const ProjMatrix* projdata,
                             const ProjMatrix* projseis,
                             const VectorDouble& propseis,
                             const VectorDouble& varseis)
{
  // Assignment of pointers
  _nprop = nprop;
  VH::fill(_meanProps, 1. / _nprop, _nprop);
  _splits = initSplit(_nprop);

  // Pass arguments to the OptimCostBinary class

  OptimCostBinary::reset(pmat, projdata, projseis, propseis, varseis);
}

/*****************************************************************************/
/*!
**  Perform the minimization
**
** \return Array of facies proportions (Dimension: [nprop][nvertex]
**
** \param[in]  facies     Array containing the Facies values (see remarks)
**                        (Dimension: npoint)
** \param[in]  splits     Array giving the facies split
**                        (Dimension: [nfacies-1][nfacies])
** \param[in]  meanprops  Array of mean of proportions (Dimension: nfacies)
** \param[in]  verbose    Verbose flag
** \param[in]  maxiter    Maximum number of iterations for Optimization algorithm
** \param[in]  eps        Tolerance for Optimization algorithm
**
** \remarks Argument 'facies' should contain values ranging from 1 to _nprop
**
*****************************************************************************/
VectorVectorDouble OptimCostColored::minimize(const VectorDouble& facies,
                                              const VectorVectorInt& splits,
                                              const VectorDouble& meanprops,
                                              bool verbose,
                                              Id maxiter,
                                              double eps)
{
  VectorDouble indic, propfac;
  VectorVectorDouble propfacs;

  // Initialization
  if (_checkFacies(facies)) return propfacs;
  // Check the split contents
  if (_checkSplits(splits)) return propfacs;
  // Check the mean proportions
  if (_checkMeanProportions(meanprops)) return propfacs;

  try
  {
    if (!isInitialized())
      my_throw("'OptimCostColored' must be initialized beforehand");
    auto npoint  = getNPoint();
    auto nvertex = getNVertex();
    Id nlevel    = _nprop - 1;
    propfacs.resize(_nprop, VectorDouble(nvertex, 0));

    // Optional printout
    if (verbose)
    {
      message("Number of points   = %d\n", npoint);
      message("Number of vertices = %d\n", nvertex);
      message("Number of facies   = %d\n", _nprop);
      message("Number of levels   = %d\n", nlevel);
      printSplits();
    }

    // Local Core allocation

    indic.resize(npoint);
    propfac.resize(nvertex);

    for (Id level = 0; level < nlevel; level++)
    {
      if (verbose) message("\nProcessing Level %d/%d\n", level + 1, nlevel);

      // Cancel the Seismic for level higher than the first one
      toggleSeismic(level == 0);

      // Get the active facies at samples
      _getFaciesToIndic(facies, _splits[level], indic);

      // Set the proportion of the active facies
      if (setMeanProportion(_getFaciesToProportion(_splits[level])))
        my_throw("Error in '_getFaciesToProportion'");

      // Perform the minimization
      propfac = OptimCostBinary::minimize(indic, verbose, maxiter, eps);
      if (propfac.empty())
        my_throw("Error in 'OptimCostBinary'");

      // Convert the results into conditional proportions
      for (Id ip = 0; ip < _nprop; ip++)
        _copyMultProportions(level, ip, propfac, propfacs);
    }
  }
  catch (const std::string& str)
  {
    // TODO : Check if std::exception can be used
    messerr("%s", str.c_str());
  }

  return propfacs;
}

/**
 * Internal function to Extract the Indicator value of a given facies
 * @param facies Array of facies values
 * @param split  Array giving the facies split (Dimension: nfacies * (nfacies-1))
 * @param indic  Array of output indicator values
 */
void OptimCostColored::_getFaciesToIndic(const VectorDouble& facies,
                                         const VectorInt& split,
                                         VectorDouble& indic) const
{
  Id facloc;
  auto npoint = getNPoint();

  for (Id i = 0; i < npoint; i++)
  {
    facloc   = static_cast<Id>(facies[i]);
    indic[i] = TEST;
    if (facloc < 1 || facloc > _nprop) continue;

    if (split[facloc - 1] == 1)
    {
      // When facies corresponds to split==1, set indic to 0
      indic[i] = 0;
    }
    else if (split[facloc - 1] == 2)
    {
      // When facies corresponds to split==2, set indic to 1
      indic[i] = 1;
    }
    else
    {
      // When facies corresponds to split==0 set indic to TEST
      indic[i] = TEST;
    }
  }
}

/*****************************************************************************/
/*!
**  Internal function to evaluate the proportion for a given split
**
** \returns The proportion of the target facies
**
** \param[in]  split      Array giving the facies split (at a given level)
**
*****************************************************************************/
double OptimCostColored::_getFaciesToProportion(const VectorInt& split) const

{
  double sum1 = 0.;
  double sum2 = 0.;
  for (Id ip = 0; ip < _nprop; ip++)
  {
    if (split[ip] == 1)
      sum1 += _meanProps[ip];
    else if (split[ip] == 2)
      sum2 += _meanProps[ip];
  }
  return (sum2 / (sum1 + sum2));
}

/*****************************************************************************/
/*!
**  Internal function to check the set of facies input values
**
** \returns Error returned code
**
** \param[in]  facies     Array containing the Facies values (see remarks)
**                        (Dimension: npoint)
**
*****************************************************************************/
Id OptimCostColored::_checkFacies(const VectorDouble& facies) const
{
  auto npoint = getNPoint();

  Id nerr = 0;
  for (Id i = 0; i < npoint; i++)
  {
    if (FFFF(facies[i])) continue;
    if (facies[i] < 1 || facies[i] > _nprop)
    {
      messerr("Error: At sample #%d - Facies (%d) should be in [1,%d]",
              i + 1, static_cast<Id>(facies[i]), _nprop);
      nerr++;
    }
  }
  return (nerr > 0);
}

/*****************************************************************************/
/*!
**  Internal function to print the Splits
**
*****************************************************************************/
void OptimCostColored::printSplits(const VectorVectorInt& splits) const
{
  Id nlevel = _nprop - 1;

  if (splits.empty())
  {
    for (Id level = 0; level < nlevel; level++)
      VH::dump(String(), _splits[level]);
  }
  else
  {
    for (Id level = 0; level < nlevel; level++)
      VH::dump(String(), splits[level]);
  }
}

/*****************************************************************************/
/*!
**  Internal function to check the validity of the Splits
**
** \returns Error returned code
**
** \param[in]  splits     Array giving the facies split
**                        (Dimension: nfacies * (nfacies-1))
**
*****************************************************************************/
Id OptimCostColored::_checkSplits(const VectorVectorInt& splits)
{
  if (splits.empty()) return 0;
  Id nlevel = _nprop - 1;

  // Check that split values are 0, 1 or 2 only
  Id nerr = 0;
  for (Id level = 0; level < nlevel; level++)
  {
    for (Id ip = 0; ip < _nprop; ip++)
    {
      if (splits[level][ip] != 0 && splits[level][ip] != 1 && splits[level][ip] != 2)
      {
        messerr("For Level=%d/%d and Facies=%d/%d, argument 'splits' is invalid (%d)",
                level + 1, nlevel, ip + 1, _nprop, splits[level][ip]);
        messerr("       It should be either 0, 1 or 2");
        nerr++;
      }
    }
  }
  if (nerr > 0) goto label_error;

  // The first level should only contain 1 and 2

  for (Id ip = 0; ip < _nprop; ip++)
    if (splits[0][ip] != 1 && splits[0][ip] != 2)
    {
      messerr("SPLIT(1,%d) is incorrect (%d)", ip + 1);
      messerr("It should either 1 or 2");
      nerr++;
    }
  if (nerr > 0) goto label_error;

  // Each level must have at least a 1 and a 2

  for (Id level = 0; level < nlevel; level++)
  {
    Id none = 0;
    Id ntwo = 0;
    for (Id ip = 0; ip < _nprop; ip++)
    {
      if (splits[level][ip] == 0) continue;
      if (splits[level][ip] == 1)
        none++;
      else
        ntwo++;
    }
    if (none <= 0 || ntwo <= 0)
    {
      messerr("At level #%d, there must be at least a 1 and a 2",
              level + 1);
      nerr++;
    }
  }
  if (nerr > 0) goto label_error;

  // All non-zero of level L should have same value at level L-1

  for (Id level = 1; level < nlevel; level++)
  {
    Id previous = -1;
    for (Id ip = 0; ip < _nprop; ip++)
    {
      if (splits[level][ip] == 0) continue;
      if (previous < 0)
        previous = splits[level - 1][ip];
      else
      {
        if (previous != splits[level - 1][ip])
        {
          messerr("Non-zero values at level #%d should share same value previous level",
                  level + 1);
          nerr++;
        }
      }
    }
  }
  if (nerr > 0) goto label_error;

  // The argument 'splits' is correct, store it

  _splits = splits;
  return 0;

  // Print the array of splits in case of error

label_error:
  printSplits();
  return 1;
}

/*****************************************************************************/
/*!
**  Internal function to check the means Proportions
**
** \returns Error returned code
**
** \param[in]  meanprops  Array of mean of proportions (Dimension: nfacies)
**
*****************************************************************************/
Id OptimCostColored::_checkMeanProportions(const VectorDouble& meanprops)
{
  if (meanprops.empty()) return 0;
  double total = 0.;
  for (Id ip = 0; ip < _nprop; ip++)
    total += meanprops[ip];
  if (ABS(total - 1.) > 1.e-4)
  {
    messerr("The Proportion Means should add up to 1.\n");
    return 1;
  }

  // The array is valid, store it
  _meanProps = meanprops;
  return 0;
}

/*****************************************************************************/
/*!
**  Internal function to convert Proportions into conditional proportions
**
** \param[in]  level      Level in the Split order
** \param[in]  ip         Rank of the reference proportion
** \param[in]  propfac    Marginal proportions
**
** \param[out] propfacs   Array of facies proportions
**                        (Dimension: nvertex * _nprop)
**
*****************************************************************************/
void OptimCostColored::_copyMultProportions(Id level,
                                            Id ip,
                                            const VectorDouble& propfac,
                                            VectorVectorDouble& propfacs)
{
  auto nvertex = getNVertex();
  Id mode      = _splits[level][ip];
  if (mode == 0) return;

  if (level == 0)
  {
    if (mode == 2)
      for (Id ivert = 0; ivert < nvertex; ivert++)
        propfacs[ip][ivert] = propfac[ivert];
    else
      for (Id ivert = 0; ivert < nvertex; ivert++)
        propfacs[ip][ivert] = 1. - propfac[ivert];
  }
  else
  {
    if (mode == 2)
      for (Id ivert = 0; ivert < nvertex; ivert++)
        propfacs[ip][ivert] *= propfac[ivert];
    else
      for (Id ivert = 0; ivert < nvertex; ivert++)
        propfacs[ip][ivert] *= 1. - propfac[ivert];
  }
}

/**
 * Provides the list of regrouped facies during Optimizaition
 * @param nfacies Number of facies
 * @param verbose Verbose flag
 * @return For each level, Vector of regrouped facies
 */
VectorVectorInt OptimCostColored::initSplit(Id nfacies, bool verbose) const
{
  Id nlevel = nfacies - 1;

  VectorVectorInt splits;
  splits.resize(nlevel, VectorInt(nfacies, 0));

  for (Id ilevel = 0; ilevel < nlevel; ilevel++)
  {
    for (Id ifacies = 0; ifacies < nfacies; ifacies++)
    {
      if (ifacies > nfacies - ilevel - 1) continue;
      if (ifacies == nfacies - ilevel - 1)
        splits[ilevel][ifacies] = 2;
      else
        splits[ilevel][ifacies] = 1;
    }
  }

  if (verbose) printSplits(splits);

  return splits;
}
} // namespace gstlrn

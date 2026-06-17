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
#include "PluriGaussian/DiscretePGS.hpp"
#include "Basic/Law.hpp"
#include "Basic/MathFunc.hpp"

namespace gstlrn
{
  DiscretePGS::DiscretePGS(
    bool flag_cumul,
    Id nconf,
    Id ndisc,
    double cmin,
    double cmax)
    : _nconf(nconf)
    , _ndisc(ndisc)
    , _flagCumul(flag_cumul)
    , _cmin(cmin)
    , _cmax(cmax)
    , _dc(0.)
    , _dp(0.)
    , _v()
    , _res()
  {
    _dc = (_cmax - _cmin) / static_cast<double>(_nconf - 1);
    _dp = 1. / static_cast<double>(_ndisc);

    _res.resize(_nconf);
    for (Id iconf = 0; iconf < _nconf; iconf++) _res[iconf].clear();

    // Define the array of thresholds

    _v.resize(_ndisc + 1);
    _v[0] = THRESH_INF;
    _v[_ndisc] = THRESH_SUP;
    for (Id idisc = 0; idisc < _ndisc; idisc++)
      _v[idisc] = law_invcdf_gaussian(static_cast<double>(idisc) * _dp);
  }

  Id DiscretePGS::_getNelem() const
  {
    return (_flagCumul) ? _ndisc + 1 : _ndisc;
  }

  double DiscretePGS::_getCoval(Id iconf) const
  {
    return (_cmin + iconf * _dc);
  }

  void DiscretePGS::printout(Id flagPrint) const
  {
    Id nelem = _getNelem();

    mestitle(0, "Precalculation of Gaussian integral");
    message("Number of Covariance Discretizations steps  = %d\n", _nconf);
    message("Lower Bound of Covariance Discretization    = %lf\n", _cmin);
    message("Upper Bound of Covariance Discretization    = %lf\n", _cmax);
    message("Covariance Discretization Interval          = %lf\n", _dc);
    if (_flagCumul)
      message("Storing the integral from -Inf to the threshold per class\n");
    else
      message("Storing the integral per discretized class\n");
    message("Covariance is discretized between %lf and %lf\n", _cmin, _cmax);

    message("\n");
    message("Number of Probability Discretizations       = %d\n", _ndisc);
    if (!_v.empty())
      printMatrix(_v, 1, _ndisc + 1, "List of Gaussian Thresholds", 0, 1);

    if (flagPrint > 0)
    {
      mestitle(2, "List of the configurations already calculated");

      for (Id iconf = 0; iconf < _nconf; iconf++)
      {
        if (_res[iconf].empty()) continue;

        if (flagPrint > 0)
          message(
            "- Configuration %d/%d (Cov=%lf)\n", iconf + 1, _nconf,
            _getCoval(iconf));

        if (flagPrint == 2)
          printMatrix(_res[iconf], nelem, nelem, String(), 0, 1);
      }
      message("\n");
    }
  }

  /****************************************************************************
   **
   ** PURPOSE: Returns the rank of the target Ctable
   **
   ** RETURNS: Rank of the Target Ctable Rank
   **
   ** IN_ARGS:  cova    : Covariance value
   **
   ** OUT_ARGS: cround  : Round discretized covariance value
   **
   *****************************************************************************/
  Id DiscretePGS::getCovrank(double cova, double* cround) const
  {
    double ecart = (cova - _cmin);
    auto placeholder = (0.5 + ecart / _dc);
    Id iconf = static_cast<Id>(placeholder);

    if (iconf < 0) iconf = 0;
    if (iconf >= _nconf) iconf = _nconf - 1;

    *cround = _getCoval(iconf);
    return (iconf);
  }

  /****************************************************************************
   **
   ** PURPOSE: Get the starting and ending ranks of the integration
   ** PURPOSE: in the array of thresholds (starting from proportions)
   **
   ** RETURNS: Index for the bound in the discretization table
   **
   ** IN_ARGS:  gaussian : Gaussian bound
   **
   *****************************************************************************/
  Id DiscretePGS::getRankFromProba(double gaussian)
  {
    Id nelem = _getNelem();
    double proba = law_cdf_gaussian(gaussian);

    Id iad = static_cast<Id>(proba / _dp);
    double vmin = _dp * iad;
    double vmax = _dp * (iad + 1);
    if (vmax - proba < proba - vmin) iad++;

    return MIN(nelem, iad);
  }

  /****************************************************************************
   **
   ** PURPOSE: Get the starting and ending ranks of the integration
   ** PURPOSE: in the array of thresholds
   **
   ** IN_ARGS:  low     : Lower bound
   ** IN_ARGS:  up      : Upper bound
   **
   ** OUT_ARGS: indmin  : index of the first discretized inside the interval
   ** OUT_ARGS: indmax  : index of the first discretized outside the interval
   **
   ** REMARKS:  The arguments (imin,imax) are returned so that we can easily
   ** REMARKS:  use them in a loop
   **
   *****************************************************************************/
  void DiscretePGS::_getRanks(double low, double up, Id* indmin, Id* indmax)
  {
    double v1;

    Id nelem = _getNelem();
    *indmin = *indmax = -1;

    // Dispatch

    if (!_flagCumul)
    {

      // Pixelated case

      for (Id idisc = 0; idisc < nelem - 1; idisc++)
      {
        v1 = (_v[idisc] + _v[idisc + 1]) / 2.;
        if (v1 < low) continue;
        if (*indmin < 0) *indmin = idisc;
        if (v1 > up)
        {
          *indmax = idisc;
          break;
        }
      }
      if (*indmax < 0) *indmax = nelem - 1;
    }
    else
    {

      // Cumulative case

      for (Id idisc = 0; idisc < nelem; idisc++)
      {
        v1 = _v[idisc];
        if (v1 < low) continue;
        if (*indmin < 0) *indmin = MAX(0, idisc - 1);
        if (v1 >= up)
        {
          *indmax = MIN(nelem - 1, idisc);
          break;
        }
      }
      if (*indmax < 0) *indmax = nelem;
    }
  }

  /****************************************************************************
   **
   ** PURPOSE: Management of one item of the structure
   **
   ** RETURNS: For Initialization: 1 if actually performed; 0 otherwise
   ** RETURNS: For Freeing: 1 if actually performed; 0 otherwise
   **
   ** IN_ARGS:  mode    : 1 for allocation; -1 for deallocation
   ** IN_ARGS:  rank    : rank within the structure
   **
   ** OUT_ARGS: nb_used : Number of defined items for this discretization level
   **                     0 if the discretization level has not been used
   ** OUT_ARGS: nb_max  : Number of pixels in the Discretized covariance array
   **
   *****************************************************************************/
  void DiscretePGS::_manage(Id mode, Id rank, Id* nb_used, Id* nb_max) const
  {
    Id nelem = _getNelem();
    Id size = nelem * nelem;
    *nb_used = 0;
    *nb_max = size;

    // Dispatch

    if (mode > 0)
    {

      // Allocation

      if (_res[rank].empty())
      {
        _res[rank].resize(size);
        for (Id i = 0; i < size; i++) _res[rank][i] = TEST;
        return;
      }
    }
    else
    {
      if (_res[rank].empty())
      {
        Id number = 0;
        for (Id i = 0; i < size; i++)
          if (FFFF(_res[rank][i])) number++;
        _res[rank].clear();
        *nb_used = number;
        return;
      }
    }
  }

  /****************************************************************************
   **
   ** PURPOSE: Get the value of the Discretized Table for given Configuration
   ** PURPOSE: If not available, the integral is calculated on the fly
   ** PURPOSE: Limitated to the 2-D case
   **
   ** IN_ARGS:  iconf0  : Rank of the Target Ctable Rank
   ** IN_ARGS:  idisc0  : Index along first dimension
   ** IN_ARGS:  jdisc0  : Index along second dimension
   **
   *****************************************************************************/
  double DiscretePGS::_integrate2D(Id iconf0, Id idisc0, Id jdisc0) const
  {
    double lower[2], upper[2], error, value, cova;
    Id infin[2], inform, nelem, iad, nb_used, nb_max;
    static double abseps = 1.e-8;
    static double releps = 0.;
    static Id maxpts = 25000;

    // Check if integral has already been defined

    if (_res[iconf0].empty()) _manage(1, iconf0, &nb_used, &nb_max);

    // Dispatch

    nelem = _getNelem();
    iad = idisc0 + nelem * jdisc0;
    cova = _getCoval(iconf0);

    // Check if the value has already been calculated

    value = _res[iconf0][iad];
    if (!FFFF(value)) return (value);

    // First call to this value, calculate it

    if (!_flagCumul)
    {

      // Pixelated case

      lower[0] = _v[idisc0];
      upper[0] = _v[idisc0 + 1];
      infin[0] = 2;
      if (lower[0] == THRESH_INF) infin[0] = 0;
      if (upper[0] == THRESH_SUP) infin[0] = 1;

      lower[1] = _v[jdisc0];
      upper[1] = _v[jdisc0 + 1];
      infin[1] = 2;
      if (lower[1] == THRESH_INF) infin[1] = 0;
      if (upper[1] == THRESH_SUP) infin[1] = 1;

      mvndst(
        2, lower, upper, infin, &cova, maxpts, abseps, releps, &error, &value,
        &inform);
      if (inform) messageAbort("Error in function 'mvndst'");
    }
    else
    {

      // Cumulative case

      lower[0] = THRESH_INF;
      lower[1] = THRESH_INF;

      upper[0] = _v[idisc0];
      infin[0] = 0;
      if (upper[0] == THRESH_SUP) infin[0] = 1;

      upper[1] = _v[jdisc0];
      infin[1] = 0;
      if (upper[1] == THRESH_SUP) infin[1] = 1;

      mvndst(
        2, lower, upper, infin, &cova, maxpts, abseps, releps, &error, &value,
        &inform);
      if (inform) messageAbort("Error in function 'mvndst'");
    }

    // Store it in the table

    _res[iconf0][iad] = value;
    return (value);
  }

  /****************************************************************************
   **
   ** PURPOSE: Get the value of the Discretized Table for given Configuration
   ** PURPOSE: If not available, the integral is calculated on the fly
   ** PURPOSE: Limitated to the 3-D case
   **
   ** IN_ARGS:  iconf0  : Rank of the Target Ctable Rank
   ** IN_ARGS:  idisc0  : Index along first dimension
   ** IN_ARGS:  jdisc0  : Index along second dimension
   ** IN_ARGS:  kdisc0  : Index along third dimension
   **
   *****************************************************************************/
  double
    DiscretePGS::_integrate3D(Id iconf0, Id idisc0, Id jdisc0, Id kdisc0) const
  {
    double lower[3], upper[3], error, value, cova;
    Id infin[3], inform, nelem, iad, nb_used, nb_max;
    static double abseps = 1.e-8;
    static double releps = 0.;
    static Id maxpts = 25000;

    // Check if integral has already been defined

    if (_res[iconf0].empty()) _manage(1, iconf0, &nb_used, &nb_max);

    // Dispatch

    nelem = _getNelem();
    iad = idisc0 + nelem * (jdisc0 + nelem * kdisc0);
    cova = _getCoval(iconf0);

    // Check if the value has already been calculated

    value = _res[iconf0][iad];
    if (!FFFF(value)) return (value);

    // First call to this value, calculate it

    if (!_flagCumul)
    {

      // Pixelated case

      lower[0] = _v[idisc0];
      upper[0] = _v[idisc0 + 1];
      infin[0] = 2;
      if (lower[0] == THRESH_INF) infin[0] = 0;
      if (upper[0] == THRESH_SUP) infin[0] = 1;

      lower[1] = _v[jdisc0];
      upper[1] = _v[jdisc0 + 1];
      infin[1] = 2;
      if (lower[1] == THRESH_INF) infin[1] = 0;
      if (upper[1] == THRESH_SUP) infin[1] = 1;

      lower[2] = _v[kdisc0];
      upper[2] = _v[kdisc0 + 1];
      infin[2] = 2;
      if (lower[2] == THRESH_INF) infin[2] = 0;
      if (upper[2] == THRESH_SUP) infin[2] = 1;

      mvndst(
        3, lower, upper, infin, &cova, maxpts, abseps, releps, &error, &value,
        &inform);
      if (inform) messageAbort("Error in function 'mvndst'");
    }
    else
    {

      // Cumulative case

      upper[0] = _v[idisc0];
      infin[0] = 0;
      if (upper[0] == THRESH_SUP) infin[0] = 1;

      upper[1] = _v[jdisc0];
      infin[1] = 0;
      if (upper[1] == THRESH_SUP) infin[1] = 1;

      upper[2] = _v[kdisc0];
      infin[2] = 0;
      if (upper[2] == THRESH_SUP) infin[2] = 1;

      mvndst(
        3, lower, upper, infin, &cova, maxpts, abseps, releps, &error, &value,
        &inform);
      if (inform) messageAbort("Error in function 'mvndst'");
    }
    // Store it Id the table

    _res[iconf0][iad] = value;
    return (value);
  }

  /****************************************************************************
  **
  ** PURPOSE: Evaluate the integral over a rectangle from the integral
  ** PURPOSE: values calculated at pixel size
  **
  ** RETURNS: Integral value
  **
  ** IN_ARGS:  iconf0  : Rank of the relevant internal structure
  ** IN_ARGS:  lows    : Array of lower values
  ** IN_ARGS:  ups     : Array of upper values
  **
  ** REMARKS: The code has been designed for boundaries of facies areas
  ** REMAKRS: parallel to the main axes
  **
  *****************************************************************************/
  double
    DiscretePGS::calculate(Id iconf0, const double lows[2], const double ups[2])
  {
    Id i1min, i1max, i2min, i2max;
    double result = 0;

    // Dispatch

    if (!_flagCumul)
    {

      // Pixelated case

      _getRanks(lows[0], ups[0], &i1min, &i1max);
      _getRanks(lows[1], ups[1], &i2min, &i2max);
      for (Id idisc = i1min; idisc < i1max; idisc++)
        for (Id jdisc = i2min; jdisc < i2max; jdisc++)
          result += _integrate2D(iconf0, idisc, jdisc);
    }
    else
    {

      // Cumulative case

      _getRanks(lows[0], ups[0], &i1min, &i1max);
      _getRanks(lows[1], ups[1], &i2min, &i2max);
      result =
        (_integrate2D(iconf0, i1max, i2max) - _integrate2D(iconf0, i1min, i2max)
         - _integrate2D(iconf0, i1max, i2min)
         + _integrate2D(iconf0, i1min, i2min));
    }
    return (result);
  }

  /****************************************************************************
  **
  ** PURPOSE: Evaluate the integral over a rectangle from the integral
  ** PURPOSE: values calculated at pixel size
  **
  ** RETURNS: Integral value
  **
  ** IN_ARGS:  rklows  : Array of indices of lower values
  ** IN_ARGS:  rkups   : Array of indices upper values
  **
  ** REMARKS: The arguments 'rklows' and ;rkups' are double although
  ** REMARKS: they designate ranks.
  **
  *****************************************************************************/
  double DiscretePGS::calculateByRank(
    Id iconf0,
    const double rklows[2],
    const double rkups[2]) const
  {
    double result = 0;

    // Dispatch

    if (!_flagCumul)
    {

      // Pixelated case

      for (Id idisc = static_cast<Id>(rklows[0]);
           idisc < static_cast<Id>(rkups[0]); idisc++)
        for (Id jdisc = static_cast<Id>(rklows[1]);
             jdisc < static_cast<Id>(rkups[1]); jdisc++)
          result += _integrate2D(iconf0, idisc, jdisc);
    }
    else
    {

      // Cumulative case

      result =
        (_integrate2D(
           iconf0, static_cast<Id>(rkups[0]), static_cast<Id>(rkups[1]))
         - _integrate2D(
           iconf0, static_cast<Id>(rklows[0]), static_cast<Id>(rkups[1]))
         - _integrate2D(
           iconf0, static_cast<Id>(rkups[0]), static_cast<Id>(rklows[1]))
         + _integrate2D(
           iconf0, static_cast<Id>(rklows[0]), static_cast<Id>(rklows[1])));
    }
    return (result);
  }

} // namespace gstlrn

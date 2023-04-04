/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "geoslib_old_f.h"

#include "Basic/NamingConvention.hpp"
#include "Estimation/CalcSimpleInterpolation.hpp"
#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"

CalcSimpleInterpolation::CalcSimpleInterpolation()
    : ACalcInterpolator(),
      _iattOut(-1),
      _flagMovAve(false),
      _flagInvDist(false),
      _flagLstSqr(false),
      _exponent(2.),
      _flagExpand(true),
      _dmax(TEST),
      _order(0)
{
}

CalcSimpleInterpolation::~CalcSimpleInterpolation()
{
}

int CalcSimpleInterpolation::_getNVar() const
{
  return getDbin()->getLocNumber(ELoc::Z);
}

bool CalcSimpleInterpolation::_check()
{
  if (! ACalcInterpolator::_check()) return false;

  if (! hasDbin()) return false;
  if (! hasDbout()) return false;

  if (_getNVar() != 1)
  {
    messerr("These methods are restricted to the Monovariate case");
    return false;
  }

  if (_flagMovAve)
  {
    if (! hasNeighParam()) return false;
  }
  if (_flagLstSqr)
  {
    if (! hasNeighParam()) return false;
  }

  return true;
}

bool CalcSimpleInterpolation::_preprocess()
{
  _iattOut = _addVariableDb(2, 1, ELoc::UNKNOWN, 0, 1, 0.);
  if (_iattOut < 0) return false;
  return true;
}

bool CalcSimpleInterpolation::_postprocess()
{
  /* Free the temporary variables */
  _cleanVariableDb(2);

  _renameVariable(2, 1, _iattOut, String(), 1);
  return true;
}

void CalcSimpleInterpolation::_rollback()
{
  _cleanVariableDb(1);
}

/****************************************************************************/
/*!
 **  Standard Kriging
 **
 ** \return  Error return code
 **
 *****************************************************************************/
bool CalcSimpleInterpolation::_run()
{
  if (_flagMovAve)
  {
    if (movave(getDbin(), getDbout(), getNeighparam(), _iattOut))
      return false;
  }

  if (_flagLstSqr)
   {
     if (lstsqr(getDbin(), getDbout(), getNeighparam(), _iattOut, _order))
       return false;
   }

  if (_flagInvDist)
  {
    if (invdist(getDbin(), getDbout(), _iattOut, _exponent, _flagExpand, _dmax))
      return false;
  }

  return true;
}

/****************************************************************************/
/*!
 **  Inverse distance estimation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  exponent    exponent of the inverse distance
 ** \param[in]  flag_expand True for expansion option (if dbin is Grid)
 ** \param[in]  dmax        Maximum search radius (if dbin is Points)
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int inverseDistance(Db *dbin,
                    Db *dbout,
                    double exponent,
                    bool flag_expand,
                    double dmax,
                    const NamingConvention &namconv)
{
  CalcSimpleInterpolation interpol;
  interpol.setDbin(dbin);
  interpol.setDbout(dbout);
  interpol.setNamingConvention(namconv);

  interpol.setFlagInvDist(true);
  interpol.setExponent(exponent);
  interpol.setFlagExpand(flag_expand);
  interpol.setDmax(dmax);

  // Run the calculator
  int error = (interpol.run()) ? 0 : 1;
  return error;
}

/****************************************************************************/
/*!
 **  Inverse distance estimation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  neighparam  ANeighParam structure
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
GSTLEARN_EXPORT int movingAverage(Db *dbin,
                                  Db *dbout,
                                  ANeighParam *neighparam,
                                  const NamingConvention &namconv)
{
  CalcSimpleInterpolation interpol;
  interpol.setDbin(dbin);
  interpol.setDbout(dbout);
  interpol.setNeighparam(neighparam);
  interpol.setNamingConvention(namconv);

  interpol.setFlagMovAve(true);

  // Run the calculator
  int error = (interpol.run()) ? 0 : 1;
  return error;
}

/****************************************************************************/
/*!
 **  Polynomial estimation using Least Squares
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  neighparam  ANeighParam structure
 ** \param[in]  order       Order of the polynomial
 ** \param[in]  namconv     Naming Convention
 **
 *****************************************************************************/
int leastSquares(Db *dbin,
                 Db *dbout,
                 ANeighParam *neighparam,
                 int order,
                 const NamingConvention &namconv)
{
  CalcSimpleInterpolation interpol;
  interpol.setDbin(dbin);
  interpol.setDbout(dbout);
  interpol.setNeighparam(neighparam);
  interpol.setNamingConvention(namconv);

  interpol.setFlagLstSqr(true);
  interpol.setOrder(order);

  // Run the calculator
  int error = (interpol.run()) ? 0 : 1;
  return error;

}


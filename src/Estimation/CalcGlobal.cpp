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
#include "geoslib_f.h"
#include "geoslib_old_f.h"

#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"
#include "Estimation/CalcGlobal.hpp"
#include "Estimation/KrigingSystem.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Basic/OptDbg.hpp"

#include <math.h>

CalcGlobal::CalcGlobal(int ivar0,
                       bool verbose)
    : ACalcInterpolator(),
      _flagArithmetic(false),
      _flagKriging(false),
      _ivar0(ivar0),
      _verbose(verbose)
{
}

CalcGlobal::~CalcGlobal()
{
}

bool CalcGlobal::_check()
{
  if (! ACalcInterpolator::_check()) return false;

  if (! hasDbin()) return false;
  if (! hasDbout()) return false;
  if (! hasModel()) return false;

  if (_flagArithmetic)
  {
    if (! getDbout()->isGrid())
    {
      messerr("'dbout'  must be a grid for Arithmetic Global estimation");
      return false;
    }
  }
  if (_ivar0 < 0 || _ivar0 >= getDbin()->getLocNumber(ELoc::Z))
  {
    messerr("The target variable (%d) must lie between 1 and the number of variables (%d)",
            _ivar0 + 1, getDbin()->getLocNumber(ELoc::Z));
    return false;
  }

  return true;
}

bool CalcGlobal::_preprocess()
{
  if (!ACalcInterpolator::_preprocess()) return false;

  return true;
}

bool CalcGlobal::_postprocess()
{
  _cleanVariableDb(2);

  return true;
}

void CalcGlobal::_rollback()
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
bool CalcGlobal::_run()
{
  if (_flagArithmetic)
  {
    if (_globalArithmetic()) return false;
  }

  if (_flagKriging)
  {
    if (_globalKriging()) return false;
  }
  return true;
}

int CalcGlobal::_globalKriging()
{
  VectorDouble rhsCum;

  // Initializations

  int ndim = getDbin()->getNDim();
  int nvar = getModel()->getVariableNumber();
  SpaceRN space(ndim);
  NeighUnique neighU = NeighUnique(false, &space);
  neighU.attach(getDbin(), getDbout());

  /* Setting options */

  KrigingSystem ksys(getDbin(), getDbout(), getModel(), &neighU);
  if (ksys.setKrigOptFlagGlobal(true)) return 1;
  if (! ksys.isReady()) return 1;

  /* Loop on the targets to be processed */

  int ng = 0;
  for (int iech_out = 0; iech_out < getDbout()->getSampleNumber(); iech_out++)
  {
    mes_process("Kriging sample", getDbout()->getSampleNumber(), iech_out);
    if (! getDbout()->isActive(iech_out)) continue;
    if (ksys.estimate(iech_out)) return 1;

    // Cumulate the R.H.S.

    VectorDouble rhs = ksys.getRHSC(_ivar0);
    if (rhsCum.empty()) rhsCum.resize(rhs.size(),0.);
    VH::addInPlace(rhsCum, rhs);
    ng++;
  }

  ksys.conclusion();

  /* Preliminary checks */

  int ntot = getDbin()->getSampleNumber(false);
  int np   = getDbin()->getSampleNumber(true);
  double cell = 1.;
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
  if (dbgrid != nullptr) cell = dbgrid->getCellSize();
  double surface = ng * cell;

  /* Average covariance over the territory */

  double cvv = model_cxx(getModel(), getDbout(), getDbout(), _ivar0, _ivar0, 0,
                         db_epsilon_distance(getDbin()));

  /* Load the scaled cumulated R.H.S. in the array rhs */

  VH::divideConstant(rhsCum, (double) ng);

  /* Derive the kriging weights */

  int nred = ksys.getNRed();
  VectorDouble lhsinv = ksys.getLHSInv();
  VectorDouble zam = ksys.getZamC();
  VectorDouble wgt(nred);
  matrix_product_safe(nred, nred, nvar, lhsinv.data(), rhsCum.data(), wgt.data());

  /* Perform the estimation */

  double estim = VH::innerProduct(rhsCum, zam);
  double stdv = cvv - VH::innerProduct(rhsCum, wgt);
  stdv = (stdv > 0) ? sqrt(stdv) : 0.;
  double cvgeo = (estim == 0. || FFFF(estim)) ? TEST : stdv / estim;

  /* Store the results in the output Global_Result struture */

  _gRes.ntot = ntot;
  _gRes.np = np;
  _gRes.ng = ng;
  _gRes.surface = surface;
  _gRes.zest = estim;
  _gRes.sse = stdv;
  _gRes.cvgeo = cvgeo;
  _gRes.cvv = cvv;
  _gRes.weights = wgt;

  /* Printout */

  if (_verbose)
  {
    mestitle(1,"Global estimation kriging");
    message("Total number of data             = %d\n", ntot);
    message("Number of active data            = %d\n", np);
    message("Number of variables              = %d\n", nvar);
    message("Cvv                              = %lf\n", cvv);
    if (FFFF(estim))
      message("Estimation by kriging            = NA\n");
    else
      message("Estimation by kriging            = %lf\n", estim);
    message("Estimation St. Dev. of the mean  = %lf\n", stdv);
    if (FFFF(cvgeo))
      message("CVgeo                            = NA\n");
    else
      message("CVgeo                            = %lf\n", cvgeo);
    message("Surface                          = %lf\n", surface);
    if (FFFF(estim))
      message("Q (Estimation * Surface)         = NA\n");
    else
      message("Q (Estimation * Surface)         = %lf\n", estim * surface);
    message("\n");
  }
  return 0;
}

int CalcGlobal::_globalArithmetic()
{
  /* Preliminary assignments */

  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
  int ntot = getDbin()->getSampleNumber(false);
  int np = getDbin()->getSampleNumber(true);
  int ng = dbgrid->getSampleNumber(true);
  double surface = ng * dbgrid->getCellSize();

  /* Average covariance over the data */

  double cxx = model_cxx(getModel(), getDbin(), getDbin(), _ivar0, _ivar0, 0, 0.);

  /* Average covariance between the data and the territory */

  double cxv = model_cxx(getModel(), getDbin(), dbgrid, _ivar0, _ivar0, 0, 0.);

  /* Average covariance over the territory */

  double cvv = model_cxx(getModel(), dbgrid, dbgrid, _ivar0, _ivar0, 0,
                         db_epsilon_distance(dbgrid));

  /* Calculating basic statistics */

  int iatt = db_attribute_identify(getDbin(), ELoc::Z, _ivar0);
  double wtot;
  double ave;
  double var;
  double mini;
  double maxi;
  db_monostat(getDbin(), iatt, &wtot, &ave, &var, &mini, &maxi);

  /* Filling the resulting structure */

  double sse = cvv - 2. * cxv + cxx;
  sse = (sse > 0) ? sqrt(sse) : 0.;
  double cvsam = (ave != 0.) ? sqrt(var) / ave : TEST;
  double cviid = cvsam / sqrt(np);
  double cvgeo = (ave != 0.) ? sse / ave : TEST;

  /* Filling the output structure */

  _gRes.ntot = ntot;
  _gRes.np = np;
  _gRes.ng = ng;
  _gRes.surface = surface;
  _gRes.zest = ave;
  _gRes.sse  = sse;
  _gRes.cvgeo = cvgeo;
  _gRes.cvv = cvv;
  _gRes.weights.resize(np, 1./np);

  if (_verbose)
  {
    mestitle(1,"Global estimation by arithmetic average");
    message("Total number of data             = %d\n", ntot);
    message("Number of active data            = %d\n", np);
    message("Sample variance                  = %lf\n", var);
    message("CVsample                         = %lf\n", cvsam);
    message("CViid                            = %lf\n", cviid);
    message("Cxx                              = %lf\n", cxx);
    message("Cxv                              = %lf\n", cxv);
    message("Cvv                              = %lf\n", cvv);
    if (FFFF(ave))
      message("Estimation by arithmetic average = NA\n");
    else
      message("Estimation by arithmetic average = %lf\n", ave);
    message("Estimation St. dev. of the mean  = %lf\n", sse);
    if (FFFF(cvgeo))
      message("CVgeo                            = NA\n");
    else
      message("CVgeo                            = %lf\n", cvgeo);
    message("Surface                          = %lf\n", surface);
    if (FFFF(ave))
      message("Q (Estimation * Surface)         = NA\n");
    else
      message("Q (Estimation * Surface)         = %lf\n", ave * surface);
    message("\n");
  }

  return 0;
}

Global_Result global_arithmetic(Db *dbin,
                                DbGrid *dbgrid,
                                Model *model,
                                int ivar0,
                                bool verbose)
{
  Global_Result gres;
  CalcGlobal global(ivar0, verbose);
  global.setDbin(dbin);
  global.setDbout(dbgrid);
  global.setModel(model);
  global.setFlagArithmetic(true);

  if (global.run())
    gres = global.getGRes();
  return gres;
}

Global_Result global_kriging(Db *dbin,
                            Db *dbout,
                            Model *model,
                            int ivar0,
                            bool verbose)
{
  Global_Result gres;
  CalcGlobal global(ivar0, verbose);
  global.setDbin(dbin);
  global.setDbout(dbout);
  global.setModel(model);
  global.setFlagKriging(true);

  if (global.run())
    gres = global.getGRes();
  return gres;
}


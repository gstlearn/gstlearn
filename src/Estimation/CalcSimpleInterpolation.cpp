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
#include "geoslib_old_f.h"

#include "Basic/NamingConvention.hpp"
#include "Basic/OptDbg.hpp"
#include "Estimation/CalcSimpleInterpolation.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Model/Model.hpp"

CalcSimpleInterpolation::CalcSimpleInterpolation()
    : ACalcInterpolator(),
      _flagEst(true),
      _flagStd(false),
      _iattEst(-1),
      _iattStd(-1),
      _flagMovAve(false),
      _flagMovMed(false),
      _flagInvDist(false),
      _flagLstSqr(false),
      _flagNearest(false),
      _exponent(2.),
      _flagExpand(true),
      _dmax(TEST),
      _order(0)
{
}

CalcSimpleInterpolation::~CalcSimpleInterpolation()
{
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

  // Some simple interpolators require a Neighborhood
  if (_flagMovAve || _flagMovMed || _flagLstSqr)
  {
    if (! hasNeigh()) return false;
  }

  // Model is required if the variance is demanded
  if (_flagStd)
  {
    if (! hasModel())
    {
      messerr("A Model is required for calculation of option 'St. Dev.'");
      return false;
    }
  }
  return true;
}

bool CalcSimpleInterpolation::_preprocess()
{
  if (!ACalcInterpolator::_preprocess()) return false;
  
  if (_flagEst)
  {
    _iattEst = _addVariableDb(2, 1, ELoc::UNKNOWN, 0, 1, 0.);
    if (_iattEst < 0) return false;
  }
  if (_flagStd)
  {
    _iattStd = _addVariableDb(2, 1, ELoc::UNKNOWN, 0, 1, 0.);
    if (_iattStd < 0) return false;
  }
  return true;
}

bool CalcSimpleInterpolation::_postprocess()
{
  /* Free the temporary variables */
  _cleanVariableDb(2);

  _renameVariable(2, VectorString(), ELoc::Z, 1, _iattEst, "estim", 1);
  _renameVariable(2, VectorString(), ELoc::Z, 1, _iattStd, "stdev", 1);
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
    if (_movave(getDbin(), getDbout(), getNeigh())) return false;
  }

  if (_flagMovMed)
  {
    if (_movmed(getDbin(), getDbout(), getNeigh())) return false;
  }

  if (_flagLstSqr)
   {
     if (_lstsqr(getDbin(), getDbout(), getNeigh())) return false;
   }

  if (_flagInvDist)
  {
    if (_invdist(getDbin(), getDbout())) return false;
  }

  if (_flagNearest)
  {
    if (_nearest(getDbin(), getDbout(), getNeigh())) return false;
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
 ** \param[in]  flag_est    True if the estimation must be calculated
 ** \param[in]  flag_std    True if the St. Dev. must be calculated
 ** \param[in]  model       Model structure (used for St. Dev.)
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int inverseDistance(Db *dbin,
                    Db *dbout,
                    double exponent,
                    bool flag_expand,
                    double dmax,
                    bool flag_est,
                    bool flag_std,
                    Model* model,
                    const NamingConvention &namconv)
{
  CalcSimpleInterpolation interpol;
  interpol.setDbin(dbin);
  interpol.setDbout(dbout);
  interpol.setFlagEst(flag_est);
  interpol.setFlagStd(flag_std);
  interpol.setModel(model);
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
 **  Moving Average estimation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  neigh       ANeigh structure
 ** \param[in]  flag_est    True if the estimation must be calculated
 ** \param[in]  flag_std    True if the St. Dev. must be calculated
 ** \param[in]  model       Model structure (used for St. Dev.)
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
GSTLEARN_EXPORT int movingAverage(Db *dbin,
                                  Db *dbout,
                                  ANeigh *neigh,
                                  bool flag_est,
                                  bool flag_std,
                                  Model *model,
                                  const NamingConvention &namconv)
{
  CalcSimpleInterpolation interpol;
  interpol.setDbin(dbin);
  interpol.setDbout(dbout);
  interpol.setNeigh(neigh);
  interpol.setFlagEst(flag_est);
  interpol.setFlagStd(flag_std);
  interpol.setModel(model);
  interpol.setNamingConvention(namconv);

  interpol.setFlagMovAve(true);

  // Run the calculator
  int error = (interpol.run()) ? 0 : 1;
  return error;
}

/****************************************************************************/
/*!
 **  Moving Median estimation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  neigh       ANeigh structure
 ** \param[in]  flag_est    True if the estimation must be calculated
 ** \param[in]  flag_std    True if the St. Dev. must be calculated
 ** \param[in]  model       Model structure (used for St. Dev.)
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
GSTLEARN_EXPORT int movingMedian(Db *dbin,
                                 Db *dbout,
                                 ANeigh *neigh,
                                 bool flag_est,
                                 bool flag_std,
                                 Model* model,
                                 const NamingConvention &namconv)
{
  CalcSimpleInterpolation interpol;
  interpol.setDbin(dbin);
  interpol.setDbout(dbout);
  interpol.setNeigh(neigh);
  interpol.setFlagEst(flag_est);
  interpol.setFlagStd(flag_std);
  interpol.setModel(model);
  interpol.setNamingConvention(namconv);

  interpol.setFlagMovMed(true);

  // Run the calculator
  int error = (interpol.run()) ? 0 : 1;
  return error;
}

/****************************************************************************/
/*!
 **  Nearest Neighbor estimation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  flag_est    True if the estimation must be calculated
 ** \param[in]  flag_std    True if the St. Dev. must be calculated
 ** \param[in]  model       Model structure (used for St. Dev.)
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
GSTLEARN_EXPORT int nearestNeighbor(Db *dbin,
                                    Db *dbout,
                                    bool flag_est,
                                    bool flag_std,
                                    Model* model,
                                    const NamingConvention &namconv)
{
  CalcSimpleInterpolation interpol;
  interpol.setDbin(dbin);
  interpol.setDbout(dbout);
  interpol.setFlagEst(flag_est);
  interpol.setFlagStd(flag_std);
  interpol.setModel(model);

  NeighMoving neighM(false, 1, 1.e6);
  interpol.setNeigh(&neighM);
  interpol.setNamingConvention(namconv);

  interpol.setFlagNearest(true);

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
 ** \param[in]  neigh       ANeigh structure
 ** \param[in]  order       Order of the polynomial
 ** \param[in]  namconv     Naming Convention
 **
 *****************************************************************************/
GSTLEARN_EXPORT int leastSquares(Db *dbin,
                                 Db *dbout,
                                 ANeigh *neigh,
                                 int order,
                                 const NamingConvention &namconv)
{
  CalcSimpleInterpolation interpol;
  interpol.setDbin(dbin);
  interpol.setDbout(dbout);
  interpol.setNeigh(neigh);
  interpol.setNamingConvention(namconv);

  interpol.setFlagLstSqr(true);
  interpol.setOrder(order);

  // Run the calculator
  int error = (interpol.run()) ? 0 : 1;
  return error;

}

/****************************************************************************/
/*!
 **  Nearest Neighbour estimation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  neigh       ANeigh structure
 **
 *****************************************************************************/
int CalcSimpleInterpolation::_nearest(Db *dbin,
                                      Db *dbout,
                                      ANeigh *neigh)
{
  VectorInt nbgh;

  /* Loop on the targets to be processed */

  for (int iech = 0; iech < dbout->getSampleNumber(); iech++)
  {
    mes_process("Estimation by Nearest Neighbor", dbout->getSampleNumber(), iech);
    if (!dbout->isActive(iech)) continue;
    if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH)
        || OptDbg::query(EDbg::RESULTS))
    {
      mestitle(1, "Target location");
      db_sample_print(dbout, iech, 1, 0, 0, 0);
    }
    VectorDouble weights;

    // Find the neighborhood
    neigh->select(iech, nbgh);

    // Perform the estimation
    if (nbgh.size() > 0)
    {
      nbgh.resize(1);
      weights.push_back(1.);
    }

    // Save the results
    _saveResults(dbin, dbout, nbgh, iech, weights);
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Moving Average estimation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  neigh       ANeigh structure
 **
 *****************************************************************************/
int CalcSimpleInterpolation::_movave(Db* dbin, Db* dbout, ANeigh* neigh)
{
  VectorInt nbgh;

  /* Loop on the targets to be processed */

  for (int iech = 0; iech < dbout->getSampleNumber(); iech++)
  {
    mes_process("Estimation by Moving Average", dbout->getSampleNumber(), iech);
    if (!dbout->isActive(iech)) continue;
    OptDbg::setCurrentIndex(iech + 1);
    if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH)
        || OptDbg::query(EDbg::RESULTS))
    {
      mestitle(1, "Target location");
      db_sample_print(dbout, iech, 1, 0, 0, 0);
    }
    VectorDouble weights;

    // Find the neighborhood
    neigh->select(iech, nbgh);

    // Perform the estimation
    for (int i = 0; i < (int) nbgh.size(); i++)
    {
      double value = dbin->getZVariable( nbgh[i], 0);
      if (FFFF(value))
      {
        nbgh.clear();
        weights.clear();
        break;
      }
      weights.push_back(1.);
    }

    // Save the results
    _saveResults(dbin, dbout, nbgh, iech, weights);
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Moving Median estimation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  neigh       ANeigh structure
 **
 *****************************************************************************/
int CalcSimpleInterpolation::_movmed(Db* dbin, Db* dbout, ANeigh* neigh)
{
  VectorInt nbgh;

  /* Loop on the targets to be processed */

   for (int iech = 0; iech < dbout->getSampleNumber(); iech++)
   {
     mes_process("Estimation by Moving Median", dbout->getSampleNumber(), iech);
     if (!dbout->isActive(iech)) continue;
     OptDbg::setCurrentIndex(iech + 1);
     if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
     {
       mestitle(1, "Target location");
       db_sample_print(dbout, iech, 1, 0, 0, 0);
     }
     VectorDouble weights;
     VectorInt nbghmed;

     // Find the neighborhood
     neigh->select(iech, nbgh);

     // Perform the estimation
     if (nbgh.size() > 0)
     {
       int rank = (int) nbgh.size() / 2;
       nbghmed.push_back(nbgh[rank]);
       weights.push_back(1.);
     }

     // Save the results
     _saveResults(dbin, dbout, nbghmed, iech, weights);
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Polynomial estimation using Least Squares
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  neigh       ANeigh structure
 **
 *****************************************************************************/
int CalcSimpleInterpolation::_lstsqr(Db* dbin, Db* dbout, ANeigh* neigh) const
{
  int ndim = dbin->getNDim();
  VectorInt nbgh;
  CovContext ctxt(1, ndim);
  const DriftList* drft = DriftFactory::createDriftListFromIRF(_order, 0, ctxt);
  int ndrift = drft->getDriftNumber();
  VectorDouble X(ndrift);
  VectorDouble B(ndrift);
  MatrixSquareSymmetric A(ndrift);

  /* Loop on the targets to be processed */

   for (int iech = 0; iech < dbout->getSampleNumber(); iech++)
   {
     mes_process("Estimation by Inverse distance", dbout->getSampleNumber(), iech);
     OptDbg::setCurrentIndex(iech + 1);
     if (!dbout->isActive(iech)) continue;
     if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
     {
       mestitle(1, "Target location");
       db_sample_print(dbout, iech, 1, 0, 0, 0);
     }

     // Find the neighborhood
     neigh->select(iech, nbgh);
     int nSize = (int) nbgh.size();
     if (nSize < ndrift)
     {
       dbout->setArray(iech, _iattEst, TEST);
       continue;
     }

     // Evaluate the least square system
     A.fill(0.);
     for (int i = 0; i < ndrift; i++) B[i] = 0.;
     for (int jech = 0; jech < nSize; jech++)
     {
       int jech1 = nbgh[jech];
       double zval = dbin->getZVariable(jech1, 0);
       if (FFFF(zval)) continue;
       VectorDouble Vdata = drft->evalDriftBySample(dbin, jech1);

       // Double loop on the drift terms
       for (int id1 = 0; id1 < ndrift; id1++)
       {
         B[id1] += zval * Vdata[id1];
         for (int id2 = 0; id2 <= id1; id2++)
         {
           A.addValue(id1,  id2, Vdata[id1] * Vdata[id2]);
         }
       }
     }

     // Solve the system
     if (A.solve(B, X) > 0) continue;

     // Evaluate the vector of drift terms at target
     VectorDouble Vtarget = drft->evalDriftBySample(dbout,  iech);

     // Perform the estimation
     double result = VH::innerProduct(X, Vtarget);

     // Assign the result
     dbout->setArray(iech, _iattEst, result);
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Inverse distance estimation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 **
 *****************************************************************************/
int CalcSimpleInterpolation::_invdist(Db *dbin, Db *dbout)
{
  if (! dbin->isGrid())
  {
    _pointInvdist(dbin, dbout);
  }
  else
  {
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(dbin);
    _gridInvdist(dbgrid, dbout);
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Inverse distance estimation when Input DB is a point file
 **
 ** \param[in]  dbin        Input Db
 ** \param[in]  dbout       Output Db
 **
 *****************************************************************************/
void CalcSimpleInterpolation::_pointInvdist(Db *dbin, Db *dbout)
{
  int ndim = dbin->getNDim();
  double dmin = dbout->getExtensionDiagonal() / 1.e5;
  VectorDouble coor(ndim);
  VectorDouble cooref(ndim);

  /* Loop on the targets to be processed */

  for (int iech = 0; iech < dbout->getSampleNumber(); iech++)
  {
    mes_process("Estimation by Inverse distance", dbout->getSampleNumber(), iech);
    OptDbg::setCurrentIndex(iech + 1);
    if (!dbout->isActive(iech)) continue;
    if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
    {
      mestitle(1, "Target location");
      db_sample_print(dbout, iech, 1, 0, 0, 0);
    }
    VectorInt nbgh;
    VectorDouble weights;
    dbout->getCoordinatesPerSampleInPlace(iech, cooref);

    /* Loop on the data points */

    for (int iech_in = 0; iech_in < dbin->getSampleNumber(); iech_in++)
    {
      if (!dbin->isActive(iech_in)) continue;
      dbin->getCoordinatesPerSampleInPlace(iech_in, coor);
      double val_neigh = dbin->getZVariable(iech_in, 0);
      if (FFFF(val_neigh)) continue;

      /* Check that the data point is a valid neighbor */

      double dist = ut_distance(ndim, coor.data(), cooref.data());
      if (!FFFF(_dmax) && dist > _dmax) continue;
      if (dist < dmin)
      {
        nbgh.clear();
        weights.clear();
        nbgh.push_back(iech_in);
        weights.push_back(1.);
        break;
      }
      double wgt = 1. / pow(dist, _exponent);
      nbgh.push_back(iech_in);
      weights.push_back(wgt);
    }

    // Save the results
    _saveResults(dbin, dbout, nbgh, iech, weights);
  }
}

/****************************************************************************/
/*!
 **  Inverse distance estimation when Input DB is a grid
 **
 ** \param[in]  dbin        Input Db
 ** \param[in]  dbout       Output Db
 **
 *****************************************************************************/
void CalcSimpleInterpolation::_gridInvdist(DbGrid *dbin, Db *dbout)
{
  int ndim = dbin->getNDim();
  int maxneigh = (int) pow(2., (double) ndim);
  double dmin = dbout->getExtensionDiagonal() / 1.e5;

  VectorDouble coor(ndim);
  VectorDouble cooref(ndim);
  VectorDouble percent(ndim);
  VectorInt indg(ndim);
  VectorInt indref(ndim);

  /* Loop on the targets to be processed */

  for (int iech = 0; iech < dbout->getSampleNumber(); iech++)
  {
    mes_process("Estimation by Inverse distance", dbout->getSampleNumber(), iech);
    OptDbg::setCurrentIndex(iech + 1);
    if (!dbout->isActive(iech)) continue;
    if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
    {
      mestitle(1, "Target location");
      db_sample_print(dbout, iech, 1, 0, 0, 0);
    }

    /* Find the grid index corresponding to the target */

    dbout->getCoordinatesPerSampleInPlace(iech, cooref);
    if (dbin->coordinateToIndicesInPlace(cooref, indref))
    {
      if (_flagEst)
        dbout->setArray(iech, _iattEst, TEST);
      if (_flagStd)
        dbout->setArray(iech, _iattStd, TEST);
      continue;
    }

    /* Loop on the neighbors */

    VectorInt nbgh;
    VectorDouble weights;
    for (int rank = 0; rank < maxneigh; rank++)
    {
      for (int idim = 0; idim < ndim; idim++) indg[idim] = indref[idim];

      /* Decompose the neighborhood rank */

      int idim = 0;
      int ind  = rank;
      while (ind > 0)
      {
        if (ind % 2 == 1) indg[idim] += 1;
        ind /= 2;
        idim++;
      }

      /* Check that the neighboring point lies within the grid */

      bool incorrect = false;
      for (idim = 0; idim < ndim && !incorrect; idim++)
      {
        if (indg[idim] >= dbin->getNX(idim))
        {
          if (_flagExpand)
            indg[idim]--;
          else
            incorrect = true;
        }
      }

      /* Process the new neighboring point */

      if (incorrect)
      {
        nbgh.clear();
        weights.clear();
        break;
      }

      /* Check the value */

      int iech_neigh   = dbin->indiceToRank(indg);
      double val_neigh = dbin->getZVariable( iech_neigh, 0);
      if (FFFF(val_neigh))
      {
        nbgh.clear();
        weights.clear();
        break;
      }

      /* Calculate the distance from neighborhood to target */

      dbin->indicesToCoordinateInPlace(indg, coor, percent);
      double dist = ut_distance(ndim, cooref.data(), coor.data());
      if (dist < dmin)
      {
        nbgh.clear();
        weights.clear();
        nbgh.push_back(iech_neigh);
        weights.push_back(1.);
        break;
      }
      double wgt = 1. / pow(dist, _exponent);
      nbgh.push_back(iech_neigh);
      weights.push_back(wgt);
    }

    // Save the results
    _saveResults(dbin, dbout, nbgh, iech, weights);
  }
}

double CalcSimpleInterpolation::_estimCalc(const Db *dbin,
                                           const VectorInt &nbgh,
                                           const VectorDouble &weights)
{

  double result = 0.;
  for (int i = 0, n = (int)nbgh.size(); i < n; i++)
  {
    int iech     = nbgh[i];
    double value = dbin->getZVariable(iech, 0);
    if (FFFF(value)) return TEST;
    result += value * weights[i];

    if (OptDbg::query(EDbg::RESULTS))
      message("Data Value = %f - Weight = %lf\n", value, weights[i]);
  }
  if (OptDbg::query(EDbg::RESULTS)) message("Estimate = %f\n", result);
  return result;
}

/**
 * Calculate the (non-optimal) variance of estimation error
 * @param dbin Db containing the data
 * @param dbout Db containing the target
 * @param nbgh List of ranks within 'dbin'
 * @param iechout Rank of the target
 * @param weights Vector of weights (same dimension as 'nbgh')
 * @return
 */
double CalcSimpleInterpolation::_stdevCalc(Db* dbin,
                                           Db* dbout,
                                           const VectorInt& nbgh,
                                           int iechout,
                                           const VectorDouble& weights) const
{
  int ndim = dbin->getNDim();
  VectorDouble coor(ndim);
  SpacePoint pout;

  dbout->getCoordinatesPerSampleInPlace(iechout, coor);
  pout.setCoords(coor);

  // Point Covariance at target
  double c00 = getModel()->eval(pout, pout);

  // Vector of Covariances between Data and Target
  VectorDouble M0x = getModel()->evalPointToDb(pout, dbin, 0, 0, true, nbgh);
  double c0x       = VH::innerProduct(M0x, weights);

  // Covariance between Data and Data
  MatrixRectangular Mxx =
    getModel()->evalCovMat(dbin, dbin, 0, 0, nbgh, nbgh);
  double cxx = Mxx.quadraticMatrix(weights, weights);

  double result = sqrt(c00 - 2. * c0x + cxx);
  if (OptDbg::query(EDbg::RESULTS)) message("St Dev = %f\n", result);

  return result;
}

void CalcSimpleInterpolation::_saveResults(Db* dbin,
                                           Db* dbout,
                                           const VectorInt& nbgh,
                                           int iech,
                                           VectorDouble& weights) const
{
  double result = TEST;
  double stdev = TEST;
  if (nbgh.size() > 0)
  {
    VH::normalize(weights, 1);
    if (_flagEst) result = _estimCalc(dbin, nbgh, weights);
    if (_flagStd) stdev = _stdevCalc(dbin, dbout, nbgh, iech, weights);
  }

  // Assign the result
  if (_flagEst) dbout->setArray(iech, _iattEst, result);
  if (_flagStd) dbout->setArray(iech, _iattStd, stdev);
}

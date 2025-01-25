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
#include "Model/ModelOptimSillsVMap.hpp"

#include "Model/Model.hpp"
#include "Variogram/Vario.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"
#include "Model/Constraints.hpp"
#include "Db/DbGrid.hpp"

#include "geoslib_old_f.h"

#define IJDIR(ijvar, ipadir) ((ijvar) * _npadir + (ipadir))
#define WT(ijvar, ipadir)       wt[IJDIR(ijvar, ipadir)]
#define GG(ijvar, ipadir)       gg[IJDIR(ijvar, ipadir)]
#define _WT(ijvar, ipadir)       _wt[IJDIR(ijvar, ipadir)]
#define _GG(ijvar, ipadir)       _gg[IJDIR(ijvar, ipadir)]
#define _WT2(ijvar, ipadir) _wt2[IJDIR(ijvar, ipadir)]
#define _GG2(ijvar, ipadir)      _gg2[IJDIR(ijvar, ipadir)]
#define TAB(ijvar, ipadir)      tabin[IJDIR(ijvar, ipadir)]
#define DD(idim, ijvar, ipadir) _dd[idim][IJDIR(ijvar, ipadir)]

#define CORRECT(idir, k)                                                       \
  (!isZero(vario->getHhByIndex(idir, k)) &&                                    \
   !FFFF(vario->getHhByIndex(idir, k)) &&                                      \
   !isZero(vario->getSwByIndex(idir, k)) &&                                    \
   !FFFF(vario->getSwByIndex(idir, k)) && !FFFF(vario->getGgByIndex(idir, k)))
#define INCORRECT(idir, k)                                                     \
  (isZero(vario->getHhByIndex(idir, k)) ||                                     \
   FFFF(vario->getHhByIndex(idir, k)) ||                                       \
   isZero(vario->getSwByIndex(idir, k)) ||                                     \
   FFFF(vario->getSwByIndex(idir, k)) || FFFF(vario->getGgByIndex(idir, k)))

ModelOptimSillsVMap::ModelOptimSillsVMap(Model* model,
                                         Constraints* constraints,
                                         const Option_AutoFit& mauto,
                                         const Option_VarioFit& optvar)
  : AModelOptimSills(model, constraints, mauto, optvar)
{
}

ModelOptimSillsVMap::ModelOptimSillsVMap(const ModelOptimSillsVMap& m)
  : AModelOptimSills(m)
{
}

ModelOptimSillsVMap& ModelOptimSillsVMap::operator=(const ModelOptimSillsVMap& m)
{
  if (this != &m)
  {
    AModelOptimSills::operator=(m);
  }
  return (*this);
}

ModelOptimSillsVMap::~ModelOptimSillsVMap()
{
}

int ModelOptimSillsVMap::loadEnvironment(const DbGrid* dbmap, bool verbose)
{
  _dbmap = dbmap;
  _modelPart._verbose = verbose;

  // Get internal dimension
  if (_getDimensions()) return 1;

  // Allocate arrays
  _allocateInternalArrays(false);

  // Initialize Model-free quantities
  _computeVMap();

  // Initialize the array of sills
  _resetSill(_ncova, _sill);

  return 0;
}

/****************************************************************************/
/*!
 **  General Routine for fitting a model using an experimental variogram
 **
 ** \return  Error return code
 **
 ** \param[in]  dbmap       Experimental Variogram Map
 ** \param[in]  verbose     Verbose flag
 **
 *****************************************************************************/
int ModelOptimSillsVMap::fit(const DbGrid* dbmap, bool verbose)
{
  // Define the environment
  if (loadEnvironment(dbmap, verbose)) return 1;

  // Initialize Model-dependent quantities
  updateFromModel();

  // Perform the sill fitting
  return fitPerform();
}

 /****************************************************************************/
/*!
 **  Fill the array of pointers on the experimental conditions
 **
 *****************************************************************************/
void ModelOptimSillsVMap::_computeVMap()
{
  const DbGrid* dbmap = _dbmap;
  int nvs2            = _nvar * (_nvar + 1) / 2;
  dbmap->rankToIndice(_nech / 2, _indg1);

  /* Load the Experimental conditions structure */

  int ipadir = 0;
  for (int iech = 0; iech < _nech; iech++)
  {
    dbmap->rankToIndice(iech, _indg2);
    double dist = distance_intra(dbmap, _nech / 2, iech, NULL);
    double wgt  = (dist > 0) ? 1. / dist : 0.;

    /* Check samples containing only undefined values */

    int ntest = 0;
    for (int ijvar = 0; ijvar < nvs2; ijvar++)
      if (!FFFF(dbmap->getZVariable(iech, ijvar))) ntest++;
    if (ntest <= 0) continue;

    for (int ijvar = 0; ijvar < nvs2; ijvar++)
    {
      _WT(ijvar, ipadir) = 0.;
      _GG(ijvar, ipadir) = 0.;

      double value = dbmap->getZVariable(iech, ijvar);
      if (FFFF(value)) continue;

      _WT(ijvar, ipadir) = wgt;
      _GG(ijvar, ipadir) = value;
    }
    ipadir++;
  }
}

/****************************************************************************/
/*!
 **  Fill the array of pointers on the moded
 **
 *****************************************************************************/
void ModelOptimSillsVMap::updateFromModel()
{
  Model* model = _modelPart._model;
  VectorDouble d0(_ndim);
  MatrixSquareGeneral tab(_nvar);
  _dbmap->rankToIndice(_nech / 2, _indg1);
  CovCalcMode mode(ECalcMember::LHS);
  mode.setAsVario(true);
  mode.setUnitary(true);

  /* Loop on the basic structures */

  for (int icov = 0; icov < _ncova; icov++)
  {
    mode.setActiveCovListFromOne(icov);

    /* Loop on the experiments */

    for (int ipadir = 0; ipadir < _npadir; ipadir++)
    {
      _dbmap->rankToIndice(ipadir, _indg2);
      for (int idim = 0; idim < _ndim; idim++)
        d0[idim] = (_indg2[idim] - _indg1[idim]) * _dbmap->getDX(idim);
      model->evaluateMatInPlace(nullptr, d0, tab, true, 1., &mode);

      /* Loop on the variables */

      int ijvar = 0;
      for (int ivar = 0; ivar < _nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
          _ge[icov].setValue(ijvar, ipadir, tab.getValue(ivar, jvar));
    }
  }
}

int ModelOptimSillsVMap::_getDimensions()
{
  int nbexp  = 0;
  int npadir = 0;
  int nvs2   = _nvar * (_nvar + 1) / 2;
  _nech      = _dbmap->getNSample();
  _nvar      = _dbmap->getLocNumber(ELoc::Z);
  _ndim      = _dbmap->getLocNumber(ELoc::X);

  /* Calculate the total number of lags */

  for (int iech = 0; iech < _nech; iech++)
  {
    int ndef = 0;
    for (int ijvar = 0; ijvar < nvs2; ijvar++)
      if (!FFFF(_dbmap->getZVariable(iech, ijvar))) ndef++;
    nbexp += ndef;
    if (ndef > 0) npadir++;
  }

  /* Setting the return arguments */

  _nbexp  = nbexp;
  _npadir = npadir;

  if (nbexp <= 0)
  {
    messerr("No active experimental variogram map samples");
    return (1);
  }
  return (0);
}

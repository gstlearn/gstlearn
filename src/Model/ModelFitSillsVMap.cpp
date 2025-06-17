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
#include "Model/ModelFitSillsVMap.hpp"

#include "Model/Model.hpp"
#include "Model/Constraints.hpp"
#include "Db/DbGrid.hpp"

#include "geoslib_old_f.h"

#define IJDIR(ijvar, ipadir) ((ijvar) * _npadir + (ipadir))
#define WT(ijvar, ipadir)       wt[IJDIR(ijvar, ipadir)]
#define GG(ijvar, ipadir)       gg[IJDIR(ijvar, ipadir)]
#define _WT(ijvar, ipadir)     _wt[IJDIR(ijvar, ipadir)]
#define _GG(ijvar, ipadir)     _gg[IJDIR(ijvar, ipadir)]

ModelFitSillsVMap::ModelFitSillsVMap(const DbGrid* dbmap,
                                     ModelCovList* model,
                                     Constraints* constraints,
                                     const ModelOptimParam& mop)
  : AModelFitSills(model, constraints, mop)
  , _dbmap(dbmap)
{
  (void) _prepare();
}

ModelFitSillsVMap::ModelFitSillsVMap(const ModelFitSillsVMap& m)
  : AModelFitSills(m)
  , _dbmap(m._dbmap)
{
  (void)_prepare();
}

ModelFitSillsVMap& ModelFitSillsVMap::operator=(const ModelFitSillsVMap& m)
{
  if (this != &m)
  {
    AModelFitSills::operator=(m);
    _dbmap = m._dbmap;
    (void)_prepare();
  }
  return (*this);
}

ModelFitSillsVMap::~ModelFitSillsVMap()
{
}

ModelFitSillsVMap* ModelFitSillsVMap::createForOptim(const DbGrid* dbmap,
                                                     ModelGeneric* model,
                                                     Constraints* constraints,
                                                     const ModelOptimParam& mop)
{
  ModelCovList* modelLocal = dynamic_cast<ModelCovList*>(model);
  if (modelLocal == nullptr)
  {
    messerr("The argument 'model' should be a 'ModelCovList'");
    return nullptr;
  }
  ModelFitSillsVMap* optim = new ModelFitSillsVMap(dbmap, modelLocal, constraints, mop);

  return optim;
}

int ModelFitSillsVMap::_prepare()
{
  // Get internal dimension
  if (_getDimensions()) return 1;

  // Allocate arrays
  _allocateInternalArrays(false);

  // Initialize Model-free quantities
  _computeVMap();

  // Initialize the array of sills
  _resetInitialSill(_sill);

  _calcmode = CovCalcMode(ECalcMember::RHS);
  _calcmode.setAsVario(true);
  _calcmode.setUnitary(true);

  return 0;
}

/****************************************************************************/
/*!
 **  General Routine for fitting a model using an experimental variogram
 **
 ** \return  Error return code
 **
 *****************************************************************************/
int ModelFitSillsVMap::fitSills()
{

  // Initialize Model-dependent quantities
  _updateFromModel();

  // In this iterative manner of Fitting Sills, the verbose flag is switched OFF
  // in order to avoid intermediate printouts
  setVerbose(false);
  int status =  _fitSills();

  return status;
}

/****************************************************************************/
/*!
 **  Fill the array of pointers on the experimental conditions
 **
 *****************************************************************************/
void ModelFitSillsVMap::_computeVMap()
{
  const DbGrid* dbmap = _dbmap;
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
    for (int ijvar = 0; ijvar < _nvs2; ijvar++)
      if (!FFFF(dbmap->getZVariable(iech, ijvar))) ntest++;
    if (ntest <= 0) continue;

    for (int ijvar = 0; ijvar < _nvs2; ijvar++)
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
void ModelFitSillsVMap::_updateFromModel()
{
  VectorDouble d0(_ndim);
  MatrixSquare tab(_nvar);
  _dbmap->rankToIndice(_nech / 2, _indg1);

  /* Loop on the basic structures */

  const CovList* cova = _model->getCovList();
  for (int icov = 0; icov < _ncova; icov++)
  {
    cova->setActiveCovListFromOne(icov);

    /* Loop on the experiments */

    for (int ipadir = 0; ipadir < _npadir; ipadir++)
    {
      _dbmap->rankToIndice(ipadir, _indg2);
      for (int idim = 0; idim < _ndim; idim++)
        d0[idim] = (_indg2[idim] - _indg1[idim]) * _dbmap->getDX(idim);
      _model->evaluateMatInPlace(nullptr, d0, tab, true, 1., &_calcmode);

      /* Loop on the variables */

      int ijvar = 0;
      for (int ivar = 0; ivar < _nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
          _ge[icov].setValue(ijvar, ipadir, tab.getValue(ivar, jvar));
    }
  }
}

int ModelFitSillsVMap::_getDimensions()
{
  int nbexp  = 0;
  int npadir = 0;
  _nech      = _dbmap->getNSample();
  _nvar      = _dbmap->getNLoc(ELoc::Z);
  _ndim      = _dbmap->getNLoc(ELoc::X);
  _nvs2      = _nvar * (_nvar + 1) / 2;
  _indg1.resize(_ndim);
  _indg2.resize(_ndim);

  /* Calculate the total number of lags */

  for (int iech = 0; iech < _nech; iech++)
  {
    int ndef = 0;
    for (int ijvar = 0; ijvar < _nvs2; ijvar++)
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

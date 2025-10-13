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
#include "Model/ModelFitSillsVario.hpp"

#include "Covariances/CovBase.hpp"
#include "Model/Constraints.hpp"
#include "Model/Model.hpp"
#include "Model/ModelOptimVario.hpp"
#include "Variogram/Vario.hpp"

#define IJDIR(ijvar, ipadir)    ((ijvar) * _npadir + (ipadir))
#define _WT(ijvar, ipadir)      _wt[IJDIR(ijvar, ipadir)]
#define _GG(ijvar, ipadir)      _gg[IJDIR(ijvar, ipadir)]
#define _WT2(ijvar, ipadir)     _wt2[IJDIR(ijvar, ipadir)]
#define _GG2(ijvar, ipadir)     _gg2[IJDIR(ijvar, ipadir)]
#define TAB(ijvar, ipadir)      tabin[IJDIR(ijvar, ipadir)]
#define DD(idim, ijvar, ipadir) _dd[idim][IJDIR(ijvar, ipadir)]

namespace gstlrn
{
ModelFitSillsVario::ModelFitSillsVario(const Vario* vario,
                                       ModelCovList* model,
                                       const Constraints* constraints,
                                       const ModelOptimParam& mop)
  : AModelFitSills(model, constraints, mop)
  , _vario(vario)
{
  (void)_prepare();
}

ModelFitSillsVario::ModelFitSillsVario(const ModelFitSillsVario& m)
  : AModelFitSills(m)
  , _vario(m._vario)
{
  (void)_prepare();
}

ModelFitSillsVario& ModelFitSillsVario::operator=(const ModelFitSillsVario& m)
{
  if (this != &m)
  {
    AModelFitSills::operator=(m);
    _vario = m._vario;
    (void)_prepare();
  }
  return (*this);
}

ModelFitSillsVario::~ModelFitSillsVario()
{
}

ModelFitSillsVario* ModelFitSillsVario::createForOptim(const Vario* vario,
                                                       ModelGeneric* model,
                                                       const Constraints* constraints,
                                                       const ModelOptimParam& mop)
{
  ModelCovList* mcv = dynamic_cast<ModelCovList*>(model);
  if (mcv == nullptr)
  {
    messerr("The argument 'model' should be a 'ModelCovList'");
    return nullptr;
  }
  auto* optim = new ModelFitSillsVario(vario, mcv, constraints, mop);

  return optim;
}

Id ModelFitSillsVario::_prepare()
{
  // Get internal dimension
  if (_getDimensions()) return 1;

  // Allocate internal arrays
  _allocateInternalArrays(true);

  // Initialize Model-free quantities
  Id wmode = _mop.getWmode();
  _wt      = _vario->computeWeightsFromVario(wmode);
  _compressArray(_wt, _wtc);
  _computeGg();
  _compressArray(_gg, _ggc);

  // Initialize the array of sills
  _resetInitialSill(_sill);

  Id norder = 0;
  if (_vario->getCalcul() == ECalcVario::GENERAL1) norder = 1;
  if (_vario->getCalcul() == ECalcVario::GENERAL2) norder = 2;
  if (_vario->getCalcul() == ECalcVario::GENERAL3) norder = 3;
  _calcmode = CovCalcMode(ECalcMember::LHS);
  _calcmode.setAsVario(true);
  _calcmode.setUnitary(true);
  _calcmode.setOrderVario(norder);
  return 0;
}

/****************************************************************************/
/*!
 **  General Routine for fitting a model using an experimental variogram
 **
 ** \return  Error return code
 **
 *****************************************************************************/
Id ModelFitSillsVario::fitSillMatrices()
{
  // Initialize Model-dependent quantities
  _updateFromModel();

  Id status = _fitSillMatrices();

  return status;
}

/****************************************************************************/
/*!
 **  Calculate the main dimensions
 **
 *****************************************************************************/
Id ModelFitSillsVario::_getDimensions()
{
  _ndim  = static_cast<Id>(_model->getNDim());
  _nvar  = _model->getNVar();
  _ncova = _model->getNCov();
  _nvs2  = _nvar * (_nvar + 1) / 2;

  Id nbexp  = 0;
  Id npadir = 0;

  /* Calculate the total number of lags */

  for (Id idir = 0; idir < _vario->getNDir(); idir++)
  {
    npadir += _vario->getNLagTotal(idir);
    for (Id ilag = 0; ilag < _vario->getNLag(idir); ilag++)
      for (Id ivar = 0; ivar < _nvar; ivar++)
        for (Id jvar = 0; jvar <= ivar; jvar++)
        {
          Id i = _vario->getDirAddress(idir, ivar, jvar, ilag, false, 1);
          if (_vario->isLagCorrect(idir, i)) nbexp++;
        }
  }

  if (nbexp <= 0)
  {
    messerr("No active experimental variogram");
    return (1);
  }

  _nbexp  = nbexp;
  _npadir = npadir;
  return (0);
}

/****************************************************************************/
/*!
 **  Fill the array of pointers on the experimental conditions
 **
 *****************************************************************************/
void ModelFitSillsVario::_computeGg()
{
  Id ipadir = 0;
  for (Id idir = 0, ndir = _vario->getNDir(); idir < ndir; idir++)
  {
    for (Id ilag = 0, nlag = _vario->getNLag(idir); ilag < nlag; ilag++, ipadir++)
    {
      Id ijvar = 0;
      for (Id ivar = ijvar = 0; ivar < _nvar; ivar++)
        for (Id jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {

          // Calculate the variogram value
          double dist        = 0.;
          _GG(ijvar, ipadir) = TEST;
          if (_vario->getFlagAsym())
          {
            Id iad     = _vario->getDirAddress(idir, ivar, jvar, ilag, false, 1);
            Id jad     = _vario->getDirAddress(idir, ivar, jvar, ilag, false, -1);
            double c00 = _vario->getC00(idir, ivar, jvar);
            double n1  = _vario->getSwByIndex(idir, iad);
            double n2  = _vario->getSwByIndex(idir, jad);
            if (n1 + n2 > 0)
            {
              double g1 = _vario->getGgByIndex(idir, iad);
              double g2 = _vario->getGgByIndex(idir, jad);
              if (_vario->isLagCorrect(idir, iad) &&
                  _vario->isLagCorrect(idir, jad))
              {
                _GG(ijvar, ipadir) = c00 - (n1 * g1 + n2 * g2) / (n1 + n2);
                dist               = (ABS(_vario->getHhByIndex(idir, iad)) +
                        ABS(_vario->getHhByIndex(idir, jad))) /
                       2.;
              }
            }
          }
          else
          {
            Id iad = _vario->getDirAddress(idir, ivar, jvar, ilag, false, 1);
            if (_vario->isLagCorrect(idir, iad))
            {
              _GG(ijvar, ipadir) = _vario->getGgByIndex(idir, iad);
              dist               = ABS(_vario->getHhByIndex(idir, iad));
            }
          }

          // Store the distances
          Id i = _vario->getDirAddress(idir, ivar, jvar, ilag, false, 1);
          for (Id idim = 0; idim < _ndim; idim++)
          {
            if (!_vario->isLagCorrect(idir, i)) continue;
            DD(idim, ijvar, ipadir) = dist * _vario->getCodir(idir, idim);
          }
        }
    }
  }
}

/*****************************************************************************/
/*!
 **  Calculates the values of a generic covariance model corresponding
 **  to the lags of an experimental variogram
 **
 *****************************************************************************/
void ModelFitSillsVario::_updateFromModel()
{
  VectorDouble d1(_ndim);

  /* Loop on the basic structures */

  for (Id icov = 0; icov < _model->getNCov(); icov++)
  {
    const CovBase* cova = _model->getCovList()->getCov(icov);
    d1.fill(0.);

    /* Loop on the experiments */

    Id ipadir = 0;
    for (Id idir = 0, ndir = _vario->getNDir(); idir < ndir; idir++)
    {
      for (Id ilag = 0, nlag = _vario->getNLag(idir); ilag < nlag; ilag++, ipadir++)
      {
        Id ijvar = 0;
        for (Id ivar = 0; ivar < _nvar; ivar++)
          for (Id jvar = 0; jvar <= ivar; jvar++, ijvar++)
          {
            Id shift = ijvar * _vario->getNLagTotal(idir);
            if (!_ge.empty()) _ge[icov].setValue(ijvar, ipadir, 0.);

            double dist = 0.;
            if (_vario->getFlagAsym())
            {
              Id iad = shift + _vario->getNLag(idir) + ilag + 1;
              Id jad = shift + _vario->getNLag(idir) - ilag - 1;
              if (!_vario->isLagCorrect(idir, iad) ||
                  !_vario->isLagCorrect(idir, jad)) continue;
              dist = (ABS(_vario->getHhByIndex(idir, iad)) +
                      ABS(_vario->getHhByIndex(idir, jad))) /
                     2.;
            }
            else
            {
              Id iad = shift + ilag;
              if (!_vario->isLagCorrect(idir, iad)) continue;
              dist = ABS(_vario->getHhByIndex(idir, iad));
            }
            for (Id idim = 0; idim < _ndim; idim++)
              d1[idim] = dist * _vario->getCodir(idir, idim);

            if (!_ge.empty())
              _ge[icov].setValue(ijvar, ipadir, cova->evalIvarIpas(1., d1, ivar, jvar, &_calcmode));

            if (!_dd.empty())
              for (Id idim = 0; idim < _ndim; idim++)
                DD(idim, ijvar, ipadir) = d1[idim];
          }
      }
    }
  }
}

/****************************************************************************/
/*!
 **  Prepare the array for Goulard's algorithm
 **  in the case of Variogram calculation
 **
 *****************************************************************************/
void ModelFitSillsVario::_prepareGoulard()
{
  VectorDouble tab(_nvar * _nvar);
  VectorDouble d0(_ndim);
  CovCalcMode mode(ECalcMember::RHS);
  mode.setAsVario(true);
  mode.setUnitary(true);

  /* Loop on the basic structures */

  const CovList* cova = _model->getCovList();
  for (Id icov = 0, ncov = _ncova; icov < ncov; icov++)
  {
    cova->setActiveCovListFromOne(icov);

    /* Loop on the experiments */

    for (Id ipadir = 0; ipadir < _npadir; ipadir++)
    {
      Id ijvar = 0;
      for (Id ivar = 0; ivar < _nvar; ivar++)
        for (Id jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {
          Id flag_test = 0;
          for (Id idim = 0; idim < _ndim && flag_test == 0; idim++)
          {
            d0[idim] = DD(idim, ijvar, ipadir);
            if (FFFF(d0[idim])) flag_test = 1;
          }
          if (flag_test)
          {
            _ge[icov].setValue(ijvar, ipadir, TEST);
          }
          else
          {
            _ge[icov].setValue(ijvar, ipadir, _model->evalIvarIpas(1., d0, ivar, jvar, &mode));
          }
        }
    }
  }
  cova->setActiveCovListFromOne(-1);
}

/****************************************************************************/
/*!
 **  Compress the weights for the experimental variograms
 **
 ** \param[in]  tabin     Uncompressed array
 **
 ** \param[out] tabout    Compressed array
 **
 *****************************************************************************/
void ModelFitSillsVario::_compressArray(const VectorDouble& tabin,
                                        VectorDouble& tabout)
{
  Id ecr    = 0;
  Id ipadir = 0;
  for (Id idir = 0, ndir = _vario->getNDir(); idir < ndir; idir++)
    for (Id ilag = 0, nlag = _vario->getNLag(idir); ilag < nlag;
         ilag++, ipadir++)
    {
      Id ijvar = 0;
      for (Id ivar = 0; ivar < _nvar; ivar++)
        for (Id jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {
          double tabval = TAB(ijvar, ipadir);
          if (!FFFF(tabval)) tabout[ecr++] = tabval;
        }
    }
}
} // namespace gstlrn
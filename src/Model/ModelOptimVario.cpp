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
#include "Model/ModelOptimVario.hpp"

#include "Model/ModelOptimSillsVario.hpp"
#include "geoslib_define.h"

#include "Model/Model.hpp"
#include "Variogram/Vario.hpp"

#define IJDIR(ijvar, ipadir) ((ijvar)*npadir + (ipadir))
#define WT(ijvar, ipadir)     wt[IJDIR(ijvar, ipadir)]

ModelOptimVario::ModelOptimVario(Model* model,
                                 Constraints* constraints,
                                 const Option_AutoFit& mauto,
                                 const Option_VarioFit& optvar)
  : AModelOptim(model, constraints, mauto, optvar)
  , _varioPart()
  , _optGoulard(model)
{
}

ModelOptimVario::ModelOptimVario(const ModelOptimVario& m)
  : AModelOptim(m)
  , _varioPart()
  , _optGoulard(m._optGoulard)
{
   _copyVarioPart(m._varioPart);
}

ModelOptimVario& ModelOptimVario::operator=(const ModelOptimVario& m)
{
  if (this != &m)
  {
    AModelOptim::operator=(m);
    _optGoulard  = m._optGoulard;
    _copyVarioPart(m._varioPart);
  }
  return (*this);
}

ModelOptimVario::~ModelOptimVario()
{
}

void ModelOptimVario::_copyVarioPart(const Vario_Part& varioPart)
{
  _varioPart._vario = varioPart._vario;
  _varioPart._wmode = varioPart._wmode;
  _varioPart._lags  = varioPart._lags;
}

bool ModelOptimVario::_checkConsistency()
{
  const Model* model = _modelPart._model;
  const Vario* vario = _varioPart._vario;

  if (vario->getDimensionNumber() != model->getDimensionNumber())
  {
    messerr("'_vario'(%d) and '_model'(%d) should have same Space Dimension",
            vario->getDimensionNumber(), model->getDimensionNumber());
    return false;
  }
  if (vario->getNVar() != model->getNVar())
  {
    messerr("'_vario'(%d) and '_model'(%d) should have same number of Variables",
      vario->getNVar(), model->getNVar());
    return false;
  }
  return true;
}

int ModelOptimVario::loadEnvironment(Vario* vario,
                                     bool flagGoulard,
                                     int wmode,
                                     bool verbose)
{
  _modelPart._verbose = verbose;
  _varioPart._vario   = vario;
  _varioPart._wmode   = wmode;
  _modelPart._optvar.setFlagGoulardUsed(flagGoulard);

  // Constitute the experimental material (using '_vario')
  if (_buildExperimental()) return 1;

  // Constitute the list of parameters
  if (_buildModelParamList()) return 1;

  // Check consistency
  if (!_checkConsistency()) return 1;

  // Instantiate Goulard algorithm (optional)
  if (flagGoulard)
  {
    _optGoulard = ModelOptimSillsVario(_modelPart._model, _constraints, _mauto,
                                       _modelPart._optvar);
    _optGoulard.loadEnvironment(vario, wmode, verbose);
  }

  return 0;
}

int ModelOptimVario::fit(Vario* vario, bool flagGoulard, int wmode, bool verbose)
{
  // Load the Environment
  if (loadEnvironment(vario, flagGoulard, wmode, verbose)) return 1;

  // Perform the optimization
  AlgorithmVario algorithm {_modelPart, _varioPart, _optGoulard};
  _performOptimization(evalCost, &algorithm, vario->getMaximumDistance(),
                       vario->getVarMatrix());

  return 0;
}

int ModelOptimVario::_buildExperimental()
{
  if (_varioPart._vario == nullptr)
  {
    messerr("Argument 'vario' must be defined beforehand");
    return 1;
  }
  const Vario* vario = _varioPart._vario;

  // Clean previous contents
  _varioPart._lags.clear();

  int nvar = vario->getNVar();
  int ndim = vario->getDimensionNumber();
  VectorDouble dd(ndim);

  for (int idir = 0, ndir = vario->getDirectionNumber(); idir < ndir; idir++)
  {
    for (int ipas = 0, npas = vario->getLagNumber(idir); ipas < npas; ipas++)
    {
      int ijvar = 0;
      for (int ivar = ijvar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {

          /* Calculate the variogram value */

          double dist = 0.;
          double gg   = TEST;
          if (vario->getFlagAsym())
          {
            int iad = vario->getDirAddress(idir, ivar, jvar, ipas, false, 1);
            int jad = vario->getDirAddress(idir, ivar, jvar, ipas, false, -1);
            double c00 = vario->getC00(idir, ivar, jvar);
            double n1  = vario->getSwByIndex(idir, iad);
            double n2  = vario->getSwByIndex(idir, jad);
            if (n1 + n2 > 0)
            {
              double g1 = vario->getGgByIndex(idir, iad);
              double g2 = vario->getGgByIndex(idir, jad);
              if (vario->isLagCorrect(idir, iad) && vario->isLagCorrect(idir, jad))
              {
                gg   = c00 - (n1 * g1 + n2 * g2) / (n1 + n2);
                dist = (ABS(vario->getHhByIndex(idir, iad)) +
                        ABS(vario->getHhByIndex(idir, jad))) / 2.;
              }
            }
          }
          else
          {
            int iad = vario->getDirAddress(idir, ivar, jvar, ipas, false, 1);
            if (vario->isLagCorrect(idir, iad))
            {
              gg   = vario->getGgByIndex(idir, iad);
              dist = ABS(vario->getHhByIndex(idir, iad));
            }
          }

          /* Define the item of the StrExp array (if defined) */

          if (FFFF(gg)) continue;
          OneLag onelag = _createOneLag(ndim, idir, ivar, jvar, gg, dist);
          _varioPart._lags.push_back(onelag);
        }
    }
  }

  // Update the weight
  VectorDouble wt = _varioPart._vario->computeWeightsFromVario(_varioPart._wmode);
  int npadir      = _varioPart._vario->getTotalLagsPerDirection();
  int ecr         = 0;
  int ipadir      = 0;

  for (int idir = 0, ndir = vario->getDirectionNumber(); idir < ndir; idir++)
  {
    for (int ipas = 0, npas = vario->getLagNumber(idir); ipas < npas; ipas++, ipadir++)
    {
      int ijvar = 0;
      for (int ivar = ijvar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
          _varioPart._lags[ecr]._weight = WT(ijvar, ipadir);
    }
  }

  return 0;
}

ModelOptimVario::OneLag ModelOptimVario::_createOneLag(int ndim,
                                                       int idir,
                                                       int ivar,
                                                       int jvar,
                                                       double gg,
                                                       double dist) const
{
  OneLag onelag;
  onelag._ivar   = ivar;
  onelag._jvar   = jvar;
  onelag._gg     = gg;
  onelag._weight = 1.;
  VectorDouble dd(ndim);
  for (int idim = 0; idim < ndim; idim++)
    dd[idim] = dist * _varioPart._vario->getCodir(idir, idim);
  onelag._P.setCoords(dd);
  return onelag;
}

double ModelOptimVario::evalCost(unsigned int nparams,
                                 const double* current,
                                 double* /*grad*/,
                                 void* my_func_data)
{
  DECLARE_UNUSED(nparams);
  AlgorithmVario* algorithm = static_cast<AlgorithmVario*>(my_func_data);
  if (algorithm == nullptr) return TEST;
  Model_Part& modelPart            = algorithm->_modelPart;
  Vario_Part& varioPart            = algorithm->_varioPart;
  ModelOptimSillsVario& optGoulard = algorithm->_goulardPart;

  // Update the Model
  _patchModel(modelPart, current);

  // Perform sill fitting using Goulard (optional)
  if (modelPart._optvar.getFlagGoulardUsed())
  {
    optGoulard.updateFromModel();
    optGoulard.fitPerform();
  }
  
  // Evaluate the Cost function
  int nlags = (int) varioPart._lags.size();
  double total = 0.;
  SpacePoint origin;
  for (int ipas = 0; ipas < nlags; ipas++)
  {
    const OneLag& lag = varioPart._lags[ipas];
    double vexp        = lag._gg;
    double vtheo = modelPart._model->eval(origin, lag._P, lag._ivar, lag._jvar, &modelPart._calcmode);
    double delta = vexp - vtheo;
    total += lag._weight * delta * delta;
  }
  _printResult("Cost Function (Variogram Fit)", modelPart, total);
  
  return total;
}

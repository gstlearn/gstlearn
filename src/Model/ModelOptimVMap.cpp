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
#include "Model/ModelOptimVMap.hpp"

#include "geoslib_define.h"

#include "Model/Model.hpp"
#include "Variogram/VMap.hpp"

#include <nlopt.h>

#define IJDIR(ijvar, ipadir) ((ijvar)*npadir + (ipadir))
#define WT(ijvar, ipadir)    _wt[IJDIR(ijvar, ipadir)]

typedef struct
{
  // Part of the structure dedicated to the Model
  ModelOptim::Model_Part& _modelPart;

  // Part relative to the Experimental variograms
  ModelOptimVmap::Vario_Part& _varioPart;

} AlgorithmVario;

ModelOptimVmap::ModelOptimVmap()
  : ModelOptim()
  , _varioPart()
  , _wt()
{
}

ModelOptimVmap::ModelOptimVmap(const ModelOptimVmap& m)
  : ModelOptim(m)
  , _varioPart()
  , _wt(m._wt)
{
   _copyVarioPart(m._varioPart);
}

void ModelOptimVmap::_copyVarioPart(const Vario_Part& varioPart)
{
  _varioPart._vario = varioPart._vario;
  _varioPart._wmode = varioPart._wmode;
  _varioPart._lags  = varioPart._lags;
}

ModelOptimVmap& ModelOptimVmap::operator=(const ModelOptimVmap& m)
{
  if (this != &m)
  {
    ModelOptim::operator=(m);
    _copyVarioPart(m._varioPart);
    _wt = m._wt;
  }
  return (*this);
}

ModelOptimVmap::~ModelOptimVmap()
{
}

bool ModelOptimVmap::_checkConsistency()
{
  const Model* model = _modelPart._model;
  const Vario* vario = _varioPart._vario;

  if (vario->getDimensionNumber() != model->getDimensionNumber())
  {
    messerr("'_vario'(%d) and '_model'(%d) should have same Space Dimension",
            vario->getDimensionNumber(), model->getDimensionNumber());
    return false;
  }
  if (vario->getVariableNumber() != model->getVariableNumber())
  {
    messerr("'_vario'(%d) and '_model'(%d) should have same number of Variables",
      vario->getVariableNumber(), model->getVariableNumber());
    return false;
  }
  return true;
}

int ModelOptimVmap::fit(Vario* vario, Model* model, int wmode, bool verbose)
{
  _modelPart._model   = model;
  _modelPart._verbose = verbose;
  _varioPart._vario   = vario;
  _varioPart._wmode   = wmode;

  // Constitute the experimental material (using '_vario')
  if (_buildExperimental()) return 1;

  // Constitute the list of parameters
  if (_buildModelParamList()) return 1;

  // Check consistency
  if (! _checkConsistency()) return 1;

  // Define the optimization criterion
  int npar    = _getParamNumber();
  nlopt_opt opt = nlopt_create(NLOPT_LN_NELDERMEAD, npar);
  nlopt_set_lower_bounds(opt, _modelPart._tablow.data());
  nlopt_set_upper_bounds(opt, _modelPart._tabupp.data());
  nlopt_srand(12345);
  nlopt_set_ftol_rel(opt, 1e-6);

  // Update the initial optimization values (due to variogram)
  updateModelParamList(vario->getMaximumDistance(), vario->getVarMatrix());

  // Define the cost function
  AlgorithmVario algorithm {_modelPart, _varioPart};
  nlopt_set_min_objective(opt, evalCost, &algorithm);

  // Perform the optimization (store the minimized value in 'minf')
  double minf;
  if (_modelPart._verbose) mestitle(1, "Model Fitting from Variogram");
  nlopt_optimize(opt, _modelPart._tabval.data(), &minf);
  nlopt_destroy(opt);
  return 0;
}

int ModelOptimVmap::_buildExperimental()
{
  if (_varioPart._vario == nullptr)
  {
    messerr("Argument 'vario' must be defined beforehand");
    return 1;
  }
  const Vario* vario = _varioPart._vario;

  // Clean previous contents
  _varioPart._lags.clear();

  int nvar = vario->getVariableNumber();
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
            double c00 = _getC00(idir, ivar, jvar);
            double n1  = vario->getSwByIndex(idir, iad);
            double n2  = vario->getSwByIndex(idir, jad);
            if (n1 + n2 > 0)
            {
              double g1 = vario->getGgByIndex(idir, iad);
              double g2 = vario->getGgByIndex(idir, jad);
              if (_isLagCorrect(idir, iad) && _isLagCorrect(idir, jad))
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
            if (_isLagCorrect(idir, iad))
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
  _computeWt();
  int npadir      = _getTotalLagsPerDirection();
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

ModelOptimVmap::OneLag ModelOptimVmap::_createOneLag(int ndim,
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

double ModelOptimVmap::_getC00(int idir, int ivar, int jvar) const
{
  const Vario* vario = _varioPart._vario;
  int iad0   = vario->getDirAddress(idir, ivar, jvar, 0, false, 0);
  int iad    = iad0;
  double c00 = vario->getSwByIndex(idir, iad);
  if (!isZero(c00) || vario->getSwByIndex(idir, iad) > 0) return c00;

  for (int ipas = 0, npas = vario->getLagNumber(idir); ipas < npas; ipas++)
  {
    iad = vario->getDirAddress(idir, ivar, jvar, ipas, false, 1);
    if (!isZero(vario->getGgByIndex(idir, iad)))
      return vario->getGgByIndex(idir, iad);
    iad = vario->getDirAddress(idir, ivar, jvar, ipas, false, -1);
    if (!isZero(vario->getGgByIndex(idir, iad)))
      return vario->getGgByIndex(idir, iad);
  }
  iad = iad0;
  return (vario->getGgByIndex(idir, iad));
}

bool ModelOptimVmap::_isLagCorrect(int idir, int k) const
{
  const Vario* vario = _varioPart._vario;
  double hh = vario->getHhByIndex(idir, k);
  if (isZero(hh) || FFFF(hh)) return false;
  double sw = vario->getSwByIndex(idir, k);
  if (isZero(sw) || FFFF(sw)) return false;
  double gg = vario->getGgByIndex(idir, k);
  return !FFFF(gg);
}

void ModelOptimVmap::_computeWt()
{
  const Vario* vario = _varioPart._vario;
  int ndir           = vario->getDirectionNumber();
  int nvar           = vario->getVariableNumber();
  int nvs2           = nvar * (nvar + 1) / 2;
  int npadir         = _getTotalLagsPerDirection();
  VectorDouble count = _computeWeightPerDirection();
  _wt.resize(npadir * nvs2);
  _wt.fill(0.);

  /* Determine the count of significant directions */

  for (int idir = 0; idir < ndir; idir++)
  {
    count[idir] = 0.;
    int npas    = vario->getLagNumber(idir);
    for (int ipas = 0; ipas < npas; ipas++)
      for (int ijvar = 0; ijvar < nvs2; ijvar++)
      {
        int shift = ijvar * vario->getLagTotalNumber(idir);
        if (vario->getFlagAsym())
        {
          int iad   = shift + vario->getLagNumber(idir) + ipas + 1;
          int jad   = shift + vario->getLagNumber(idir) - ipas - 1;
          double n1 = vario->getSwByIndex(idir, iad);
          double n2 = vario->getSwByIndex(idir, jad);
          if (_isLagCorrect(idir, iad)) count[idir] += n1;
          if (_isLagCorrect(idir, jad)) count[idir] += n2;
        }
        else
        {
          int iad   = shift + ipas;
          double nn = vario->getSwByIndex(idir, iad);
          if (_isLagCorrect(idir, iad)) count[idir] += nn;
        }
      }
  }

  int ipadir = 0;
  switch (_varioPart._wmode)
  {
    case 1:
      ipadir = 0;
      for (int idir = 0; idir < ndir; idir++)
      {
        int npas = vario->getLagNumber(idir);
        for (int ipas = 0; ipas < npas; ipas++, ipadir++)
        {
          if (isZero(count[idir])) continue;
          for (int ijvar = 0; ijvar < nvs2; ijvar++)
          {
            int shift = ijvar * vario->getLagTotalNumber(idir);
            if (vario->getFlagAsym())
            {
              int iad = shift + vario->getLagNumber(idir) + ipas + 1;
              int jad = shift + vario->getLagNumber(idir) - ipas - 1;
              if (_isLagCorrect(idir, iad) && _isLagCorrect(idir, jad))
                WT(ijvar, ipadir) = count[idir];
            }
            else
            {
              int iad = shift + ipas;
              if (_isLagCorrect(idir, iad)) WT(ijvar, ipadir) = count[idir];
            }
          }
        }
      }
      break;

    case 2:
      ipadir = 0;
      for (int idir = 0; idir < ndir; idir++)
      {
        int npas = vario->getLagNumber(idir);
        for (int ipas = 0; ipas < npas; ipas++, ipadir++)
        {
          if (isZero(count[idir])) continue;
          for (int ijvar = 0; ijvar < nvs2; ijvar++)
          {
            int shift = ijvar * vario->getLagTotalNumber(idir);
            if (vario->getFlagAsym())
            {
              int iad = shift + vario->getLagNumber(idir) + ipas + 1;
              int jad = shift + vario->getLagNumber(idir) - ipas - 1;
              if (_isLagCorrect(idir, iad) || _isLagCorrect(idir, jad))
                continue;
              double n1 = vario->getSwByIndex(idir, iad);
              double n2 = vario->getSwByIndex(idir, jad);
              double d1 = ABS(vario->getHhByIndex(idir, iad));
              double d2 = ABS(vario->getHhByIndex(idir, jad));
              if (d1 > 0 && d2 > 0)
                WT(ijvar, ipadir) =
                  sqrt((n1 + n2) * (n1 + n2) / (n1 * d1 + n2 * d2) / 2.);
            }
            else
            {
              int iad = shift + ipas;
              if (_isLagCorrect(idir, iad)) continue;
              double nn = vario->getSwByIndex(idir, iad);
              double dd = ABS(vario->getHhByIndex(idir, iad));
              if (dd > 0) WT(ijvar, ipadir) = nn / dd;
            }
          }
        }
      }
      break;

    case 3:
      ipadir = 0;
      for (int idir = 0; idir < ndir; idir++)
      {
        int npas = vario->getLagNumber(idir);
        for (int ipas = 0; ipas < npas; ipas++, ipadir++)
        {
          if (isZero(count[idir])) continue;
          for (int ijvar = 0; ijvar < nvs2; ijvar++)
          {
            int shift = ijvar * vario->getLagTotalNumber(idir);
            if (vario->getFlagAsym())
            {
              int iad = shift + vario->getLagNumber(idir) + ipas + 1;
              int jad = shift + vario->getLagNumber(idir) - ipas - 1;
              if (_isLagCorrect(idir, iad) && _isLagCorrect(idir, jad))
                WT(ijvar, ipadir) = 1. / vario->getLagNumber(idir);
            }
            else
            {
              int iad = shift + ipas;
              if (_isLagCorrect(idir, iad))
                WT(ijvar, ipadir) = 1. / vario->getLagNumber(idir);
            }
          }
        }
      }
      break;

    default:
      ipadir = 0;
      for (int idir = 0; idir < ndir; idir++)
      {
        int npas = vario->getLagNumber(idir);
        for (int ipas = 0; ipas < npas; ipas++, ipadir++)
        {
          if (isZero(count[idir])) continue;
          for (int ijvar = 0; ijvar < nvs2; ijvar++)
          {
            int shift = ijvar * vario->getLagTotalNumber(idir);
            if (vario->getFlagAsym())
            {
              int iad = shift + vario->getLagNumber(idir) + ipas + 1;
              int jad = shift + vario->getLagNumber(idir) - ipas - 1;
              if (_isLagCorrect(idir, iad) && _isLagCorrect(idir, jad))
                WT(ijvar, ipadir) = 1.;
            }
            else
            {
              int iad = shift + ipas;
              if (_isLagCorrect(idir, iad)) WT(ijvar, ipadir) = 1.;
            }
          }
        }
      }
      break;
  }

  /* Scaling by direction and by variable */

  for (int ijvar = 0; ijvar < nvs2; ijvar++)
  {
    ipadir = 0;
    for (int idir = 0; idir < ndir; idir++)
    {
      double total = 0.;
      int npas     = vario->getLagNumber(idir);
      for (int ipas = 0; ipas < npas; ipas++, ipadir++)
      {
        if (isZero(count[idir])) continue;
        if (WT(ijvar, ipadir) > 0 && !FFFF(WT(ijvar, ipadir)))
          total += WT(ijvar, ipadir);
      }
      if (isZero(total)) continue;
      ipadir -= vario->getLagNumber(idir);
      for (int ipas = 0, npas = vario->getLagNumber(idir); ipas < npas;
           ipas++, ipadir++)
      {
        if (isZero(count[idir])) continue;
        if (WT(ijvar, ipadir) > 0 && !FFFF(WT(ijvar, ipadir)))
          WT(ijvar, ipadir) /= total;
      }
    }
  }

  /* Scaling by variable variances */

  int ijvar0 = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++, ijvar0++)
    {
      double ratio =
        (vario->getVar(ivar, jvar) > 0 && vario->getVar(jvar, ivar) > 0)
          ? sqrt(vario->getVar(ivar, jvar) * vario->getVar(jvar, ivar)) : 1.;
      ipadir = 0;
      for (int idir = 0; idir < ndir; idir++)
      {
        int npas = vario->getLagNumber(idir);
        for (int ipas = 0; ipas < npas; ipas++, ipadir++)
          if (!FFFF(WT(ijvar0, ipadir))) WT(ijvar0, ipadir) /= ratio;
      }
    }
}

VectorDouble ModelOptimVmap::_computeWeightPerDirection()
{
  const Vario* vario = _varioPart._vario;
  int ndir = vario->getDirectionNumber();
  int nvar = vario->getVariableNumber();
  int nvs2 = nvar * (nvar + 1) / 2;
  VectorDouble count(ndir);

  /* Determine the count of significant directions */

  for (int idir = 0; idir < ndir; idir++)
  {
    count[idir] = 0.;
    for (int ipas = 0, npas = vario->getLagNumber(idir); ipas < npas; ipas++)
      for (int ijvar = 0; ijvar < nvs2; ijvar++)
      {
        int shift = ijvar * vario->getLagTotalNumber(idir);
        if (vario->getFlagAsym())
        {
          int iad   = shift + vario->getLagNumber(idir) + ipas + 1;
          int jad   = shift + vario->getLagNumber(idir) - ipas - 1;
          double n1 = vario->getSwByIndex(idir, iad);
          double n2 = vario->getSwByIndex(idir, jad);
          if (_isLagCorrect(idir, iad)) count[idir] += n1;
          if (_isLagCorrect(idir, jad)) count[idir] += n2;
        }
        else
        {
          int iad   = shift + ipas;
          double nn = vario->getSwByIndex(idir, iad);
          if (_isLagCorrect(idir, iad)) count[idir] += nn;
        }
      }
  }
  return count;
}

int ModelOptimVmap::_getTotalLagsPerDirection() const
{
  const Vario* vario = _varioPart._vario;
  int npatot = 0;
  int ndir   = vario->getDirectionNumber();
  for (int idir = 0; idir < ndir; idir++)
    npatot += vario->getLagTotalNumber(idir);
  return npatot;
}

double ModelOptimVmap::evalCost(unsigned int nparams,
                                 const double* current,
                                 double* /*grad*/,
                                 void* my_func_data)
{
  DECLARE_UNUSED(nparams);
  AlgorithmVario* algorithm = static_cast<AlgorithmVario*>(my_func_data);
  if (algorithm == nullptr) return TEST;
  Model_Part& modelPart = algorithm->_modelPart;
  Vario_Part& varioPart = algorithm->_varioPart;
  
  // Update the Model
  _patchModel(modelPart, current);

  // Evaluate the Cost function
  int nlags = (int) varioPart._lags.size();
  double total = 0.;
  SpacePoint origin;
  CovCalcMode calcmod;
  calcmod.setAsVario(true);
  for (int ipas = 0; ipas < nlags; ipas++)
  {
    const OneLag& lag = varioPart._lags[ipas];
    double vexp        = lag._gg;
    double vtheo = modelPart._model->eval(origin, lag._P, lag._ivar, lag._jvar, &calcmod);
    double delta = vexp - vtheo;
    total += lag._weight * delta * delta;
  }
  _printResult("Cost Function (Variogram Fit)", modelPart, total);
  
  return total;
}

int vmap_auto_fit(const DbGrid* dbmap,
                  Model* model,
                  bool verbose,
                  const Option_AutoFit& mauto_arg,
                  const Constraints& cons_arg,
                  const Option_VarioFit& optvar_arg)
{
  int npar0, npar, flag_reduce;
  VectorDouble varchol, scale, param, lower, upper;

  // Copy of const referencse into local classes

  Option_AutoFit mauto    = mauto_arg;
  Constraints constraints = cons_arg;
  Option_VarioFit optvar  = optvar_arg;

  /* Initializations */

  int error      = 1;
  int nbexp      = 0;
  int status     = 0;
  int norder     = 0;
  int npadir     = 0;
  int ncova      = model->getCovaNumber();
  int nvar       = model->getVariableNumber();
  int ndim       = model->getDimensionNumber();
  double hmax    = 0.;
  double gmax    = 0.;
  StrMod* strmod = nullptr;
  VectorDouble angles(ndim, 0.);
  mauto.setVerbose(verbose);

  /* Preliminary checks */

  if (nvar != dbmap->getLocNumber(ELoc::Z))
  {
    messerr("Number of variables in Db (%d) must match the one in Model (%d)",
            model->getVariableNumber(), dbmap->getLocNumber(ELoc::Z));
    goto label_end;
  }
  if (constraints.isConstraintSillDefined())
  {
    if (!optvar.getFlagGoulardUsed())
    {
      messerr("When Constraints on the sum of Sills are defined");
      messerr("The Goulard option must be switched ON");
      goto label_end;
    }
    if (!FFFF(constraints.getConstantSillValue()))
      constraints.expandConstantSill(nvar);
  }
  if (st_get_vmap_dimension(dbmap, nvar, &npadir, &nbexp)) goto label_end;
  hmax = dbmap->getExtensionDiagonal();
  st_vmap_varchol_manage(dbmap, varchol);

  /* Scale the parameters in the Option_AutoFit structure */

  st_mauto_rescale(nvar, varchol, mauto);

  /* Free the keypair mechanism strings */

  st_keypair_sill(-1, model);
  st_keypair_results(-1, 0, 0, NULL, NULL);

  /* Create the experimental structures */

  hmax /= 2.;
  DBMAP = dbmap;
  INDG1.fill(0., ndim);
  INDG2.resize(ndim);

  /* Core allocation */

  if (st_manage_recint(optvar, 0, ndim, nvar, 0, ncova, npadir)) goto label_end;
  st_load_vmap(npadir, RECINT.gg, RECINT.wt);

  /* Generate the default values */

  npar0 =
    st_vmap_auto_count(dbmap, model, constraints, optvar, param, lower, upper);

  /* Create the Model structures */

  strmod = st_model_auto_strmod_alloc(model, NULL, npar0, norder, hmax, angles,
                                      optvar, &npar);
  if (strmod == nullptr) goto label_end;
  if (npar == 0)
  {
    messerr("The VMAP Automatic Fitting procedure");
    messerr("does not allow using Goulard algorithm only");
    goto label_end;
  }
  scale.resize(npar);

  /* Set the default values and bounds */

  st_model_auto_pardef(strmod, npar, hmax, varchol, angles, param, lower,
                       upper);
  st_model_auto_scldef(strmod, npar, hmax, varchol, scale);
  st_model_auto_constraints_apply(strmod, npar, constraints, param, lower,
                                  upper);

  /* Minimization algorithm */

  STRMOD            = strmod;
  MAUTO             = mauto;
  CONSTRAINTS       = constraints;
  ST_PREPAR_GOULARD = st_prepar_goulard_vmap;
  do
  {
    st_model_auto_strmod_print(1, strmod, constraints, mauto, param, lower,
                               upper, npar, nbexp);
    status =
      foxleg_f(nbexp, npar, 0, MatrixRectangular(), param, lower, upper, scale,
               mauto, 0, st_strmod_vmap_evaluate, RECINT.gg, RECINT.wt);
    if (status > 0) goto label_end;

    st_model_auto_strmod_print(0, strmod, constraints, mauto, param, lower,
                               upper, npar, nbexp);
    flag_reduce = st_model_auto_strmod_reduce(strmod, &npar, hmax, gmax, param,
                                              lower, upper, constraints, mauto);

    /* Reset the parameters to their default values */

    if (status < 0)
    {
      for (int i = 0; i < npar; i++) param[i] = lower[i] = upper[i] = TEST;
      st_model_auto_pardef(strmod, npar, hmax, varchol, angles, param, lower,
                           upper);
    }
    st_model_auto_scldef(strmod, npar, hmax, varchol, scale);
  } while (flag_reduce && npar > 0);

  /* Perform the last cosmetic updates */

  st_model_post_update(strmod, optvar);

  /* Set the returned error code */

  if (mauto.getVerbose() && status < 0)
    messerr("Convergence not reached after %d iterations (%d parameters)",
            mauto.getMaxiter(), npar);

  /* Store the sills in the keypair mechanism */

  st_keypair_sill(1, model);

  error = 0;

label_end:
  st_model_auto_strmod_free(strmod);
  return (error);
}

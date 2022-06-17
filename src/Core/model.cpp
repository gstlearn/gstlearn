/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_f.h"
#include "geoslib_enum.h"
#include "geoslib_old_f.h"
#include "Drifts/DriftFactory.hpp"
#include "Drifts/DriftList.hpp"
#include "Drifts/EDrift.hpp"
#include "Drifts/ADrift.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/File.hpp"
#include "Covariances/CovLMC.hpp"
#include "Covariances/CovLMGradient.hpp"
#include "Covariances/CovLMCTapering.hpp"
#include "Covariances/CovLMCConvolution.hpp"
#include "Covariances/CovLMCAnamorphosis.hpp"
#include "Covariances/ACovAnisoList.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovGradientNumerical.hpp"
#include "Model/CovInternal.hpp"
#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Anamorphosis/AnamDiscreteDD.hpp"
#include "Anamorphosis/AnamDiscreteIR.hpp"
#include "Variogram/Vario.hpp"
#include "Space/SpaceRN.hpp"
#include "Basic/Law.hpp"
#include "Basic/String.hpp"
#include "Db/Db.hpp"
#include "csparse_f.h"

#include <math.h>

/*! \cond */
#define AD(ivar,jvar)          (ivar) + nvar * (jvar)
#define COVTAB(ivar,jvar)       covtab[AD(ivar,jvar)]
#define C00TAB(ivar,jvar)       c00tab[AD(ivar,jvar)]
#define AIC(icov,ivar,jvar)     aic[(icov)*nvar*nvar + AD(ivar,jvar)]
#define VALPRO(ivar)            valpro[(ivar)]
#define VECPRO(ivar,jvar)       vecpro[AD(ivar,jvar)]
#define CC(ivar,jvar)           cc[AD(ivar,jvar)]
#define DISC1(i,idim)          (koption->disc1[(idim) * koption->ntot + (i)])
#define DISC2(i,idim)          (koption->disc2[(idim) * koption->ntot + (i)])
#define G(i,j)                 (G[(i) * nech + j])
#define Gmatrix(i,j)           (Gmatrix[(j) * nech + i])
/*! \endcond */

static CovInternal *COVINT = nullptr;
int NDIM_LOCAL = 0;
VectorDouble X1_LOCAL = VectorDouble();
VectorDouble X2_LOCAL = VectorDouble();

/*****************************************************************************/
/*!
 **  Update the Model in the case of Non-stationary parameters
 **  This requires the knowledge of the two end-points
 **
 ** \param[in]  covint       Internal structure for non-stationarity
 **                          or NULL (for stationary case)
 ** \param[in]  model        Model structure
 **
 *****************************************************************************/
void model_nostat_update(CovInternal *covint, Model *model)
{
  if (!model->isNoStat()) return;
  if (covint == NULL) return;
  COVINT = covint;

  const ANoStat *nostat = model->getNoStat();
  nostat->updateModel(model, covint->getIcas1(), covint->getIech1(),
                      covint->getIcas2(), covint->getIech2());
}

/****************************************************************************/
/*!
 **  Check the compatibility between the Model and the Db
 **
 ** \return  Error returned code
 **
 ** \param[in]  model Model structure
 ** \param[in]  db    Db structure
 **
 *****************************************************************************/
static int st_check_environ(const Model *model, const Db *db)
{
  if (model->getDimensionNumber() == db->getNDim()) return (0);
  messerr("Dimension of the Db (%d) does not match dimension of the Model (%d)",
          db->getNDim(), model->getDimensionNumber());
  return (1);
}

/****************************************************************************/
/*!
 **  Check if the model is defined
 **
 ** \return  Error returned code
 **
 ** \param[in]  model Model structure
 **
 *****************************************************************************/
static int st_check_model(const Model *model)
{
  if (model != nullptr) return (0);
  messerr("No Model is defined");
  return (1);
}

/****************************************************************************/
/*!
 **  Check if the variable rank is correct
 **
 ** \return  Error return code
 **
 ** \param[in]  nvar Total number of variables
 ** \param[in]  ivar Rank of the variable
 **
 *****************************************************************************/
static int st_check_variable(int nvar, int ivar)
{
  if (ivar >= 0 && ivar < nvar) return (0);
  messerr("Error for the variable rank (%d)", ivar);
  messerr("It should lie within [0,%d[", nvar);
  return (1);
}

/****************************************************************************/
/*!
 **  Initializes the covtab array
 ***
 ** \param[in]  flag_init 1 If covtab() must be initialized; 0 otherwise
 ** \param[in]  model       Model structure
 ** \param[in]  covtab      Array to be initialized
 **
 *****************************************************************************/
void model_covtab_init(int flag_init, Model *model, double *covtab)
{
  int nvar = model->getVariableNumber();
  if (flag_init)
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
        COVTAB(ivar,jvar)= 0.;
}

/****************************************************************************/
/*!
 **  Scale the array COVTAB
 **
 ** \param[in]  nvar      Number of variables
 ** \param[in]  norme     Number of values used for scaling
 ** \param[in,out] covtab Input/output array covtab
 **
 *****************************************************************************/
static void st_covtab_rescale(int nvar, double norme, double *covtab)
{
  int ivar, jvar;

  if (norme != 0.)
    for (ivar = 0; ivar < nvar; ivar++)
      for (jvar = 0; jvar < nvar; jvar++)
        COVTAB(ivar,jvar) /= norme;
  return;
}

/*****************************************************************************/
/*!
 **  Evaluates a basic covariance structure for a vector of increments
 **  This basic structure is normalized
 **
 ** \return  Value of the basic structure
 **
 ** \param[in]  model       Structure containing the model
 ** \param[in]  icov        Rank of the Basic structure
 ** \param[in]  member      Member of the Kriging System (ECalcMember)
 ** \param[in]  d1          vector of increment (or NULL)
 **
 *****************************************************************************/
double model_calcul_basic(Model *model,
                          int icov,
                          const ECalcMember &member,
                          const VectorDouble &d1)
{
  // const CovAniso* cova = model->getCova(icov);
  // TODO: Why having to use ACov rather than CovAniso?
  const ACov *cova = dynamic_cast<const ACov*>(model->getCova(icov));

  if (member != ECalcMember::LHS && model->isCovaFiltered(icov))
    return (0.);
  else
    return cova->evalIvarIpas(0, 0, 1., d1);
}

/*****************************************************************************/
/*!
 **  Returns the covariances for an increment
 **  This is the generic internal function
 **  It can be called for stationary or non-stationary case
 **
 ** \param[in]  covint       Internal structure for non-stationarityAddress for the next term after the drift
 **                          or NULL (for stationary case)
 ** \param[in]  model        Model structure
 ** \param[in]  mode         CovCalcMode structure
 ** \param[in]  flag_init    Initialize the array beforehand
 **
 ** \param[out] d1          Working array (dimension = ndim) or NULL
 ** \param[out] covtab      output covariance (dimension = nvar * nvar)
 **
 *****************************************************************************/
void model_calcul_cov(CovInternal *covint,
                      Model *model,
                      const CovCalcMode &mode,
                      int flag_init,
                      double /*weight*/,
                      VectorDouble d1,
                      double *covtab)
{
  // Load the non-stationary parameters if needed

  if (model->isNoStat()) model_nostat_update(covint, model);

  // Evaluate the Model

  MatrixSquareGeneral mat = model->evalNvarIpas(1., d1, VectorDouble(), mode);

  int nvar = model->getVariableNumber();
  if (mat.getNTotal() != nvar * nvar)
    my_throw("Error in loading Covariance calculation into COVTAB");
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
    {
      double value = mat.getValue(ivar, jvar);
      if (flag_init)
        COVTAB(ivar,jvar)= value;
      else
        COVTAB(ivar,jvar) += value;
      }
  return;
}

/*****************************************************************************/
/*!
 **  Returns the covariances for an increment
 **
 ** \param[in]  model       Structure containing the model
 ** \param[in]  mode        CovCalcMode structure
 ** \param[in]  ivar        Rank of the first variable
 ** \param[in]  jvar        Rank of the second variable
 ** \param[in]  d1          vector of increment (or NULL)
 **
 *****************************************************************************/
double model_calcul_cov_ij(Model *model,
                           const CovCalcMode &mode,
                           int ivar,
                           int jvar,
                           const VectorDouble &d1)
{

  /* Modify the member in case of properties */

  if (model->getCovMode() != EModelProperty::NONE)
    my_throw("Model transformation is not programmed for this entry");

  /* Call the generic model calculation module */

  // TODO Correct this which has something to do with pure virtual eval although implemented in Acov
  // compared to eval0 which is not implemented with such arguments.
  double value = model->evalIvarIpas(ivar, jvar, 1., d1,VectorDouble(), mode);

  return value;
}

/****************************************************************************/
/*!
 **  Identify the coordinates of the two end-points
 **  (used by the external covariance function)
 **
 ** \return A (protected) pointer on the Covariance Internal class
 **
 *****************************************************************************/
const CovInternal* get_external_covariance()
{
  return COVINT;
}

/****************************************************************************/
/*!
 **  Returns the drift vector for a point
 **
 ** \param[in]  model  Model structure
 ** \param[in]  member Member of the Kriging System (ECalcMember)
 ** \param[in]  db     Db structure
 ** \param[in]  iech   Rank of the sample
 **
 ** \param[out] drftab  Array of drift values
 **
 *****************************************************************************/
void model_calcul_drift(Model *model,
                        const ECalcMember &member,
                        const Db *db,
                        int iech,
                        double *drftab)
{
  VectorDouble drft = model->evalDriftVec(db, iech, member);
  for (int il = 0; il < (int) drft.size(); il++)
    drftab[il] = drft[il];
  return;
}

/****************************************************************************/
/*!
 **  Deallocate the Model structure
 **
 ** \return  Pointer to the freed structure
 **
 ** \param[in]  model Model to be freed
 **
 *****************************************************************************/
Model* model_free(Model *model)

{
  /* Initializations */

  if (model == nullptr) return (model);
  delete model;
  return (NULL);
}

/****************************************************************************/
/*!
 **  Check if the non-stationary Model has a given non-stationary parameter
 **
 ** \return  1 if the given non-stationary parameter is defined; 0 otherwise
 **
 ** \param[in]   model    Model structure
 ** \param[in]   type0    Requested type (EConsElem)
 **
 *****************************************************************************/
int is_model_nostat_param(Model *model, const EConsElem &type0)
{
  if (!model->isNoStat()) return 1;
  const ANoStat *nostat = model->getNoStat();
  if (nostat->isDefinedByType(-1, type0)) return 1;

  return (0);
}

/****************************************************************************/
/*!
 **  Create a generic Model with nugget effect only
 **
 ** \return  Pointer to the Model structure newly allocated
 **
 ** \param[in]  ndim Space dimension
 ** \param[in]  nvar Number of variables
 **
 *****************************************************************************/
Model* model_default(int ndim, int nvar)
{
  /* Create the empty Model */

  CovContext ctxt = CovContext(nvar, ndim, 1.);
  Model* model = new Model(ctxt);

  /* Add the nugget effect variogram model */

  CovLMC covs(ctxt.getSpace());
  CovAniso cov(ECov::NUGGET, ctxt);
  covs.addCov(&cov);
  model->setCovList(&covs);
  return model;
}

/****************************************************************************/
/*!
 **  Add a basic covariance
 **
 ** \return  Error return code
 **
 ** \param[in]  model           Pointer to the Model structure
 ** \param[in]  type            Type of the basic structure
 ** \param[in]  flag_aniso      1 if the basic structure is anisotropic
 ** \param[in]  flag_rotation   1 if the basic structure is rotated
 **                             (only when anisotropic)
 ** \param[in]  range           Isotropic range of the basic structure
 ** \param[in]  param           Auxiliary parameter
 ** \param[in]  aniso_ranges    Array giving the anisotropy ranges
 **                             Only used when flag_aniso
 **                             (Dimension = ndim)
 ** \param[in]  aniso_rotmat    Anisotropy rotation matrix
 **                             Only used when flag_aniso && flag_rotation
 **                             (Dimension = ndim  * ndim)
 ** \param[in]  sill            Sill matrix (optional)
 **                             (Dimension = nvar * nvar)
 **
 *****************************************************************************/
int model_add_cova(Model *model,
                   const ECov &type,
                   int flag_aniso,
                   int flag_rotation,
                   double range,
                   double param,
                   const VectorDouble &aniso_ranges,
                   const VectorDouble &aniso_rotmat,
                   const VectorDouble &sill,
                   double /*ball_radius*/)
{
  CovAniso* cova = nullptr;
  if (st_check_model(model)) return 1;

  // Add a new element in the structure

  cova = new CovAniso(type, model->getContext());
  if (model->isFlagGradientFunctional() && ! cova->hasCovDerivative())
  {
    messerr("This covariance is not compatible with Functional Gradient Calculation");
    return 1;
  }

  cova->setParam(param);
  if (flag_aniso)
  {
    cova->setRanges(aniso_ranges);
    if (flag_rotation) cova->setAnisoRotation(aniso_rotmat);
  }
  else
    cova->setRange(range);

  if (static_cast<int>(sill.size()) > 0) cova->setSill(sill);
  model->addCov(cova);

  return 0;
}

/****************************************************************************/
/*!
 **  Add a basic drift
 **
 ** \return  Error return code
 **
 ** \param[in]  model      Pointer to the Model structure
 ** \param[in]  type       Type of the basic drift function (EDrift)
 ** \param[in]  rank_fex   Rank of the external drift (starting from 0)
 **
 ** \remark  If the variables are NOT linked, the number of drift equations
 ** \remark  is equal to: Number of variables * Number of drift functions
 ** \remark  If the variables are linked, the number of drift equations
 ** \remark  is equal to the number of drift functions.
 **
 *****************************************************************************/
int model_add_drift(Model *model, const EDrift &type, int rank_fex)
{
  ADriftElem *drift;
  if (st_check_model(model)) return 1;
  drift = DriftFactory::createDriftFunc(type, model->getContext());
  drift->setRankFex(rank_fex);
  model->addDrift(drift);
  delete drift;
  return 0;
}

/****************************************************************************/
/*!
 **  For a given basic structure, get the reduction factor to convert the
 **  theoretical to practical scale
 **
 ** \return  Conversion factor
 **
 ** \param[in]  type      Type of the basic structure
 ** \param[in]  param     Value of the third parameter
 **
 *****************************************************************************/
double cova_get_scale_factor(const ECov &type, double param)
{
  CovContext ctxt = CovContext(1, 1);
  ACovFunc *cova = CovFactory::createCovFunc(type, ctxt);
  cova->setParam(param);
  return cova->getScadef();
}

/*****************************************************************************/
/*!
 **  Calculate the linear model of coregionalization starting from the
 **  coregionalization matrix
 **
 ** \return  Error return code.
 **
 ** \param[in]  model    Model structure
 **
 ** \param[out]  aic     array of 'aic' values
 ** \param[out]  valpro  array of eigen values
 ** \param[out]  vecpro  array of eigen vectors
 **
 ** \remark  In case of error, the message is printed by the routine
 ** \remark  Warning: in the case of linked drift, the test of definite
 ** \remark  positiveness is bypassed as we are not in the scope of the
 ** \remark  linear model of coregionalization anymore.
 ** \remark  As a consequence the array "aic()" is not evaluated
 **
 *****************************************************************************/
int model_update_coreg(Model *model,
                       double *aic,
                       double *valpro,
                       double *vecpro)
{
  int ncova = model->getCovaNumber();
  int nvar  = model->getVariableNumber();

  /* Calculate the eigen values and vectors of the coregionalization matrix */

  for (int icov = 0; icov < ncova; icov++)
  {
    if (!is_matrix_definite_positive(
        nvar, model->getCova(icov)->getSill().getValues().data(), valpro,
        vecpro, 0))
    {
      messerr("Warning: the model is not authorized");
      messerr("The coregionalization matrix for the structure %d is not definite positive",
          icov + 1);
      return 1;
    }

    /* Calculate the factor matrix */

    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
        AIC(icov,ivar,jvar)= VECPRO(ivar,jvar) * sqrt(VALPRO(jvar));
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Calculate the value of the model for a set of distances
 **
 ** \return  Error return code
 **
 ** \param[in]  model      Model structure
 ** \param[in]  ivar       Rank of the first variable
 ** \param[in]  jvar       Rank of the second variable
 ** \param[in]  rank_sel   Rank of the basic structure (optional)
 ** \param[in]  flag_norm  1 if the model is normalized
 ** \param[in]  flag_cov   1 if the result must be given in covariance
 ** \param[in]  nugget_opt option for the nugget effect basic structure
 ** \li                     0 : no particular option
 ** \li                     1 : discard the nugget effect
 ** \li                    -1 : only consider the nugget effect
 ** \param[in]  nostd      0 standard; +-1 special; ITEST normalized
 ** \param[in]  norder     Order of the Generalized Variogram
 ** \param[in]  member     Member of the Kriging System (ECalcMember)
 ** \param[in]  nh         Number of increments
 ** \param[in]  codir      Array giving the direction coefficients
 ** \param[in]  h          Vector of increments
 **
 ** \param[out] g          Array containing the model values
 **
 ** \remark  When rank_sel is positive, it indicates the rank of the only
 ** \remark  basic structure to be accounted for
 **
 *****************************************************************************/
int model_evaluate(Model *model,
                   int ivar,
                   int jvar,
                   int rank_sel,
                   int flag_norm,
                   int flag_cov,
                   int nugget_opt,
                   int nostd,
                   int norder,
                   const ECalcMember &member,
                   int nh,
                   VectorDouble &codir,
                   double *h,
                   double *g)
{
  CovCalcMode mode;
  mode.update(member, nugget_opt, nostd, rank_sel, flag_norm, flag_cov);
  if (norder > 0) mode.setOrderVario(norder);

  /* Preliminary checks */

  if (st_check_model(model)) return 1;
  int ndim = model->getDimensionNumber();
  int nvar = model->getVariableNumber();
  if (st_check_variable(nvar, ivar)) return 1;
  if (st_check_variable(nvar, jvar)) return 1;

  /* Core allocation */

  VectorDouble d1(ndim);
  VectorDouble covtab(nvar * nvar);

  /* Normalize the direction vector codir */

  vario_fix_codir(ndim, codir);

  /* Loop on the lags */

  for (int ih = 0; ih < nh; ih++)
  {
    for (int idim = 0; idim < ndim; idim++)
      d1[idim] = h[ih] * codir[idim];
    model_calcul_cov(NULL,model, mode, 1, 1., d1, covtab.data());
    g[ih] = COVTAB(ivar, jvar);
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Calculate the value of the model for a set of distances (non stationary)
 **
 ** \return  Error return code
 **
 ** \param[in]  model      Model structure
 ** \param[in]  ivar       Rank of the first variable
 ** \param[in]  jvar       Rank of the second variable
 ** \param[in]  rank_sel   Rank of the basic structure (optional)
 ** \param[in]  flag_norm  1 if the model is normalized
 ** \param[in]  flag_cov   1 if the result must be given in covariance
 ** \param[in]  nugget_opt option for the nugget effect basic structure
 ** \li                     0 : no particular option
 ** \li                     1 : discard the nugget effect
 ** \li                    -1 : only consider the nugget effect
 ** \param[in]  nostd      0 standard; +-1 special; ITEST normalized
 ** \param[in]  norder     Order of the Generalized Variogram
 ** \param[in]  member     Member of the Kriging System (ECalcMember)
 ** \param[in]  db1        First Db structure
 ** \param[in]  iech1      First sample
 ** \param[in]  db2        Second Db structure
 ** \param[in]  iech2      Second sample
 ** \param[in]  nh         Number of increments
 ** \param[in]  codir      Array giving the direction coefficients
 ** \param[in]  h          Vector of increments
 **
 ** \param[out] g          Array containing the model values
 **
 ** \remark  When rank_sel is positive, it indicates the rank of the only
 ** \remark  basic structure to be accounted for
 **
 *****************************************************************************/
int model_evaluate_nostat(Model *model,
                          int ivar,
                          int jvar,
                          int rank_sel,
                          int flag_norm,
                          int flag_cov,
                          int nugget_opt,
                          int nostd,
                          int norder,
                          const ECalcMember &member,
                          Db *db1,
                          int iech1,
                          Db *db2,
                          int iech2,
                          int nh,
                          VectorDouble &codir,
                          double *h,
                          double *g)
{
  CovCalcMode mode;
  mode.update(member, nugget_opt, nostd, rank_sel, flag_norm, flag_cov);
  if (norder > 0) mode.setOrderVario(norder);

  /* Preliminary checks */

  if (st_check_model(model)) return 1;
  int ndim = model->getDimensionNumber();
  int nvar = model->getVariableNumber();
  if (st_check_variable(nvar, ivar)) return 1;
  if (st_check_variable(nvar, jvar)) return 1;
  CovInternal covint(1, iech1, 2, iech2, ndim, db1, db2);

  /* Core allocation */

  VectorDouble d1(ndim, 0.);
  VectorDouble covtab(nvar * nvar);

  /* Normalize the direction vector codir */

  vario_fix_codir(ndim, codir);

  /* Calculate the C(0) term (used only for covariance or covariogram) */

  double c00 = model->getContext().getCovar0(ivar, jvar);

  model_calcul_cov(&covint, model, mode, 1, 1., d1, covtab.data());
  double var0 = COVTAB(ivar, jvar);
  if (c00 <= 0. || FFFF(c00)) c00 = var0;

  /* Loop on the lags */

  for (int ih = 0; ih < nh; ih++)
  {
    for (int idim = 0; idim < ndim; idim++)
      d1[idim] = h[ih] * codir[idim];
    model_calcul_cov(&covint, model, mode, 1, 1., d1, covtab.data());
    g[ih] = COVTAB(ivar, jvar);
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Evaluate the model on a regular grid
 **
 ** \param[in]  model      Model structure
 ** \param[in]  db         Db structure
 ** \param[in]  ivar       Rank of the first variable
 ** \param[in]  jvar       Rank of the second variable
 ** \param[in]  flag_norm  1 if the model is normalized
 ** \param[in]  flag_cov   1 if the result must be given in covariance
 **
 ** \param[out] g          Array containing the model values
 **
 *****************************************************************************/
int model_grid(Model *model,
               Db *db,
               int ivar,
               int jvar,
               int flag_norm,
               int flag_cov,
               double *g)
{
  CovCalcMode mode;
  mode.update(ECalcMember::LHS, 0, 0, -1, flag_norm, flag_cov);

  /* Preliminary checks */

  if (st_check_model(model)) return 1;
  if (st_check_environ(model, db)) return 1;
  int ndim = model->getDimensionNumber();
  int nvar = model->getVariableNumber();

  /* Core allocation */

  VectorDouble d1(ndim,0.);
  VectorDouble covtab(nvar * nvar);

  /* Initialization */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
    g[iech] = TEST;

  /* Loop on the lags */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    db_sample_load(db, ELoc::X, iech, d1.data());
    model_calcul_cov(NULL,model, mode, 1, 1., d1, covtab.data());
    g[iech] = COVTAB(ivar, jvar);
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Returns the number of external drift functions
 **
 ** \return  Number of external drift functions
 **
 ** \param[in]  model Model structure
 **
 *****************************************************************************/
int model_nfex(Model *model)

{
  return model->getExternalDriftNumber();
}

/****************************************************************************/
/*!
 **  Filter a basic drift function
 **
 ** \param[in]  model   Model structure
 ** \param[in]  rank    Rank of the basic drift structure to be filtered
 **                     (numbered starting from 0)
 ** \param[in]  filter  1 to filter; 0 to unfilter
 **
 *****************************************************************************/
void model_drift_filter(Model *model, int rank, int filter)
{
  model->setDriftFiltered(rank, filter);
}

/****************************************************************************/
/*!
 **  Evaluate the average model between two Dbs
 **
 ** \return  Average model value
 **
 ** \param[in]  model Model structure
 ** \param[in]  db1   First Db
 ** \param[in]  db2   Second Db
 ** \param[in]  ivar  Rank of the first variable
 ** \param[in]  jvar  Rank of the second variable
 ** \param[in]  seed  Seed for the random number generator
 ** \param[in]  eps   Epsilon used for randomization in calculation of CVV
 **
 *****************************************************************************/
double model_cxx(Model *model,
                 Db *db1,
                 Db *db2,
                 int ivar,
                 int jvar,
                 int seed,
                 double eps)
{
  CovCalcMode mode;

  /* Initializations */

  if (st_check_model(model)) return TEST;
  if (st_check_environ(model, db1)) return TEST;
  if (st_check_environ(model, db2)) return TEST;
  int ndim = model->getDimensionNumber();
  int nvar = model->getVariableNumber();
  if (st_check_variable(nvar, ivar)) return TEST;
  if (st_check_variable(nvar, jvar)) return TEST;
  if (seed != 0) law_set_random_seed(seed);

  /* Core allocation */

  VectorDouble d1(ndim, 0.);
  VectorDouble covtab(nvar * nvar, 0.);

  /* Loop on the first sample */

  double norme = 0.;
  for (int iech1 = 0; iech1 < db1->getSampleNumber(); iech1++)
  {
    if (!db1->isActive(iech1)) continue;
    double w1 = db1->getWeight(iech1);
    if (w1 == 0.) continue;

    /* Loop on the second sample */

    for (int iech2 = 0; iech2 < db2->getSampleNumber(); iech2++)
    {
      if (!db2->isActive(iech2)) continue;
      double w2 = db2->getWeight(iech2);
      if (w2 == 0.) continue;

      /* Loop on the dimension of the space */

      int skip = 0;
      for (int i = 0; i < ndim && skip == 0; i++)
      {
        double v1 = db1->getCoordinate(iech1, i);
        double v2 = db2->getCoordinate(iech2, i);
        if (eps != 0.) v2 += eps * law_uniform(-0.5, 0.5);
        if (FFFF(v1) || FFFF(v2)) skip = 1;
        d1[i] = v1 - v2;
      }
      if (skip) continue;

      model_calcul_cov(NULL,model, mode, 0, w1 * w2, d1, covtab.data());
      norme += w1 * w2;
    }
  }

  /* Scaling */

  st_covtab_rescale(nvar, norme, covtab.data());
  return COVTAB(ivar, jvar);
}

/****************************************************************************/
/*!
 **  Establish the covariance matrix between two Dbs
 **  where samples are selected by ranks
 **
 ** \return Array containing the covariance matrix
 **
 ** \param[in]  model  Model structure
 ** \param[in]  db1    First Db
 ** \param[in]  nsize1 Number of selected samples
 ** \param[in]  ranks1 Array giving ranks of selected samples (optional)
 ** \param[in]  db2    Second Db
 ** \param[in]  nsize2 Number of selected samples
 ** \param[in]  ranks2 Array giving ranks of selected samples (optional)
 ** \param[in]  ivar0  Rank of the first variable (-1: all variables)
 ** \param[in]  jvar0  Rank of the second variable (-1: all variables)
 ** \param[in]  flag_norm  1 if the model is normalized
 ** \param[in]  flag_cov   1 if the result must be given in covariance
 **
 ** \remarks The covariance matrix (returned) must be freed by calling routine
 ** \remarks Its dimension is nsize1 * nsize2 * nvar * nvar
 ** \remarks where 'nvar' is the number of active variables (1 or nvar)
 ** \remarks The covariance matrix is established for the first variable
 ** \remarks and returned as a covariance
 ** \remarks As the ranks are used, no test is performed on any selection
 ** \remarks but only ranks positive or null are considered
 **
 *****************************************************************************/
double* model_covmat_by_ranks(Model *model,
                              Db *db1,
                              int nsize1,
                              const int *ranks1,
                              Db *db2,
                              int nsize2,
                              const int *ranks2,
                              int ivar0,
                              int jvar0,
                              int flag_norm,
                              int flag_cov)
{
  CovCalcMode mode;
  mode.update(ECalcMember::LHS, 0, 0, -1, flag_norm, flag_cov);
  if (st_check_model(model)) return nullptr;
  if (st_check_environ(model, db1)) return nullptr;
  if (st_check_environ(model, db2)) return nullptr;
  int ndim = model->getDimensionNumber();
  int nvar = model->getVariableNumber();
  int nvar1 = nvar;
  if (ivar0 >= 0)
  {
    nvar1 = 1;
    if (st_check_variable(nvar, ivar0)) return nullptr;
  }
  int nvar2 = nvar;
  if (jvar0 >= 0)
  {
    nvar2 = 1;
    if (st_check_variable(nvar, jvar0)) return nullptr;
  }
  if (ranks1 == nullptr) nvar1 = db1->getSampleNumber();
  if (ranks2 == nullptr) nvar2 = db1->getSampleNumber();

  /* Core allocation */

  VectorDouble d1(ndim,0.);
  VectorDouble covtab(nvar * nvar,0.);
  double* covmat = (double*) mem_alloc(sizeof(double) * nsize1 * nsize2, 0);
  if (covmat != nullptr)
  {

    /* Loop on the number of variables */

    int ecr = 0;
    for (int ivar = 0; ivar < nvar1; ivar++)
    {
      if (ivar0 >= 0) ivar = ivar0;

      /* Loop on the first sample */

      for (int i1 = 0; i1 < nsize1; i1++)
      {
        int iech1 = (ranks1 != nullptr) ? ranks1[i1] : i1;
        if (iech1 < 0) continue;

        /* Loop on the second variable */

        for (int jvar = 0; jvar < nvar2; jvar++)
        {
          if (jvar0 >= 0) jvar = jvar0;

          /* Loop on the second sample */

          for (int i2 = 0; i2 < nsize2; i2++)
          {
            int iech2 = (ranks2 != nullptr) ? ranks2[i2] : i2;
            if (iech2 < 0) continue;

            /* Loop on the dimension of the space */

            double value = TEST;
            int skip = 0;
            for (int i = 0; i < ndim && skip == 0; i++)
            {
              double v1 = db1->getCoordinate(iech1, i);
              double v2 = db2->getCoordinate(iech2, i);
              if (FFFF(v1) || FFFF(v2)) skip = 1;
              d1[i] = v1 - v2;
            }
            if (!skip)
            {
              model_calcul_cov(NULL, model, mode, 1, 1., d1, covtab.data());
              value = COVTAB(ivar, jvar);
            }
            covmat[ecr++] = value;
          }
        }
      }
    }
  }
  return (covmat);
}

/****************************************************************************/
/*!
 **  Establish the drift rectangular matrix for a given Db
 **
 ** \return Error return code
 **
 ** \param[in]  model  Model structure
 ** \param[in]  member Member of the Kriging System (ECalcMember)
 ** \param[in]  db     Db structure
 **
 ** \param[out] drfmat The drift matrix
 **                    (Dimension = nech * nvar * nfeq * nvar)
 **
 *****************************************************************************/
int model_drift_mat(Model *model,
                     const ECalcMember &member,
                     Db *db,
                     double *drfmat)
{
  if (st_check_model(model)) return 1;
  if (st_check_environ(model, db)) return 1;
  int nvar = model->getVariableNumber();
  int nbfl = model->getDriftNumber();
  int nfeq = model->getDriftEquationNumber();
  int nech = db->getSampleNumber();
  VectorDouble drftab(nbfl,0.);

  /* Loop on the variables */

  int ecr = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
  {

    /* Loop on the samples */

    for (int iech = 0; iech < nech; iech++)
    {
      if (!db->isActive(iech)) continue;
      model_calcul_drift(model, member, db, iech, drftab.data());

      /* Loop on the drift functions */

      if (model->isFlagLinked())
      {
        for (int ib = 0; ib < nfeq; ib++)
        {
          double value = 0.;
          for (int il = 0; il < nbfl; il++)
            value += drftab[il] * model->getCoefDrift(ivar, il, ib);
          drfmat[ecr++] = value;
        }
      }
      else
      {
        for (int jvar = 0; jvar < nvar; jvar++)
          for (int jl = 0; jl < nbfl; jl++)
          {
            int jb = jvar + nvar * jl;
            double value = 0.;
            for (int il = 0; il < nbfl; il++)
              value += drftab[il] * model->getCoefDrift(ivar, il, jb);
            drfmat[ecr++] = value;
          }
      }
    }
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Establish the drift vector for a given sample of the Db
 **
 ** \param[in]  model  Model structure
 ** \param[in]  member Member of the Kriging System (ECalcMember)
 ** \param[in]  db     Db structure
 ** \param[in]  iech   Rank of the particular sample
 **
 ** \param[out] vector Returned vector
 **                    (Dimension = nvar * nfeq)
 **
 *****************************************************************************/
int model_drift_vector(Model *model,
                        const ECalcMember &member,
                        Db *db,
                        int iech,
                        double *vector)
{
  if (st_check_model(model)) return 1;
  if (st_check_environ(model, db)) return 1;
  int nvar = model->getVariableNumber();
  int nbfl = model->getDriftNumber();
  int nfeq = model->getDriftEquationNumber();
  VectorDouble drftab(nbfl, 0.);

  /* Initialize the covariance matrix */

  for (int i = 0; i < nvar * nfeq; i++) vector[i] = TEST;

  model_calcul_drift(model, member, db, iech, drftab.data());

  int ecr = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int ib = 0; ib < nfeq; ib++)
    {
      double value = 0.;
      for (int il = 0; il < nbfl; il++)
        value += drftab[il] * model->getCoefDrift(ivar, il, ib);
      vector[ecr++] = value;
    }
  return 0;
}

/*****************************************************************************/
/*!
 **  Patches the value of the drift coefficient in the model for the
 **  rank of variable 'iv' and of the equation 'ib'. The rank of the
 **  drift function is found by matching the type in the basis of the
 **  drift functions available.
 **  Note : this function is only used when the model has linked
 **  drift functions
 **
 ** \param[in]  model Model structure
 ** \param[in]  iv    rank of the variable
 ** \param[in]  ib    rank of the equation
 ** \param[in]  type  type of the drift function (EDrift)
 ** \param[in]  rank  rank of the external drift
 ** \param[in]  value value to be added to the drift coefficient
 **
 *****************************************************************************/
static void st_drift_modify(Model *model,
                            int iv,
                            int ib,
                            const EDrift &type,
                            int rank,
                            double value)
{
  /* Look for the drift function */

  int il = -1;
  for (int i = 0; i < model->getDriftNumber() && il < 0; i++)
    if (model->getDriftType(i) == type &&
        model->getDrift(i)->getRankFex() == rank) il = i;
  if (il < 0) messageAbort("st_drift_modify");

  /* Patch the drift coefficient */

  model->setCoefDrift(iv, il, ib, model->getCoefDrift(iv, il, ib) + value);
}

/*****************************************************************************/
/*!
 **  Updates the drift component of the new_model for the variable
 **  'iv' considered as the derivative along 'mode' of the
 **  first variable in the (old) model
 **
 ** \param[in]  iv    rank of the variable
 ** \param[in]  mode  type of the derivative (::ENUM_MODEL_DERIVATIVES)
 ** \param[in]  model Model structure
 **
 ** \param[out] new_model Model structure
 **
 ** \remark  The new_model must have been dimensioned beforehand
 **
 *****************************************************************************/
static void st_drift_derivative(int iv,
                                int mode,
                                const Model *model,
                                Model *new_model)

{
  for (int ib = 0; ib < model->getDriftEquationNumber(); ib++)
    for (int il = 0; il < model->getDriftNumber(); il++)
    {
      double value = model->getCoefDrift(0, il, ib);
      if (value == 0) continue;

      EDrift type = model->getDriftType(il);
      int rank = model->getRankFext(il);
      switch (mode)
      {
        case MODEL_DERIVATIVE_NONE: /* Simple copy */
          st_drift_modify(new_model, iv, ib, type, rank, value);
          break;

        case MODEL_DERIVATIVE_X: /* Derivative along X */
          switch (type.toEnum())
          {
            case EDrift::E_X:
              st_drift_modify(new_model, iv, ib, EDrift::UC, 0, value);
              break;
            case EDrift::E_X2:
              st_drift_modify(new_model, iv, ib, EDrift::X, 0, 2. * value);
              break;
            case EDrift::E_XY:
              st_drift_modify(new_model, iv, ib, EDrift::Y, 0, value);
              break;
            case EDrift::E_XZ:
              st_drift_modify(new_model, iv, ib, EDrift::Z, 0, value);
              break;
            case EDrift::E_X3:
              st_drift_modify(new_model, iv, ib, EDrift::X2, 0, 3. * value);
              break;
            case EDrift::E_X2Y:
              st_drift_modify(new_model, iv, ib, EDrift::XY, 0, 2. * value);
              break;
            case EDrift::E_XY2:
              st_drift_modify(new_model, iv, ib, EDrift::Y2, 0, value);
              break;
            default:
              break;
          }
          break;

        case MODEL_DERIVATIVE_Y: /* Derivative along Y */
          switch (type.toEnum())
          {
            case EDrift::E_Y:
              st_drift_modify(new_model, iv, ib, EDrift::UC, 0, value);
              break;
            case EDrift::E_Y2:
              st_drift_modify(new_model, iv, ib, EDrift::Y, 0, 2. * value);
              break;
            case EDrift::E_XY:
              st_drift_modify(new_model, iv, ib, EDrift::X, 0, value);
              break;
            case EDrift::E_YZ:
              st_drift_modify(new_model, iv, ib, EDrift::Z, 0, value);
              break;
            case EDrift::E_Y3:
              st_drift_modify(new_model, iv, ib, EDrift::Y2, 0, 3. * value);
              break;
            case EDrift::E_XY2:
              st_drift_modify(new_model, iv, ib, EDrift::XY, 0, 2. * value);
              break;
            case EDrift::E_X2Y:
              st_drift_modify(new_model, iv, ib, EDrift::X2, 0, value);
              break;
            default:
              break;
          }
          break;

        default:
          break;
      }
    }
}

/****************************************************************************/
/*!
 **  Duplicates a Model from another Model
 **
 ** \return  The modified Model structure
 **
 ** \param[in]  model       Input Model
 ** \param[in]  ball_radius Radius for Gradient calculation
 ** \param[in]  mode        Type of transformation
 ** \li                     -1 for Data (SK)
 ** \li                      0 for Data
 ** \li                      1 for Data - Gradient
 **
 *****************************************************************************/
Model* model_duplicate(const Model *model, double ball_radius, int mode)

{
  Model *new_model;
  const CovAniso *cova;
  const ADriftElem *drft;
  int new_nvar, nfact;
  double sill;
  bool flag_linked, flag_gradient;

  // Preliminary checks

  new_model = nullptr;
  int nvar = model->getVariableNumber();
  int ndim = model->getDimensionNumber();
  int ncova = model->getCovaNumber();
  int nbfl = model->getDriftNumber();
  nfact = new_nvar = 0;
  flag_linked = false;
  flag_gradient = false;

  // Create the new model (linked drift functions)

  switch (mode)
  {
    case -1:                    // Data (SK)
    case 0:                     // Data
      new_nvar = nvar;
      nfact = 1;
      flag_linked = false;
      flag_gradient = false;
      break;

    case 1:                     // Data - Gradient
      if (nvar != 1 || ndim != 2)
      {
        messerr("This procedure is limited to a single variable in 2-D");
        return new_model;
      }
      new_nvar = 3;
      nfact = 6;
      flag_linked = true;
      flag_gradient = true;
      break;
  }
  CovContext ctxt(model->getContext());
  ctxt.setNVar(new_nvar);
  new_model = new Model(ctxt);
  if (new_model == nullptr) return new_model;

  // **************************************
  // Create the basic covariance structures
  // **************************************

  ACovAnisoList* covs = nullptr;
  if (flag_gradient)
    covs = new CovLMGradient();
  else
    covs = new CovLMC();

  int lec = 0;
  for (int icov = 0; icov < ncova; icov++)
  {
    cova = model->getCova(icov);
    sill = model->getSill(icov, 0, 0);
    for (int ifact = 0; ifact < nfact; ifact++, lec++)
    {
      CovAniso* covnew = nullptr;
      if (flag_gradient)
        covnew = new CovGradientNumerical(cova->getType(),ball_radius,ctxt);
      else
        covnew = new CovAniso(cova->getType(), ctxt);
      covnew->setParam(cova->getParam());
      if (cova->getFlagAniso())
      {
        covnew->setRanges(cova->getRanges());
        if (cova->getFlagRotation())
          covnew->setAnisoRotation(cova->getAnisoRotation());
      }
      else
        covnew->setRange(cova->getRange());

      /* Modify the Sill */;

      switch (mode)
      {
        case 0:
        case -1:
          for (int ivar = 0; ivar < new_nvar; ivar++)
            for (int jvar = 0; jvar < new_nvar; jvar++)
              covnew->setSill(ivar, jvar, cova->getSill(ivar,jvar));
          break;

        case 1:                   // Data - Gradient
          covnew->initSill(0.);
          if (ifact == 0)
          {
            covnew->setSill(0, 0,  sill);
          }
          else if (ifact == 1)
          {
            covnew->setSill(0, 1, -sill);
            covnew->setSill(1, 0,  sill);
          }
          else if (ifact == 2)
          {
            covnew->setSill(1, 1,  sill);
          }
          else if (ifact == 3)
          {
            covnew->setSill(0, 2, -sill);
            covnew->setSill(2, 0,  sill);
          }
          else if (ifact == 4)
          {
            covnew->setSill(1, 2, -sill);
            covnew->setSill(2, 1, -sill);
          }
          else if (ifact == 5)
          {
            covnew->setSill(2, 2, sill);
          }
          else
          {
            my_throw("Argument 'ifact' invalid");
          }
          break;
      }
      covs->addCov(covnew);
      delete covnew;
    }
  }
  new_model->setCovList(covs);
  delete covs;

  // *********************************
  // Create the basic drift structures
  // *********************************

  if (mode >= 0)
  {
    DriftList drifts = DriftList(ctxt.getSpace());
    drifts.setFlagLinked(flag_linked);
    for (int il = 0; il < nbfl; il++)
    {
      drft = model->getDrift(il);
      ADriftElem *newdrft = DriftFactory::createDriftFunc(drft->getType(), ctxt);
      newdrft->setRankFex(drft->getRankFex());
      drifts.addDrift(newdrft);
      delete newdrft;
      drifts.setFiltered(il, model->isDriftFiltered(il));
    }
    new_model->setDriftList(&drifts);

    // Update the drift for the derivatives

    if (mode == 1)
    {
      int nval = new_nvar * new_model->getDriftEquationNumber()
                 * new_model->getDriftNumber();
      for (int i = 0; i < nval; i++)
        new_model->setCoefDriftByRank(i, 0.);
      st_drift_derivative(0, MODEL_DERIVATIVE_NONE, model, new_model);
      st_drift_derivative(1, MODEL_DERIVATIVE_X, model, new_model);
      st_drift_derivative(2, MODEL_DERIVATIVE_Y, model, new_model);
    }
  }

  return (new_model);
}

/****************************************************************************/
/*!
 **  Normalize the model
 **
 ** \param[in]  model         Model structure
 ** \param[in]  flag_verbose  1 for a verbose output
 **
 *****************************************************************************/
int model_normalize(Model *model, int flag_verbose)

{
  int nvar = model->getVariableNumber();
  int ncov = model->getCovaNumber();
  VectorDouble total(nvar,0.);

  /* Calculate the total sills for each variable */

  int flag_norm = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    total[ivar] = model->getTotalSill(ivar, ivar);
    if (total[ivar] == 0.) return 1;
    total[ivar] = sqrt(total[ivar]);
    if (ABS(total[ivar] - 1.) > EPSILON6) flag_norm = 1;
  }

  /* Scale the different sills for the different variables */

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      for (int icov = 0; icov < ncov; icov++)
      {
        double sill = model->getCova(icov)->getSill(ivar, jvar);
        sill /= total[ivar] * total[jvar];
        model->getCova(icov)->setSill(ivar, jvar, sill);
      }

  /* Printout */

  if (flag_verbose && flag_norm)
  {
    message("The model has been normalized\n");
    for (int ivar = 0; ivar < nvar; ivar++)
      message("- Variable %d : Scaling factor = %lf\n", ivar + 1,
              total[ivar] * total[ivar]);
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Stabilize the model (in the monovariate case)
 **
 ** \return  Error returned code
 **
 ** \param[in]  model         Model structure
 ** \param[in]  flag_verbose  1 for a verbose output
 ** \param[in]  percent       Percentage of nugget effect added
 **
 ** \remark  If the model only contains GAUSSIAN structures, add
 ** \remark  a NUGGET EFFECT structure with a sill equal to a percentage
 ** \remark  of the total sill of the GAUSSIAN component(s)
 **
 ** \remark  This function does not do anything in the multivariate case
 **
 *****************************************************************************/
int model_stabilize(Model *model, int flag_verbose, double percent)
{
  int nvar = model->getVariableNumber();
  if (nvar > 1) return 0;
  if (percent <= 0.) return 0;
  int ncov = model->getCovaNumber();

  /* Check if the model only contains GAUSSIAN components */

  double total = 0.;
  for (int icov = 0; icov < ncov; icov++)
  {
    CovAniso* cova = model->getCova(icov);
    if (cova->getType() != ECov::GAUSSIAN) return (0);
    total += model->getSill(icov, 0, 0);
  }
  total = total * percent / 100.;

  /* Update each Gaussian component */

  for (int icov = 0; icov < ncov; icov++)
    model->getCova(icov)->setSill(0, 0, 1. - total);

  /* Add a NUGGET EFFECT component */

  if (model_add_cova(model, ECov::NUGGET, 0, 0, 0., 0., VectorDouble(),
                     VectorDouble(), VectorDouble(1, total), 0.)) return 1;

  /* Printout */

  if (flag_verbose)
  {
    message("The model which only contains Gaussian components\n");
    message("has been stabilized by adding a small Nugget Effect\n");
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Update the model for fitting Covariance or Covariogram
 **
 ** \param[in]  model         Model structure
 ** \param[in]  c0            Array of variance values at the origin
 ** \param[in]  flag_verbose  1 for verbose output
 **
 ** \param[out] flag_nugget  1 if a nugget component must be added
 ** \param[out] nugget       Array of sills for the nugget component
 **
 *****************************************************************************/
void model_covupdt(Model *model,
                   double *c0,
                   int flag_verbose,
                   int *flag_nugget,
                   double *nugget)
{
  /// TODO : dead code ?
  CovAniso *cova;
  double *silltot, *range, diff;
  int i, icov, jcov, nvar, ncova, rank_nugget, rank_exceed, ivar, jvar;
  int *rank, flag_update, flag_rescale;

  /* Initializations */

  silltot = range = nullptr;
  rank = nullptr;
  nvar = model->getVariableNumber();
  ncova = model->getCovaNumber();
  flag_update = flag_rescale = 0;

  /* Core allocation */

  rank = (int*) mem_alloc(sizeof(int) * ncova, 1);
  range = (double*) mem_alloc(sizeof(double) * ncova, 1);
  silltot = (double*) mem_alloc(sizeof(double) * nvar * nvar, 1);
  for (i = 0; i < nvar * nvar; i++)
    silltot[i] = 0.;

  /* Sort the basic structures by increasing range */
  rank_nugget = -1;
  for (icov = 0; icov < ncova; icov++)
  {
    cova = model->getCova(icov);
    if (cova->getType() == ECov::NUGGET) rank_nugget = icov;
    rank[icov] = icov;
    range[icov] = cova->getRange();
  }
  ut_sort_double(0, ncova, rank, range);

  /* Loop on the basic structures, in order to : */
  /* - cumulate the sills (excluding the nugget effect component) */
  /* - find the rank of the structure which exceeds the total variance */

  rank_exceed = -1;
  for (jcov = 0; jcov < ncova && rank_exceed < 0; jcov++)
  {
    icov = rank[ncova - 1 - jcov];
    cova = model->getCova(icov);
    if (cova->getType() == ECov::NUGGET) continue;
    for (ivar = 0; ivar < nvar; ivar++)
    {
      silltot[AD(ivar, ivar)] += model->getSill(icov, ivar, ivar);
      if (silltot[AD(ivar, ivar)] > c0[AD(ivar, ivar)]) rank_exceed = icov;
    }
  }

  if (rank_exceed >= 0)
  {
    flag_rescale = (rank_exceed == 0);
    if (flag_rescale) rank_nugget = rank_exceed;
    if (flag_verbose)
    {
      message("Error in the Covariance or Covariogram Model\n");
      message("The cumulated sill exceeds the experimental C(0)\n");

      if (rank_exceed > 0)
      {
        message("The following basic structures are discarded : ");
        for (jcov = rank_exceed; jcov < ncova; jcov++)
        {
          icov = rank[ncova - 1 - jcov];
          message(" #%d", icov + 1);
        }
        message("\n");
      }
      else
      {
        message("All the structures are discarded\n");
        message("except the structure #%d which is rescaled\n",
                rank[ncova - 1 - rank_exceed] + 1);
      }
    }

    /* Discard the exceeded basic structures */

    for (jcov = rank_exceed; jcov < ncova; jcov++)
    {
      icov = rank[ncova - 1 - jcov];
      cova = model->getCova(icov);
      if (cova->getType() == ECov::NUGGET) continue;
      for (ivar = 0; ivar < nvar; ivar++)
        for (jvar = 0; jvar < nvar; jvar++)
          model->getCova(icov)->setSill(ivar, jvar, 0.);
    }

    /* Update the cumulated sill */

    for (i = 0; i < nvar * nvar; i++)
      silltot[i] = 0.;
    for (jcov = 0; jcov < ncova; jcov++)
    {
      icov = rank[ncova - 1 - jcov];
      cova = model->getCova(icov);
      if (cova->getType() == ECov::NUGGET) continue;
      for (ivar = 0; ivar < nvar; ivar++)
        silltot[AD(ivar, ivar)] += model->getSill(icov, ivar, ivar);
    }
  }

  /* Calculate the additional nugget effect */
  for (ivar = 0; ivar < nvar; ivar++)
  {
    diff = c0[AD(ivar, ivar)] - silltot[AD(ivar, ivar)];
    if (diff > 0) flag_update = 1;
    for (jvar = 0; jvar < nvar; jvar++)
    {
      if (rank_nugget >= 0)
        model->setSill(rank_nugget, ivar, jvar, (ivar == jvar) ? diff : 0.);
      else
        nugget[AD(ivar, jvar)] = (ivar == jvar) ? diff : 0.;
    }
  }

  /* Returning arguments */

  rank = (int*) mem_free((char* ) rank);
  range = (double*) mem_free((char* ) range);
  silltot = (double*) mem_free((char* ) silltot);
  *flag_nugget = flag_update && (rank_nugget < 0);
  if (flag_verbose && (*flag_nugget))
  {
    message(
        "A Nugget Effect component is added so as to match the experimental variance\n");
  }
  return;
}

/****************************************************************************/
/*!
 **  Evaluate the drift with a given set of coefficients
 **
 ** \param[in]  model   Model structure
 ** \param[in]  db      Db structure
 ** \param[in]  iech    Rank of the sample
 ** \param[in]  ivar    Rank of the variable
 ** \param[in]  coef    Array of coefficients (optional)
 **
 ** \param[out] drftab  Working array
 **
 *****************************************************************************/
double model_drift_evaluate(int /*verbose*/,
                            Model *model,
                            const Db *db,
                            int iech,
                            int ivar,
                            double *coef,
                            double *drftab)
{
  if (st_check_environ(model, db)) return TEST;

  model_calcul_drift(model, ECalcMember::LHS, db, iech, drftab);

  /* Check if all the drift terms are defined */

  for (int il = 0; il < model->getDriftNumber(); il++)
    if (FFFF(drftab[il])) return (TEST);

  /* Perform the correction */

  double drift = 0.;
  for (int ib = 0; ib < model->getDriftEquationNumber(); ib++)
  {
    double value = 0.;
    for (int il = 0; il < model->getDriftNumber(); il++)
      value += drftab[il] * model->getCoefDrift(ivar, il, ib);
    drift += value * coef[ib];
  }
  return (drift);
}

/****************************************************************************/
/*!
 **  Returns the characteristics of the covariance
 **
 ** \param[in]  type           Type of the covariance
 **
 ** \param[out] cov_name       Name of the covariance
 ** \param[out] flag_range     range definition
 ** \li                         +1 if the range is defined
 ** \li                         -1 if the range is redundant with the sill
 ** \param[out] flag_param     1 if the third parameter is defined
 ** \param[out] min_order      Minimum IRF order for validity
 ** \param[out] max_ndim       Maximum dimension for validity
 ** \param[out] flag_int_1d    Integral range in 1-D
 ** \param[out] flag_int_2d    Integral range in 2-D
 ** \param[out] flag_aniso     1 if anisotropy is meaningful
 ** \param[out] flag_rotation  1 if an anisotropy rotation is meaningful
 ** \param[out] scale          Scaling parameter
 ** \param[out] parmax         Maximum value for the third parameter
 **
 *****************************************************************************/
void model_cova_characteristics(const ECov &type,
                                char cov_name[STRING_LENGTH],
                                int *flag_range,
                                int *flag_param,
                                int *min_order,
                                int *max_ndim,
                                int *flag_int_1d,
                                int *flag_int_2d,
                                int *flag_aniso,
                                int *flag_rotation,
                                double *scale,
                                double *parmax)
{
  SpaceRN space(1); // Use 1-D in order to retrieve all covariances
  CovContext ctxt = CovContext(1, 1);
  ACovFunc *cov = CovFactory::createCovFunc(type, ctxt);
  (void) gslStrcpy((char*) cov_name, cov->getCovName().c_str());
  *flag_range = cov->hasRange();
  *flag_param = cov->hasParam();
  *min_order = cov->getMinOrder();
  *max_ndim = cov->getMaxNDim();
  *flag_int_1d = cov->hasInt1D();
  *flag_int_2d = cov->hasInt2D();
  *flag_aniso = (((*flag_range) != 0) && (*max_ndim < 0 || *max_ndim > 1));
  *flag_rotation = ((*flag_aniso) && (*max_ndim < 0 || *max_ndim > 1));
  *scale = cov->getScadef();
  *parmax = cov->getParMax();
  return;
}

/*****************************************************************************/
/*!
 **  Calculates variogram values by sampling a model
 **
 ** \return  Error return code
 **
 ** \param[in]  vario     Vario structure
 ** \param[in]  model     Model structure
 ** \param[in]  flag_norm  1 if the model is normalized
 ** \param[in]  flag_cov  1 if the result must be given in covariance
 **
 *****************************************************************************/
int model_sample(Vario *vario, Model *model, int flag_norm, int flag_cov)
{
  int ndim = vario->getDimensionNumber();
  int ndir = vario->getDirectionNumber();
  int nvar = model->getVariableNumber();
  CovCalcMode mode;
  mode.update(ECalcMember::LHS, 0, 0, -1, flag_norm, flag_cov);

  /* Core allocation */

  VectorDouble d1(ndim,0.);
  VectorDouble covtab(nvar * nvar, 0.);
  vario->setNVar(nvar);

  // Internal dimensioning

  vario->internalVariableResize();
  vario->internalDirectionResize();

  /* Calculate the C(0) constant term */

  model_calcul_cov(NULL,model, mode, 1, 1., VectorDouble(), covtab.data());
  for (int i = 0; i < nvar * nvar; i++)
    vario->setVarIndex(i, covtab[i]);

  /* Loop on the directions */

  for (int idir = 0; idir < ndir; idir++)
  {

    /* Loop on the variogram lags */

    int npas = vario->getLagNumber(idir);
    for (int ipas = 0; ipas < npas; ipas++)
    {

      /* Loop on the variables */

      int ijvar = 0;
      for (int ivar = 0; ivar < vario->getVariableNumber(); ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {
          int i = vario->getDirAddress(idir, ivar, jvar, ipas, false, 0);
          vario->setSwByIndex(idir, i, 1.);
          vario->setHhByIndex(idir, i, ipas * vario->getDPas(idir));
          for (int idim = 0; idim < ndim; idim++)
            d1[idim] = vario->getHhByIndex(idir, i) * vario->getCodir(idir, idim);
          model_calcul_cov(NULL,model, mode, 1, 1., d1, covtab.data());
          vario->setGgByIndex(idir, i, COVTAB(ivar, jvar));
        }
    }
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Combine two monovariate models into a bivariate model (residuals model)
 **
 ** \return  Pointer to the newly created Model structure
 **
 ** \param[in]  model1      First input Model
 ** \param[in]  model2      Second input Model
 ** \param[in]  r           Correlation coefficient
 **
 ** \remarks: The drift is not copied into the new model
 ** \remarks: It has been extended to the case where only one model is defined
 **
 *****************************************************************************/
Model* model_combine(const Model *model1, const Model *model2, double r)
{
  Model *model;

  if (model1 == nullptr)
  {
    messerr("This function requires at least one model defined");
    return nullptr;
  }
  if (model1 != nullptr && model1->getVariableNumber() != 1)
  {
    messerr("This function can only combine monovariate models");
    return nullptr;
  }
  if (model2 != nullptr && model2->getVariableNumber() != 1)
  {
    messerr("This function can only combine monovariate models");
    return nullptr;
  }
  if (model1 == nullptr)
  {
    model = model2->duplicate();
    return model;
  }
  if (model2 == nullptr)
  {
    model = model1->duplicate();
    return model;
  }
  if (model1->getDimensionNumber() != model2->getDimensionNumber())
  {
    messerr("The two models to be combined must share the space dimension");
    return nullptr;
  }
  if (model1->isFlagLinked() || model2->isFlagLinked())
  {
    messerr("This function cannot combine models with linked drifts");
    return nullptr;
  }

  /* Create the output model */

  double field = MAX(model1->getField(), model2->getField());

  VectorDouble mean(2);
  VectorDouble cova0(4);
  VectorDouble sill(4);
  mean[0] = model1->getContext().getMean(0);
  mean[1] = model2->getContext().getMean(0);
  cova0[0] = 1.;
  cova0[1] = r;
  cova0[2] = r;
  cova0[3] = 1.;

  // Creating the context
  CovContext ctxt = CovContext(2, model1->getDimensionNumber(),
                               field, mean, cova0);

  // Creating the new Model
  model = new Model(ctxt);

  /* Add the covariance of the first Model */

  int ncov = 0;
  for (int i = 0; i < model1->getCovaNumber(); i++)
  {
    const CovAniso* cova = model1->getCova(i);
    sill[0] = cova->getSill(0, 0);
    sill[1] = sill[2] = r * cova->getSill(0, 0);
    sill[3] = r * r * cova->getSill(0, 0);
    if (model_add_cova(model, cova->getType(), cova->getFlagAniso(),
                       cova->getFlagRotation(), cova->getRange(),
                       cova->getParam(), cova->getRanges(),
                       cova->getAnisoRotMatVec(), sill, 0.))
    {
      delete model;
      return nullptr;
    }
    ncov++;
  }

  /* Add the covariance of the second Model */

  for (int i = 0; i < model2->getCovaNumber(); i++)
  {
    const CovAniso* cova = model2->getCova(i);
    sill[0] = 0.;
    sill[1] = sill[2] = 0.;
    sill[3] = (1. - r * r) * cova->getSill(0, 0);
    if (model_add_cova(model, cova->getType(), cova->getFlagAniso(),
                       cova->getFlagRotation(), cova->getRange(),
                       cova->getParam(), cova->getRanges(),
                       cova->getAnisoRotMatVec(), sill,0.))
    {
      delete model;
      return nullptr;
    }

    ncov++;
  }
  return model;
}

/****************************************************************************/
/*!
 **  Returns the number of structures (different from Bugget)
 **
 ** \return  Number of structures (Nugget excluded)
 **
 ** \param[in]  model   Model structure
 **
 *****************************************************************************/
int model_get_nonugget_cova(Model *model)

{
  int nstruc = 0;
  for (int icov = 0; icov < model->getCovaNumber(); icov++)
  {
    CovAniso* cova = model->getCova(icov);
    if (cova->getType() != ECov::NUGGET) nstruc++;
  }
  return nstruc;
}

/****************************************************************************/
/*!
 **  Calculate the regularized model as an experimental variogram
 **
 ** \return  Error return code
 **
 ** \param[in]  model     Model structure
 ** \param[in]  vario     Vario structure
 ** \param[in]  db        Db discretization grid structure
 **
 *****************************************************************************/
int model_regularize(Model *model,
                     Vario *vario,
                     Db *db,
                     int /*opt_norm*/,
                     double /*nug_ratio*/)
{
  CovCalcMode mode;
  if (st_check_model(model)) return 1;
  if (st_check_environ(model, db)) return 1;
  int ndim = model->getDimensionNumber();
  int nvar = model->getVariableNumber();

  /* Preliminary checks */

  if (!is_grid(db))
  {
    messerr("This calculation facility is dedicated to grid architecture");
    return 1;
  }
  int nech = db->getSampleNumber();
  int norme = nech * nech;
  vario->setNVar(nvar);
  vario->internalVariableResize();
  vario->internalDirectionResize();

  /* Core allocation */

  VectorDouble dd(ndim, 0.);
  VectorDouble c00tab(nvar * nvar, 0);
  VectorDouble covtab(nvar * nvar, 0);

  /* Calculate the Cvv (for a zero-shift) */

  model_covtab_init(1, model, c00tab.data());
  for (int iech = 0; iech < nech; iech++)
    for (int jech = 0; jech < nech; jech++)
    {
      for (int idim = 0; idim < ndim; idim++)
        dd[idim] = db->getDistance1D(iech, jech, idim);
      model_calcul_cov(NULL,model, mode, 0, 1, dd, c00tab.data());
    }
  st_covtab_rescale(nvar, norme, c00tab.data());

  /* Initialize the variance array */

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      vario->setVar(ivar, jvar, C00TAB(ivar, jvar));

  /* Loop on the directions */

  for (int idir = 0; idir < vario->getDirectionNumber(); idir++)
  {

    /* Loop on the number of lags */

    for (int ipas = 0; ipas < vario->getLagNumber(idir); ipas++)
    {
      model_covtab_init(1, model, covtab.data());
      double dist = ipas * vario->getDPas(idir);

      for (int iech = 0; iech < nech; iech++)
        for (int jech = 0; jech < nech; jech++)
        {
          for (int idim = 0; idim < ndim; idim++)
          {
            double v1 = db->getCoordinate(iech, idim);
            double v2 = db->getCoordinate(jech, idim)
                + dist * vario->getCodir(idir, idim);
            dd[idim] = v1 - v2;
          }
          model_calcul_cov(NULL,model, mode, 0, 1, dd, covtab.data());
        }
      st_covtab_rescale(nvar, norme, covtab.data());

      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++)
        {
          int iad = vario->getDirAddress(idir, ivar, jvar, ipas, false, 0);
          vario->setGgByIndex(idir, iad, C00TAB(ivar,jvar)- COVTAB(ivar,jvar));
          vario->setHhByIndex(idir, iad, dist);
          vario->setSwByIndex(idir, iad, 1);
        }
    }
  }
  return 0;
}

/*****************************************************************************/
/*!
 **  Establish and invert a covariance matrix using Incomplete Cholesky method
 **
 ** \param[in]  verbose    Verbose option
 ** \param[in]  db         Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  npivot_max Maximum number of pivots (or 0)
 ** \param[in]  eta        Precision (or TEST)
 ** \param[in]  nsize1     Number of pivots already selected
 ** \param[in]  ranks1     Ranks of pivots already selected
 ** \param[in]  center     Ooptioanl Centering point (for increments)
 ** \param[in]  flag_sort  Reordering flag (see remarks)
 **
 ** \param[out] npivot_arg Number of pivots
 ** \param[out] Pret       Array of indices of the retained samples (from 1)
 **                        Dimension: nech
 ** \param[out] Gret       Rectangular matrix
 **                        Dimension: nech * npivot_arg
 **
 ** \remark The output arrays Pret and Gret should be freed by calling function
 **
 ** \remark The array G contains as many lines as there are samples
 ** \remark If flag_sort = FALSE, the first lines concentrate on pivots,
 ** \remark   and the other points are located afterwards
 ** \remark If flag_sort = TRUE, the lines are sorted in the same order as the
 ** \remark   initial set of samples
 **
 ** \remark The incomplete Cholsky algorithm stops when either the next pivot
 ** \remark value is below 'eta' or when maximum number of pivots 'npivot_max'
 ** \remark has been reached
 **
 ** \remark If the center point is provided in 'center', the calculations
 ** \remark of covariance of increments are calculated instead. Then 'center'
 ** \remark must provide the coordinates of the origin point.
 **
 *****************************************************************************/
int model_covmat_inchol(int verbose,
                        Db *db,
                        Model *model,
                        double eta,
                        int npivot_max,
                        int nsize1,
                        int *ranks1,
                        double *center,
                        int flag_sort,
                        int *npivot_arg,
                        int **Pret,
                        double **Gret)
{
  int *pvec, i, j, npivot, jstar, nech, error, flag_incr;
  double *G, *Gmatrix, *diag, *crit, g, residual, maxdiag, tol, b, c00;
  VectorDouble d1;
  CovCalcMode mode;

  error = 1;
  nech = db->getSampleNumber();
  pvec = nullptr;
  diag = crit = G = Gmatrix = nullptr;
  flag_incr = (center != nullptr);

  if (npivot_max <= 0) npivot_max = nech;
  npivot_max = MIN(npivot_max, nech);
  d1.resize(db->getNDim());
  diag = (double*) mem_alloc(sizeof(double) * nech, 0);
  if (diag == nullptr) goto label_end;
  crit = (double*) mem_alloc(sizeof(double) * (1 + nech), 0);
  if (crit == nullptr) goto label_end;
  pvec = (int*) mem_alloc(sizeof(int) * nech, 0);
  if (pvec == nullptr) goto label_end;
  model_calcul_cov(NULL,model, mode, 1, 1., VectorDouble(), &c00);
  for (i = 0; i < nech; i++)
    pvec[i] = i;

  residual = 0.;
  for (i = 0; i < nech; i++)
  {
    if (flag_incr)
    {
      double covar2;

      for (int idim = 0; idim < 3; idim++)
        d1[idim] = db->getCoordinate(pvec[i], idim) - center[idim];
      model_calcul_cov(NULL,model, mode, 1, 1., d1, &covar2);
      diag[i] = 2. * (c00 - covar2);
    }
    else
    {
      diag[i] = c00;
    }
    residual += diag[i];
  }
  tol = (!FFFF(eta)) ? eta * residual : 0.;
  jstar = npivot = 0;

  // Main loop

  while ((residual > tol) && (npivot < npivot_max))
  {
    // Initialize and add a new zeros column to matrix G[]
    G = (double*) mem_realloc((char* ) G, (npivot + 1) * nech * sizeof(double),
                              0);
    if (G == nullptr) goto label_end;
    for (i = 0; i < nech; i++)
      G(npivot,i) = 0.;

    // Find best new element jstar (index of maximum along diagonal)
    jstar = 0;
    if (npivot < nsize1)
    {
      jstar = ranks1[npivot];
    }
    else if (npivot != 0)
    {
      maxdiag = 0.0;
      for (i = npivot; i < nech; i++)
      {
        if (diag[i] > maxdiag)
        {
          jstar = i;
          maxdiag = diag[i];
        }
      }
    }

    // Update permutation pvec (not necessary if jstar = npivot)
    if (jstar != npivot)
    {
      i = pvec[jstar];
      pvec[jstar] = pvec[npivot];
      pvec[npivot] = i;
      diag[npivot] = diag[jstar];

      // Update rows elements on G
      for (j = 0; j <= npivot; j++)
      {
        g = G(j, jstar);
        G(j,jstar) = G(j, npivot);
        G(j,npivot) = g;
      }
    }

    // Calculate the diagonal element of G
    G(npivot,npivot) = sqrt(diag[jstar]);

    // Calculate the new column of G
    for (i = npivot + 1; i < nech; i++)
    {
      if (flag_incr)
      {
        double covar1, covar2, covar3;

        (void) distance_intra(db, pvec[i], pvec[npivot], d1.data());
        model_calcul_cov(NULL,model, mode, 1, 1., d1, &covar1);

        for (int idim = 0; idim < 3; idim++)
          d1[idim] = db->getCoordinate(pvec[npivot], idim) - center[idim];
        model_calcul_cov(NULL,model, mode, 1, 1., d1, &covar2);

        for (int idim = 0; idim < 3; idim++)
          d1[idim] = db->getCoordinate(pvec[i], idim) - center[idim];
        model_calcul_cov(NULL,model, mode, 1, 1., d1, &covar3);

        G(npivot,i) = covar1 - covar2 - covar3 + c00;
      }
      else
      {
        // Calculate the covariance column C(:, npivot)
        (void) distance_intra(db, pvec[i], pvec[npivot], d1.data());
        model_calcul_cov(NULL,model, mode, 1, 1., d1, &G(npivot, i));
      }
    }
    if (npivot != 0)
    {
      for (i = npivot + 1; i < nech; i++)
        for (j = 0; j < npivot; j++)
          G(npivot,i) -= G(j,i) * G(j, npivot);
    }
    for (i = npivot + 1; i < nech; i++)
      G(npivot,i) /= G(npivot, npivot);

    // Updates diagonal elements
    for (i = npivot + 1; i < nech; i++)
    {
      if (flag_incr)
      {
        double covar2;

        for (int idim = 0; idim < 3; idim++)
          d1[idim] = db->getCoordinate(pvec[i], idim) - center[idim];
        model_calcul_cov(NULL,model, mode, 1, 1., d1, &covar2);

        b = 2. * (c00 - covar2);
      }
      else
      {
        b = c00;
      }
      for (j = 0; j <= npivot; j++)
        b -= G(j,i) * G(j, i);
      diag[i] = b;
    }

    // Save the new residual element
    residual = 0.;
    for (i = npivot + 1; i < nech; i++)
      residual += diag[i];
    crit[npivot] = diag[npivot] + residual;
    npivot++;
  }

  // Last column
  if (npivot == nech - 1)
  {
    G = (double*) mem_realloc((char* ) G, (npivot + 1) * nech * sizeof(double),
                              0);
    if (G == nullptr) goto label_end;
    for (i = 0; i < nech; i++)
      G(npivot,i) = 0.;
    G(npivot,npivot) = sqrt(diag[npivot]);
    crit[npivot] = 0.;
    npivot++;
  }

  // Return arguments
  *npivot_arg = npivot;

  // Normalize the criterion
  for (i = 0; i < npivot; i++)
    crit[i] /= (double) nech;

  // Reorder the output G matrix
  Gmatrix = (double*) mem_alloc(npivot * nech * sizeof(double), 0);
  if (Gmatrix == nullptr) goto label_end;
  for (j = 0; j < npivot; j++)
    for (i = 0; i < nech; i++)
    {
      if (flag_sort)
        Gmatrix(pvec[i],j) = G(j, i);
      else
        Gmatrix(i,j) = G(j, i);
    }
  *Gret = Gmatrix;

  // Renumber starting from 1
  for (i = 0; i < nech; i++)
    pvec[i]++;
  *Pret = pvec;

  // Printout of the order of the retained samples
  if (verbose)
  {
    message("Number of pivots = %d\n", npivot);
    print_imatrix("Order", 0, 1, 1, npivot, NULL, pvec);
    print_matrix("Criterion", 0, 1, 1, npivot, NULL, crit);
  }

  /* Set the error return code */

  error = 0;

  /* Core deallocation */

  label_end: diag = (double*) mem_free((char* ) diag);
  crit = (double*) mem_free((char* ) crit);
  G = (double*) mem_free((char* ) G);
  return (error);
}

/*****************************************************************************/
/*!
 **  Returns the st. dev. at a given increment for a given model
 **
 ** \param[in]  model       Structure containing the model
 ** \param[in]  db1         First Db
 ** \param[in]  iech1       Rank in the first Db
 ** \param[in]  iech2       Rank in the second Db
 ** \param[in]  verbose     Verbose flag
 ** \param[in]  factor      Multiplicative factor for st. deviation
 **
 *****************************************************************************/
double model_calcul_stdev(Model* model,
                          Db* db1,
                          int iech1,
                          Db* /*db2*/,
                          int iech2,
                          int verbose,
                          double factor)
{
  double c00, cov;
  CovCalcMode mode;

  /* Covariance at origin */

  int ndim = db1->getNDim();
  VectorDouble d1(ndim, 0.);
  model_calcul_cov(NULL,model, mode, 1, 1., d1, &c00);

  /* Covariance at increment */

  for (int idim = 0; idim < ndim; idim++)
    d1[idim] = db1->getDistance1D(iech1, iech2, idim);
  model_calcul_cov(NULL,model, mode, 1, 1., d1, &cov);
  double stdev = factor * sqrt(c00 - cov);

  if (verbose)
  {
    message("Db1(%d) - Db2(%d)", iech1 + 1, iech2 + 1);
    message(" - Incr=");
    for (int idim = 0; idim < ndim; idim++)
      message(" %lf", d1[idim]);
    message(" - c(0)=%lf cov=%lf stdev=%lf\n", c00, cov, stdev);
  }
  return (stdev);
}

/****************************************************************************/
/*!
 **  Establish the covariance matrix between two Dbs
 **  where samples are selected by ranks
 **  The output is stored in a Sparse Matrix
 **
 ** \return Array containing the covariance matrix
 **
 ** \param[in]  model  Model structure
 ** \param[in]  db1    First Db
 ** \param[in]  nsize1 Number of selected samples
 ** \param[in]  ranks1 Array giving ranks of selected samples (optional)
 ** \param[in]  db2    Second Db
 ** \param[in]  nsize2 Number of selected samples
 ** \param[in]  ranks2 Array giving ranks of selected samples (optional)
 ** \param[in]  ivar0  Rank of the first variable (-1: all variables)
 ** \param[in]  jvar0  Rank of the second variable (-1: all variables)
 ** \param[in]  flag_norm  1 if the model is normalized
 ** \param[in]  flag_cov   1 if the result must be given in covariance
 **
 ** \remarks The covariance matrix (returned) must be freed by calling routine
 ** \remarks The covariance matrix is established for the first variable
 ** \remarks and returned as a covariance
 ** \remarks As the ranks are used, no test is performed on any selection
 ** \remarks but only ranks positive or null are considered
 **
 *****************************************************************************/
cs* model_covmat_by_ranks_cs(Model *model,
                             Db *db1,
                             int nsize1,
                             const int *ranks1,
                             Db *db2,
                             int nsize2,
                             const int *ranks2,
                             int ivar0,
                             int jvar0,
                             int flag_norm,
                             int flag_cov)
{
  CovCalcMode mode;
  mode.update(ECalcMember::LHS, 0, 0, -1, flag_norm, flag_cov);
  if (st_check_model(model)) return nullptr;
  if (st_check_environ(model, db1)) return nullptr;
  if (st_check_environ(model, db2)) return nullptr;
  int ndim  = model->getDimensionNumber();
  int nvar  = model->getVariableNumber();
  int nvar1 = nvar;
  if (ivar0 >= 0)
  {
    nvar1 = 1;
    if (st_check_variable(nvar, ivar0)) return nullptr;
  }
  int nvar2 = nvar;
  if (jvar0 >= 0)
  {
    nvar2 = 1;
    if (st_check_variable(nvar, jvar0)) return nullptr;
  }
  if (ranks1 == nullptr) nsize1 = db1->getSampleNumber();
  if (ranks2 == nullptr) nsize2 = db1->getSampleNumber();

  /* Core allocation */

  VectorDouble d1(ndim, 0.);
  VectorDouble covtab(nvar * nvar, 0.);

  // Constitute the triplet

  cs* T = cs_spalloc(0, 0, 1, 1, 1);
  if (T == nullptr) return nullptr;

  /* Loop on the number of variables */

  for (int ivar = 0; ivar < nvar1; ivar++)
  {
    if (ivar0 >= 0) ivar = ivar0;

    /* Loop on the first sample */

    for (int i1 = 0; i1 < nsize1; i1++)
    {
      int iech1 = (ranks1 != nullptr) ? ranks1[i1] : i1;
      if (iech1 < 0) continue;

      /* Loop on the second variable */

      for (int jvar = 0; jvar < nvar2; jvar++)
      {
        if (jvar0 >= 0) jvar = jvar0;

        /* Loop on the second sample */

        for (int i2 = 0; i2 < nsize2; i2++)
        {
          int iech2 = (ranks2 != nullptr) ? ranks2[i2] : i2;
          if (iech2 < 0) continue;

          /* Determine the indices */

          int ecr1 = ivar * nsize1 + i1;
          int ecr2 = jvar * nsize2 + i2;

          // Save calculations due to matrix symmetry

          if (ecr2 > ecr1) continue;

          /* Loop on the dimension of the space */

          double value = TEST;
          int skip = 0;
          for (int i = 0; i < ndim && skip == 0; i++)
          {
            double v1 = db1->getCoordinate(iech1, i);
            double v2 = db2->getCoordinate(iech2, i);
            if (FFFF(v1) || FFFF(v2)) skip = 1;
            d1[i] = v1 - v2;
          }
          if (skip) continue;

          model_calcul_cov(NULL,model, mode, 1, 1., d1, covtab.data());
          value = COVTAB(ivar, jvar);
          if (ABS(value) < EPSILON10) continue;
          if (! cs_entry(T, ecr1, ecr2, value))
          {
            T = cs_spfree(T);
            return nullptr;
          }
          if (ecr1 != ecr2)
          {
            if (! cs_entry(T, ecr2, ecr1, value))
            {
              T = cs_spfree(T);
              return nullptr;
            }
          }
        }
      }
    }
  }

  // Convert from triplet to sparse matrix

  cs* covmat = cs_triplet(T);
  T = cs_spfree(T);
  return (covmat);
}

/****************************************************************************/
/*!
 **  Establish the covariance matrix between two Dbs
 **
 ** \param[in]  model Model structure
 ** \param[in]  db1   First Db
 ** \param[in]  db2   Second Db
 ** \param[in]  ivar0 Rank of the first variable (-1: all variables)
 ** \param[in]  jvar0 Rank of the second variable (-1: all variables)
 ** \param[in]  flag_norm  1 if the model is normalized
 ** \param[in]  flag_cov   1 if the result must be given in covariance
 **
 ** \param[out] covmat The covariance matrix
 **                    (Dimension = (nactive * nvar) [squared])
 **                    nactive: Number of samples active
 **                    nvar   : Number of selected variables (1 or nvar)
 **
 ** \remarks: If db2 is not provided, it is set to db1
 **
 *****************************************************************************/
int model_covmat(Model *model,
                 Db *db1,
                 Db *db2,
                 int ivar0,
                 int jvar0,
                 int flag_norm,
                 int flag_cov,
                 double *covmat)
{
  CovCalcMode mode;
  mode.update(ECalcMember::LHS, 0, 0, -1, flag_norm, flag_cov);
  if (db2 == nullptr) db2 = db1;
  if (st_check_model(model)) return 1;
  if (st_check_environ(model, db1)) return 1;
  if (st_check_environ(model, db2)) return 1;
  int ndim = model->getDimensionNumber();
  int nvar = model->getVariableNumber();
  int nech1 = db1->getSampleNumber();
  int nech2 = db2->getSampleNumber();
  int nvar1 = nvar;
  if (ivar0 >= 0)
  {
    nvar1 = 1;
    if (st_check_variable(nvar, ivar0)) return 1;
  }
  int nvar2 = nvar;
  if (jvar0 >= 0)
  {
    nvar2 = 1;
    if (st_check_variable(nvar, jvar0)) return 1;
  }

  /* Core allocation */

  VectorDouble d1(ndim, 0);
  VectorDouble covtab(nvar * nvar,0.);

  /* Loop on the first variable */

  int ecr = 0;
  for (int ivar = 0; ivar < nvar1; ivar++)
  {
    if (ivar0 >= 0) ivar = ivar0;

    /* Loop on the first sample */

    for (int iech1 = 0; iech1 < nech1; iech1++)
    {
      if (!db1->isActive(iech1)) continue;

      /* Loop on the second variable */

      for (int jvar = 0; jvar < nvar2; jvar++)
      {
        if (jvar0 >= 0) jvar = jvar0;

        /* Loop on the second sample */

        for (int iech2 = 0; iech2 < nech2; iech2++)
        {
          if (!db2->isActive(iech2)) continue;

          /* Loop on the dimension of the space */

          double value = TEST;
          int skip = 0;
          for (int i = 0; i < ndim && skip == 0; i++)
          {
            double v1 = db1->getCoordinate(iech1, i);
            double v2 = db2->getCoordinate(iech2, i);
            if (FFFF(v1) || FFFF(v2)) skip = 1;
            d1[i] = v1 - v2;
          }
          if (!skip)
          {
            if (model->isNoStat())
            {
              CovInternal covint(1, iech1, 2, iech2, ndim, db1, db2);
              model_calcul_cov(&covint, model, mode, 1, 1., d1, covtab.data());
            }
            else
            {
              model_calcul_cov(NULL,model, mode, 1, 1., d1, covtab.data());
            }
            value = COVTAB(ivar, jvar);
          }
          covmat[ecr++] = value;
        }
      }
    }
  }
  return 0;
}

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
#include "geoslib_enum.h"
#include "geoslib_old_f.h"
#include "geoslib_f.h"
#include "Drifts/DriftFactory.hpp"
#include "Drifts/EDrift.hpp"
#include "Drifts/ADrift.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/File.hpp"
#include "Covariances/ACovAnisoList.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovGradientNumerical.hpp"
#include "Model/CovInternal.hpp"
#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"
#include "Model/ModTrans.hpp"
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
  if (flag_init) for (int ivar = 0; ivar < nvar; ivar++)
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

  if (norme != 0.) for (ivar = 0; ivar < nvar; ivar++)
    for (jvar = 0; jvar < nvar; jvar++)
      COVTAB(ivar,jvar)/= norme;
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
  //  const CovAniso* cova = model->getCova(icov);
  // TODO: Why having to use ACov rather than CovAniso?
  const ACov *cova = dynamic_cast<const ACov*>(model->getCova(icov));

  if (member != ECalcMember::LHS && model->isCovaFiltered(icov))
    return (0.);
  else
    return cova->eval(0, 0, 1., d1);
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
 ** \param[in]  weight       Weight attached to this calculation
 **
 ** \param[out] d1          Working array (dimension = ndim) or NULL
 ** \param[out] covtab      output covariance (dimension = nvar * nvar)
 **
 *****************************************************************************/
void model_calcul_cov_direct(CovInternal *covint,
                             Model *model,
                             const CovCalcMode &mode,
                             int flag_init,
                             double weight,
                             VectorDouble d1,
                             double *covtab)
{
  // Load the non-stationary parameters if needed

  if (model->isNoStat()) model_nostat_update(covint, model);

  // Evaluate the Model

  MatrixSquareGeneral mat;
  if (d1.empty())
    mat = model->getCovAnisoList()->ACov::eval0(mode);
  else
    mat = model->getCovAnisoList()->ACov::eval(1., d1, VectorDouble(), mode);

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
 **  This is the function for convolution case
 **
 ** \param[in]  cov_nostat   Internal structure (should be NULL)
 ** \param[in]  model        Model structure
 ** \param[in]  mode         CovCalcMode structure
 ** \param[in]  flag_init    Initialize the array beforehand
 ** \param[in]  weight       Weight attached to this calculation
 **
 ** \param[out] d1          Working array (dimension = ndim) or NULL
 ** \param[out] covtab      output covariance (dimension = nvar * nvar)
 **
 *****************************************************************************/
static void st_model_calcul_cov_convolution(CovInternal *cov_nostat,
                                            Model *model,
                                            const CovCalcMode &mode,
                                            int flag_init,
                                            double weight,
                                            VectorDouble d1,
                                            double *covtab)
{

  /* This function is not programmed yet */

  messerr("The convolution covariance calculation function");
  messerr("is not programmed yet");

  return;
}

/*****************************************************************************/
/*!
 **  Returns the covariances for an increment
 **  This is the function for Hermitian case
 **
 ** \param[in]  cov_nostat   Internal structure (should be NULL)
 ** \param[in]  model        Model structure
 ** \param[in]  mode         CovCalcMode structure
 ** \param[in]  flag_init    Initialize the array beforehand
 ** \param[in]  weight       Weight attached to this calcAddress for the next term after the driftulation
 **
 ** \param[out] d1          Working array (dimension = ndim) or NULL
 ** \param[out] covtab      output covariance (dimension = nvar * nvar)
 **
 ** \remark This function is only programmed in the monovariate case
 **
 *****************************************************************************/
static void st_model_calcul_cov_anam_hermitian(CovInternal *cov_nostat,
                                               Model *model,
                                               const CovCalcMode &mode,
                                               int flag_init,
                                               double weight,
                                               VectorDouble d1,
                                               double *covtab)
{
  double rho, covfin, dist2, coeff, psin2, rn, rhon;
  int iclass;

  /* Initializations */

  const ModTrans &modtrs = model->getModTrans();
  AnamHermite *anam_hermite = dynamic_cast<AnamHermite*>(modtrs.getAnam());

  /* Check if the distance is zero */

  dist2 = rho = covfin = coeff = 0.;
  if (!d1.empty()) dist2 = matrix_norm(d1.data(), model->getDimensionNumber());
  if (dist2 > 0)
  {

    /* Calculate the generic variogram value */

    model_calcul_cov_direct(NULL, model, mode, 1, 1, d1, &rho);
  }

  /* Update the covariance */

  if (modtrs.getAnamIClass() == 0)
  {

    /*********************************************/
    /* Structure for the whole discretized grade */
    /*********************************************/

    if (dist2 <= 0.)
    {
      rn = 1.;
      for (iclass = 1; iclass < modtrs.getAnamNClass(); iclass++)
      {
        // TODO verifier qu'il faut bien mettre la moyenne par classe ici
        psin2 = modtrs.getAnamMeans(iclass) * modtrs.getAnamMeans(iclass);
        rn *= anam_hermite->getRCoef();
        switch (mode.getMember().toEnum())
        {
          case ECalcMember::E_LHS:
            coeff = 1. / (rn * rn);
            break;
          case ECalcMember::E_RHS:
            coeff = 1. / rn;
            break;
          case ECalcMember::E_VAR:
            coeff = 1.;
            break;
        }
        covfin += coeff * psin2;
      }
    }
    else
    {
      rn = 1.;
      rhon = 1.;
      for (iclass = 1; iclass < modtrs.getAnamNClass(); iclass++)
      {
        // TODO verifier qu'il faut bien mettre la moyenne par classe ici
        psin2 = modtrs.getAnamMeans(iclass) * modtrs.getAnamMeans(iclass);
        rn *= anam_hermite->getRCoef();
        rhon *= rho;
        switch (mode.getMember().toEnum())
        {
          case ECalcMember::E_LHS:
            coeff = 1.;
            break;

          case ECalcMember::E_RHS:
            coeff = 1. / rn;
            break;

          case ECalcMember::E_VAR:
            coeff = 1.;
            break;
        }
        covfin += coeff * psin2 * rhon;
      }
    }
  }
  else
  {

    /**************************************************/
    /* Structure for the factor 'modtrs.anam_iclass' */
    /**************************************************/

    if (dist2 <= 0.)
    {
      rn = pow(anam_hermite->getRCoef(), (double) modtrs.getAnamIClass());
      switch (mode.getMember().toEnum())
      {
        case ECalcMember::E_LHS:
          coeff = 1.;
          break;

        case ECalcMember::E_RHS:
          coeff = rn;
          break;

        case ECalcMember::E_VAR:
          coeff = 1.;
          break;
      }
      covfin = coeff;
    }
    else
    {
      rn = pow(anam_hermite->getRCoef(), (double) modtrs.getAnamIClass());
      rhon = pow(rho, (double) modtrs.getAnamIClass());
      switch (mode.getMember().toEnum())
      {
        case ECalcMember::E_LHS:
          coeff = rn * rn;
          break;
        case ECalcMember::E_RHS:
          coeff = rn;
          break;
        case ECalcMember::E_VAR:
          coeff = 1.;
          break;
      }
      covfin = coeff * rhon;
    }
  }

  if (flag_init)
    covtab[0] = covfin * weight;
  else
    covtab[0] += covfin * weight;

  return;
}

/*****************************************************************************/
/*!
 **  Returns the covariances for an increment
 **  This is the function for Discrete Diffusion case
 **
 ** \param[in]  cov_nostat   Internal structure (should be NULL)
 ** \param[in]  model        Model structure
 ** \param[in]  mode         CovCalcMode structure
 ** \param[in]  flag_init    Initialize the array beforehand
 ** \param[in]  weight       Weight attached to this calculation
 **
 ** \param[out] d1          Working array (dimension = ndim) or NULL
 ** \param[out] covtab      output covariance (dimension = nvar * nvar)
 **
 ** \remark This function is only programmed in the monovariate case
 **
 *****************************************************************************/
static void st_model_calcul_cov_anam_DD(CovInternal *cov_nostat,
                                        Model *model,
                                        const CovCalcMode &mode,
                                        int flag_init,
                                        double weight,
                                        VectorDouble d1,
                                        double *covtab)
{
  double gamref, covfin, csi, li, mui, dist2, coeff;
  int iclass, ndim;

  /* Initializations */

  const ModTrans &modtrs = model->getModTrans();
  AnamDiscreteDD *anam_discrete_DD = dynamic_cast<AnamDiscreteDD*>(modtrs.getAnam());
  ndim = model->getDimensionNumber();

  /* Check if the distance is zero */

  dist2 = gamref = covfin = coeff = 0.;
  if (!d1.empty()) dist2 = matrix_norm(d1.data(), model->getDimensionNumber());
  if (dist2 > 0)
  {
    VectorDouble covint(ndim * ndim);

    /* Calculate the generic variogram value */

    model_calcul_cov_direct(NULL, model, mode, 1, 1., VectorDouble(),
                            covint.data());
    model_calcul_cov_direct(NULL, model, mode, 0, -1, d1, covint.data());
    gamref = covint[0];
  }

  /* Update the covariance */

  if (modtrs.getAnamIClass() == 0)
  {

    /*********************************************/
    /* Structure for the whole discretized grade */
    /*********************************************/

    if (dist2 <= 0.)
      for (iclass = 1; iclass < modtrs.getAnamNClass(); iclass++)
      {
        csi = anam_discrete_DD->getDDStatCnorm(iclass);
        mui = anam_discrete_DD->getDDStatMul(iclass);
        switch (mode.getMember().toEnum())
        {
          case ECalcMember::E_LHS:
            coeff = csi * csi;
            break;
          case ECalcMember::E_RHS:
            coeff = csi * csi / mui;
            break;
          case ECalcMember::E_VAR:
            coeff = csi * csi;
            break;
        }
        covfin += coeff;
      }
    else
      for (iclass = 1; iclass < modtrs.getAnamNClass(); iclass++)
      {
        li = anam_discrete_DD->getDDStatLambda(iclass);
        csi = anam_discrete_DD->getDDStatCnorm(iclass);
        mui = anam_discrete_DD->getDDStatMul(iclass);
        switch (mode.getMember().toEnum())
        {
          case ECalcMember::E_LHS:
            coeff = csi * csi;
            break;

          case ECalcMember::E_RHS:
            coeff = csi * csi / mui;
            break;

          case ECalcMember::E_VAR:
            coeff = csi * csi;
            break;
        }
        covfin += coeff * exp(-li * gamref);
      }
  }
  else
  {

    /**************************************************/
    /* Structure for the factor 'modtrs.anam_iclass' */
    /**************************************************/

    if (dist2 <= 0.)
    {
      mui = anam_discrete_DD->getDDStatMul(modtrs.getAnamIClass());
      switch (mode.getMember().toEnum())
      {
        case ECalcMember::E_LHS:
          coeff = 1.;
          break;

        case ECalcMember::E_RHS:
          coeff = mui;
          break;

        case ECalcMember::E_VAR:
          coeff = 1.;
          break;
      }
      covfin = coeff;
    }
    else
    {
      mui = anam_discrete_DD->getDDStatMul(modtrs.getAnamIClass());
      li = anam_discrete_DD->getDDStatLambda(modtrs.getAnamIClass());
      switch (mode.getMember().toEnum())
      {
        case ECalcMember::E_LHS:
          coeff = mui * mui;
          break;
        case ECalcMember::E_RHS:
          coeff = mui;
          break;
        case ECalcMember::E_VAR:
          coeff = 1.;
          break;
      }
      covfin = coeff * exp(-li * gamref);
    }
  }

  if (flag_init)
    covtab[0] = covfin * weight;
  else
    covtab[0] += covfin * weight;

  return;
}

/*****************************************************************************/
/*!
 **  Calculate the residual covariance
 **
 ** \param[in]  model  : Model structure
 ** \param[in]  icut0  : Rank of the reference class
 ** \param[in]  covwrk : Working array
 **
 *****************************************************************************/
static double st_cov_residual(Model *model,
                              int icut0,
                              const VectorDouble &covwrk)
{
  ModTrans &modtrs = model->getModTrans();
  AnamDiscreteIR *anam_discrete_IR = dynamic_cast<AnamDiscreteIR*>(modtrs.getAnam());

  /* Get the pointer of the basic structure for the current model */

  int icov = 0;
  for (int icut = 0; icut < icut0; icut++)
    icov += (int) modtrs.getAnamStrCount()[icut];

  /* Loop of the covariance basic structures */

  double covloc = 0.;
  int number = (int) modtrs.getAnamStrCount()[icut0];
  for (int i = 0; i < number; i++, icov++)
    covloc += covwrk[icov] * model->getSill(icov, 0, 0);
  covloc *= anam_discrete_IR->getIRStatR(icut0 + 1);

  return (covloc);
}

/*****************************************************************************/
/*!
 **  Calculate the sum of the residual covariances
 **
 ** \param[in]  model  : Model structure
 ** \param[in]  icut0  : Rank of the reference class
 ** \param[in]  covwrk : Working array
 **
 *****************************************************************************/
static double st_covsum_residual(Model *model,
                                 int icut0,
                                 const VectorDouble &covwrk)
{
  double covsum = 0.;
  for (int icut = 0; icut <= icut0; icut++)
    covsum += st_cov_residual(model, icut, covwrk);
  return (1. + covsum);
}

/*****************************************************************************/
/*!
 **  Returns the covariances for an increment
 **  This is the function for Discrete Indicator Residuals case
 **
 ** \param[in]  cov_nostat   Internal structure (should be NULL)
 ** \param[in]  model        Model structure
 ** \param[in]  mode         CovCalcMode structure
 ** \param[in]  flag_init    Initialize the array beforehand
 ** \param[in]  weight       Weight attached to this calculation
 **
 ** \param[out] d1          Working array (dimension = ndim) or NULL
 ** \param[out] covtab      output covariance (dimension = nvar * nvar)
 **
 ** \remark This function is only programmed in the monovariate case
 **
 *****************************************************************************/
static void st_model_calcul_cov_anam_IR(CovInternal *cov_nostat,
                                        Model *model,
                                        const CovCalcMode &mode,
                                        int flag_init,
                                        double weight,
                                        VectorDouble d1,
                                        double *covtab)
{
  double covfin, bi, r, dist2, cov1, cov2, dcov;
  int icut, nclass, icut0, ncut, ndim, ncova;

  /* Initializations */

  ModTrans &modtrs = model->getModTrans();
  AnamDiscreteIR *anam_discrete_IR = dynamic_cast<AnamDiscreteIR*>(modtrs.getAnam());
  nclass = modtrs.getAnamNClass();
  ncut = nclass - 1;
  ncova = model->getCovaNumber();
  ;
  ndim = model->getDimensionNumber();
  r = anam_discrete_IR->getRCoef();
  VectorDouble covwrk(ncova, 0.);
  VectorDouble covint(ndim * ndim, 0.);

  /* Calculate the generic variogram value */

  model_calcul_cov_direct(NULL, model, mode, flag_init, 1., d1, covint.data());

  /* Modification of the covariance */

  covfin = 0.;
  icut0 = modtrs.getAnamIClass() - 1;
  if (icut0 < 0)
  {

    /* Structure for the whole discretized grade */

    cov2 = 1.;
    for (icut = 0; icut < ncut; icut++)
    {
      cov1 = cov2;
      bi = anam_discrete_IR->getIRStatB(icut);
      cov2 = pow(st_covsum_residual(model, icut, covwrk), r);
      dcov = cov2 - cov1;
      covfin += bi * bi * dcov;
    }
  }
  else
  {
    /* Check if the distance is zero */

    dist2 = 0.;
    if (!d1.empty())
      dist2 = matrix_norm(d1.data(), model->getDimensionNumber());

    /* Structure for the factor 'modtrs.anam_iclass' */

    if (dist2 <= 0)
    {
      covfin = anam_discrete_IR->getIRStatR(icut0 + 1);
    }
    else
    {
      cov1 = pow(st_covsum_residual(model, icut0 - 1, covwrk), r);
      cov2 = pow(st_covsum_residual(model, icut0, covwrk), r);
      covfin = cov2 - cov1;
    }
  }

  if (flag_init)
    covtab[0] = covfin * weight;
  else
    covtab[0] += covfin * weight;

  return;
}

/*****************************************************************************/
/*!
 **  Returns the covariances for an increment
 **  This is the function for tapering case
 **
 ** \param[in]  cov_nostat   Internal structure (should be NULL)
 ** \param[in]  model        Model structure
 ** \param[in]  mode         CovCalcMode structure
 ** \param[in]  flag_init    Initialize the array beforehand
 ** \param[in]  weight       Weight attached to this calculation
 **
 ** \param[out] d1          Working array (dimension = ndim) or NULL
 ** \param[out] covtab      output covariance (dimension = nvar * nvar)
 **
 *****************************************************************************/
static void st_model_calcul_cov_tapering(CovInternal *cov_nostat,
                                         Model *model,
                                         const CovCalcMode &mode,
                                         int flag_init,
                                         double weight,
                                         VectorDouble d1,
                                         double *covtab)
{
  double h;
  int nvar = model->getVariableNumber();

  /* Calculate the generic covariance value */

  model_calcul_cov_direct(NULL, model, mode, flag_init, weight, d1, covtab);

  /* Calculate the tapering effect */

  ModTrans &modtrs = model->getModTrans();
  h = 0.;
  if (!d1.empty())
    h = sqrt(matrix_norm(d1.data(), model->getDimensionNumber())) / modtrs.getTape()->getRange();
  double cov_tape = modtrs.getTape()->evaluate(h);

  /* Update all covariance values */

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      COVTAB(ivar,jvar)*= cov_tape;

  return;
}

/*****************************************************************************/
/*!
 **  Returns the covariances for an increment
 **
 ** \param[in]  model       Structure containing the model
 ** \param[in]  mode        CovCalcMode structure
 ** \param[in]  flag_init   initialize the array beforehand
 ** \param[in]  weight      Weight attached to this calculation
 ** \param[in]  d1          vector of increment (or NULL)
 **
 ** \param[out] covtab      output covariance (dimension = nvar * nvar)
 **
 *****************************************************************************/
void model_calcul_cov(Model *model,
                      CovCalcMode &mode,
                      int flag_init,
                      double weight,
                      const VectorDouble &d1,
                      double *covtab)
{
  /* Modify the member in case of properties */

  if (model->getModTransMode() != EModelProperty::NONE)
  {
    int anam_var = model->getModTrans().getAnamPointBlock();
    /* 'anam_var' is negative if model evaluation is called from dk() */
    /* This modification is performed in model_anamorphosis_set_factor() */
    // TODO : this must be checked! anam_vario always equal to LHS, RHS or VAR ???
    if (anam_var >= 0) mode.setMember(ECalcMember::fromValue(anam_var));
  }

  /* Call the generic model calculation module */

  model->generic_cov_function(NULL, model, mode, flag_init, weight, d1, covtab);

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

  if (model->getModTransMode() != EModelProperty::NONE)
  my_throw("Model transformation is not programmed for this entry");

  /* Call the generic model calculation module */

  // TODO Correct this which has something to do with pure virtual eval although implemented in Acov
  // compared to eval0 which is not implemented with such arguments.
  double value;
  if (d1.empty())
    value = model->getCovAnisoList()->eval0(ivar, jvar, mode);
  else
    value = model->getCovAnisoList()->ACov::eval(ivar, jvar, 1., d1,
                                                 VectorDouble(), mode);

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

/*****************************************************************************/
/*!
 **  Returns the covariances for an increment
 **  The increment is calculated between two samples of two Dbs
 **  This function is available in the non-stationary case
 **
 ** \param[in]  model       Structure containing the model
 ** \param[in]  mode        CovCalcMode structure
 ** \param[in]  covint      Pointer to the Internal Covariance structure
 ** \param[in]  flag_init   initialize the array beforehand
 ** \param[in]  weight      Weight attached to this calculation
 **
 ** \param[out] d1          Working array (dimension = ndim) or NULL
 ** \param[out] covtab      output covariance (dimension = nvar * nvar)
 **
 *****************************************************************************/
void model_calcul_cov_nostat(Model *model,
                             CovCalcMode &mode,
                             CovInternal *covint,
                             int flag_init,
                             double weight,
                             VectorDouble &d1,
                             double *covtab)
{
  /* Modify the member in case of properties */

  if (model->getModTransMode() != EModelProperty::NONE)
  {
    my_throw("ModTrans not yet implemented");
    int anam_var = model->getModTrans().getAnamPointBlock();

    /* 'anam_var' is negative if model evaluation is called from df() */
    /* This modification is performed in model_anamorphosis_set_factor() */
    // TODO : this must be checked! anam_vario always equal to LHS, RHS or VAR ???
    if (anam_var >= 0) mode.setMember(ECalcMember::fromValue(anam_var));
  }

  /* Call the generic model calculation module */

  model->generic_cov_function(covint, model, mode, flag_init, weight, d1,
                              covtab);

  return;
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
  for (int il = 0; il < model->getDriftNumber(); il++)
    drftab[il] = model->evaluateDrift(db, iech, il, member);
  return;
}

/****************************************************************************/
/*!
 **  Calculate the C0 terms
 **
 ** \param[in]  model     Model structure
 ** \param[in]  koption   Koption structure
 ** \param[in]  covtab    array of cumulated covariances
 **
 ** \param[out]  var0     array for C0[] (Dimension = nvar * nvar)
 **
 *****************************************************************************/
void model_variance0(Model *model,
                     Koption *koption,
                     double *covtab,
                     double *var0)
{
  int i, j, ecr, ivar, jvar, idim, nscale;
  CovCalcMode mode;

  /* Initializations */

  nscale = 1;
  int nvar = model->getVariableNumber();
  int ndim = model->getDimensionNumber();
  VectorDouble d1(ndim, 0.);
  mode.setMember(ECalcMember::VAR);

  switch (koption->calcul.toEnum())
  {
    case EKrigOpt::E_PONCTUAL:
      nscale = 1;
      model_calcul_cov(model, mode, 1, 1., d1, covtab);
      break;

    case EKrigOpt::E_BLOCK:
      nscale = koption->ntot;
      model_covtab_init(1, model, covtab);
      for (i = 0; i < nscale; i++)
        for (j = 0; j < nscale; j++)
        {
          for (idim = 0; idim < model->getDimensionNumber(); idim++)
            d1[idim] = DISC1(i,idim) - DISC2(j, idim);
          model_calcul_cov(model, mode, 0, 1., d1, covtab);
        }
      nscale = nscale * nscale;
      break;

    case EKrigOpt::E_DRIFT:
      nscale = 1;
      model_covtab_init(1, model, covtab);
      break;
  }

  /* Normalization */

  for (ivar = 0; ivar < nvar; ivar++)
    for (jvar = 0; jvar < nvar; jvar++)
      COVTAB(ivar,jvar)/= (double) nscale;

      /* Returning arguments */

  for (ivar = ecr = 0; ivar < nvar; ivar++)
    for (jvar = 0; jvar < nvar; jvar++, ecr++)
      var0[ecr] = COVTAB(ivar, jvar);

  return;
}

/****************************************************************************/
/*!
 **  Calculate the C0 terms
 **  This function is available in the non-stationary case
 **
 ** \param[in]  model     Model structure
 ** \param[in]  koption   Koption structure
 ** \param[in]  covint    Covariance Internal structure
 ** \param[in]  covtab    array of cumulated covariances
 **
 ** \param[out]  var0     array for C0[] (Dimension = nvar * nvar)
 **
 *****************************************************************************/
void model_variance0_nostat(Model *model,
                            Koption *koption,
                            CovInternal *covint,
                            double *covtab,
                            double *var0)
{
  int i, j, ecr, ivar, jvar, idim, nscale;
  CovCalcMode mode;

  /* Initializations */

  nscale = 1;
  int nvar = model->getVariableNumber();
  int ndim = model->getDimensionNumber();
  VectorDouble d1(ndim, 0.);

  mode.setMember(ECalcMember::VAR);
  switch (koption->calcul.toEnum())
  {
    case EKrigOpt::E_PONCTUAL:
      nscale = 1;
      model_calcul_cov_nostat(model, mode, covint, 1, 1., d1, covtab);
      break;

    case EKrigOpt::E_BLOCK:
      nscale = koption->ntot;
      model_covtab_init(1, model, covtab);
      for (i = 0; i < nscale; i++)
        for (j = 0; j < nscale; j++)
        {
          for (idim = 0; idim < model->getDimensionNumber(); idim++)
            d1[idim] = DISC1(i,idim) - DISC2(j, idim);
          model_calcul_cov_nostat(model, mode, covint, 0, 1., d1, covtab);
        }
      nscale = nscale * nscale;
      break;

    case EKrigOpt::E_DRIFT:
      nscale = 1;
      model_covtab_init(1, model, covtab);
      break;
  }

  for (ivar = ecr = 0; ivar < nvar; ivar++)
    for (jvar = 0; jvar < nvar; jvar++, ecr++)
      var0[ecr] = COVTAB(ivar,jvar)/ (double) nscale;

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
 **  Allocate the Model structure
 **
 ** \return  Pointer to the Model structure newly allocated
 **
 ** \param[in]  ndim          Space dimension
 ** \param[in]  nvar          Number of variables
 ** \param[in]  field         Field extension
 ** \param[in]  flag_linked   1 if the variables are linked; 0 otherwise
 ** \param[in]  flag_gradient 1 if the model is used for Gradient calculation
 ** \param[in]  ball_radius   Radius for Gradient calculation
 ** \param[in]  mean          Array for the mean (only used for KS)
 **                           (dimension: nvar)
 ** \param[in]  covar0        Array of variance-covariance at origin
 **                           (only used for covariance or covariogram)
 **                           (dimension: nvar*nvar)
 **
 *****************************************************************************/
Model* model_init(int ndim,
                  int nvar,
                  double field,
                  int flag_linked,
                  double ball_radius,
                  bool flag_gradient,
                  const VectorDouble &mean,
                  const VectorDouble &covar0)
{
  Model *model = nullptr;

  /// TODO : Force SpaceRN creation (modèle poubelle)
  CovContext ctxt = CovContext(nvar, ndim);
  ctxt.setField(field);
  ctxt.setBallRadius(ball_radius);
  if (static_cast<int>(mean.size()) > 0) ctxt.setMean(mean);
  if (static_cast<int>(covar0.size()) > 0) ctxt.setCovar0(covar0);

  model = new Model(ctxt, flag_gradient, flag_linked);

  /* Set the error return flag */

  model_setup(model);

  return (model);
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
  Model *model;
  double sill;
  int error;

  /* Initializations */

  error = 1;

  /* Create the empty Model */

  model = model_init(ndim, nvar, 1.);
  if (model == nullptr) return model;

  /* Add the nugget effect variogram model */

  sill = 1.;
  if (model_add_cova(model, ECov::NUGGET, 0, 0, 0., 0., VectorDouble(),
                     VectorDouble(), VectorDouble(1, sill))) goto label_end;

  /* Set the error return flag */

  model_setup(model);
  error = 0;

  label_end: if (error) model = model_free(model);
  return (model);
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
                   const VectorDouble &sill)
{
  if (st_check_model(model)) return 1;

  // Add a new element in the structure

  if (model->isFlagGradient())
  {
    CovGradientNumerical covgrad(type, model->getContext());
//    if (! covgrad.isGradientCompatible()) return 1; TODO incorporte this type of selection
    covgrad.setParam(param);
    if (flag_aniso)
    {
      covgrad.setRanges(aniso_ranges);
      if (flag_rotation) covgrad.setAnisoRotation(aniso_rotmat);
    }
    else
      covgrad.setRange(range);

    if (static_cast<int>(sill.size()) > 0) covgrad.setSill(sill);
    model->addCova(&covgrad);
  }
  else
  {
    CovAniso cova(type, model->getContext());
    cova.setParam(param);
    if (flag_aniso)
    {
      cova.setRanges(aniso_ranges);
      if (flag_rotation) cova.setAnisoRotation(aniso_rotmat);
    }
    else
      cova.setRange(range);

    if (static_cast<int>(sill.size()) > 0) cova.setSill(sill);
    model->addCova(&cova);
  }

  /* Set the error return code */

  model_setup(model);

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
  int error;

  /* Initializations */

  error = 1;
  if (st_check_model(model)) goto label_end;

  // Allocate the new element

  drift = DriftFactory::createDriftFunc(type, model->getContext());
  drift->setRankFex(rank_fex);
  model->addDrift(drift);

  /* Set the error return code */

  model_setup(model);
  error = 0;

  label_end: return (error);
}

/****************************************************************************/
/*!
 **  Cancel any property
 **
 ** \return  Error return code
 **
 ** \param[in]  model      Pointer to the Model structure
 **
 *****************************************************************************/
int model_add_no_property(Model *model)

{
  int error;

  /* Initializations */

  error = 1;

  ModTrans &modtrs = model->getModTrans();
  if (st_check_model(model)) goto label_end;

  // Cancel the property

  modtrs.cancelProperty();

  /* Set the calling function */

  model_setup(model);

  /* Set the error return code */

  error = 0;

  label_end: return (error);
}

/****************************************************************************/
/*!
 **  Add the Model convolution parameters
 **
 ** \return  Error return code
 **
 ** \param[in]  model      Pointer to the Model structure
 ** \param[in]  conv_type  Type of the convolution function (starting from 1)
 ** \param[in]  conv_idir  Orientation of the convolution (starting from 1)
 ** \param[in]  conv_ndisc Number of discretization points per direction
 ** \param[in]  conv_range Convolution parameter
 **
 *****************************************************************************/
int model_add_convolution(Model *model,
                          int conv_type,
                          int conv_idir,
                          int conv_ndisc,
                          double conv_range)
{
  int error;

  /* Initializations */

  error = 1;
  ModTrans &modtrs = model->getModTrans();
  if (st_check_model(model)) goto label_end;

  /* Load the Convolution parameters */

  if (modtrs.addConvolution(conv_type, conv_idir, conv_ndisc, conv_range))
    goto label_end;

  /* Set the calling function */

  model_setup(model);

  /* Set the error return code */

  error = 0;

  label_end: return (error);
}

/****************************************************************************/
/*!
 **  Define the ranks of the factors of interest
 **
 ** \return  Error return code
 **
 ** \param[in]  model       Pointer to the Model structure
 ** \param[in]  anam_iclass Rank of the target factor (starting from 1)
 **                         0 for the whole discretized grade variable
 **
 ** \remark  This function overwrites some items of properties such as:
 ** \remark  - the rank of anam_iclass (passed as argument)
 ** \remark  - the anam_var element is set to -1 in order to avoid choosing
 ** \remark    the covariance mode as a function of anam_var rather than member
 **
 *****************************************************************************/
int model_anamorphosis_set_factor(Model *model, int anam_iclass)
{
  int error;

  /* Initializations */

  error = 1;
  ModTrans &modtrs = model->getModTrans();
  if (st_check_model(model)) goto label_end;

  /* Preliminary checks */

  if (!(anam_iclass == 0 || anam_iclass < modtrs.getAnamNClass()))
  {
    messerr("The rank of the active factor (%d) is incorrect", anam_iclass);
    messerr("It should lie between 1 and the number of factors (%d)",
            modtrs.getAnamNClass() - 1);
    messerr("or be set to 0 to estimate the whole discretized grade");
    goto label_end;
  }

  modtrs.setAnamIClass(anam_iclass);
  modtrs.setAnamVar(-1);

  /* Set the error return code */

  error = 0;

  label_end: return (error);
}

/****************************************************************************/
/*!
 **  Add the Model parameters accounting for Anamorphosis transform
 **
 ** \return  Error return code
 **
 ** \param[in]  model       Pointer to the Model structure
 ** \param[in]  anam_type   Type of the anamorphosis (EAnam)
 ** \param[in]  anam_nclass Number of classes
 ** \param[in]  anam_iclass Rank of the target factor (starting from 1)
 **                         0 for the whole discretized grade variable
 ** \param[in]  anam_var    Type of calculation
 **                         1 Covariance or variogram of punctual variable
 **                         2 Cross-Covariance or cross-variogram between
 **                           punctual and block variables
 **                         3 Covariance or variogram of block variable
 ** \param[in]  anam_coefr  Change of support coefficient (1 for point)
 ** \param[in]  anam_coefs  Information Effect coefficient (1 for point)
 ** \param[in]  anam_strcnt Array giving the number of structures per model
 **                         (used for Discrete Indicator Residuals case)
 ** \param[in]  anam_stats  Array of statistics
 **
 *****************************************************************************/
int model_add_anamorphosis(Model *model,
                           const EAnam &anam_type,
                           int anam_nclass,
                           int anam_iclass,
                           int anam_var,
                           double anam_coefr,
                           double anam_coefs,
                           VectorDouble &anam_strcnt,
                           VectorDouble &anam_stats)
{
  int error;

  /* Initializations */

  error = 1;
  ModTrans &modtrs = model->getModTrans();
  if (st_check_model(model)) goto label_end;

  /* Preliminary checks */

  if (model->getVariableNumber() > 1)
  {
    messerr("The Anamorphosis Mode is only programmed in Monovariate case");
    goto label_end;
  }

  /* Load the exponentiation parameters */

  if (modtrs.addAnamorphosis(anam_type, anam_nclass, anam_iclass, anam_var,
                             anam_coefr, anam_coefs, anam_strcnt, anam_stats))
    goto label_end;

  /* Set the calling function */

  model_setup(model);

  /* Set the error return code */

  error = 0;

  label_end: return (error);
}

/****************************************************************************/
/*!
 **  Add the Model Tapering parameters
 **
 ** \return  Error return code
 **
 ** \param[in]  model      Pointer to the Model structure
 ** \param[in]  tape_type  Type of the tapering function (starting from 1)
 ** \param[in]  tape_range Range of the tapering function
 **
 *****************************************************************************/
int model_add_tapering(Model *model, int tape_type, double tape_range)
{
  int error;

  /* Initializations */

  error = 1;
  ModTrans &modtrs = model->getModTrans();
  if (st_check_model(model)) goto label_end;

  /* Preliminary check */

  if (model->getDimensionNumber() > modtrs.getTape()->getMaxNDim())
  {
    messerr("The selected tapering function is not compatible with");
    messerr("the space dimension used for the Model (%d)",
            model->getDimensionNumber());
    goto label_end;
  }

  /* Load the tapering parameters */

  if (modtrs.addTapering(tape_type, tape_range)) goto label_end;

  /* Set the calling function */

  model_setup(model);

  /* Set the error return code */

  error = 0;

  label_end: return (error);
}

/****************************************************************************/
/*!
 **  For a given basic structure, get the reduction factor to convert the
 **  theoretical to practical scale
 **
 ** \return  Convertion factor
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

/****************************************************************************/
/*!
 **  Sets the functions for covariance and drift
 **
 ** \param[in,out]  model Model structure
 **
 *****************************************************************************/
void model_setup(Model *model)

{
  if (model == nullptr) return;

  /**********************************/
  /* Allocation of auxiliary arrays */
  /**********************************/

  /* Define the generic covariance function */

  switch (model->getModTransMode().toEnum())
  {
    case EModelProperty::E_NONE:
      model->generic_cov_function = model_calcul_cov_direct;
      break;

    case EModelProperty::E_CONV:
      model->generic_cov_function = st_model_calcul_cov_convolution;
      break;

    case EModelProperty::E_ANAM:
      switch (model->getModTrans().getAnam()->getType().toEnum())
      {
        case EAnam::E_HERMITIAN:
          model->generic_cov_function = st_model_calcul_cov_anam_hermitian;
          break;

        case EAnam::E_DISCRETE_DD:
          model->generic_cov_function = st_model_calcul_cov_anam_DD;
          break;

        case EAnam::E_DISCRETE_IR:
          model->generic_cov_function = st_model_calcul_cov_anam_IR;
          break;

        default:
          messerr("The Model modified by Properties is not available");
          messerr("For the following Anamorphosis type (%d)",
                  model->getModTrans().getAnam()->getType().getValue());
          break;
      }
      break;

    case EModelProperty::E_TAPE:
      model->generic_cov_function = st_model_calcul_cov_tapering;
      break;

    default:
      break;
  }
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
  int ivar, jvar, icov, ncova, nvar, error;

  /* Initializations */

  error = 1;
  ncova = model->getCovaNumber();
  nvar = model->getVariableNumber();

  /* Calculate the eigen values and vectors of the coregionalization matrix */

  for (icov = 0; icov < ncova; icov++)
  {
    if (!is_matrix_definite_positive(
        nvar, model->getCova(icov)->getSill().getValues().data(), valpro,
        vecpro, 0))
    {
      messerr("Warning: the model is not authorized");
      messerr(
          "The coregionalization matrix for the structure %d is not definite positive",
          icov + 1);
      goto label_end;
    }

    /* Calculate the factor matrix */

    for (ivar = 0; ivar < nvar; ivar++)
      for (jvar = 0; jvar < nvar; jvar++)
        AIC(icov,ivar,jvar)= VECPRO(ivar,jvar) * sqrt(VALPRO(jvar));
      }

      /* Set the error return code */

  error = 0;

  label_end: return (error);
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
  int error = 1;
  double *covtab = nullptr;
  CovCalcMode mode;
  mode.update(nugget_opt, nostd, member, rank_sel, flag_norm, flag_cov);
  if (norder > 0) mode.setOrderVario(norder);

  /* Preliminary checks */

  if (st_check_model(model)) return 1;
  int ndim = model->getDimensionNumber();
  int nvar = model->getVariableNumber();
  if (st_check_variable(nvar, ivar)) return 1;
  if (st_check_variable(nvar, jvar)) return 1;

  /* Core allocation */

  VectorDouble d1(ndim);
  covtab = (double*) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == nullptr) goto label_end;

  /* Normalize the direction vector codir */

  vario_fix_codir(ndim, codir);

  /* Loop on the lags */

  for (int ih = 0; ih < nh; ih++)
  {
    for (int idim = 0; idim < ndim; idim++)
      d1[idim] = h[ih] * codir[idim];
    model_calcul_cov(model, mode, 1, 1., d1, covtab);
    g[ih] = COVTAB(ivar, jvar);
  }

  /* Set the error return code */

  error = 0;

  label_end: covtab = (double*) mem_free((char* ) covtab);
  return (error);
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
  double *covtab, c00, var0;
  VectorDouble d1;

  /* Initializations */

  int error = 1;
  covtab = nullptr;
  CovCalcMode mode;
  mode.update(nugget_opt, nostd, member, rank_sel, flag_norm, flag_cov);
  if (norder > 0) mode.setOrderVario(norder);

  /* Preliminary checks */

  if (st_check_model(model)) return 1;
  int ndim = model->getDimensionNumber();
  int nvar = model->getVariableNumber();
  if (st_check_variable(nvar, ivar)) return 1;
  if (st_check_variable(nvar, jvar)) return 1;
  CovInternal covint(1, iech1, 2, iech2, ndim, db1, db2);

  /* Core allocation */

  covtab = (double*) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == nullptr) goto label_end;

  /* Normalize the direction vector codir */

  vario_fix_codir(ndim, codir);

  /* Calculate the C(0) term (used only for covariance or covariogram) */

  c00 = model->getContext().getCovar0(ivar, jvar);
  d1.resize(ndim, 0.);

  model_calcul_cov_nostat(model, mode, &covint, 1, 1., d1, covtab);
  var0 = COVTAB(ivar, jvar);
  if (c00 <= 0. || FFFF(c00)) c00 = var0;

  /* Loop on the lags */

  for (int ih = 0; ih < nh; ih++)
  {
    for (int idim = 0; idim < ndim; idim++)
      d1[idim] = h[ih] * codir[idim];
    model_calcul_cov_nostat(model, mode, &covint, 1, 1., d1, covtab);
    g[ih] = COVTAB(ivar, jvar);
  }

  /* Set the error return code */

  error = 0;

  label_end: covtab = (double*) mem_free((char* ) covtab);
  return (error);
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
  double *covtab;
  int iech, nvar, ndim, error;
  VectorDouble d1;

  /* Initializations */

  error = 1;
  covtab = nullptr;
  CovCalcMode mode;
  mode.update(0, 0, ECalcMember::LHS, -1, flag_norm, flag_cov);

  /* Preliminary checks */

  if (st_check_model(model)) goto label_end;
  if (st_check_environ(model, db)) goto label_end;
  ndim = model->getDimensionNumber();
  nvar = model->getVariableNumber();

  /* Core allocation */

  d1.resize(ndim);
  covtab = (double*) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == nullptr) goto label_end;

  /* Initialization */

  for (iech = 0; iech < db->getSampleNumber(); iech++)
    g[iech] = TEST;

  /* Loop on the lags */

  for (iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    db_sample_load(db, ELoc::X, iech, d1.data());
    model_calcul_cov(model, mode, 1, 1., d1, covtab);
    g[iech] = COVTAB(ivar, jvar);
  }

  /* Set the error return code */

  error = 0;

  label_end: covtab = (double*) mem_free((char* ) covtab);
  return (error);
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
  int il, nfex;

  /* Initializations */

  nfex = 0;
  if (model->getDriftNumber() <= 0) return (nfex);

  /* Loop on the drift functions */

  for (il = 0; il < model->getDriftNumber(); il++)
  {
    if (model->getDrift(il)->getType() == EDrift::F) nfex++;
  }
  return (nfex);
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
  if (rank < 0 || rank >= model->getDriftNumber()) return;
  model->setDriftFiltered(rank, filter);
  return;
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
  double *covtab, cxx, v1, v2, w1, w2, norme;
  int ndim, nvar, iech1, iech2, i, skip;
  VectorDouble d1;
  CovCalcMode mode;

  /* Initializations */

  cxx = TEST;
  covtab = nullptr;
  if (st_check_model(model)) goto label_end;
  if (st_check_environ(model, db1)) goto label_end;
  if (st_check_environ(model, db2)) goto label_end;
  ndim = model->getDimensionNumber();
  nvar = model->getVariableNumber();
  if (st_check_variable(nvar, ivar)) goto label_end;
  if (st_check_variable(nvar, jvar)) goto label_end;
  if (seed != 0) law_set_random_seed(seed);

  /* Core allocation */

  d1.resize(ndim, 0.);
  covtab = (double*) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == nullptr) goto label_end;
  model_covtab_init(1, model, covtab);

  /* Loop on the first sample */

  norme = 0.;
  for (iech1 = 0; iech1 < db1->getSampleNumber(); iech1++)
  {
    if (!db1->isActive(iech1)) continue;
    w1 = db1->getWeight(iech1);
    if (w1 == 0.) continue;

    /* Loop on the second sample */

    for (iech2 = 0; iech2 < db2->getSampleNumber(); iech2++)
    {
      if (!db2->isActive(iech2)) continue;
      w2 = db2->getWeight(iech2);
      if (w2 == 0.) continue;

      /* Loop on the dimension of the space */

      for (i = skip = 0; i < ndim && skip == 0; i++)
      {
        v1 = db1->getCoordinate(iech1, i);
        v2 = db2->getCoordinate(iech2, i);
        if (eps != 0.) v2 += eps * law_uniform(-0.5, 0.5);
        if (FFFF(v1) || FFFF(v2)) skip = 1;
        d1[i] = v1 - v2;
      }
      if (skip) continue;

      model_calcul_cov(model, mode, 0, w1 * w2, d1, covtab);
      norme += w1 * w2;
    }
  }

  /* Scaling */

  st_covtab_rescale(nvar, norme, covtab);
  cxx = COVTAB(ivar, jvar);

  /* Free memory */

  label_end: covtab = (double*) mem_free((char* ) covtab);
  return (cxx);
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
  double *covmat, *covtab, v1, v2, value;
  int ndim, nvar, nvar1, nvar2, iech1, iech2, i, skip, ecr, error, i1, i2;
  VectorDouble d1;

  /* Initializations */

  error = 1;
  covtab = covmat = nullptr;
  CovCalcMode mode;
  mode.update(0, 0, ECalcMember::LHS, -1, flag_norm, flag_cov);
  if (st_check_model(model)) goto label_end;
  if (st_check_environ(model, db1)) goto label_end;
  if (st_check_environ(model, db2)) goto label_end;
  ndim = model->getDimensionNumber();
  nvar = nvar1 = nvar2 = model->getVariableNumber();
  if (ivar0 >= 0)
  {
    nvar1 = 1;
    if (st_check_variable(nvar, ivar0)) goto label_end;
  }
  if (jvar0 >= 0)
  {
    nvar2 = 1;
    if (st_check_variable(nvar, jvar0)) goto label_end;
  }
  if (ranks1 == nullptr) nvar1 = db1->getSampleNumber();
  if (ranks2 == nullptr) nvar2 = db1->getSampleNumber();

  /* Core allocation */

  d1.resize(ndim, 0.);
  covtab = (double*) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == nullptr) goto label_end;
  covmat = (double*) mem_alloc(sizeof(double) * nsize1 * nsize2, 0);
  if (covmat == nullptr) goto label_end;

  /* Loop on the number of variables */

  ecr = 0;
  for (int ivar = 0; ivar < nvar1; ivar++)
  {
    if (ivar0 >= 0) ivar = ivar0;

    /* Loop on the first sample */

    for (i1 = 0; i1 < nsize1; i1++)
    {
      iech1 = (ranks1 != nullptr) ? ranks1[i1] :
                                    i1;
      if (iech1 < 0) continue;

      /* Loop on the second variable */

      for (int jvar = 0; jvar < nvar2; jvar++)
      {
        if (jvar0 >= 0) jvar = jvar0;

        /* Loop on the second sample */

        for (i2 = 0; i2 < nsize2; i2++)
        {
          iech2 = (ranks2 != nullptr) ? ranks2[i2] :
                                        i2;
          if (iech2 < 0) continue;

          /* Loop on the dimension of the space */

          value = TEST;
          for (i = skip = 0; i < ndim && skip == 0; i++)
          {
            v1 = db1->getCoordinate(iech1, i);
            v2 = db2->getCoordinate(iech2, i);
            if (FFFF(v1) || FFFF(v2)) skip = 1;
            d1[i] = v1 - v2;
          }
          if (!skip)
          {
            model_calcul_cov(model, mode, 1, 1., d1, covtab);
            value = COVTAB(ivar, jvar);
          }
          covmat[ecr++] = value;
        }
      }
    }
  }

  /* Set the error returned code */

  error = 0;

  /* Free memory */

  label_end: covtab = (double*) mem_free((char* ) covtab);
  if (error) covmat = (double*) mem_free((char* ) covmat);
  return (covmat);
}

/****************************************************************************/
/*!
 **  Establish the drift rectangular matrix for a given Db
 **
 ** \param[in]  model  Model structure
 ** \param[in]  member Member of the Kriging System (ECalcMember)
 ** \param[in]  db     Db structure
 **
 ** \param[out] drfmat The drift matrix
 **                    (Dimension = nech * nvar * nfeq * nvar)
 **
 *****************************************************************************/
void model_drift_mat(Model *model,
                     const ECalcMember &member,
                     Db *db,
                     double *drfmat)
{
  int nech, nvar, nbfl, nfeq, ecr, jb;
  double *drftab, value;

  /* Initializations */

  drftab = nullptr;
  if (st_check_model(model)) goto label_end;
  if (st_check_environ(model, db)) goto label_end;
  nvar = model->getVariableNumber();
  nbfl = model->getDriftNumber();
  nfeq = model->getDriftEquationNumber();
  nech = db->getSampleNumber();

  /* Core allocation */

  drftab = (double*) mem_alloc(sizeof(double) * nbfl, 0);
  if (drftab == nullptr) goto label_end;

  ecr = 0;

  /* Loop on the variables */

  for (int ivar = 0; ivar < nvar; ivar++)
  {

    /* Loop on the samples */

    for (int iech = 0; iech < nech; iech++)
    {
      if (!db->isActive(iech)) continue;
      model_calcul_drift(model, member, db, iech, drftab);

      /* Loop on the drift functions */

      if (model->isFlagLinked())
      {
        for (int ib = 0; ib < nfeq; ib++)
        {
          value = 0.;
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
            jb = jvar + nvar * jl;
            value = 0.;
            for (int il = 0; il < nbfl; il++)
              value += drftab[il] * model->getCoefDrift(ivar, il, jb);
            drfmat[ecr++] = value;
          }
      }
    }
  }

  label_end: drftab = (double*) mem_free((char* ) drftab);
  return;
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
void model_drift_vector(Model *model,
                        const ECalcMember &member,
                        Db *db,
                        int iech,
                        double *vector)
{
  int nvar, nbfl, nfeq, ivar, ib, il, ecr, i;
  double *drftab, value;

  /* Initializations */

  drftab = nullptr;
  if (st_check_model(model)) goto label_end;
  if (st_check_environ(model, db)) goto label_end;
  nvar = model->getVariableNumber();
  nbfl = model->getDriftNumber();
  nfeq = model->getDriftEquationNumber();

  /* Core allocation */

  drftab = (double*) mem_alloc(sizeof(double) * nbfl, 0);
  if (drftab == nullptr) goto label_end;

  /* Initialize the covariance matrix */

  for (i = 0; i < nvar * nfeq; i++)
    vector[i] = TEST;

  model_calcul_drift(model, member, db, iech, drftab);

  ecr = 0;
  for (ivar = 0; ivar < nvar; ivar++)
    for (ib = 0; ib < nfeq; ib++)
    {
      value = 0.;
      for (il = 0; il < nbfl; il++)
        value += drftab[il] * model->getCoefDrift(ivar, il, ib);
      vector[ecr++] = value;
    }

  label_end: drftab = (double*) mem_free((char* ) drftab);
  return;
}

/****************************************************************************/
/*!
 **  Returns the matrix of covariance at the origin
 **
 ** \param[in]  model Model structure
 ** \param[in]  mode  CovCalcMode structure
 **
 ** \param[out] covtab Working array (Dimension: nvar * nvar)
 ** \param[out] d1     Working array (Dimension: ndim)
 ** \param[out] c00tab The covariance matrix (Dimension = nvar * nvar)
 **
 *****************************************************************************/
static void st_matrix_c00(Model *model,
                          CovCalcMode &mode,
                          double *covtab,
                          VectorDouble d1,
                          double *c00tab)
{
  int ivar1, ivar2, nvar;
  double c00, var0;

  nvar = model->getVariableNumber();
  for (ivar1 = 0; ivar1 < nvar; ivar1++)
    for (ivar2 = 0; ivar2 < nvar; ivar2++)
    {
      c00 = model->getContext().getCovar0(ivar1, ivar2);
      model_calcul_cov(model, mode, 1, 1., VectorDouble(), covtab);
      var0 = COVTAB(ivar1, ivar2);
      if (c00 <= 0. || FFFF(c00)) c00 = var0;
      C00TAB(ivar1,ivar2)= c00;
    }
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
  int i, il;

  /* Look for the drift function */

  for (i = 0, il = -1; i < model->getDriftNumber() && il < 0; i++)
    if (model->getDriftType(i) == type && model->getDrift(i)->getRankFex()
        == rank) il = i;
  if (il < 0) messageAbort("st_drift_modify");

  /* Patch the drift coefficient */

  model->setCoefDrift(iv, il, ib, model->getCoefDrift(iv, il, ib) + value);

  return;
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
  int il, ib, rank;
  double value;

  for (ib = 0; ib < model->getDriftEquationNumber(); ib++)
    for (il = 0; il < model->getDriftNumber(); il++)
    {
      value = model->getCoefDrift(0, il, ib);
      if (value == 0) continue;

      EDrift type = model->getDriftType(il);
      rank = model->getRankFext(il);
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
  return;
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
  int flag_linked, new_nvar, nfact;
  double sill;
  bool flag_gradient;

  // Preliminary checks

  new_model = nullptr;
  int nvar = model->getVariableNumber();
  int ndim = model->getDimensionNumber();
  int ncova = model->getCovaNumber();
  int nbfl = model->getDriftNumber();
  flag_linked = nfact = new_nvar = 0;
  flag_gradient = false;

  // Create the new model (linked drift functions)

  switch (mode)
  {
    case -1:                    // Data (SK)
    case 0:                     // Data
      new_nvar = nvar;
      nfact = 1;
      flag_linked = 0;
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
      flag_linked = 1;
      flag_gradient = true;
      break;
  }
  new_model = model_init(ndim, new_nvar, model->getField(), flag_linked,
                         ball_radius, flag_gradient,
                         model->getContext().getMean(),
                         model->getContext().getCovar0());
  if (new_model == nullptr) return new_model;

  // ****************************************
  // Create the basic covariance structures
  // ****************************************

  int lec = 0;
  for (int icov = 0; icov < ncova; icov++)
  {
    cova = model->getCova(icov);
    sill = model->getSill(icov, 0, 0);
    for (int ifact = 0; ifact < nfact; ifact++, lec++)
    {
      if (model_add_cova(new_model, cova->getType(), cova->getFlagAniso(),
                         cova->getFlagRotation(), cova->getRange(),
                         cova->getParam(), cova->getRanges(),
                         cova->getAnisoRotMatVec(), VectorDouble()))
        return new_model;

      /* Modify the Covariance calculation type */;

      switch (mode)
      {
        case 0:
        case -1:
          for (int ivar = 0; ivar < new_nvar; ivar++)
            for (int jvar = 0; jvar < new_nvar; jvar++)
              new_model->setSill(lec, ivar, jvar,
                                 model->getSill(icov, ivar, jvar));
          break;

        case 1:                   // Data - Gradient
          for (int ivar = 0; ivar < new_nvar; ivar++)
            for (int jvar = 0; jvar < new_nvar; jvar++)
              new_model->setSill(lec, ivar, jvar, 0.);
          if (ifact == 0)
          {
            new_model->setSill(lec, 0, 0, sill);
          }
          else if (ifact == 1)
          {
            new_model->setSill(lec, 0, 1, -sill);
            new_model->setSill(lec, 1, 0, sill);
          }
          else if (ifact == 2)
          {
            new_model->setSill(lec, 1, 1, sill);
          }
          else if (ifact == 3)
          {
            new_model->setSill(lec, 0, 2, -sill);
            new_model->setSill(lec, 2, 0, sill);
          }
          else if (ifact == 4)
          {
            new_model->setSill(lec, 1, 2, -sill);
            new_model->setSill(lec, 2, 1, -sill);
          }
          else if (ifact == 5)
          {
            new_model->setSill(lec, 2, 2, sill);
          }
          else
          {
            my_throw("Argument 'ifact' invalid");
          }
          break;
      }
    }
  }

  // *********************************
  // Create the basic drift structures
  // *********************************

  if (mode >= 0)
  {
    for (int il = 0; il < nbfl; il++)
    {
      drft = model->getDrift(il);
      ADriftElem *newdrft = DriftFactory::createDriftFunc(
          drft->getType(), new_model->getContext());
      newdrft->setRankFex(drft->getRankFex());
      new_model->addDrift(newdrft);
      new_model->setDriftFiltered(il, model->isDriftFiltered(il));
    }

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

  // Set the error return code

  model_setup(new_model);

  return (new_model);
}

///****************************************************************************/
///*!
// **  Modifies a monovariate Model into a multivariate Model
// **
// ** \return  The modified Model structure
// **
// ** \param[in]  model       Input Model
// ** \param[in]  new_nvar    New number of variables
// ** \param[in]  mean        Array for means (optional)
// ** \param[in]  vars        Array for variances (optional)
// ** \param[in]  corr        Array for correlations (optional)
// **
// *****************************************************************************/
//Model *model_modify(Model  *model,
//                                int     new_nvar,
//                                double *mean,
//                                double *vars,
//                                double *corr)
//{
/// TODO [Cova] : to be restored ?
//  Model  *new_model;
//  int     ivar,jvar,nvar,icov,ncova,il,nbfl,error,ndim;
//  double  sill;
//  Cova   *cova,*cova_new;
//  Drift   *drft;
//  VectorDouble ranges;
//
//  /* Initializations */
//
//  new_model = nullptr;
//  error  = 1;
//  nvar   = model->getNVar();
//  ndim   = model->getNDim();
//  ncova  = model->getNCova();
//  nbfl   = model->getNDrift();
//  ranges.resize(ndim);
//
//  /* Preliminary checks */
//
//  if (nvar != 1)
//  {
//    messerr("This procedure is limited to a monovariate input model");
//    goto label_end;
//  }
//  if (new_nvar <= 1)
//  {
//    messerr("This procedure must only be used when new_nvar(%d) is larger than 1",new_nvar);
//    goto label_end;
//  }
//  if (vars != nullptr && ! is_matrix_non_negative(1,new_nvar,vars,0))
//  {
//    messerr("You provided vars[]. It must be non negative");
//    goto label_end;
//  }
//  if (vars != nullptr && vars[0] == 0.)
//  {
//    messerr("You provided vars[]. It must have vars[0] != 0");
//    goto label_end;
//  }
//  if (corr != nullptr && ! is_matrix_correlation(new_nvar,corr))
//  {
//    messerr("You provided corr[]. It must be a correlation matrix");
//    goto label_end;
//  }
//
//  /* Create the new model */
//
//  new_model = model_init(model->getNDim(),new_nvar,model->getFlagLinked(),
//                         model->getField(),0.,VectorDouble(),VectorDouble());
//
//  /* Set the mean (if provided) */
//
//  for (ivar = 1; ivar < new_nvar; ivar++)
//    new_model->setMean(ivar, (mean != nullptr) ? mean[ivar] :
//                                                         model->getMean(0));
//
//  /* Set the variance-covariance at the origin (if provided) */
//
//  new_model->setCovar0(0,0,model->getCovar0(0,0));
//  for (ivar=0; ivar<new_nvar; ivar++)
//    for (jvar=0; jvar<new_nvar; jvar++)
//    {
//      if (vars != nullptr)
//        new_model->setCovar0(ivar,jvar,vars[AD(ivar,jvar)]);
//      else
//        new_model->setCovar0(ivar,jvar,(ivar == 0 && jvar == 0) ?
//          model->getCovar0(0,0) : 0.);
//    }
//
//  /******************************************/
//  /* Create the basic covariance structures */
//  /******************************************/
//
//  for (icov=0; icov<ncova; icov++)
//  {
//    cova = model->getCova(icov);
//    sill = model->getSill(icov,0,0);
//    model_convert_cova_to_ranges(ndim,cova,ranges);
//
//    /* Copy the basic structure */
//
//    if (model_add_cova(new_model,
//                       cova->getType(),cova->getFlagAniso(),
//                       cova->getFlagRotation(),cova->getRange(),
//                       cova->getParam(),
//                       ranges,cova->getAnisoRotMat(),
//                       VectorDouble())) goto label_end;
//    cova_new = new_model->getCova(icov);
//    cova_new->setFlagFilter(cova->getFlagFilter());
//
//    /* Scale the sills */
//
//    for (ivar=0; ivar<new_nvar; ivar++)
//      for (jvar=0; jvar<new_nvar; jvar++)
//      {
//        double value =
//            (vars == nullptr) ?
//            sill : (sill * sqrt(vars[ivar] * vars[jvar]) *
//                    corr[ivar * new_nvar + jvar] / vars[0]);
//        new_model->getCova(icov)->setSill(ivar,jvar,value);
//      }
//  }
//
//  /*************************************/
//  /* Create the basic drift structures */
//  /*************************************/
//
//  for (il=0; il<nbfl; il++)
//  {
//    drft = model->getDrift(il);
//    if (model_add_drift(new_model,
//                        drft->getDriftType(),
//                        drft->getRankFex())) goto label_end;
//    new_model->getDrift(il)->setFlagFilter(drft->getFlagFilter());
//  }
//
//  /* Set the error return code */
//
//  if (model_setup(new_model)) goto label_end;
//  error = 0;
//
//label_end:
//  if (error) new_model = model_free(new_model);
//
//  return(new_model);
//}
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
  double *total;
  int nvar, ncov, error, ivar, jvar, icov, flag_norm;

  /* Initializations */

  error = 1;
  nvar = model->getVariableNumber();
  ncov = model->getCovaNumber();
  total = nullptr;

  /* Core allocation */

  total = (double*) mem_alloc(sizeof(double) * nvar, 1);

  /* Calculate the total sills for each variable */

  for (ivar = flag_norm = 0; ivar < nvar; ivar++)
  {
    total[ivar] = model->getCovAnisoList()->getTotalSill(ivar, ivar);
    if (total[ivar] == 0.) goto label_end;
    total[ivar] = sqrt(total[ivar]);
    if (ABS(total[ivar] - 1.) > EPSILON6) flag_norm = 1;
  }

  /* Scale the different sills for the different variables */

  for (ivar = 0; ivar < nvar; ivar++)
    for (jvar = 0; jvar < nvar; jvar++)
      for (icov = 0; icov < ncov; icov++)
      {
        double sill = model->getCova(icov)->getSill(ivar, jvar);
        sill /= total[ivar] * total[jvar];
        model->getCova(icov)->setSill(ivar, jvar, sill);
      }

  /* Printout */

  if (flag_verbose && flag_norm)
  {
    message("The model has been normalized\n");
    for (ivar = 0; ivar < nvar; ivar++)
      message("- Variable %d : Scaling factor = %lf\n", ivar + 1,
              total[ivar] * total[ivar]);
  }

  /* Error return code */

  error = 0;

  label_end: total = (double*) mem_free((char* ) total);
  return (error);
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
  CovAniso *cova;
  double total;
  int nvar, ncov, error, icov;

  /* Initializations */

  error = 1;
  nvar = model->getVariableNumber();
  if (nvar > 1) return (0);
  if (percent <= 0.) return (0);
  ncov = model->getCovaNumber();

  /* Check if the model only contains GAUSSIAN components */

  total = 0.;
  for (icov = 0; icov < ncov; icov++)
  {
    cova = model->getCova(icov);
    if (cova->getType() != ECov::GAUSSIAN) return (0);
    total += model->getSill(icov, 0, 0);
  }
  total = total * percent / 100.;

  /* Update each Gaussian component */

  for (icov = 0; icov < ncov; icov++)
  {
    cova = model->getCova(icov);
    model->getCova(icov)->setSill(0, 0, 1. - total);
  }

  /* Add a NUGGET EFFECT component */

  if (model_add_cova(model, ECov::NUGGET, 0, 0, 0., 0., VectorDouble(),
                     VectorDouble(), VectorDouble(1, total))) goto label_end;

  /* Printout */

  if (flag_verbose)
  {
    message("The model which only contains Gaussian components\n");
    message("has been stabilized by adding a small Nugget Effect\n");
  }

  /* Error return code */

  error = 0;

  label_end: return (error);
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
        model->getCova(rank_nugget)->setSill(ivar, jvar, (ivar == jvar) ? diff :
                                                                          0.);
      else
        nugget[AD(ivar, jvar)] = (ivar == jvar) ? diff :
                                                  0.;
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
 ** \param[in]  verbose Verbose option
 ** \param[in]  model   Model structure
 ** \param[in]  db      Db structure
 ** \param[in]  iech    Rank of the sample
 ** \param[in]  ivar    Rank of the variable
 ** \param[in]  coef    Array of coefficients (optional)
 **
 ** \param[out] drftab  Working array
 **
 *****************************************************************************/
double model_drift_evaluate(int verbose,
                            Model *model,
                            const Db *db,
                            int iech,
                            int ivar,
                            double *coef,
                            double *drftab)
{
  double drift, value;
  int il, ib;

  /* Initializations */

  drift = TEST;
  if (st_check_environ(model, db)) goto label_end;

  model_calcul_drift(model, ECalcMember::LHS, db, iech, drftab);

  /* Check if all the drift terms are defined */

  for (il = 0; il < model->getDriftNumber(); il++)
    if (FFFF(drftab[il])) return (TEST);

  /* Perform the correction */

  drift = 0.;
  for (ib = 0; ib < model->getDriftEquationNumber(); ib++)
  {
    value = 0.;
    for (il = 0; il < model->getDriftNumber(); il++)
      value += drftab[il] * model->getCoefDrift(ivar, il, ib);
    drift += value * coef[ib];
  }
  label_end: return (drift);
}

/****************************************************************************/
/*!
 **  Ask the characteristics of the Model structure
 **
 ** \return  Pointer to the newly allocated model
 ** \return  or NULL if a problem occured
 **
 ** \param[in]  ndim      Space dimension
 ** \param[in]  nvar      Number of variables
 ** \param[in]  order     Order of the IRF
 ** \param[in]  flag_sill 1 if the sill must be defined; 0 otherwise
 ** \param[in]  flag_norm 1 if the total sill must be normalize to 1
 **                       (only valid in the monovariate case)
 ** \param[in]  model_in  Input Model structure
 **
 *****************************************************************************/
Model* input_model(int ndim,
                   int nvar,
                   int order,
                   int flag_sill,
                   int flag_norm,
                   Model *model_in)
{
  /// TODO [Cova] : to be restored ?
//  int    i,flag_def,error,ncova;
//  Model *model;
//  Cova  *cova,*cova_in;
//
//  /* Initializations */
//
//  error    = 1;
//  flag_def = (model_in != nullptr);
//  model    = nullptr;
//
//  /* Core allocation */
//
//  model = model_init(ndim,nvar,0,0.,0.,VectorDouble(),VectorDouble());
//  if (model == nullptr) goto label_end;
//
//  /* Number of Basic structures */
//
//  ncova = _lire_int("Count of Basic structures",flag_def,
//                    (flag_def) ? model_in->getNCova() : ITEST,
//                    1,ITEST);
//
//  /* Loop on the basic structures */
//
//  for (i=0; i<ncova; i++)
//  {
//    // Add the new covariance to the Model
//
//    cova = model->addCova();
//
//    /* Ask for the parameters of the basic structure */
//
//    cova_in = nullptr;
//    if (flag_def && i < model_in->getNCova())
//      cova_in = model_in->getCova(i);
//
//    // Define the covariance interactively
//
//    cova->input(order,flag_sill,cova_in);
//  }
//
//  /* Normalization */
//
//  if (flag_norm && nvar == 1)
//  {
//    if (model_normalize(model,1)) goto label_end;
//  }
//
//  /* Set the error returned code */
//
//  error = 0;
//
//label_end:
//  if (error) model = model_free(model);
//  return(model);
  return nullptr;
}

/****************************************************************************/
/*!
 **  Returns the number of basic structures in the Model
 **
 ** \return  Number of basic structures
 **
 ** \param[in]  model     Model structure
 **
 *****************************************************************************/
int model_dimension(Model *model)
{
  return (model->getCovaNumber());
}

/****************************************************************************/
/*!
 **  Ask the characteristics of one Covariance structure
 **
 ** \return  Error returned code
 **
 ** \param[in]  model     Model structure
 ** \param[in]  icov      Rank of the Covariance structure (from 0)
 **
 ** \param[out]  cov_type      Type of the covariance (enum of ECov)
 ** \param[out]  flag_aniso    1 for anisotropy and 0 otherwise
 ** \param[out]  param         Parameter
 ** \param[out]  sill          Array of sills (Dimension = nvar * nvar)
 ** \param[out]  aniso_rotmat  Rotation matrix (Dimension = ndim * ndim)
 ** \param[out]  aniso_ranges  Rotation ranges (Dimension = ndim)
 **
 *****************************************************************************/
int model_extract_cova(Model *model,
                       int icov,
                       ECov *cov_type,
                       int *flag_aniso,
                       double *param,
                       VectorDouble &sill,
                       VectorDouble &aniso_rotmat,
                       VectorDouble &aniso_ranges)
{
  CovAniso *cova;
  int ndim;

  /* Initializations */

  if (icov < 0 || icov >= model->getCovaNumber()) return (1);
  cova = model->getCova(icov);
  ndim = model->getDimensionNumber();

  /* Returning arguments */

  *cov_type = cova->getType();
  *flag_aniso = !cova->isIsotrop();
  *param = cova->getParam();
  sill = cova->getSill().getValues();

  if (!cova->getAnisoRotMatVec().empty())
    aniso_rotmat = cova->getAnisoRotMatVec();
  else
  {
    aniso_rotmat.resize(ndim * ndim);
    ut_rotation_init(ndim, aniso_rotmat.data());
  }
  aniso_ranges = cova->getRanges();

  return (0);
}

/****************************************************************************/
/*!
 **  Ask the characteristics of one Property
 **
 ** \param[in]  model     Model structure
 **
 ** \param[out]  tape_range    Range of the tapering
 **
 ** \remark This function is only a template: at the current stage, it only
 ** \remark returns the values that may have been updated by a model fitting
 ** \remark procedure.
 **
 *****************************************************************************/
void model_extract_properties(Model *model, double *tape_range)
{
  ModTrans &modtrs = model->getModTrans();

  *tape_range = modtrs.getTape()->getRange();
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
 ** \param[out] flag_aniso     1 if anisotropy is meaningfull
 ** \param[out] flag_rotation  1 if an anisotropy rotation is meaningfull
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
  SpaceRN space(1); // Retrieve all covariances
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
  double *covtab;
  int i, idir, ndir, ipas, npas, idim, ndim, error, nvar, ivar, jvar, ijvar;
  VectorDouble d1;

  error = 1;
  ndim = vario->getDimensionNumber();
  ndir = vario->getDirectionNumber();
  nvar = model->getVariableNumber();
  covtab = nullptr;
  CovCalcMode mode;
  mode.update(0, 0, ECalcMember::LHS, -1, flag_norm, flag_cov);

  /* Core allocation */

  d1.resize(ndim);
  covtab = (double*) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == nullptr) goto label_end;
  vario->setNVar(nvar);

  // Internal redimensioning

  vario->internalVariableResize();
  vario->internalDirectionResize();

  /* Calculate the C(0) constant term */

  model_calcul_cov(model, mode, 1, 1., VectorDouble(), covtab);
  for (i = 0; i < nvar * nvar; i++)
    vario->setVars(i, covtab[i]);

  /* Loop on the directions */

  for (idir = 0; idir < ndir; idir++)
  {
    npas = vario->getLagNumber(idir);

    /* Loop on the variogram lags */

    for (ipas = 0; ipas < npas; ipas++)
    {

      /* Loop on the variables */

      for (ivar = ijvar = 0; ivar < vario->getVariableNumber(); ivar++)
        for (jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {
          i = vario->getDirAddress(idir, ivar, jvar, ipas, false, 0);
          vario->setSw(idir, i, 1.);
          vario->setHh(idir, i, ipas * vario->getDPas(idir));
          for (idim = 0; idim < ndim; idim++)
            d1[idim] = vario->getHh(idir, i) * vario->getCodir(idir, idim);
          model_calcul_cov(model, mode, 1, 1., d1, covtab);
          vario->setGg(idir, i, COVTAB(ivar, jvar));
        }
    }
  }

  /* Set the error returned code */

  error = 0;

  label_end: covtab = (double*) mem_free((char* ) covtab);
  return (error);
}

/****************************************************************************/
/*!
 **  Establish the covariance multivariate vector between a given sample
 **  and a given variable and all the samples and variables of a Db
 **
 ** \param[in]  model      Model structure
 ** \param[in]  db         Db structure
 ** \param[in]  ivar       Rank of the target variable
 ** \param[in]  iech       Rank of the target sample
 ** \param[in]  flag_norm  1 if the model is normalized
 ** \param[in]  flag_cov   1 if the result must be given in covariance
 **
 ** \param[out] vector     Returned vector
 **                    (Dimension = neq) where neq = nactive * nvar
 **
 *****************************************************************************/
void model_vector_multivar(Model *model,
                           Db *db,
                           int ivar,
                           int iech,
                           int flag_norm,
                           int flag_cov,
                           double *vector)
{
  double *covtab, *c00tab;
  int ndim, nvar, jech, i, skip, nech, ecr, jvar;
  VectorDouble d1;

  /* Initializations */

  covtab = c00tab = nullptr;
  CovCalcMode mode;
  mode.update(0, 0, ECalcMember::LHS, -1, flag_norm, flag_cov);
  if (st_check_model(model)) goto label_end;
  if (st_check_environ(model, db)) goto label_end;
  ndim = model->getDimensionNumber();
  nvar = model->getVariableNumber();
  nech = db->getSampleNumber();
  if (st_check_variable(nvar, ivar)) goto label_end;

  /* Core allocation */

  d1.resize(ndim, 0.);
  covtab = (double*) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == nullptr) goto label_end;
  c00tab = (double*) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (c00tab == nullptr) goto label_end;

  /* Calculate the C(0) term */

  st_matrix_c00(model, mode, covtab, d1, c00tab);

  /* Loop on the sample */

  ecr = 0;
  for (jvar = 0; jvar < nvar; jvar++)
    for (jech = 0; jech < nech; jech++)
    {
      if (!db->isActive(jech)) continue;

      /* Loop on the dimension of the space */

      for (i = skip = 0; i < ndim && skip == 0; i++)
      {
        d1[i] = db->getDistance1D(iech, jech, i);
        if (FFFF(d1[i])) skip = 1;
      }
      if (skip) continue;

      model_calcul_cov(model, mode, 1, 1., d1, covtab);
      vector[ecr++] = COVTAB(ivar, jvar);
    }

  /* Free memory */

  label_end: covtab = (double*) mem_free((char* ) covtab);
  c00tab = (double*) mem_free((char* ) c00tab);
  return;
}

/****************************************************************************/
/*!
 **  Establish the covariance vector between a given sample
 **  and all the samples of a Db
 **
 ** \param[in]  model      Model structure
 ** \param[in]  db1        Data Db
 ** \param[in]  db2        Target Db
 ** \param[in]  ivar       Rank of the first variable
 ** \param[in]  jvar       Rank of the second variable
 ** \param[in]  jech       Rank of the particular sample (in db2)
 ** \param[in]  flag_norm  1 if the model is normalized
 ** \param[in]  flag_cov   1 if the result must be given in covariance
 **
 ** \param[out] vector     Returned vector
 **
 *****************************************************************************/
void model_vector(Model *model,
                  Db *db1,
                  Db *db2,
                  int ivar,
                  int jvar,
                  int jech,
                  int flag_norm,
                  int flag_cov,
                  double *vector)
{
  double *covtab, v1, v2;
  int ndim, nvar, iech, i, skip, nech;
  VectorDouble d1;

  /* Initializations */

  covtab = nullptr;
  CovCalcMode mode;
  mode.update(0, 0, ECalcMember::LHS, -1, flag_norm, flag_cov);
  if (st_check_model(model)) goto label_end;
  if (st_check_environ(model, db1)) goto label_end;
  if (st_check_environ(model, db2)) goto label_end;
  ndim = model->getDimensionNumber();
  nvar = model->getVariableNumber();
  nech = db1->getSampleNumber();
  if (st_check_variable(nvar, ivar)) goto label_end;
  if (st_check_variable(nvar, jvar)) goto label_end;

  /* Core allocation */

  d1.resize(ndim, 0.);
  covtab = (double*) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == nullptr) goto label_end;

  /* Initialize the covariance matrix */

  for (i = 0; i < nech; i++)
    vector[i] = TEST;

  /* Loop on the sample */

  for (iech = 0; iech < nech; iech++)
  {
    if (!db1->isActive(iech)) continue;

    /* Loop on the dimension of the space */

    for (i = skip = 0; i < ndim && skip == 0; i++)
    {
      v1 = db1->getCoordinate(iech, i);
      v2 = db2->getCoordinate(jech, i);
      if (FFFF(v1) || FFFF(v2)) skip = 1;
      d1[i] = v1 - v2;
    }
    if (skip) continue;

    model_calcul_cov(model, mode, 1, 1., d1, covtab);
    vector[iech] = COVTAB(ivar, jvar);
  }

  /* Free memory */

  label_end: covtab = (double*) mem_free((char* ) covtab);
  return;
}

/****************************************************************************/
/*!
 **  Establish the covariance vector between a given sample
 **  and all the samples of a Db (Non-stationary case)
 **
 ** \param[in]  model      Model structure
 ** \param[in]  db         Db structure
 ** \param[in]  ivar       Rank of the first variable
 ** \param[in]  jvar       Rank of the second variable
 ** \param[in]  iech       Rank of the particular sample
 **
 ** \param[out] vector     Returned vector
 **
 *****************************************************************************/
void model_vector_nostat(Model *model,
                         Db *db,
                         int ivar,
                         int jvar,
                         int iech,
                         double *vector)
{
  double *covtab, value;
  int ndim, nvar, jech, i, skip, nech;
  VectorDouble d1;
  CovCalcMode mode;

  /* Initializations */

  covtab = nullptr;
  if (st_check_model(model)) goto label_end;
  if (st_check_environ(model, db)) goto label_end;
  ndim = model->getDimensionNumber();
  nvar = model->getVariableNumber();
  nech = db->getSampleNumber();
  if (st_check_variable(nvar, ivar)) goto label_end;
  if (st_check_variable(nvar, jvar)) goto label_end;

  /* Core allocation */

  d1.resize(ndim, 0.);
  covtab = (double*) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == nullptr) goto label_end;

  /* Initialize the covariance matrix */

  for (i = 0; i < nech; i++)
    vector[i] = TEST;

  /* Loop on the sample */

  for (jech = 0; jech < nech; jech++)
  {
    if (!db->isActive(jech)) continue;

    /* Loop on the dimension of the space */

    for (i = skip = 0; i < ndim && skip == 0; i++)
    {
      d1[i] = db->getDistance1D(iech, jech, i);
      if (FFFF(d1[i])) skip = 1;
    }
    if (skip) continue;

    CovInternal covint(1, iech, 1, jech, ndim, db, db);
    model_calcul_cov_nostat(model, mode, &covint, 1, 1., d1, covtab);
    value = COVTAB(ivar, jvar);
    vector[jech] = value;
  }

  /* Free memory */

  label_end: covtab = (double*) mem_free((char* ) covtab);
  return;
}

/****************************************************************************/
/*!
 **  Calculate the maximum distance for a model
 **
 ** \return  Maximum distance
 **
 ** \param[in]  model      Model structure
 **
 *****************************************************************************/
double model_maximum_distance(Model *model)

{
  return model->getCovAnisoList()->getMaximumDistance();
}

/****************************************************************************/
/*!
 **  For a given basic structure, convert the theoretical range (scale) into
 **  the practical range (which is the one actually stored in gstlearn)
 **
 ** \return  Range
 **
 ** \param[in]  type      Type of the basic structure
 ** \param[in]  scale     Theoretical range
 ** \param[in]  param     Third parameter
 **
 *****************************************************************************/
double model_scale2range(const ECov &type, double scale, double param)
{
  double factor, range;

  factor = cova_get_scale_factor(type, param);
  range = scale * factor;
  return (range);
}

/****************************************************************************/
/*!
 **  For a given basic structure, convert the practical range into
 **  the theoretical range (scale)
 **
 ** \return  Scale
 **
 ** \param[in]  type      Type of the basic structure
 ** \param[in]  range     Practical range
 ** \param[in]  param     Third parameter
 **
 *****************************************************************************/
double model_range2scale(const ECov &type, double range, double param)
{
  double factor, scale;

  factor = cova_get_scale_factor(type, param);
  scale = range / factor;
  return (scale);
}

/****************************************************************************/
/*!
 **  Get the field parameter from a Model
 **
 ** \return  Field value
 **
 ** \param[in]  model     Model structure
 **
 *****************************************************************************/
double model_get_field(Model *model)
{
  return (model->getField());
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
 ** \remarks: It has been exptended to the case where only one model is defined
 **
 *****************************************************************************/
Model* model_combine(const Model *model1, const Model *model2, double r)
{
  Model *model;
  const CovAniso *cova;
  double field;
  int error, i, ncov;
  VectorDouble sill, mean, cova0;

  /* Initializations */

  error = 1;
  model = nullptr;
  sill.resize(4);
  mean.resize(2);
  cova0.resize(4);
  if (model1 == nullptr)
  {
    messerr("This function requires at least one model defined");
    return (model);
  }
  if (model1 != nullptr && model1->getVariableNumber() != 1)
  {
    messerr("This function can only combine monovariate models");
    return (model);
  }
  if (model2 != nullptr && model2->getVariableNumber() != 1)
  {
    messerr("This function can only combine monovariate models");
    return (model);
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
    return (model);
  }
  if (model1->isFlagLinked() || model2->isFlagLinked())
  {
    messerr("This function cannot combine models with linked drifts");
    return (model);
  }

  /* Create the output model */

  field = MAX(model1->getField(), model2->getField());
  mean[0] = model1->getContext().getMean(0);
  mean[1] = model2->getContext().getMean(0);
  cova0[0] = 1.;
  cova0[1] = r;
  cova0[2] = r;
  cova0[3] = 1.;
  double radius = MAX(model1->getContext().getBallRadius(),
                      model2->getContext().getBallRadius());
  model = model_init(model1->getDimensionNumber(), 2, field, 0, radius, false,
                     mean, cova0);
  if (model == nullptr) return model;

  ncov = 0;

  /* Add the covariance of the first Model */

  for (i = 0; i < model1->getCovaNumber(); i++)
  {
    cova = model1->getCova(i);
    sill[0] = cova->getSill(0, 0);
    sill[1] = sill[2] = r * cova->getSill(0, 0);
    sill[3] = r * r * cova->getSill(0, 0);
    if (model_add_cova(model, cova->getType(), cova->getFlagAniso(),
                       cova->getFlagRotation(), cova->getRange(),
                       cova->getParam(), cova->getRanges(),
                       cova->getAnisoRotMatVec(), sill)) goto label_end;
    ncov++;
  }

  /* Add the covariance of the second Model */

  for (i = 0; i < model2->getCovaNumber(); i++)
  {
    cova = model2->getCova(i);
    sill[0] = 0.;
    sill[1] = sill[2] = 0.;
    sill[3] = (1. - r * r) * cova->getSill(0, 0);
    if (model_add_cova(model, cova->getType(), cova->getFlagAniso(),
                       cova->getFlagRotation(), cova->getRange(),
                       cova->getParam(), cova->getRanges(),
                       cova->getAnisoRotMatVec(), sill)) goto label_end;
    ncov++;
  }

  /* Set the error return code */

  error = 0;

  label_end: if (error) model = model_free(model);
  return (model);
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
  CovAniso *cova;
  int nstruc, icov;

  /* Loop on the basic structures */

  nstruc = 0;
  for (icov = 0; icov < model->getCovaNumber(); icov++)
  {
    cova = model->getCova(icov);
    if (cova->getType() != ECov::NUGGET) nstruc++;
  }
  return (nstruc);
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
 ** \param[in]  opt_norm  Option for normalization
 ** \param[in]  nug_ratio Ratio of the nugget effect
 **
 *****************************************************************************/
int model_regularize(Model *model,
                     Vario *vario,
                     Db *db,
                     int opt_norm,
                     double nug_ratio)
{
  double *covtab, *c00tab, v1, v2, norme, dist;
  int idim, ndim, nvar, idir, nech, ipas, iech, jech, ivar, jvar, iad, error;
  VectorDouble dd;
  CovCalcMode mode;

  /* Initializations */

  error = 1;
  c00tab = covtab = nullptr;
  if (st_check_model(model)) goto label_end;
  if (st_check_environ(model, db)) goto label_end;
  ndim = model->getDimensionNumber();
  nvar = model->getVariableNumber();

  /* Preliminary checks */

  if (!is_grid(db))
  {
    messerr("This calculation facility is dedicated to grid architecture");
    goto label_end;
  }
  nech = db->getSampleNumber();
  norme = nech * nech;
  vario->setNVar(nvar);
  vario->internalVariableResize();
  vario->internalDirectionResize();

  /* Core allocation */

  dd.resize(ndim, 0.);
  c00tab = (double*) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (c00tab == nullptr) goto label_end;
  covtab = (double*) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == nullptr) goto label_end;

  /* Calculate the Cvv (for a zero-shift) */

  model_covtab_init(1, model, c00tab);
  for (iech = 0; iech < nech; iech++)
    for (jech = 0; jech < nech; jech++)
    {
      for (idim = 0; idim < ndim; idim++)
      {
        dd[idim] = db->getDistance1D(iech, jech, idim);
      }
      model_calcul_cov(model, mode, 0, 1, dd, c00tab);
    }
  st_covtab_rescale(nvar, norme, c00tab);

  /* Initialize the variance array */

  for (ivar = 0; ivar < nvar; ivar++)
    for (jvar = 0; jvar < nvar; jvar++)
      vario->setVars(ivar, jvar, C00TAB(ivar, jvar));

  /* Loop on the directions */

  for (idir = 0; idir < vario->getDirectionNumber(); idir++)
  {

    /* Loop on the number of lags */

    for (ipas = 0; ipas < vario->getLagNumber(idir); ipas++)
    {
      model_covtab_init(1, model, covtab);
      dist = ipas * vario->getDPas(idir);

      for (iech = 0; iech < nech; iech++)
        for (jech = 0; jech < nech; jech++)
        {
          for (idim = 0; idim < ndim; idim++)
          {
            v1 = db->getCoordinate(iech, idim);
            v2 = db->getCoordinate(jech, idim)
                + dist * vario->getCodir(idir, idim);
            dd[idim] = v1 - v2;
          }
          model_calcul_cov(model, mode, 0, 1, dd, covtab);
        }
      st_covtab_rescale(nvar, norme, covtab);

      for (ivar = 0; ivar < nvar; ivar++)
        for (jvar = 0; jvar <= ivar; jvar++)
        {
          iad = vario->getDirAddress(idir, ivar, jvar, ipas, false, 0);
          vario->setGg(idir, iad, C00TAB(ivar,jvar)- COVTAB(ivar,jvar));
          vario->setHh(idir, iad, dist);
          vario->setSw(idir, iad, 1);
        }
    }
  }

  /* Set the error return code */

  error = 0;

  label_end: c00tab = (double*) mem_free((char* ) c00tab);
  covtab = (double*) mem_free((char* ) covtab);
  return (error);
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
  model_calcul_cov(model, mode, 1, 1., VectorDouble(), &c00);
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
      model_calcul_cov(model, mode, 1, 1., d1, &covar2);
      diag[i] = 2. * (c00 - covar2);
    }
    else
    {
      diag[i] = c00;
    }
    residual += diag[i];
  }
  tol = (!FFFF(eta)) ? eta * residual :
                       0.;
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
        model_calcul_cov(model, mode, 1, 1., d1, &covar1);

        for (int idim = 0; idim < 3; idim++)
          d1[idim] = db->getCoordinate(pvec[npivot], idim) - center[idim];
        model_calcul_cov(model, mode, 1, 1., d1, &covar2);

        for (int idim = 0; idim < 3; idim++)
          d1[idim] = db->getCoordinate(pvec[i], idim) - center[idim];
        model_calcul_cov(model, mode, 1, 1., d1, &covar3);

        G(npivot,i) = covar1 - covar2 - covar3 + c00;
      }
      else
      {
        // Calculate the covariance column C(:, npivot)
        (void) distance_intra(db, pvec[i], pvec[npivot], d1.data());
        model_calcul_cov(model, mode, 1, 1., d1, &G(npivot, i));
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
        model_calcul_cov(model, mode, 1, 1., d1, &covar2);

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

/****************************************************************************/
/*!
 **  Calculate the maximum order of the GRF
 **
 ** \return  Maximum order
 **
 ** \param[in]  model      Model structure
 **
 *****************************************************************************/
int model_maximum_order(Model *model)

{
  int order, max_order;

  if (model == nullptr) return (-1);

  max_order = 0;
  for (int il = 0; il < model->getDriftNumber(); il++)
  {
    ADriftElem *drft = model->getDrift(il);
    order = drft->getOrderIRF();
    if (order > max_order) max_order = order;
  }
  return (max_order);
}

/****************************************************************************/
/*!
 **  Find if a given drift function is defined in the current model
 **
 ** \return  1 if the drift function is used; 0 otherwise
 **
 ** \param[in]  model      Model structure
 ** \param[in]  type0      Drift function to be found (EDrift)
 **
 *****************************************************************************/
int model_is_drift_defined(Model *model, const EDrift &type0)
{
  if (model == nullptr) return (0);
  for (int il = 0; il < model->getDriftNumber(); il++)
  {
    if (model->getDriftType(il) == type0) return (1);
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Returns the st. dev. at a given increment for a given model
 **
 ** \param[in]  model       Structure containing the model
 ** \param[in]  db1         First Db
 ** \param[in]  iech1       Rank in the first Db
 ** \param[in]  db2         Second Db
 ** \param[in]  iech2       Rank in the second Db
 ** \param[in]  verbose     Verbose flag
 ** \param[in]  factor      Multiplicative factor for st. deviation
 **
 *****************************************************************************/
double model_calcul_stdev(Model *model,
                          Db *db1,
                          int iech1,
                          Db *db2,
                          int iech2,
                          int verbose,
                          double factor)
{
  int ndim;
  double c00, cov, stdev;
  CovCalcMode mode;

  /* Initializations */

  ndim = db1->getNDim();

  /* Covariance at origin */

  VectorDouble d1(ndim, 0.);
  model_calcul_cov(model, mode, 1, 1., d1, &c00);

  /* Covariance at increment */

  for (int idim = 0; idim < ndim; idim++)
    d1[idim] = db1->getDistance1D(iech1, iech2, idim);
  model_calcul_cov(model, mode, 1, 1., d1, &cov);

  stdev = factor * sqrt(c00 - cov);

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
  double *covtab, v1, v2, value;
  int ndim, nvar, nvar1, nvar2, iech1, iech2, i, skip, error, i1, i2;
  VectorDouble d1;
  cs *T = nullptr;
  cs *covmat = nullptr;

  /* Initializations */

  error = 1;
  covtab = nullptr;
  CovCalcMode mode;
  mode.update(0, 0, ECalcMember::LHS, -1, flag_norm, flag_cov);
  if (st_check_model(model)) goto label_end;
  if (st_check_environ(model, db1)) goto label_end;
  if (st_check_environ(model, db2)) goto label_end;
  ndim = model->getDimensionNumber();
  nvar = nvar1 = nvar2 = model->getVariableNumber();
  if (ivar0 >= 0)
  {
    nvar1 = 1;
    if (st_check_variable(nvar, ivar0)) goto label_end;
  }
  if (jvar0 >= 0)
  {
    nvar2 = 1;
    if (st_check_variable(nvar, jvar0)) goto label_end;
  }
  if (ranks1 == nullptr) nsize1 = db1->getSampleNumber();
  if (ranks2 == nullptr) nsize2 = db1->getSampleNumber();

  /* Core allocation */

  d1.resize(ndim, 0.);
  covtab = (double*) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == nullptr) goto label_end;

  // Constitute the triplet

  T = cs_spalloc(0, 0, 1, 1, 1);
  if (T == nullptr) goto label_end;

  /* Loop on the number of variables */

  for (int ivar = 0; ivar < nvar1; ivar++)
  {
    if (ivar0 >= 0) ivar = ivar0;

    /* Loop on the first sample */

    for (i1 = 0; i1 < nsize1; i1++)
    {
      iech1 = (ranks1 != nullptr) ? ranks1[i1] :
                                    i1;
      if (iech1 < 0) continue;

      /* Loop on the second variable */

      for (int jvar = 0; jvar < nvar2; jvar++)
      {
        if (jvar0 >= 0) jvar = jvar0;

        /* Loop on the second sample */

        for (i2 = 0; i2 < nsize2; i2++)
        {
          iech2 = (ranks2 != nullptr) ? ranks2[i2] :
                                        i2;
          if (iech2 < 0) continue;

          /* Determine the indices */

          int ecr1 = ivar * nsize1 + i1;
          int ecr2 = jvar * nsize2 + i2;

          // Save calculations due to matrix symmetry

          if (ecr2 > ecr1) continue;

          /* Loop on the dimension of the space */

          value = TEST;
          for (i = skip = 0; i < ndim && skip == 0; i++)
          {
            v1 = db1->getCoordinate(iech1, i);
            v2 = db2->getCoordinate(iech2, i);
            if (FFFF(v1) || FFFF(v2)) skip = 1;
            d1[i] = v1 - v2;
          }
          if (skip) continue;

          model_calcul_cov(model, mode, 1, 1., d1, covtab);
          value = COVTAB(ivar, jvar);
          if (ABS(value) < EPSILON10) continue;
          if (! cs_entry(T, ecr1, ecr2, value)) goto label_end;
          if (ecr1 != ecr2)
          {
            if (! cs_entry(T, ecr2, ecr1, value)) goto label_end;
          }
        }
      }
    }
  }

  // Convert from triplet to sparse matrix

  covmat = cs_triplet(T);
  T = cs_spfree(T);
  error = 0;

  /* Free memory */

  label_end: T = cs_spfree(T);
  covtab = (double*) mem_free((char* ) covtab);
  if (error) covmat = cs_spfree(covmat);
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
void model_covmat(Model *model,
                  Db *db1,
                  Db *db2,
                  int ivar0,
                  int jvar0,
                  int flag_norm,
                  int flag_cov,
                  double *covmat)
{
  double *covtab, v1, v2, value;
  int ndim, nvar, nvar1, nvar2, iech1, iech2, i, skip, nech1, nech2, ecr;
  VectorDouble d1;

  /* Initializations */

  covtab = nullptr;
  CovCalcMode mode;
  mode.update(0, 0, ECalcMember::LHS, -1, flag_norm, flag_cov);
  if (db2 == nullptr) db2 = db1;
  if (st_check_model(model)) goto label_end;
  if (st_check_environ(model, db1)) goto label_end;
  if (st_check_environ(model, db2)) goto label_end;
  ndim = model->getDimensionNumber();
  nvar = model->getVariableNumber();
  nech1 = db1->getSampleNumber();
  nech2 = db2->getSampleNumber();
  nvar1 = nvar2 = nvar;
  if (ivar0 >= 0)
  {
    nvar1 = 1;
    if (st_check_variable(nvar, ivar0)) goto label_end;
  }
  if (jvar0 >= 0)
  {
    nvar2 = 1;
    if (st_check_variable(nvar, jvar0)) goto label_end;
  }

  /* Core allocation */

  d1.resize(ndim, 0);
  covtab = (double*) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == nullptr) goto label_end;

  /* Loop on the first variable */

  ecr = 0;
  for (int ivar = 0; ivar < nvar1; ivar++)
  {
    if (ivar0 >= 0) ivar = ivar0;

    /* Loop on the first sample */

    for (iech1 = 0; iech1 < nech1; iech1++)
    {
      if (!db1->isActive(iech1)) continue;

      /* Loop on the second variable */

      for (int jvar = 0; jvar < nvar2; jvar++)
      {
        if (jvar0 >= 0) jvar = jvar0;

        /* Loop on the second sample */

        for (iech2 = 0; iech2 < nech2; iech2++)
        {
          if (!db2->isActive(iech2)) continue;

          /* Loop on the dimension of the space */

          value = TEST;
          for (i = skip = 0; i < ndim && skip == 0; i++)
          {
            v1 = db1->getCoordinate(iech1, i);
            v2 = db2->getCoordinate(iech2, i);
            if (FFFF(v1) || FFFF(v2)) skip = 1;
            d1[i] = v1 - v2;
          }
          if (!skip)
          {
            if (model->isNoStat())
            {
              CovInternal covint(1, iech1, 2, iech2, ndim, db1, db2);
              model_calcul_cov_nostat(model, mode, &covint, 1, 1., d1, covtab);
            }
            else
            {
              model_calcul_cov(model, mode, 1, 1., d1, covtab);
            }
            value = COVTAB(ivar, jvar);
          }
          covmat[ecr++] = value;
        }
      }
    }
  }

  /* Free memory */

  label_end: covtab = (double*) mem_free((char* ) covtab);
  return;
}

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
#include "geoslib_e.h"
#include "Drifts/DriftFactory.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Covariances/ACovAnisoList.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovGradientNumerical.hpp"
#include "Space/SpaceRN.hpp"

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

static double EPS = 1.e-6;

int NDIM_LOCAL = 0;
VectorDouble X1_LOCAL = VectorDouble();
VectorDouble X2_LOCAL = VectorDouble();

/*****************************************************************************/
/*!
 **  Update the Model in the case of Non-stationary parameters
 **  This requires the knowledge of the two end-points
 **
 ** \param[in]  cov_nostat   Internal structure for non-stationarity
 **                          or NULL (for stationary case)
 ** \param[in]  model        Model structure
 **
 *****************************************************************************/
GEOSLIB_API void model_nostat_update(CovNostatInternal *cov_nostat,
                                     Model* model)
{
  if (!model->isNoStat()) return;
  if (cov_nostat == NULL) return;

  int iech1 = cov_nostat->iech1;
  int iech2 = cov_nostat->iech2;

  model->getNoStat().updateModel(model, iech1, iech2);
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
static int st_check_environ(Model *model, Db *db)
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
  if (model != (Model *) NULL) return (0);
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
GEOSLIB_API void model_covtab_init(int flag_init, Model *model, double *covtab)
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
 ** \param[in]  nvar   Number of variables
 ** \param[in]  norme  Number of values used for scaling
 ** \param[in]  covtab Input array covtab
 **
 ** \param[out]  covtab Output array covtab
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
 ** \param[in]  member      Member of the Kriging System (::ENUM_MEMBERS)
 ** \param[in]  d1          vector of increment (or NULL)
 **
 *****************************************************************************/
GEOSLIB_API double model_calcul_basic(Model *model,
                                      int icov,
                                      int member,
                                      const VectorDouble& d1)
{
  //  const CovAniso* cova = model->getCova(icov);
  // TODO: Why having to use ACov rather than CvoAniso
  const ACov* cova = dynamic_cast<const ACov*>(model->getCova(icov));

  if (member != MEMBER_LHS && model->isCovaFiltered(icov))
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
 ** \param[in]  cov_nostat   Internal structure for non-stationarity
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
static void st_model_calcul_cov_direct(CovNostatInternal *cov_nostat,
                                       Model* model,
                                       const CovCalcMode& mode,
                                       int flag_init,
                                       double weight,
                                       VectorDouble d1,
                                       double *covtab)
{
  // Load the non-stationary parameters if needed

  if (model->isNoStat()) model_nostat_update(cov_nostat, model);

  // Evaluate the Model

  MatrixCSGeneral mat;
  if (d1.empty())
    mat = model->getCovAnisoList()->ACov::eval0(mode);
  else
    mat = model->getCovAnisoList()->ACov::eval(1., d1, SpacePoint(), mode);

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
static void st_model_calcul_cov_convolution(CovNostatInternal *cov_nostat,
                                            Model *model,
                                            const CovCalcMode& mode,
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
 ** \param[in]  weight       Weight attached to this calculation
 **
 ** \param[out] d1          Working array (dimension = ndim) or NULL
 ** \param[out] covtab      output covariance (dimension = nvar * nvar)
 **
 ** \remark This function is only programmed in the monovariate case
 **
 *****************************************************************************/
static void st_model_calcul_cov_anam_hermitian(CovNostatInternal *cov_nostat,
                                               Model *model,
                                               const CovCalcMode& mode,
                                               int flag_init,
                                               double weight,
                                               VectorDouble d1,
                                               double *covtab)
{
  double rho, covfin, dist2, coeff, psin2, rn, rhon;
  int iclass;

  /* Initializations */

  ModTrans& modtrs = model->getModTrans();
  AnamHermite* anam_hermite = dynamic_cast<AnamHermite*>(modtrs.getAnam());

  /* Check if the distance is zero */

  dist2 = rho = covfin = coeff = 0.;
  if (!d1.empty()) dist2 = matrix_norm(d1.data(), model->getDimensionNumber());
  if (dist2 > 0)
  {

    /* Calculate the generic variogram value */

    st_model_calcul_cov_direct(NULL, model, mode, 1, 1, d1, &rho);
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
        switch (mode.getMember())
        {
          case MEMBER_LHS:
            coeff = 1. / (rn * rn);
            break;
          case MEMBER_RHS:
            coeff = 1. / rn;
            break;
          case MEMBER_VAR:
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
        switch (mode.getMember())
        {
          case MEMBER_LHS:
            coeff = 1.;
            break;

          case MEMBER_RHS:
            coeff = 1. / rn;
            break;

          case MEMBER_VAR:
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
      switch (mode.getMember())
      {
        case MEMBER_LHS:
          coeff = 1.;
          break;

        case MEMBER_RHS:
          coeff = rn;
          break;

        case MEMBER_VAR:
          coeff = 1.;
          break;
      }
      covfin = coeff;
    }
    else
    {
      rn = pow(anam_hermite->getRCoef(), (double) modtrs.getAnamIClass());
      rhon = pow(rho, (double) modtrs.getAnamIClass());
      switch (mode.getMember())
      {
        case MEMBER_LHS:
          coeff = rn * rn;
          break;
        case MEMBER_RHS:
          coeff = rn;
          break;
        case MEMBER_VAR:
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
static void st_model_calcul_cov_anam_DD(CovNostatInternal *cov_nostat,
                                        Model *model,
                                        const CovCalcMode& mode,
                                        int flag_init,
                                        double weight,
                                        VectorDouble d1,
                                        double *covtab)
{
  double gamref, covfin, csi, li, mui, dist2, coeff;
  int iclass, ndim;

  /* Initializations */

  ModTrans& modtrs = model->getModTrans();
  AnamDiscreteDD* anam_discrete_DD = dynamic_cast<AnamDiscreteDD*>(modtrs.getAnam());
  ndim = model->getDimensionNumber();

  /* Check if the distance is zero */

  dist2 = gamref = covfin = coeff = 0.;
  if (!d1.empty()) dist2 = matrix_norm(d1.data(), model->getDimensionNumber());
  if (dist2 > 0)
  {
    VectorDouble covint(ndim * ndim);

    /* Calculate the generic variogram value */

    st_model_calcul_cov_direct(NULL, model, mode, 1, 1., VectorDouble(),
                               covint.data());
    st_model_calcul_cov_direct(NULL, model, mode, 0, -1, d1, covint.data());
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
        switch (mode.getMember())
        {
          case MEMBER_LHS:
            coeff = csi * csi;
            break;
          case MEMBER_RHS:
            coeff = csi * csi / mui;
            break;
          case MEMBER_VAR:
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
        switch (mode.getMember())
        {
          case MEMBER_LHS:
            coeff = csi * csi;
            break;

          case MEMBER_RHS:
            coeff = csi * csi / mui;
            break;

          case MEMBER_VAR:
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
      switch (mode.getMember())
      {
        case MEMBER_LHS:
          coeff = 1.;
          break;

        case MEMBER_RHS:
          coeff = mui;
          break;

        case MEMBER_VAR:
          coeff = 1.;
          break;
      }
      covfin = coeff;
    }
    else
    {
      mui = anam_discrete_DD->getDDStatMul(modtrs.getAnamIClass());
      li = anam_discrete_DD->getDDStatLambda(modtrs.getAnamIClass());
      switch (mode.getMember())
      {
        case MEMBER_LHS:
          coeff = mui * mui;
          break;
        case MEMBER_RHS:
          coeff = mui;
          break;
        case MEMBER_VAR:
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
                              const VectorDouble& covwrk)
{
  ModTrans& modtrs = model->getModTrans();
  AnamDiscreteIR* anam_discrete_IR = dynamic_cast<AnamDiscreteIR*>(modtrs.getAnam());

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
                                 const VectorDouble& covwrk)
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
static void st_model_calcul_cov_anam_IR(CovNostatInternal *cov_nostat,
                                        Model *model,
                                        const CovCalcMode& mode,
                                        int flag_init,
                                        double weight,
                                        VectorDouble d1,
                                        double *covtab)
{
  double covfin, bi, r, dist2, cov1, cov2, dcov;
  int icut, nclass, icut0, ncut, ndim, ncova;

  /* Initializations */

  ModTrans& modtrs = model->getModTrans();
  AnamDiscreteIR* anam_discrete_IR = dynamic_cast<AnamDiscreteIR*>(modtrs.getAnam());
  nclass = modtrs.getAnamNClass();
  ncut = nclass - 1;
  ncova = model->getCovaNumber();
  ;
  ndim = model->getDimensionNumber();
  r = anam_discrete_IR->getRCoef();
  VectorDouble covwrk(ncova, 0.);
  VectorDouble covint(ndim * ndim, 0.);

  /* Calculate the generic variogram value */

  st_model_calcul_cov_direct(NULL, model, mode, flag_init, 1., d1,
                             covint.data());

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
static void st_model_calcul_cov_tapering(CovNostatInternal *cov_nostat,
                                         Model *model,
                                         const CovCalcMode& mode,
                                         int flag_init,
                                         double weight,
                                         VectorDouble d1,
                                         double *covtab)
{
  double h;
  int nvar = model->getVariableNumber();

  /* Calculate the generic covariance value */

  st_model_calcul_cov_direct(NULL, model, mode, flag_init, weight, d1, covtab);

  /* Calculate the tapering effect */

  ModTrans& modtrs = model->getModTrans();
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
GEOSLIB_API void model_calcul_cov(Model *model,
                                  CovCalcMode& mode,
                                  int flag_init,
                                  double weight,
                                  const VectorDouble& d1,
                                  double *covtab)
{
  /* Modify the member in case of properties */

  if (model->getModTransMode() != MODEL_PROPERTY_NONE)
  {
    int anam_var = model->getModTrans().getAnamPointBlock();
    /* 'anam_var' is negative if model evaluation is called from dk() */
    /* This modification is performed in model_anamorphosis_set_factor() */
    if (anam_var >= 0) mode.setMember((ENUM_MEMBERS) anam_var);
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
GEOSLIB_API double model_calcul_cov_ij(Model *model,
                                       const CovCalcMode& mode,
                                       int ivar,
                                       int jvar,
                                       const VectorDouble& d1)
{

  /* Modify the member in case of properties */

  if (model->getModTransMode() != MODEL_PROPERTY_NONE)
  my_throw("Model transformation is not programmed for this entry");

  /* Call the generic model calculation module */

  // TODO Correct this which has something to do with pure virtual eval although implemented in Acov
  // compared to eval0 which is not implemented with such arguments.
  double value;
  if (d1.empty())
    value = model->getCovAnisoList()->eval0(ivar, jvar, mode);
  else
    value = model->getCovAnisoList()->ACov::eval(ivar, jvar, 1., d1,
                                                 SpacePoint(), mode);

  return value;
}

///*****************************************************************************/
///*!
// **  In the non-stationary case and if an external covariance function is
// **  present in the Model, allocate the relevant arrays and store
// **  the coordinates of the two end-points
// **
// ** \param[in]  cov_nostat   Internal structure for non-stationarity
// **                          or NULL (for stationary case)
// ** \param[in]  model        Model structure
// **
// *****************************************************************************/
//static void st_model_nostat_for_external(CovNostatInternal *cov_nostat,
//                                         Model *model)
//{
//  /// TODO [Cova] : Lea to be restored
//  int is_external,ndim;
//  Cova *cova;
//
//  // Check if an external structure is required
//
//  is_external = 0;
//  for (int icov=0; icov<model->getNCova(); icov++)
//  {
//    cova = model->getCova(icov);
//    if (cova->getType() < 0) is_external = 1;
//  }
//  if (! is_external) return;
//
//  // Check the dimension of already existing arrays
//
//  if (cov_nostat->db1 == (Db *) NULL) return;
//  if (cov_nostat->db2 == (Db *) NULL) return;
//  ndim = get_NDIM(cov_nostat->db1);
//
//  // Load the coordinates
//
//  if (ndim != (int) X1_LOCAL.size()) X1_LOCAL.resize(ndim);
//  if (ndim != (int) X2_LOCAL.size()) X2_LOCAL.resize(ndim);
//  for (int idim=0; idim<ndim; idim++)
//  {
//    X1_LOCAL[idim] = get_IDIM(cov_nostat->db1,cov_nostat->iech1,idim);
//    X2_LOCAL[idim] = get_IDIM(cov_nostat->db2,cov_nostat->iech2,idim);
//  }
//}
/****************************************************************************/
/*!
 **  Identify the coordinates of the two end-points
 **  (used by the external covariance function)
 **
 ** \param[out] E_Cov    The External_Cov structure
 **
 ** \remarks The arguments 'x1' and 'x2' are assigned
 **
 *****************************************************************************/
GEOSLIB_API void fill_external_cov_model(External_Cov& E_Cov)
{
  E_Cov.x1 = X1_LOCAL;
  E_Cov.x2 = X2_LOCAL;
}

/*****************************************************************************/
/*!
 **  Returns the covariances for an increment
 **  The increment is calculated between two samples of two Dbs
 **  This function is available in the non-stationary case
 **
 ** \param[in]  model       Structure containing the model
 ** \param[in]  mode        CovCalcMode structure
 ** \param[in]  flag_init   initialize the array beforehand
 ** \param[in]  weight      Weight attached to this calculation
 ** \param[in]  db1         First Db structure
 ** \param[in]  iech1       Rank of the sample in the First Db
 ** \param[in]  db2         Second Db structure
 ** \param[in]  iech2       Rank of the sample in the Second Db
 **
 ** \param[out] d1          Working array (dimension = ndim) or NULL
 ** \param[out] covtab      output covariance (dimension = nvar * nvar)
 **
 *****************************************************************************/
GEOSLIB_API void model_calcul_cov_nostat(Model *model,
                                         CovCalcMode& mode,
                                         int flag_init,
                                         double weight,
                                         Db *db1,
                                         int iech1,
                                         Db *db2,
                                         int iech2,
                                         VectorDouble& d1,
                                         double *covtab)
{
  CovNostatInternal cov_nostat;

  /* Load the non_stationary parameters */

  cov_nostat.iech1 = iech1;
  cov_nostat.iech2 = iech2;
  cov_nostat.db1 = db1;
  cov_nostat.db2 = db2;

  /* Modify the member in case of properties */

  if (model->getModTransMode() != MODEL_PROPERTY_NONE)
  {
    my_throw("ModTrans not yet implemented");
    int anam_var = model->getModTrans().getAnamPointBlock();

    /* 'anam_var' is negative if model evaluation is called from df() */
    /* This modification is performed in model_anamorphosis_set_factor() */
    if (anam_var >= 0) mode.setMember((ENUM_MEMBERS) anam_var);
  }

  /* Store the information from the data base if needed by an external */
  /* covariance function */

  /// TODO [Cova] : Lea to be restored
  // st_model_nostat_for_external(&cov_nostat,model);
  /* Call the generic model calculation module */

  model->generic_cov_function(&cov_nostat, model, mode, flag_init, weight, d1,
                              covtab);

  return;
}

/****************************************************************************/
/*!
 **  Returns the drift vector for a point
 **
 ** \param[in]  model  Model structure
 ** \param[in]  member Member of the Kriging System (::ENUM_MEMBERS)
 ** \param[in]  db     Db structure
 ** \param[in]  iech   Rank of the sample
 **
 ** \param[out] drftab  Array of drift values
 **
 *****************************************************************************/
GEOSLIB_API void model_calcul_drift(Model *model,
                                    int member,
                                    Db *db,
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
GEOSLIB_API void model_variance0(Model *model,
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
  mode.setMember(MEMBER_VAR);

  switch (koption->calcul)
  {
    case KOPTION_PONCTUAL:
      nscale = 1;
      model_calcul_cov(model, mode, 1, 1., d1, covtab);
      break;

    case KOPTION_BLOCK:
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

    case KOPTION_DRIFT:
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
 ** \param[in]  db0       Db structure
 ** \param[in]  iech0     Rank of the sample in the Db
 ** \param[in]  covtab    array of cumulated covariances
 **
 ** \param[out]  var0     array for C0[] (Dimension = nvar * nvar)
 **
 *****************************************************************************/
GEOSLIB_API void model_variance0_nostat(Model *model,
                                        Koption *koption,
                                        Db *db0,
                                        int iech0,
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

  mode.setMember(MEMBER_VAR);
  switch (koption->calcul)
  {
    case KOPTION_PONCTUAL:
      nscale = 1;
      model_calcul_cov_nostat(model, mode, 1, 1., db0, iech0, db0, iech0, d1,
                              covtab);
      break;

    case KOPTION_BLOCK:
      nscale = koption->ntot;
      model_covtab_init(1, model, covtab);
      for (i = 0; i < nscale; i++)
        for (j = 0; j < nscale; j++)
        {
          for (idim = 0; idim < model->getDimensionNumber(); idim++)
            d1[idim] = DISC1(i,idim) - DISC2(j, idim);
          model_calcul_cov_nostat(model, mode, 0, 1., db0, iech0, db0, iech0,
                                  d1, covtab);
        }
      nscale = nscale * nscale;
      break;

    case KOPTION_DRIFT:
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
GEOSLIB_API Model *model_free(Model *model)

{
  /* Initializations */

  if (model == (Model *) NULL) return (model);
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
 ** \param[in]   type0    Requested type (ENUM_CONS)
 **
 *****************************************************************************/
GEOSLIB_API int is_model_nostat_param(Model *model, ENUM_CONS type0)
{
  NoStatArray& nostat = model->getNoStat();
  if (nostat.isDefinedByType(-1, type0)) return 1;

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
GEOSLIB_API Model *model_init(int ndim,
                              int nvar,
                              double field,
                              int flag_linked,
                              double ball_radius,
                              bool flag_gradient,
                              const VectorDouble& mean,
                              const VectorDouble& covar0)
{
  int error = 1;
  Model* model = (Model *) NULL;

  ASpaceObject::createGlobalSpace(SPACE_RN, ndim); // TODO Avoid this artificial setting
  CovContext ctxt = CovContext(nvar, 2, field);
  ctxt.setBallRadius(ball_radius);
  if (mean.size() > 0)   ctxt.setMean(mean);
  if (covar0.size() > 0) ctxt.setCovar0(covar0);

  model = new Model(ctxt, flag_gradient, flag_linked);

  /* Set the error return flag */

  if (model_setup(model)) goto label_end;

  /* Set the error return code */

  error = 0;

  label_end: if (error) model = model_free(model);
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
GEOSLIB_API Model *model_default(int ndim, int nvar)
{
  Model *model;
  double sill;
  int error;

  /* Initializations */

  error = 1;

  /* Create the empty Model */

  model = model_init(ndim, nvar, 1.);

  /* Add the nugget effect variogram model */

  sill = 1.;
  if (model_add_cova(model, COV_NUGGET, 0, 0, 0., 0., VectorDouble(),
                     VectorDouble(), VectorDouble(sill))) goto label_end;

  /* Set the error return flag */

  if (model_setup(model)) goto label_end;
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
 ** \param[in]  type            Type of the basic structure (::ENUM_COVS)
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
GEOSLIB_API int model_add_cova(Model *model,
                               int type,
                               int flag_aniso,
                               int flag_rotation,
                               double range,
                               double param,
                               const VectorDouble& aniso_ranges,
                               const VectorDouble& aniso_rotmat,
                               const VectorDouble& sill)
{
  /// TODO : use ENUM_COVS
  if (st_check_model(model)) return 1;

  // Add a new element in the structure

  if (model->isFlagGradient())
  {
    CovGradientNumerical covgrad((ENUM_COVS) type, model->getContext());
    covgrad.setParam(param);
    if (flag_aniso)
    {
      covgrad.setRanges(aniso_ranges);
      if (flag_rotation) covgrad.setAnisoRotation(aniso_rotmat);
    }
    else
      covgrad.setRange(range);

    if (sill.size() > 0) covgrad.setSill(sill);
    model->addCova(&covgrad);
  }
  else
  {
    CovAniso cova((ENUM_COVS) type, model->getContext());
    cova.setParam(param);
    if (flag_aniso)
    {
      cova.setRanges(aniso_ranges);
      if (flag_rotation) cova.setAnisoRotation(aniso_rotmat);
    }
    else
      cova.setRange(range);

    if (sill.size() > 0) cova.setSill(sill);
    model->addCova(&cova);
  }

  /* Set the error return code */

  if (model_setup(model)) return 1;

  return 0;
}

/****************************************************************************/
/*!
 **  Add a basic drift
 **
 ** \return  Error return code
 **
 ** \param[in]  model      Pointer to the Model structure
 ** \param[in]  type       Type of the basic drift function (::ENUM_DRIFTS)
 ** \param[in]  rank_fex   Rank of the external drift (starting from 0)
 **
 ** \remark  If the variables are NOT linked, the number of drift equations
 ** \remark  is equal to: Number of variables * Number of drift functions
 ** \remark  If the variables are linked, the number of drift equations
 ** \remark  is equal to the number of drift functions.
 **
 *****************************************************************************/
GEOSLIB_API int model_add_drift(Model *model, int type, int rank_fex)
{
  ADriftElem *drift;
  int error;

  /* Initializations */

  error = 1;
  if (st_check_model(model)) goto label_end;

  // Allocate the new element

  drift = DriftFactory::createDriftFunc((ENUM_DRIFTS) type, model->getContext());
  drift->setRankFex(rank_fex);
  model->addDrift(drift);

  /* Set the error return code */

  if (model_setup(model)) goto label_end;
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
GEOSLIB_API int model_add_no_property(Model *model)

{
  int error;

  /* Initializations */

  error = 1;

  ModTrans& modtrs = model->getModTrans();
  if (st_check_model(model)) goto label_end;

  // Cancel the property

  modtrs.cancelProperty();

  /* Set the calling function */

  if (model_setup(model)) goto label_end;

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
GEOSLIB_API int model_add_convolution(Model *model,
                                      int conv_type,
                                      int conv_idir,
                                      int conv_ndisc,
                                      double conv_range)
{
  int error;

  /* Initializations */

  error = 1;
  ModTrans& modtrs = model->getModTrans();
  if (st_check_model(model)) goto label_end;

  /* Load the Convolution parameters */

  if (modtrs.addConvolution(conv_type, conv_idir, conv_ndisc, conv_range))
    goto label_end;

  /* Set the calling function */

  if (model_setup(model)) goto label_end;

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
GEOSLIB_API int model_anamorphosis_set_factor(Model *model, int anam_iclass)
{
  int error;

  /* Initializations */

  error = 1;
  ModTrans& modtrs = model->getModTrans();
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
 ** \param[in]  anam_type   Type of the anamorphosis
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
GEOSLIB_API int model_add_anamorphosis(Model *model,
                                       int anam_type,
                                       int anam_nclass,
                                       int anam_iclass,
                                       int anam_var,
                                       double anam_coefr,
                                       double anam_coefs,
                                       VectorDouble& anam_strcnt,
                                       VectorDouble& anam_stats)
{
  int error;

  /* Initializations */

  error = 1;
  ModTrans& modtrs = model->getModTrans();
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

  if (model_setup(model)) goto label_end;

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
GEOSLIB_API int model_add_tapering(Model *model,
                                   int tape_type,
                                   double tape_range)
{
  int error;

  /* Initializations */

  error = 1;
  ModTrans& modtrs = model->getModTrans();
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

  if (model_setup(model)) goto label_end;

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
GEOSLIB_API double cova_get_scale_factor(int type, double param)
{
  ACovFunc* cova = CovFactory::createCovFunc((ENUM_COVS) type, CovContext());
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
GEOSLIB_API int model_setup(Model* model)

{
  int error = 1;
  if (st_check_model(model)) goto label_end;

  /**********************************/
  /* Allocation of auxiliary arrays */
  /**********************************/

  /* Define the generic covariance function */

  switch (model->getModTransMode())
  {
    case MODEL_PROPERTY_NONE:
      model->generic_cov_function = st_model_calcul_cov_direct;
      break;

    case MODEL_PROPERTY_CONV:
      model->generic_cov_function = st_model_calcul_cov_convolution;
      break;

    case MODEL_PROPERTY_ANAM:
      switch (model->getModTrans().getAnam()->getType())
      {
        case ANAM_HERMITIAN:
          model->generic_cov_function = st_model_calcul_cov_anam_hermitian;
          break;

        case ANAM_DISCRETE_DD:
          model->generic_cov_function = st_model_calcul_cov_anam_DD;
          break;

        case ANAM_DISCRETE_IR:
          model->generic_cov_function = st_model_calcul_cov_anam_IR;
          break;

        default:
          messerr("The Model modified by Properties is not available");
          messerr("For the following Anamorphosis type (%d)",
                  model->getModTrans().getAnam()->getType());
          goto label_end;
      }
      break;

    case MODEL_PROPERTY_TAPE:
      model->generic_cov_function = st_model_calcul_cov_tapering;
      break;
  }

  /* Set the error return code */

  error = 0;

  label_end: return (error);
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
GEOSLIB_API int model_update_coreg(Model *model,
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
 ** \param[in]  member     Member of the Kriging System (::ENUM_MEMBERS)
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
GEOSLIB_API int model_evaluate(Model *model,
                               int ivar,
                               int jvar,
                               int rank_sel,
                               int flag_norm,
                               int flag_cov,
                               int nugget_opt,
                               int nostd,
                               int norder,
                               int member,
                               int nh,
                               VectorDouble& codir,
                               double *h,
                               double *g)
{
  int error = 1;
  double* covtab = (double *) NULL;
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
  covtab = (double *) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == (double *) NULL) goto label_end;

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

  label_end: covtab = (double *) mem_free((char * ) covtab);
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
 ** \param[in]  member     Member of the Kriging System (::ENUM_MEMBERS)
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
GEOSLIB_API int model_evaluate_nostat(Model *model,
                                      int ivar,
                                      int jvar,
                                      int rank_sel,
                                      int flag_norm,
                                      int flag_cov,
                                      int nugget_opt,
                                      int nostd,
                                      int norder,
                                      int member,
                                      Db *db1,
                                      int iech1,
                                      Db *db2,
                                      int iech2,
                                      int nh,
                                      VectorDouble& codir,
                                      double *h,
                                      double *g)
{
  double *covtab, var0, c00;
  int ih, nvar, idim, ndim, error;
  VectorDouble d1;

  /* Initializations */

  error = 1;
  covtab = (double *) NULL;
  CovCalcMode mode;
  mode.update(nugget_opt, nostd, member, rank_sel, flag_norm, flag_cov);
  if (norder > 0) mode.setOrderVario(norder);

  /* Preliminary checks */

  if (st_check_model(model)) goto label_end;
  ndim = model->getDimensionNumber();
  nvar = model->getVariableNumber();
  if (st_check_variable(nvar, ivar)) goto label_end;
  if (st_check_variable(nvar, jvar)) goto label_end;

  /* Core allocation */

  covtab = (double *) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == (double *) NULL) goto label_end;

  /* Normalize the direction vector codir */

  vario_fix_codir(ndim, codir);

  /* Calculate the C(0) term (used only for covariance or covariogram) */

  c00 = model->getContext().getCovar0(ivar, jvar);
  d1.resize(ndim, 0.);
  model_calcul_cov_nostat(model, mode, 1, 1., db1, iech1, db2, iech2, d1,
                          covtab);
  var0 = COVTAB(ivar, jvar);
  if (c00 <= 0. || FFFF(c00)) c00 = var0;

  /* Loop on the lags */

  for (ih = 0; ih < nh; ih++)
  {
    for (idim = 0; idim < ndim; idim++)
      d1[idim] = h[ih] * codir[idim];
    model_calcul_cov_nostat(model, mode, 1, 1., db1, iech1, db2, iech2, d1,
                            covtab);
    g[ih] = COVTAB(ivar, jvar);
  }

  /* Set the error return code */

  error = 0;

  label_end: covtab = (double *) mem_free((char * ) covtab);
  return (error);
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
GEOSLIB_API int model_grid(Model *model,
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
  covtab = (double *) NULL;
  CovCalcMode mode;
  mode.update(0, 0, MEMBER_LHS, -1, flag_norm, flag_cov);

  /* Preliminary checks */

  if (st_check_model(model)) goto label_end;
  if (st_check_environ(model, db)) goto label_end;
  ndim = model->getDimensionNumber();
  nvar = model->getVariableNumber();

  /* Core allocation */

  d1.resize(ndim);
  covtab = (double *) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == (double *) NULL) goto label_end;

  /* Initialization */

  for (iech = 0; iech < get_NECH(db); iech++)
    g[iech] = TEST;

  /* Loop on the lags */

  for (iech = 0; iech < get_NECH(db); iech++)
  {
    if (!db->isActive(iech)) continue;
    db_sample_load(db, LOC_X, iech, d1.data());
    model_calcul_cov(model, mode, 1, 1., d1, covtab);
    g[iech] = COVTAB(ivar, jvar);
  }

  /* Set the error return code */

  error = 0;

  label_end: covtab = (double *) mem_free((char * ) covtab);
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
GEOSLIB_API int model_nfex(Model *model)

{
  int il, nfex;

  /* Initializations */

  nfex = 0;
  if (model->getDriftNumber() <= 0) return (nfex);

  /* Loop on the drift functions */

  for (il = 0; il < model->getDriftNumber(); il++)
  {
    if (model->getDrift(il)->getType() == DRIFT_F) nfex++;
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
GEOSLIB_API void model_drift_filter(Model *model, int rank, int filter)
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
GEOSLIB_API double model_cxx(Model *model,
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
  covtab = (double *) NULL;
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
  covtab = (double *) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == (double *) NULL) goto label_end;
  model_covtab_init(1, model, covtab);

  /* Loop on the first sample */

  norme = 0.;
  for (iech1 = 0; iech1 < get_NECH(db1); iech1++)
  {
    if (!db1->isActive(iech1)) continue;
    w1 = db1->getWeight(iech1);
    if (w1 == 0.) continue;

    /* Loop on the second sample */

    for (iech2 = 0; iech2 < get_NECH(db2); iech2++)
    {
      if (!db2->isActive(iech2)) continue;
      w2 = db2->getWeight(iech2);
      if (w2 == 0.) continue;

      /* Loop on the dimension of the space */

      for (i = skip = 0; i < ndim && skip == 0; i++)
      {
        v1 = get_IDIM(db1, iech1, i);
        v2 = get_IDIM(db2, iech2, i);
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

  label_end: covtab = (double *) mem_free((char * ) covtab);
  return (cxx);
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
 *****************************************************************************/
GEOSLIB_API void model_covmat(Model *model,
                              Db *db1,
                              Db *db2,
                              int ivar0,
                              int jvar0,
                              int flag_norm,
                              int flag_cov,
                              double *covmat)
{
  double *covtab, *covtab0, v1, v2, value;
  int ndim, nvar, nvar1, nvar2, iech1, iech2, i, skip, nech1, nech2, ecr;
  VectorDouble d1;

  /* Initializations */

  covtab = covtab0 = (double *) NULL;
  CovCalcMode mode;
  mode.update(0, 0, MEMBER_LHS, -1, flag_norm, flag_cov);
  if (st_check_model(model)) goto label_end;
  if (st_check_environ(model, db1)) goto label_end;
  if (st_check_environ(model, db2)) goto label_end;
  ndim = model->getDimensionNumber();
  nvar = model->getVariableNumber();
  nech1 = get_NECH(db1);
  nech2 = get_NECH(db2);
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

  d1.resize(ndim, 0.);
  covtab0 = (double *) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab0 == (double *) NULL) goto label_end;
  covtab = (double *) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == (double *) NULL) goto label_end;

  /* Calculate the C(0) term */

  model_calcul_cov(model, mode, 1, 1., VectorDouble(), covtab0);

  /* Loop on the first sample */

  ecr = 0;

  /* Loop on the first variable */

  for (int ivar = 0; ivar < nvar1; ivar++)
  {
    if (ivar0 >= 0) ivar = ivar0;

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
            v1 = get_IDIM(db1, iech1, i);
            v2 = get_IDIM(db2, iech2, i);
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

  /* Free memory */

  label_end: covtab = (double *) mem_free((char * ) covtab);
  covtab0 = (double *) mem_free((char * ) covtab0);
  return;
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
 ** \param[in]  ranks1 Array giving ranks of selected samples
 ** \param[in]  db2    Second Db
 ** \param[in]  nsize2 Number of selected samples
 ** \param[in]  ranks2 Array giving ranks of selected samples
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
GEOSLIB_API double *model_covmat_by_ranks(Model *model,
                                          Db *db1,
                                          int nsize1,
                                          int *ranks1,
                                          Db *db2,
                                          int nsize2,
                                          int *ranks2,
                                          int ivar0,
                                          int jvar0,
                                          int flag_norm,
                                          int flag_cov)
{
  double *covmat, *covtab, *covtab0, v1, v2, value;
  int ndim, nvar, nvar1, nvar2, iech1, iech2, i, skip, ecr, error, i1, i2;
  VectorDouble d1;

  /* Initializations */

  error = 1;
  covtab = covtab0 = covmat = (double *) NULL;
  CovCalcMode mode;
  mode.update(0, 0, MEMBER_LHS, -1, flag_norm, flag_cov);
  if (st_check_model(model)) goto label_end;
  if (st_check_environ(model, db1)) goto label_end;
  if (st_check_environ(model, db2)) goto label_end;
  ndim = model->getDimensionNumber();
  nvar = model->getVariableNumber();
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

  d1.resize(ndim, 0.);
  covtab = (double *) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == (double *) NULL) goto label_end;
  covtab0 = (double *) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab0 == (double *) NULL) goto label_end;
  covmat = (double *) mem_alloc(sizeof(double) * nsize1 * nsize2, 0);
  if (covmat == (double *) NULL) goto label_end;

  /* Calculate the C(0) term */

  model_calcul_cov(model, mode, 1, 1., VectorDouble(), covtab0);
  ecr = 0;

  /* Loop on the number of variables */

  for (int ivar = 0; ivar < nvar1; ivar++)
  {
    if (ivar0 >= 0) ivar = ivar0;

    /* Loop on the first sample */

    for (i1 = 0; i1 < nsize1; i1++)
    {
      iech1 = (ranks1 != (int *) NULL) ? ranks1[i1] :
                                         i1;
      if (iech1 < 0) continue;

      /* Loop on the second variable */

      for (int jvar = 0; jvar < nvar2; jvar++)
      {
        if (jvar0 >= 0) jvar = jvar0;

        /* Loop on the second sample */

        for (i2 = 0; i2 < nsize2; i2++)
        {
          iech2 = (ranks2 != (int *) NULL) ? ranks2[i2] :
                                             i2;
          if (iech2 < 0) continue;

          /* Loop on the dimension of the space */

          value = TEST;
          for (i = skip = 0; i < ndim && skip == 0; i++)
          {
            v1 = get_IDIM(db1, iech1, i);
            v2 = get_IDIM(db2, iech2, i);
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

  label_end: covtab = (double *) mem_free((char * ) covtab);
  covtab0 = (double *) mem_free((char * ) covtab0);
  if (error) covmat = (double *) mem_free((char * ) covmat);
  return (covmat);
}

/****************************************************************************/
/*!
 **  Establish the drift rectangular matrix for a given Db
 **
 ** \param[in]  model  Model structure
 ** \param[in]  member Member of the Kriging System (::ENUM_MEMBERS)
 ** \param[in]  db     Db structure
 **
 ** \param[out] drfmat The drift matrix
 **                    (Dimension = nech * nvar * nfeq * nvar)
 **
 *****************************************************************************/
GEOSLIB_API void model_drift_mat(Model *model,
                                 int member,
                                 Db *db,
                                 double *drfmat)
{
  int nech, nvar, nbfl, nfeq, ecr, jb;
  double *drftab, value;

  /* Initializations */

  drftab = (double *) NULL;
  if (st_check_model(model)) goto label_end;
  if (st_check_environ(model, db)) goto label_end;
  nvar = model->getVariableNumber();
  nbfl = model->getDriftNumber();
  nfeq = model->getDriftEquationNumber();
  nech = get_NECH(db);

  /* Core allocation */

  drftab = (double *) mem_alloc(sizeof(double) * nbfl, 0);
  if (drftab == (double *) NULL) goto label_end;

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

  label_end: drftab = (double *) mem_free((char * ) drftab);
  return;
}

/****************************************************************************/
/*!
 **  Establish the drift vector for a given sample of the Db
 **
 ** \param[in]  model  Model structure
 ** \param[in]  member Member of the Kriging System (::ENUM_MEMBERS)
 ** \param[in]  db     Db structure
 ** \param[in]  iech   Rank of the particular sample
 **
 ** \param[out] vector Returned vector
 **                    (Dimension = nvar * nfeq)
 **
 *****************************************************************************/
GEOSLIB_API void model_drift_vector(Model *model,
                                    int member,
                                    Db *db,
                                    int iech,
                                    double *vector)
{
  int nvar, nbfl, nfeq, ivar, ib, il, ecr, i;
  double *drftab, value;

  /* Initializations */

  drftab = (double *) NULL;
  if (st_check_model(model)) goto label_end;
  if (st_check_environ(model, db)) goto label_end;
  nvar = model->getVariableNumber();
  nbfl = model->getDriftNumber();
  nfeq = model->getDriftEquationNumber();

  /* Core allocation */

  drftab = (double *) mem_alloc(sizeof(double) * nbfl, 0);
  if (drftab == (double *) NULL) goto label_end;

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

  label_end: drftab = (double *) mem_free((char * ) drftab);
  return;
}

/****************************************************************************/
/*!
 **  Establish the covariance matrix between two Dbs (Non stationary case)
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
 *****************************************************************************/
GEOSLIB_API void model_covmat_nostat(Model *model,
                                     Db *db1,
                                     Db *db2,
                                     int ivar0,
                                     int jvar0,
                                     int flag_norm,
                                     int flag_cov,
                                     double *covmat)
{
  double *covtab, *covtab0, v1, v2, value;
  int ndim, nvar, nvar1, nvar2, iech1, iech2, i, skip, nech1, nech2, ecr;
  VectorDouble d1;

  /* Initializations */

  covtab = covtab0 = (double *) NULL;
  CovCalcMode mode;
  mode.update(0, 0, MEMBER_LHS, -1, flag_norm, flag_cov);
  if (st_check_model(model)) goto label_end;
  if (st_check_environ(model, db1)) goto label_end;
  if (st_check_environ(model, db2)) goto label_end;
  ndim = model->getDimensionNumber();
  nvar = model->getVariableNumber();
  nech1 = get_NECH(db1);
  nech2 = get_NECH(db2);
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
  covtab = (double *) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == (double *) NULL) goto label_end;
  covtab0 = (double *) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab0 == (double *) NULL) goto label_end;

  /* Calculate the C(0) term */

  model_calcul_cov(model, mode, 1, 1., VectorDouble(), covtab0);
  ecr = 0;

  /* Loop on the first variable */

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
            v1 = get_IDIM(db1, iech1, i);
            v2 = get_IDIM(db2, iech2, i);
            if (FFFF(v1) || FFFF(v2)) skip = 1;
            d1[i] = v1 - v2;
          }
          if (!skip)
          {
            model_calcul_cov_nostat(model, mode, 1, 1., db1, iech1, db2, iech2,
                                    d1, covtab);
            value = COVTAB(ivar, jvar);
          }
          covmat[ecr++] = value;
        }
      }
    }
  }

  /* Free memory */

  label_end: covtab = (double *) mem_free((char * ) covtab);
  covtab0 = (double *) mem_free((char * ) covtab0);
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
                          CovCalcMode& mode,
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

  /****************************************************************************/
  /*!
   **  Establish the multivariate covariance matrix for one Db
   **
   ** \param[in]  model Model structure
   ** \param[in]  db    Db structure
   ** \param[in]  flag_norm  1 if the model is normalized
   ** \param[in]  flag_cov   1 if the result must be given in covariance
   **
   ** \param[out] covmat The covariance matrix
   **                    (Dimension = neq * neq) where neq = nactive * nvar
   **
   *****************************************************************************/
GEOSLIB_API void model_covmat_multivar(Model *model,
                                       Db *db,
                                       int flag_norm,
                                       int flag_cov,
                                       double *covmat)
{
  double *covtab, *c00tab, v1, v2, value;
  int ndim, nvar, nech, ivar1, ivar2, iech1, iech2, i, skip, ecr;
  VectorDouble d1;

  /* Initializations */

  covtab = c00tab = (double *) NULL;
  CovCalcMode mode;
  mode.update(0, 0, MEMBER_LHS, -1, flag_norm, flag_cov);
  if (st_check_model(model)) goto label_end;
  if (st_check_environ(model, db)) goto label_end;
  ndim = model->getDimensionNumber();
  nvar = model->getVariableNumber();
  nech = get_NECH(db);

  /* Core allocation */

  d1.resize(ndim, 0.);
  covtab = (double *) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == (double *) NULL) goto label_end;
  c00tab = (double *) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (c00tab == (double *) NULL) goto label_end;

  /* Calculate the C(0) term */

  st_matrix_c00(model, mode, covtab, d1, c00tab);

  /* Loop on the first sample */

  ecr = 0;
  for (ivar1 = 0; ivar1 < nvar; ivar1++)
    for (iech1 = 0; iech1 < nech; iech1++)
    {
      if (!db->isActive(iech1)) continue;

      /* Loop on the second sample */

      for (ivar2 = 0; ivar2 < nvar; ivar2++)
        for (iech2 = 0; iech2 < nech; iech2++)
        {
          if (!db->isActive(iech2)) continue;

          /* Loop on the dimension of the space */

          value = TEST;
          for (i = skip = 0; i < ndim && skip == 0; i++)
          {
            v1 = get_IDIM(db, iech1, i);
            v2 = get_IDIM(db, iech2, i);
            if (FFFF(v1) || FFFF(v2)) skip = 1;
            d1[i] = v1 - v2;
          }
          if (!skip)
          {
            model_calcul_cov(model, mode, 1, 1., d1, covtab);
            value = COVTAB(ivar1, ivar2);
            if (flag_norm) value /= sqrt(C00TAB(ivar1,ivar1)* C00TAB(ivar2,ivar2));
          }
          covmat[ecr++] = value;
        }
    }

  /* Free memory */

  label_end: covtab = (double *) mem_free((char * ) covtab);
  c00tab = (double *) mem_free((char * ) c00tab);
  return;
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
 ** \param[in]  type  type of the drift function
 ** \param[in]  rank  rank of the external drift
 ** \param[in]  value value to be added to the drift coefficient
 **
 *****************************************************************************/
static void st_drift_modify(Model *model,
                            int iv,
                            int ib,
                            int type,
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
                                Model *model,
                                Model *new_model)

{
  int il, ib, type, rank;
  double value;

  for (ib = 0; ib < model->getDriftEquationNumber(); ib++)
    for (il = 0; il < model->getDriftNumber(); il++)
    {
      value = model->getCoefDrift(0, il, ib);
      if (value == 0) continue;

      type = model->getDriftType(il);
      rank = model->getRankFext(il);
      switch (mode)
      {
        case MODEL_DERIVATIVE_NONE: /* Simple copy */
          st_drift_modify(new_model, iv, ib, type, rank, value);
          break;

        case MODEL_DERIVATIVE_X: /* Derivative along X */
          switch (type)
          {
            case DRIFT_X:
              st_drift_modify(new_model, iv, ib, DRIFT_1, 0, value);
              break;
            case DRIFT_X2:
              st_drift_modify(new_model, iv, ib, DRIFT_X, 0, 2. * value);
              break;
            case DRIFT_XY:
              st_drift_modify(new_model, iv, ib, DRIFT_Y, 0, value);
              break;
            case DRIFT_XZ:
              st_drift_modify(new_model, iv, ib, DRIFT_Z, 0, value);
              break;
            case DRIFT_X3:
              st_drift_modify(new_model, iv, ib, DRIFT_X2, 0, 3. * value);
              break;
            case DRIFT_X2Y:
              st_drift_modify(new_model, iv, ib, DRIFT_XY, 0, 2. * value);
              break;
            case DRIFT_XY2:
              st_drift_modify(new_model, iv, ib, DRIFT_Y2, 0, value);
              break;
            default:
              break;
          }
          break;

        case MODEL_DERIVATIVE_Y: /* Derivative along Y */
          switch (type)
          {
            case DRIFT_Y:
              st_drift_modify(new_model, iv, ib, DRIFT_1, 0, value);
              break;
            case DRIFT_Y2:
              st_drift_modify(new_model, iv, ib, DRIFT_Y, 0, 2. * value);
              break;
            case DRIFT_XY:
              st_drift_modify(new_model, iv, ib, DRIFT_X, 0, value);
              break;
            case DRIFT_YZ:
              st_drift_modify(new_model, iv, ib, DRIFT_Z, 0, value);
              break;
            case DRIFT_Y3:
              st_drift_modify(new_model, iv, ib, DRIFT_Y2, 0, 3. * value);
              break;
            case DRIFT_XY2:
              st_drift_modify(new_model, iv, ib, DRIFT_XY, 0, 2. * value);
              break;
            case DRIFT_X2Y:
              st_drift_modify(new_model, iv, ib, DRIFT_X2, 0, value);
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
 **  Duplicates a Model from another Model (1 variable in 2-D)
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
GEOSLIB_API Model *model_duplicate(Model *model, double ball_radius, int mode)

{
  Model *new_model;
  CovAniso *cova;
  ADriftElem *drft;
  int flag_linked, new_nvar, nfact;
  double sill;
  bool flag_gradient;

  // Preliminary checks

  new_model = (Model *) NULL;
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
                         ball_radius, flag_gradient, model->getContext().getMean(),
                         model->getContext().getCovar0());

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
      ADriftElem* newdrft = DriftFactory::createDriftFunc(drft->getType(),
                                                          new_model->getContext());
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
        new_model->setCoefDrift(i, 0.);
      st_drift_derivative(0, MODEL_DERIVATIVE_NONE, model, new_model);
      st_drift_derivative(1, MODEL_DERIVATIVE_X, model, new_model);
      st_drift_derivative(2, MODEL_DERIVATIVE_Y, model, new_model);
    }
  }

  // Set the error return code

  if (model_setup(new_model)) return new_model;

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
//GEOSLIB_API Model *model_modify(Model  *model,
//                                int     new_nvar,
//                                double *mean,
//                                double *vars,
//                                double *corr)
//{
//  TODO : Dead code ?
//  Model  *new_model;
//  int     ivar,jvar,nvar,icov,ncova,il,nbfl,error,ndim;
//  double  sill;
//  Cova   *cova,*cova_new;
//  Drift   *drft;
//  VectorDouble ranges;
//
//  /* Initializations */
//
//  new_model = (Model *) NULL;
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
//  if (vars != (double *) NULL && ! is_matrix_non_negative(1,new_nvar,vars,0))
//  {
//    messerr("You provided vars[]. It must be non negative");
//    goto label_end;
//  }
//  if (vars != (double *) NULL && vars[0] == 0.)
//  {
//    messerr("You provided vars[]. It must have vars[0] != 0");
//    goto label_end;
//  }
//  if (corr != (double *) NULL && ! is_matrix_correlation(new_nvar,corr))
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
//    new_model->setMean(ivar, (mean != (double *) NULL) ? mean[ivar] :
//                                                         model->getMean(0));
//
//  /* Set the variance-covariance at the origin (if provided) */
//
//  new_model->setCovar0(0,0,model->getCovar0(0,0));
//  for (ivar=0; ivar<new_nvar; ivar++)
//    for (jvar=0; jvar<new_nvar; jvar++)
//    {
//      if (vars != (double *) NULL)
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
//            (vars == (double *) NULL || corr == (double *) NULL) ?
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
GEOSLIB_API int model_normalize(Model *model, int flag_verbose)

{
  double *total;
  int nvar, ncov, error, ivar, jvar, icov, flag_norm;

  /* Initializations */

  error = 1;
  nvar = model->getVariableNumber();
  ncov = model->getCovaNumber();
  total = (double *) NULL;

  /* Core allocation */

  total = (double *) mem_alloc(sizeof(double) * nvar, 1);

  /* Calculate the total sills for each variable */

  for (ivar = flag_norm = 0; ivar < nvar; ivar++)
  {
    total[ivar] = model->getCovAnisoList()->getTotalSill(ivar, ivar);
    if (total[ivar] == 0.) goto label_end;
    total[ivar] = sqrt(total[ivar]);
    if (ABS(total[ivar] - 1.) > EPS) flag_norm = 1;
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

  label_end: total = (double *) mem_free((char * ) total);
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
GEOSLIB_API int model_stabilize(Model *model, int flag_verbose, double percent)
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
    if (cova->getType() != COV_GAUSSIAN) return (0);
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

  if (model_add_cova(model, COV_NUGGET, 0, 0, 0., 0., VectorDouble(),
                     VectorDouble(), VectorDouble(total))) goto label_end;

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
GEOSLIB_API void model_covupdt(Model *model,
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

  silltot = range = (double *) NULL;
  rank = (int *) NULL;
  nvar = model->getVariableNumber();
  ncova = model->getCovaNumber();
  flag_update = flag_rescale = 0;

  /* Core allocation */

  rank = (int *) mem_alloc(sizeof(int) * ncova, 1);
  range = (double *) mem_alloc(sizeof(double) * ncova, 1);
  silltot = (double *) mem_alloc(sizeof(double) * nvar * nvar, 1);
  for (i = 0; i < nvar * nvar; i++)
    silltot[i] = 0.;

  /* Sort the basic structures by increasing range */
  rank_nugget = -1;
  for (icov = 0; icov < ncova; icov++)
  {
    cova = model->getCova(icov);
    if (cova->getType() == COV_NUGGET) rank_nugget = icov;
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
    if (cova->getType() == COV_NUGGET) continue;
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
      if (cova->getType() == COV_NUGGET) continue;
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
      if (cova->getType() == COV_NUGGET) continue;
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

  rank = (int *) mem_free((char * ) rank);
  range = (double *) mem_free((char * ) range);
  silltot = (double *) mem_free((char * ) silltot);
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
GEOSLIB_API double model_drift_evaluate(int verbose,
                                        Model *model,
                                        Db *db,
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

  model_calcul_drift(model, MEMBER_LHS, db, iech, drftab);

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
GEOSLIB_API Model *input_model(int ndim,
                               int nvar,
                               int order,
                               int flag_sill,
                               int flag_norm,
                               Model *model_in)
{
  /// TODO [Cova] : to be restored
//  int    i,flag_def,error,ncova;
//  Model *model;
//  Cova  *cova,*cova_in;
//
//  /* Initializations */
//
//  error    = 1;
//  flag_def = (model_in != (Model *) NULL);
//  model    = (Model *) NULL;
//
//  /* Core allocation */
//
//  model = model_init(ndim,nvar,0,0.,0.,VectorDouble(),VectorDouble());
//  if (model == (Model *) NULL) goto label_end;
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
//    cova_in = (Cova *) NULL;
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
GEOSLIB_API int model_dimension(Model *model)
{
  return (model->getCovaNumber());
}

/****************************************************************************/
/*!
 **  Ask the characteristics of one Cova structure
 **
 ** \return  Error returned code
 **
 ** \param[in]  model     Model structure
 ** \param[in]  icov      Rank of the Covariance structure (from 0)
 **
 ** \param[out]  cov_type      Type of the covariance
 **                            (Starting from 1)
 ** \param[out]  flag_aniso    1 for anisotropy and 0 otherwise
 ** \param[out]  param         Parameter
 ** \param[out]  sill          Array of sills (Dimension = nvar * nvar)
 ** \param[out]  aniso_rotmat  Rotation matrix (Dimension = ndim * ndim)
 ** \param[out]  aniso_ranges  Rotation ranges (Dimension = ndim)
 **
 *****************************************************************************/
GEOSLIB_API int model_extract_cova(Model *model,
                                   int icov,
                                   int *cov_type,
                                   int *flag_aniso,
                                   double *param,
                                   VectorDouble& sill,
                                   VectorDouble& aniso_rotmat,
                                   VectorDouble& aniso_ranges)
{
  CovAniso *cova;
  int ndim;

  /* Initializations */

  if (icov < 0 || icov >= model->getCovaNumber()) return (1);
  cova = model->getCova(icov);
  ndim = model->getDimensionNumber();

  /* Returning arguments */

  *cov_type = cova->getType() + 1;
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
GEOSLIB_API void model_extract_properties(Model *model, double *tape_range)
{
  ModTrans& modtrs = model->getModTrans();

  *tape_range = modtrs.getTape()->getRange();
}

/****************************************************************************/
/*!
 **  Returns the characteristics of the covariance
 **
 ** \param[in]  rank   Rank of the covariance
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
GEOSLIB_API void model_cova_characteristics(int rank,
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
  CovContext ctxt = CovContext(1, 2, 0.);
  ACovFunc* cov = CovFactory::createCovFunc((ENUM_COVS) rank, ctxt);
  (void) strcpy((char *) cov_name, cov->getCovName().c_str());
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
GEOSLIB_API int model_sample(Vario *vario,
                             Model *model,
                             int flag_norm,
                             int flag_cov)
{
  double *covtab;
  int i, idir, ndir, ipas, npas, idim, ndim, error, nvar, ivar, jvar, ijvar;
  VectorDouble d1;

  error = 1;
  ndim = vario->getDimensionNumber();
  ndir = vario->getDirectionNumber();
  nvar = model->getVariableNumber();
  covtab = (double *) NULL;
  CovCalcMode mode;
  mode.update(0, 0, MEMBER_LHS, -1, flag_norm, flag_cov);

  /* Core allocation */

  d1.resize(ndim);
  covtab = (double *) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == (double *) NULL) goto label_end;
  vario->resize(ndim, nvar);

  /* Calculate the C(0) constant term */

  model_calcul_cov(model, mode, 1, 1., VectorDouble(), covtab);
  for (int i = 0; i < nvar * nvar; i++)
    vario->setVars(i, covtab[i]);

  /* Loop on the directions */

  for (idir = 0; idir < ndir; idir++)
  {
    Dir& dir = vario->getDirs(idir);
    npas = dir.getNPas();

    /* Loop on the variogram lags */

    for (ipas = 0; ipas < npas; ipas++)
    {

      /* Loop on the variables */

      for (ivar = ijvar = 0; ivar < vario->getVariableNumber(); ivar++)
        for (jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {
          i = dir.getAddress(ivar, jvar, ipas, false, 0);
          dir.setSw(i, 1.);
          dir.setHh(i, ipas * dir.getDPas());
          for (idim = 0; idim < ndim; idim++)
            d1[idim] = dir.getHh(i) * dir.getCodir(idim);
          model_calcul_cov(model, mode, 1, 1., d1, covtab);
          dir.setGg(i, COVTAB(ivar, jvar));
        }
    }
  }

  /* Set the error returned code */

  error = 0;

  label_end: covtab = (double *) mem_free((char * ) covtab);
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
GEOSLIB_API void model_vector_multivar(Model *model,
                                       Db *db,
                                       int ivar,
                                       int iech,
                                       int flag_norm,
                                       int flag_cov,
                                       double *vector)
{
  double *covtab, *c00tab, v1, v2;
  int ndim, nvar, jech, i, skip, nech, ecr, jvar;
  VectorDouble d1;

  /* Initializations */

  covtab = c00tab = (double *) NULL;
  CovCalcMode mode;
  mode.update(0, 0, MEMBER_LHS, -1, flag_norm, flag_cov);
  if (st_check_model(model)) goto label_end;
  if (st_check_environ(model, db)) goto label_end;
  ndim = model->getDimensionNumber();
  nvar = model->getVariableNumber();
  nech = get_NECH(db);
  if (st_check_variable(nvar, ivar)) goto label_end;

  /* Core allocation */

  d1.resize(ndim, 0.);
  covtab = (double *) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == (double *) NULL) goto label_end;
  c00tab = (double *) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (c00tab == (double *) NULL) goto label_end;

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
        v1 = get_IDIM(db, iech, i);
        v2 = get_IDIM(db, jech, i);
        if (FFFF(v1) || FFFF(v2)) skip = 1;
        d1[i] = v1 - v2;
      }
      if (skip) continue;

      model_calcul_cov(model, mode, 1, 1., d1, covtab);
      vector[ecr++] = COVTAB(ivar, jvar);
    }

  /* Free memory */

  label_end: covtab = (double *) mem_free((char * ) covtab);
  c00tab = (double *) mem_free((char * ) c00tab);
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
GEOSLIB_API void model_vector(Model *model,
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

  covtab = (double *) NULL;
  CovCalcMode mode;
  mode.update(0, 0, MEMBER_LHS, -1, flag_norm, flag_cov);
  if (st_check_model(model)) goto label_end;
  if (st_check_environ(model, db1)) goto label_end;
  if (st_check_environ(model, db2)) goto label_end;
  ndim = model->getDimensionNumber();
  nvar = model->getVariableNumber();
  nech = get_NECH(db1);
  if (st_check_variable(nvar, ivar)) goto label_end;
  if (st_check_variable(nvar, jvar)) goto label_end;

  /* Core allocation */

  d1.resize(ndim, 0.);
  covtab = (double *) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == (double *) NULL) goto label_end;

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
      v1 = get_IDIM(db1, iech, i);
      v2 = get_IDIM(db2, jech, i);
      if (FFFF(v1) || FFFF(v2)) skip = 1;
      d1[i] = v1 - v2;
    }
    if (skip) continue;

    model_calcul_cov(model, mode, 1, 1., d1, covtab);
    vector[iech] = COVTAB(ivar, jvar);
  }

  /* Free memory */

  label_end: covtab = (double *) mem_free((char * ) covtab);
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
GEOSLIB_API void model_vector_nostat(Model *model,
                                     Db *db,
                                     int ivar,
                                     int jvar,
                                     int iech,
                                     double *vector)
{
  double *covtab, v1, v2, value;
  int ndim, nvar, jech, i, skip, nech;
  VectorDouble d1;
  CovCalcMode mode;

  /* Initializations */

  covtab = (double *) NULL;
  if (st_check_model(model)) goto label_end;
  if (st_check_environ(model, db)) goto label_end;
  ndim = model->getDimensionNumber();
  nvar = model->getVariableNumber();
  nech = get_NECH(db);
  if (st_check_variable(nvar, ivar)) goto label_end;
  if (st_check_variable(nvar, jvar)) goto label_end;

  /* Core allocation */

  d1.resize(ndim, 0.);
  covtab = (double *) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == (double *) NULL) goto label_end;

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
      v1 = get_IDIM(db, iech, i);
      v2 = get_IDIM(db, jech, i);
      if (FFFF(v1) || FFFF(v2)) skip = 1;
      d1[i] = v1 - v2;
    }
    if (skip) continue;

    model_calcul_cov_nostat(model, mode, 1, 1., db, iech, db, jech, d1, covtab);
    value = COVTAB(ivar, jvar);
    vector[jech] = value;
  }

  /* Free memory */

  label_end: covtab = (double *) mem_free((char * ) covtab);
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
GEOSLIB_API double model_maximum_distance(Model *model)

{
  return model->getCovAnisoList()->getMaximumDistance();
}

/****************************************************************************/
/*!
 **  For a given basic structure, convert the theoretical range (scale) into
 **  the practical range (which is the one actually stored in Geoslib)
 **
 ** \return  Range
 **
 ** \param[in]  type      Type of the basic structure
 ** \param[in]  scale     Theoretical range
 ** \param[in]  param     Third parameter
 **
 *****************************************************************************/
GEOSLIB_API double model_scale2range(int type, double scale, double param)
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
GEOSLIB_API double model_range2scale(int type, double range, double param)
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
GEOSLIB_API double model_get_field(Model *model)
{
  return (model->getField());
}

/****************************************************************************/
/*!
 **  Combine two basic models into a bivariate model (residuals model)
 **
 ** \return  Pointer to the newly created Model structure
 **
 ** \param[in]  model1      First input Model
 ** \param[in]  model2      Second input Model
 ** \param[in]  r           Correlation coefficient
 **
 ** \remarks: The drift is not copied into the new model
 **
 *****************************************************************************/
GEOSLIB_API Model *model_combine(Model *model1, Model *model2, double r)
{
  Model *model;
  CovAniso *cova;
  double field;
  int error, i, ncov;
  VectorDouble sill, mean, cova0;

  /* Initializations */

  error = 1;
  model = (Model *) NULL;
  sill.resize(4);
  mean.resize(2);
  cova0.resize(4);
  if (model1 == (Model *) NULL || model2 == (Model *) NULL)
  {
    messerr("This function requires two defined models");
    return (model);
  }
  if (model1->getVariableNumber() != 1 || model2->getVariableNumber() != 1)
  {
    messerr("This function can only combine monovariate models");
    return (model);
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
GEOSLIB_API int model_get_nonugget_cova(Model *model)

{
  CovAniso *cova;
  int nstruc, icov;

  /* Loop on the basic structures */

  nstruc = 0;
  for (icov = 0; icov < model->getCovaNumber(); icov++)
  {
    cova = model->getCova(icov);
    if (cova->getType() != COV_NUGGET) nstruc++;
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
GEOSLIB_API int model_regularize(Model *model,
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
  c00tab = covtab = (double *) NULL;
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
  nech = get_NECH(db);
  norme = nech * nech;
  vario->resize(ndim, nvar);

  /* Core allocation */

  dd.resize(ndim, 0.);
  c00tab = (double *) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (c00tab == (double *) NULL) goto label_end;
  covtab = (double *) mem_alloc(sizeof(double) * nvar * nvar, 0);
  if (covtab == (double *) NULL) goto label_end;

  /* Calculate the Cvv (for a zero-shift) */

  model_covtab_init(1, model, c00tab);
  for (iech = 0; iech < nech; iech++)
    for (jech = 0; jech < nech; jech++)
    {
      for (idim = 0; idim < ndim; idim++)
      {
        v1 = get_IDIM(db, iech, idim);
        v2 = get_IDIM(db, jech, idim);
        dd[idim] = v1 - v2;
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
    Dir& dir = vario->getDirs(idir);

    /* Loop on the number of lags */

    for (ipas = 0; ipas < dir.getNPas(); ipas++)
    {
      model_covtab_init(1, model, covtab);
      dist = ipas * dir.getDPas();

      for (iech = 0; iech < nech; iech++)
        for (jech = 0; jech < nech; jech++)
        {
          for (idim = 0; idim < ndim; idim++)
          {
            v1 = get_IDIM(db, iech, idim);
            v2 = get_IDIM(db, jech, idim) + dist * dir.getCodir(idim);
            dd[idim] = v1 - v2;
          }
          model_calcul_cov(model, mode, 0, 1, dd, covtab);
        }
      st_covtab_rescale(nvar, norme, covtab);

      for (ivar = 0; ivar < nvar; ivar++)
        for (jvar = 0; jvar <= ivar; jvar++)
        {
          iad = dir.getAddress(ivar, jvar, ipas, false, 0);
          dir.setGg(iad, C00TAB(ivar,jvar)- COVTAB(ivar,jvar));
          dir.setHh(iad, dist);
          dir.setSw(iad, 1);
        }
    }
  }

  /* Set the error return code */

  error = 0;

  label_end: c00tab = (double *) mem_free((char * ) c00tab);
  covtab = (double *) mem_free((char * ) covtab);
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
GEOSLIB_API int model_covmat_inchol(int verbose,
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
  nech = get_NECH(db);
  pvec = (int *) NULL;
  diag = crit = G = Gmatrix = (double *) NULL;
  flag_incr = (center != (double *) NULL);

  if (npivot_max <= 0) npivot_max = nech;
  npivot_max = MIN(npivot_max, nech);
  d1.resize(db->getNDim());
  diag = (double *) mem_alloc(sizeof(double) * nech, 0);
  if (diag == (double *) NULL) goto label_end;
  crit = (double *) mem_alloc(sizeof(double) * (1 + nech), 0);
  if (crit == (double *) NULL) goto label_end;
  pvec = (int *) mem_alloc(sizeof(int) * nech, 0);
  if (pvec == (int *) NULL) goto label_end;
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
        d1[idim] = get_IDIM(db, pvec[i], idim) - center[idim];
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
    G = (double *) mem_realloc((char * ) G,
                               (npivot + 1) * nech * sizeof(double), 0);
    if (G == (double *) NULL) goto label_end;
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
          d1[idim] = get_IDIM(db, pvec[npivot], idim) - center[idim];
        model_calcul_cov(model, mode, 1, 1., d1, &covar2);

        for (int idim = 0; idim < 3; idim++)
          d1[idim] = get_IDIM(db, pvec[i], idim) - center[idim];
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
          d1[idim] = get_IDIM(db, pvec[i], idim) - center[idim];
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
    G = (double *) mem_realloc((char * ) G,
                               (npivot + 1) * nech * sizeof(double), 0);
    if (G == (double *) NULL) goto label_end;
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
  Gmatrix = (double *) mem_alloc(npivot * nech * sizeof(double), 0);
  if (Gmatrix == (double *) NULL) goto label_end;
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

  label_end: diag = (double *) mem_free((char * ) diag);
  crit = (double *) mem_free((char * ) crit);
  G = (double *) mem_free((char * ) G);
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
GEOSLIB_API int model_maximum_order(Model *model)

{
  int order, max_order;

  if (model == (Model *) NULL) return (-1);

  max_order = 0;
  for (int il = 0; il < model->getDriftNumber(); il++)
  {
    ADriftElem* drft = model->getDrift(il);
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
 ** \param[in]  type0      Drift function to be found
 **
 *****************************************************************************/
GEOSLIB_API int model_is_drift_defined(Model *model, int type0)
{
  if (model == (Model *) NULL) return (0);
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
GEOSLIB_API double model_calcul_stdev(Model *model,
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
    d1[idim] = get_IDIM(db1, iech1, idim) - get_IDIM(db2, iech2, idim);
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


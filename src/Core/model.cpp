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

#include "Drifts/DriftFactory.hpp"
#include "Drifts/DriftList.hpp"
#include "Drifts/ADrift.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Covariances/CovAnisoList.hpp"
#include "Covariances/CovLMGradient.hpp"
#include "Covariances/CovLMCTapering.hpp"
#include "Covariances/CovLMCConvolution.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovGradientNumerical.hpp"
#include "Model/CovInternal.hpp"
#include "Model/Model.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Anamorphosis/AnamDiscreteIR.hpp"
#include "Variogram/Vario.hpp"
#include "Space/SpaceRN.hpp"
#include "Basic/String.hpp"
#include "Db/Db.hpp"
#include "Basic/Memory.hpp"

#include <math.h>

/*! \cond */
#define AD(ivar,jvar)          (ivar) + nvar * (jvar)
#define AIC(icov,ivar,jvar)     aic[(icov)*nvar*nvar + AD(ivar,jvar)]
#define VALPRO(ivar)            valpro[(ivar)]
#define VECPRO(ivar,jvar)       vecpro[AD(ivar,jvar)]
#define CC(ivar,jvar)           cc[AD(ivar,jvar)]
#define DISC1(i,idim)          (koption->disc1[(idim) * koption->ntot + (i)])
#define DISC2(i,idim)          (koption->disc2[(idim) * koption->ntot + (i)])
#define G(i,j)                 (G[(i) * nech + j])
#define Gmatrix(i,j)           (Gmatrix[(j) * nech + i])
/*! \endcond */

int NDIM_LOCAL = 0;
VectorDouble X1_LOCAL = VectorDouble();
VectorDouble X2_LOCAL = VectorDouble();

/****************************************************************************/
/*!
 **  Duplicates a Model from another Model for Gradients
 **
 ** \return  The modified Model structure
 **
 ** \param[in]  model       Input Model
 ** \param[in]  ball_radius Radius for Gradient calculation
 **
 *****************************************************************************/
Model* model_duplicate_for_gradient(const Model *model, double ball_radius)

{
  Model *new_model;
  const CovAniso *cova;
  int new_nvar, nfact;
  double sill;

  // Preliminary checks

  new_model = nullptr;
  int nvar  = model->getNVar();
  int ndim  = model->getDimensionNumber();
  int ncova = model->getCovaNumber();

  // Create the new model (linked drift functions)

  if (nvar != 1 || ndim != 2)
  {
    messerr("This procedure is limited to a single variable in 2-D");
    return new_model;
  }

  new_nvar = 3;
  nfact = 6;
  CovContext ctxt(model->getContext());
  ctxt.setNVar(new_nvar);
  new_model = new Model(ctxt);
  if (new_model == nullptr) return new_model;

  // **************************************
  // Create the basic covariance structures
  // **************************************

  CovAnisoList* covs = new CovLMGradient(ctxt);

  int lec = 0;
  for (int icov = 0; icov < ncova; icov++)
  {
    cova = model->getCova(icov);
    sill = model->getSill(icov, 0, 0);
    for (int ifact = 0; ifact < nfact; ifact++, lec++)
    {
      CovAniso* covnew = nullptr;
      covnew = new CovGradientNumerical(cova->getType(),ball_radius,ctxt);
      covnew->setParam(cova->getParam());
      if (cova->getFlagAniso())
      {
        covnew->setRanges(cova->getRanges());
        if (cova->getFlagRotation())
          covnew->setAnisoRotation(cova->getAnisoRotation());
      }
      else
        covnew->setRangeIsotropic(cova->getRange());

      /* Modify the Sill */;

      covnew->initSill(0.);
      if (ifact == 0)
      {
        covnew->setSill(0, 0, sill);
      }
      else if (ifact == 1)
      {
        covnew->setSill(0, 1, -sill);
        covnew->setSill(1, 0, sill);
      }
      else if (ifact == 2)
      {
        covnew->setSill(1, 1, sill);
      }
      else if (ifact == 3)
      {
        covnew->setSill(0, 2, -sill);
        covnew->setSill(2, 0, sill);
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
      covs->addCov(covnew);
      delete covnew;
    }
  }
  new_model->setCovAnisoList(covs);
  delete covs;

  // *********************************
  // Create the basic drift structures
  // *********************************

  DriftList* drifts = DriftFactory::createDriftListForGradients(model->getDriftList(), ctxt);
  new_model->setDriftList(drifts);
  delete drifts;
  return (new_model);
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
                   const double *c0,
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
  nvar = model->getNVar();
  ncova = model->getCovaNumber();
  flag_update = 0;

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
          model->setSill(icov, ivar, jvar, 0.);
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

  mem_free((char* ) rank);
  mem_free((char* ) range);
  mem_free((char* ) silltot);
  *flag_nugget = flag_update && (rank_nugget < 0);
  if (flag_verbose && (*flag_nugget))
  {
    message("A Nugget Effect component is added so as to match the experimental variance\n");
  }
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
  auto space = SpaceRN::create(1); // Use 1-D in order to retrieve all covariances
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
  delete cov;
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
  if (model1 != nullptr && model1->getNVar() != 1)
  {
    messerr("This function can only combine monovariate models");
    return nullptr;
  }
  if (model2 != nullptr && model2->getNVar() != 1)
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

  VectorDouble mean(2);
  VectorDouble cova0(4);
  MatrixSquareSymmetric sill(2);
  mean[0] = model1->getContext().getMean(0);
  mean[1] = model2->getContext().getMean(0);
  cova0[0] = 1.;
  cova0[1] = r;
  cova0[2] = r;
  cova0[3] = 1.;

  // Creating the context
  CovContext ctxt = CovContext(2, model1->getDimensionNumber(), mean, cova0);

  // Creating the new Model
  model = new Model(ctxt);

  /* Add the covariance of the first Model */

  for (int i = 0; i < model1->getCovaNumber(); i++)
  {
    const CovAniso* cova = model1->getCova(i);
    sill.setValue(0, 0, cova->getSill(0, 0));
    sill.setValue(1, 0, r * cova->getSill(0, 0));
    sill.setValue(1, 1, r * r * cova->getSill(0, 0));
    model->addCovFromParam(cova->getType(), cova->getRange(), 0., cova->getParam(),
                           cova->getRanges(), sill, cova->getAnisoAngles());
  }

  /* Add the covariance of the second Model */

  for (int i = 0; i < model2->getCovaNumber(); i++)
  {
    const CovAniso* cova = model2->getCova(i);
    sill.setValue(0,0, 0.);
    sill.setValue(0,1, 0.);
    sill.setValue(1,1, (1. - r * r) * cova->getSill(0, 0));
    model->addCovFromParam(cova->getType(), cova->getRange(), 0., cova->getParam(),
                           cova->getRanges(), sill, cova->getAnisoAngles());
  }
  return model;
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
 ** \param[in]  center     Optional Centering point (for increments)
 ** \param[in]  flag_sort  Reordering flag (see remarks)
 ** \param[in]  mode       CovCalcMode structure
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
                        const int *ranks1,
                        const double *center,
                        int flag_sort,
                        int *npivot_arg,
                        int **Pret,
                        double **Gret,
                        const CovCalcMode*  mode)
{
  int *pvec, i, j, npivot, jstar, nech, error, flag_incr;
  double *G, *Gmatrix, *diag, *crit, g, residual, maxdiag, tol, b, c00;
  VectorDouble d1;

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
  c00 = model->evaluateOneGeneric(nullptr, VectorDouble(), 1., mode);
  for (i = 0; i < nech; i++)
    pvec[i] = i;

  residual = 0.;
  for (i = 0; i < nech; i++)
  {
    if (flag_incr)
    {
      double covar2 = 0.;

      for (int idim = 0; idim < 3; idim++)
        d1[idim] = db->getCoordinate(pvec[i], idim) - center[idim];
      covar2 = model->evaluateOneGeneric(nullptr, d1, 1., mode);
      diag[i] = 2. * (c00 - covar2);
    }
    else
    {
      diag[i] = c00;
    }
    residual += diag[i];
  }
  tol = (!FFFF(eta)) ? eta * residual : 0.;
  npivot = 0;

  // Main loop

  while ((residual > tol) && (npivot < npivot_max))
  {
    // Initialize and add a new zeros column to matrix G[]
    G = (double*) mem_realloc((char* ) G, (npivot + 1) * nech * sizeof(double), 0);
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
        double covar1 = 0.;
        double covar2 = 0.;
        double covar3 = 0.;

        (void) distance_intra(db, pvec[i], pvec[npivot], d1.data());
        covar1 = model->evaluateOneGeneric(nullptr, d1, 1., mode);

        for (int idim = 0; idim < 3; idim++)
          d1[idim] = db->getCoordinate(pvec[npivot], idim) - center[idim];
        covar2 = model->evaluateOneGeneric(nullptr, d1, 1., mode);

        for (int idim = 0; idim < 3; idim++)
          d1[idim] = db->getCoordinate(pvec[i], idim) - center[idim];
        covar3 = model->evaluateOneGeneric(nullptr, d1, 1., mode);

        G(npivot,i) = covar1 - covar2 - covar3 + c00;
      }
      else
      {
        // Calculate the covariance column C(:, npivot)
        (void) distance_intra(db, pvec[i], pvec[npivot], d1.data());
        G(npivot, i) = model->evaluateOneGeneric(nullptr, d1, 1., mode);
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
        double covar2 = 0.;

        for (int idim = 0; idim < 3; idim++)
          d1[idim] = db->getCoordinate(pvec[i], idim) - center[idim];
        covar2 = model->evaluateOneGeneric(nullptr, d1, 1., mode);

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
    G = (double*) mem_realloc((char* ) G, (npivot + 1) * nech * sizeof(double), 0);
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

  label_end:
  mem_free((char* ) diag);
  mem_free((char* ) crit);
  mem_free((char* ) G);
  return (error);
}


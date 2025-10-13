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
#include "Model/Model.hpp"
#include "Anamorphosis/AnamDiscreteIR.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Basic/String.hpp"
#include "Basic/Utilities.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovAnisoList.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovGradientGeneric.hpp"
#include "Covariances/CovLMCConvolution.hpp"
#include "Covariances/CovLMCTapering.hpp"
#include "Db/Db.hpp"
#include "Drifts/ADrift.hpp"
#include "Drifts/DriftList.hpp"
#include "Model/CovInternal.hpp"
#include "Space/SpaceRN.hpp"
#include "Variogram/Vario.hpp"
#include "geoslib_old_f.h"

#include <cmath>

/*! \cond */
#define AD(ivar, jvar)        (ivar) + nvar*(jvar)
#define AIC(icov, ivar, jvar) aic[(icov) * nvar * nvar + AD(ivar, jvar)]
#define VALPRO(ivar)          valpro[(ivar)]
#define VECPRO(ivar, jvar)    vecpro[AD(ivar, jvar)]
#define CC(ivar, jvar)        cc[AD(ivar, jvar)]
#define DISC1(i, idim)        (koption->disc1[(idim) * koption->ntot + (i)])
#define DISC2(i, idim)        (koption->disc2[(idim) * koption->ntot + (i)])
#define G(i, j)               (G[(i) * nech + j])
#define Gmatrix(i, j)         (Gmatrix[(j) * nech + i])
/*! \endcond */

namespace gstlrn
{
Id NDIM_LOCAL = 0;
VectorDouble X1_LOCAL;
VectorDouble X2_LOCAL;

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
void model_cova_characteristics(const ECov& type,
                                char cov_name[STRING_LENGTH],
                                Id* flag_range,
                                Id* flag_param,
                                Id* min_order,
                                Id* max_ndim,
                                Id* flag_int_1d,
                                Id* flag_int_2d,
                                Id* flag_aniso,
                                Id* flag_rotation,
                                double* scale,
                                double* parmax)
{
  auto space = SpaceRN::create(1); // Use 1-D in order to retrieve all covariances
  CovContext ctxt(1, 1);
  ACovFunc* cov = CovFactory::createCovFunc(type, ctxt);
  (void)gslStrcpy(static_cast<char*>(cov_name), STRING_LENGTH, cov->getCovName().c_str());
  *flag_range    = cov->hasRange();
  *flag_param    = cov->hasParam();
  *min_order     = cov->getMinOrder();
  *max_ndim      = static_cast<Id>(cov->getMaxNDim());
  *flag_int_1d   = cov->hasInt1D();
  *flag_int_2d   = cov->hasInt2D();
  *flag_aniso    = (((*flag_range) != 0) && (*max_ndim < 0 || *max_ndim > 1));
  *flag_rotation = ((*flag_aniso) && (*max_ndim < 0 || *max_ndim > 1));
  *scale         = cov->getScadef();
  *parmax        = cov->getParMax();
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
Model* model_combine(const Model* model1, const Model* model2, double r)
{
  Model* model;

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
  if (model1->getNDim() != model2->getNDim())
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
  MatrixSymmetric sill(2);
  mean[0]  = model1->getMean(0);
  mean[1]  = model2->getMean(0);
  cova0[0] = 1.;
  cova0[1] = r;
  cova0[2] = r;
  cova0[3] = 1.;

  // Creating the context
  CovContext ctxt(2, static_cast<Id>(model1->getNDim()), cova0);

  // Creating the new Model
  model = new Model(ctxt);
  model->setMeans(mean);
  /* Add the covariance of the first Model */

  for (Id i = 0; i < model1->getNCov(); i++)
  {
    const CovAniso* cova = model1->getCovAniso(i);
    sill.setValue(0, 0, cova->getSill(0, 0));
    sill.setValue(1, 0, r * cova->getSill(0, 0));
    sill.setValue(1, 1, r * r * cova->getSill(0, 0));
    model->addCovFromParam(cova->getType(), cova->getRangeIso(), 0., cova->getParam(),
                           cova->getRanges(), sill, cova->getAnisoAngles());
  }

  /* Add the covariance of the second Model */

  for (Id i = 0; i < model2->getNCov(); i++)
  {
    const CovAniso* cova = model2->getCovAniso(i);
    sill.setValue(0, 0, 0.);
    sill.setValue(0, 1, 0.);
    sill.setValue(1, 1, (1. - r * r) * cova->getSill(0, 0));
    model->addCovFromParam(cova->getType(), cova->getRangeIso(), 0., cova->getParam(),
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
 ** \param[out] pvec       Array of indices of the retained samples (from 1)
 **                        Dimension: nech
 ** \param[out] Gmatrix    Rectangular matrix
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
Id model_covmat_inchol(Id verbose,
                       Db* db,
                       Model* model,
                       double eta,
                       Id npivot_max,
                       Id nsize1,
                       const Id* ranks1,
                       const double* center,
                       Id flag_sort,
                       Id* npivot_arg,
                       VectorInt& pvec,
                       VectorDouble& Gmatrix,
                       const CovCalcMode* mode)
{
  Id i, j, npivot, jstar, nech, flag_incr;
  double g, residual, maxdiag, tol, b, c00;
  VectorDouble d1;
  VectorDouble diag;
  VectorDouble crit;
  VectorDouble G;

  nech      = db->getNSample();
  flag_incr = (center != nullptr);

  if (npivot_max <= 0) npivot_max = nech;
  npivot_max = MIN(npivot_max, nech);
  d1.resize(db->getNDim());
  diag.resize(nech);
  crit.resize(1 + nech);
  pvec.resize(nech);
  for (i = 0; i < nech; i++) pvec[i] = i;
  c00 = model->evaluateOneGeneric(nullptr, VectorDouble(), 1., mode);

  residual = 0.;
  for (i = 0; i < nech; i++)
  {
    if (flag_incr)
    {
      double covar2 = 0.;

      for (Id idim = 0; idim < 3; idim++)
        d1[idim] = db->getCoordinate(pvec[i], idim) - center[idim];
      covar2  = model->evaluateOneGeneric(nullptr, d1, 1., mode);
      diag[i] = 2. * (c00 - covar2);
    }
    else
    {
      diag[i] = c00;
    }
    residual += diag[i];
  }
  tol    = (!FFFF(eta)) ? eta * residual : 0.;
  npivot = 0;

  // Main loop

  while ((residual > tol) && (npivot < npivot_max))
  {
    // Initialize and add a new zeros column to matrix G[]
    G.resize((npivot + 1) * nech);
    for (i = 0; i < nech; i++)
      G(npivot, i) = 0.;

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
          jstar   = i;
          maxdiag = diag[i];
        }
      }
    }

    // Update permutation pvec (not necessary if jstar = npivot)
    if (jstar != npivot)
    {
      i            = pvec[jstar];
      pvec[jstar]  = pvec[npivot];
      pvec[npivot] = i;
      diag[npivot] = diag[jstar];

      // Update rows elements on G
      for (j = 0; j <= npivot; j++)
      {
        g            = G(j, jstar);
        G(j, jstar)  = G(j, npivot);
        G(j, npivot) = g;
      }
    }

    // Calculate the diagonal element of G
    G(npivot, npivot) = sqrt(diag[jstar]);

    // Calculate the new column of G
    for (i = npivot + 1; i < nech; i++)
    {
      if (flag_incr)
      {
        double covar1 = 0.;
        double covar2 = 0.;
        double covar3 = 0.;

        (void)distance_intra(db, pvec[i], pvec[npivot], d1.data());
        covar1 = model->evaluateOneGeneric(nullptr, d1, 1., mode);

        for (Id idim = 0; idim < 3; idim++)
          d1[idim] = db->getCoordinate(pvec[npivot], idim) - center[idim];
        covar2 = model->evaluateOneGeneric(nullptr, d1, 1., mode);

        for (Id idim = 0; idim < 3; idim++)
          d1[idim] = db->getCoordinate(pvec[i], idim) - center[idim];
        covar3 = model->evaluateOneGeneric(nullptr, d1, 1., mode);

        G(npivot, i) = covar1 - covar2 - covar3 + c00;
      }
      else
      {
        // Calculate the covariance column C(:, npivot)
        (void)distance_intra(db, pvec[i], pvec[npivot], d1.data());
        G(npivot, i) = model->evaluateOneGeneric(nullptr, d1, 1., mode);
      }
    }
    if (npivot != 0)
    {
      for (i = npivot + 1; i < nech; i++)
        for (j = 0; j < npivot; j++)
          G(npivot, i) -= G(j, i) * G(j, npivot);
    }
    for (i = npivot + 1; i < nech; i++)
      G(npivot, i) /= G(npivot, npivot);

    // Updates diagonal elements
    for (i = npivot + 1; i < nech; i++)
    {
      if (flag_incr)
      {
        double covar2 = 0.;

        for (Id idim = 0; idim < 3; idim++)
          d1[idim] = db->getCoordinate(pvec[i], idim) - center[idim];
        covar2 = model->evaluateOneGeneric(nullptr, d1, 1., mode);

        b = 2. * (c00 - covar2);
      }
      else
      {
        b = c00;
      }
      for (j = 0; j <= npivot; j++)
        b -= G(j, i) * G(j, i);
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
    G.resize((npivot + 1) * nech);
    for (i = 0; i < nech; i++)
      G(npivot, i) = 0.;
    G(npivot, npivot) = sqrt(diag[npivot]);
    crit[npivot]      = 0.;
    npivot++;
  }

  // Return arguments
  *npivot_arg = npivot;

  // Normalize the criterion
  for (i = 0; i < npivot; i++)
    crit[i] /= static_cast<double>(nech);

  // Reorder the output G matrix
  Gmatrix.resize(npivot * nech);
  for (j = 0; j < npivot; j++)
    for (i = 0; i < nech; i++)
    {
      if (flag_sort)
        Gmatrix(pvec[i], j) = G(j, i);
      else
        Gmatrix(i, j) = G(j, i);
    }

  // Renumber starting from 1
  for (i = 0; i < nech; i++)
    pvec[i]++;

  // Printout of the order of the retained samples
  if (verbose)
  {
    message("Number of pivots = %d\n", npivot);
    print_imatrix("Order", 0, 1, 1, npivot, NULL, pvec.data());
    print_matrix("Criterion", 0, 1, 1, npivot, NULL, crit.data());
  }

  return 0;
}
} // namespace gstlrn

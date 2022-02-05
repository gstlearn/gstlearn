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
#include "geoslib_old_f.h"
#include "Basic/Utilities.hpp"
#include "Db/Db.hpp"

/*! \cond */
#define F_IATT(iatt)  ((iatt) < 0 || (iatt) >= db->getNMax())
#define MM(idim,jdim) (mm[(idim) * ndim + jdim])
/*! \endcond */

/****************************************************************************/
/*!
 **  Load the valid data and calculate the its weight
 **
 ** \return  Error returned code
 **
 ** \param[in]  db     Db structure
 ** \param[in]  flag_z 1 if the variable is defined
 ** \param[in]  flag_w 1 if the weight is defined
 ** \param[in]  iech   Rank of the sample
 ** \param[in]  iatt   Rank of the attribute (used if flag_z)
 **
 ** \param[out]  coor   Array of coordinates
 ** \param[out]  wz     Weighted value
 **
 *****************************************************************************/
static int st_cgi_data(Db *db,
                       int flag_z,
                       int flag_w,
                       int iech,
                       int iatt,
                       double *coor,
                       double *wz)
{
  double value, weight;
  int idim;

  /* Check if the variable is defined */

  value = 1.;
  if (flag_z)
  {
    value = db->getArray(iech, iatt);
    if (FFFF(value)) return (1);
    if (value < 0.)
    {
      messerr("The variable cannot be negative (Sample %d = %lf)", iech + 1,
              value);
      messerr("Procedure is interrupted");
      return (1);
    }
  }

  /* Check if the weight is defined */

  weight = 1.;
  if (flag_w)
  {
    weight = db->getWeight(iech);
    if (FFFF(weight)) return (1);
    if (weight < 0.)
    {
      messerr("The weight cannot be negative (Sample %d = %lf)", iech + 1,
              weight);
      messerr("Procedure is interrupted");
      return (1);
    }
  }

  /* Check if the sample has defined coordinates */

  db_sample_load(db, ELoc::X, iech, coor);
  for (idim = 0; idim < db->getNDim(); idim++)
    if (FFFF(coor[idim])) return (1);

  /* Returning argument */

  *wz = value * weight;

  return (0);
}

/****************************************************************************/
/*!
 **  Calculate the Center of Gravity
 **
 ** \return  Error returned code
 **
 ** \param[in]  db   Db structure
 ** \param[in]  iatt Rank of the optional attribute
 **
 ** \param[out] center   Array for the center of gravity
 ** \param[out] mvalue   Array of eigen values (normalized)
 ** \param[out] mvector  Array of eigen vectors
 ** \param[out] inertia  Value of the inertia
 ** \param[out] wztot    Sum of weights
 **
 ** \remark The array mvalue must be dimensionned to ndim
 ** \remark The array mvector must be dimensionned to ndim * ndim
 **
 *****************************************************************************/
int cgi(Db *db,
                        int iatt,
                        double *center,
                        double *mvalue,
                        double *mvector,
                        double *inertia,
                        double *wztot)
{
  int iech, nech, idim, jdim, ndim, error, flag_z, flag_w;
  double *coor, *mm, wz, sum;

  /* Initializations */

  error = 1;
  nech = db->getSampleNumber();
  ndim = db->getNDim();
  coor = mm = nullptr;
  flag_z = db->isUIDDefined(iatt);
  flag_w = db->hasWeight();

  /* Core allocation */

  coor = db_sample_alloc(db, ELoc::X);
  if (coor == nullptr) goto label_end;
  mm = (double*) mem_alloc(sizeof(double) * ndim * ndim, 0);
  if (mm == nullptr) goto label_end;

  /* Initialize the arrays to zero */

  for (idim = 0; idim < ndim; idim++)
  {
    coor[idim] = center[idim] = 0.;
    for (jdim = 0; jdim < ndim; jdim++)
      MM(idim,jdim) = 0.;
  }

  /* Calculate the Center of Gravity */

  (*wztot) = 0.;
  for (iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    if (st_cgi_data(db, flag_z, flag_w, iech, iatt, coor, &wz)) continue;
    for (idim = 0; idim < ndim; idim++)
      center[idim] += wz * coor[idim];
    (*wztot) += wz;
  }
  if ((*wztot) <= 0.)
  {
    messerr("The sum of the weights must be positive : %lf", (*wztot));
    goto label_end;
  }
  for (idim = 0; idim < ndim; idim++)
    center[idim] /= (*wztot);

  /* Calculate the inertia and the weighted PCA */

  (*inertia) = 0.;
  for (iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    if (st_cgi_data(db, flag_z, flag_w, iech, iatt, coor, &wz)) continue;
    for (idim = 0; idim < ndim; idim++)
      coor[idim] -= center[idim];
    for (idim = 0; idim < ndim; idim++)
    {
      (*inertia) += wz * coor[idim] * coor[idim];
      for (jdim = 0; jdim < ndim; jdim++)
        MM(idim,jdim) += wz * coor[idim] * coor[jdim];
    }
  }

  /* Normation */

  (*inertia) /= (*wztot);
  for (idim = 0; idim < ndim; idim++)
    for (jdim = 0; jdim < ndim; jdim++)
      MM(idim,jdim) /= (*wztot);

  /* Calculate the eigen values and vectors */

  if (matrix_eigen(mm, ndim, mvalue, mvector)) goto label_end;

  /* Normation */

  sum = 0.;
  for (idim = 0; idim < ndim; idim++)
    sum += mvalue[idim];

  /* Set the error return code */

  error = 0;

  label_end: coor = db_sample_free(coor);
  mm = (double*) mem_free((char* ) mm);
  return (error);
}

/****************************************************************************/
/*!
 **  Calculate several Spatial indices
 **
 ** \return  Error returned code
 **
 ** \param[in]  db   Db structure
 **
 ** \param[out] totab    Total Abundance
 ** \param[out] parea    Positive area
 ** \param[out] eqarea   Equivalent area
 **
 ** \remark This functions have been developped in the scope of the UE
 ** \remark program Fisboat, DG-Fish, STREP #502572
 **
 *****************************************************************************/
int spatial(Db *db,
                            double *totab,
                            double *parea,
                            double *eqarea)
{
  double w, z, sum, top, bot, maille;
  int iech;

  /* Initializations */

  top = bot = sum = 0.;
  maille = (is_grid(db)) ? db_grid_maille(db) :
                           1.;

  /* Loop on the samples */

  for (iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    z = db->getVariable(iech, 0);
    if (FFFF(z)) continue;
    w = db->getWeight(iech);
    if (FFFF(w)) continue;
    if (z > 0) sum += w;
    top += w * z;
    bot += w * z * z;
  }
  top *= maille;
  bot *= maille;

  /* Returning arguments */

  *totab = top;
  *parea = sum;
  *eqarea = (bot == 0.) ? TEST :
                          top * top / bot;
  return (0);
}


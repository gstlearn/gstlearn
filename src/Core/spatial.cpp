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
 ** \param[in]  flag_z True if the variable is defined
 ** \param[in]  flag_w True if the weight is defined
 ** \param[in]  iech   Rank of the sample
 ** \param[in]  iatt   Rank of the attribute (used if flag_z)
 **
 ** \param[out]  coor   Array of coordinates
 ** \param[out]  wz     Weighted value
 **
 *****************************************************************************/
static int st_cgi_data(Db *db,
                       bool flag_z,
                       bool flag_w,
                       int iech,
                       int iatt,
                       VectorDouble& coor,
                       double *wz)
{
  /* Check if the variable is defined */

  double value = 1.;
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

  double weight = 1.;
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

  db_sample_load(db, ELoc::X, iech, coor.data());
  for (int idim = 0; idim < db->getNDim(); idim++)
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
 ** \remark 'mvalue' must be dimensioned to ndim
 ** \remark 'mvector' must be dimensioned to ndim * ndim
 **
 *****************************************************************************/
int cgi(Db *db,
        int iatt,
        VectorDouble& center,
        VectorDouble& mvalue,
        MatrixRectangular& mvector,
        double *inertia,
        double *wztot)
{
  double wz;
  int nech = db->getSampleNumber();
  int ndim = db->getNDim();
  bool flag_z = db->isUIDDefined(iatt);
  bool flag_w = db->hasLocVariable(ELoc::W);
  VectorDouble coor(ndim, 0.);
  MatrixSquareSymmetric mm(ndim);
  center.fill(0.);

  /* Calculate the Center of Gravity */

  (*wztot) = 0.;
  for (int iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    if (st_cgi_data(db, flag_z, flag_w, iech, iatt, coor, &wz)) continue;
    for (int idim = 0; idim < ndim; idim++)
      center[idim] += wz * coor[idim];
    (*wztot) += wz;
  }
  if ((*wztot) <= 0.)
  {
    messerr("The sum of the weights must be positive : %lf", (*wztot));
    return 1;
  }
  for (int idim = 0; idim < ndim; idim++)
    center[idim] /= (*wztot);

  /* Calculate the inertia and the weighted PCA */

  (*inertia) = 0.;
  for (int iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    if (st_cgi_data(db, flag_z, flag_w, iech, iatt, coor, &wz)) continue;
    for (int idim = 0; idim < ndim; idim++)
      coor[idim] -= center[idim];
    for (int idim = 0; idim < ndim; idim++)
    {
      (*inertia) += wz * coor[idim] * coor[idim];
      for (int jdim = 0; jdim <= idim; jdim++)
        mm.updValue(idim,jdim, EOperator::ADD,wz * coor[idim] * coor[jdim]);
    }
  }

  /* Normation */

  (*inertia) /= (*wztot);
  for (int idim = 0; idim < ndim; idim++)
    for (int jdim = 0; jdim <= idim; jdim++)
      mm.updValue(idim,jdim, EOperator::DIVOPP,*wztot);

  /* Calculate the eigen values and vectors */

  if (mm.computeEigen()) return 1;
  mvalue = mm.getEigenValues();
  mvector = *mm.getEigenVectors();

  /* Normation */

  double sum = 0.;
  for (int idim = 0; idim < ndim; idim++)
    sum += mvalue[idim];

  return 0;
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
int spatial(Db *db, double *totab, double *parea, double *eqarea)
{
  double top = 0.;
  double bot = 0.;
  double sum = 0.;
  double maille = (db->isGrid()) ? db_grid_maille(db) : 1.;

  /* Loop on the samples */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    double z = db->getLocVariable(ELoc::Z,iech, 0);
    if (FFFF(z)) continue;
    double w = db->getWeight(iech);
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
  *eqarea = (bot == 0.) ? TEST : top * top / bot;
  return (0);
}


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
#include <Spatial/Indices.hpp>
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"

/****************************************************************************/
/*!
 **  Load the valid data and calculate its weight
 **
 ** \return  Error returned code
 **
 ** \param[in]  db     Db structure
 ** \param[in]  flag_w True if the weight is defined
 ** \param[in]  iech   Rank of the sample
 ** \param[in]  name   Name of the attribute (used if flag_z)
 **
 ** \param[out]  coor   Array of coordinates
 ** \param[out]  wz     Weighted value
 **
 *****************************************************************************/
static int st_cgi_data(Db *db,
                       bool flag_w,
                       int iech,
                       const String& name,
                       VectorDouble& coor,
                       double *wz)
{
  // Check if the sample is masked off

  if (! db->isActive(iech)) return 1;

  /* Check if the variable is defined */

  double value = 1.;
  if (! name.empty())
  {
    value = db->getValue(name, iech);
    if (FFFF(value)) return (1);
    if (value < 0.)
    {
      messerr("The variable cannot be negative (Sample %d = %lf)", iech + 1, value);
      messerr("Procedure is interrupted");
      return 1;
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
      messerr("The weight cannot be negative (Sample %d = %lf)", iech + 1, weight);
      messerr("Procedure is interrupted");
      return 1;
    }
  }

  /* Check if the sample has defined coordinates */

  db->getCoordinatesPerSampleInPlace(iech, coor);
  for (int idim = 0; idim < db->getNDim(); idim++)
    if (FFFF(coor[idim])) return 1;

  /* Returning argument */

  *wz = value * weight;

  return 0;
}

/****************************************************************************/
/*!
 **  Calculate the Center of Gravity
 **
 ** \return  Error returned code
 **
 ** \param[in]  db   Db structure
 ** \param[in]  name Name of the optional attribute
 **
 *****************************************************************************/
cgiOutput cgi(Db *db, const String &name)
{
  // Initializations
  double wz;
  int nech = db->getSampleNumber();
  int ndim = db->getNDim();
  bool flag_w = db->hasLocVariable(ELoc::W);

  // Prepare the output structure
  cgiOutput res;

  /* Calculate the Center of Gravity */

  res.wztot = 0.;
  res.nvalid = 0;
  res.center.resize(ndim, 0.);
  VectorDouble coor(ndim, 0.);
  for (int iech = 0; iech < nech; iech++)
  {
    if (st_cgi_data(db, flag_w, iech, name, coor, &wz)) continue;
    for (int idim = 0; idim < ndim; idim++)
      res.center[idim] += wz * coor[idim];
    res.wztot += wz;
    res.nvalid++;
  }
  if (res.wztot <= 0.)
  {
    messerr("The sum of the weights must be positive : %lf", res.wztot);
    return res;
  }
  for (int idim = 0; idim < ndim; idim++)
    res.center[idim] /= res.wztot;

  /* Calculate the inertia and the weighted PCA */

  res.inertia = 0.;
  MatrixSquareSymmetric mm(ndim);
  for (int iech = 0; iech < nech; iech++)
  {
    if (st_cgi_data(db, flag_w, iech, name, coor, &wz)) continue;
    for (int idim = 0; idim < ndim; idim++)
      coor[idim] -= res.center[idim];
    for (int idim = 0; idim < ndim; idim++)
    {
      res.inertia += wz * coor[idim] * coor[idim];
      for (int jdim = 0; jdim <= idim; jdim++)
        mm.updValue(idim, jdim, EOperator::ADD, wz * coor[idim] * coor[jdim]);
    }
  }

  /* Normation */

  res.inertia /= res.wztot;
  for (int idim = 0; idim < ndim; idim++)
    for (int jdim = 0; jdim <= idim; jdim++)
      mm.updValue(idim, jdim, EOperator::DIVIDE, res.wztot);

  /* Calculate the eigen values and vectors */

  if (mm.computeEigen()) return res;
  res.mvalues = mm.getEigenValues();
  res.mvectors = *mm.getEigenVectors();
  double r  = res.mvalues[0] / res.mvalues[1];
  double e2 = (res.mvectors.getValue(1,0) / res.mvectors.getValue(0,0));
  e2 = e2 * e2;
  res.iso = 1. / sqrt(r);

  // Calculate the axes of the interia axes, so that:
  // i)  the ratio of the radii of the ellipse is proportional to the
  //     ratio of the eigen values
  // ii) the inertia of the (full) ellipse is equal to the inertia
  //     of the input variable

  double sx1 = res.mvectors.getValue(0,0) / ABS(res.mvectors.getValue(0,0));
  double sy1 = res.mvectors.getValue(1,0) / ABS(res.mvectors.getValue(1,0));
  double sx2 = res.mvectors.getValue(0,1) / ABS(res.mvectors.getValue(0,1));
  double sy2 = res.mvectors.getValue(1,1) / ABS(res.mvectors.getValue(1,1));

  // Calculate values of the structure #
  MatrixRectangular axes(4, 2);
  double K = sqrt(r) * (3. + r) / (2. * sqrt(2.) * pow(1. + r, 1.5));
  axes.setValue(0,0, res.center[0] + sx1 * sqrt(res.inertia / (K*(1. + e2))));
  axes.setValue(0,1, res.center[1] + sy1 * sqrt(res.inertia / (K*(1. + (1./e2)))));
  axes.setValue(1,0, 2. * res.center[0] - axes.getValue(0,0));
  axes.setValue(1,1, 2. * res.center[1] - axes.getValue(0,1));
  axes.setValue(2,0, res.center[0] + sx2 * sqrt(res.inertia / (K*r*(1. + (1./e2)))));
  axes.setValue(2,1, res.center[1] + sy2 * sqrt(res.inertia / (K*r*(1. + e2))));
  axes.setValue(3,0, 2. * res.center[0] - axes.getValue(2,0));
  axes.setValue(3,1, 2. * res.center[1] - axes.getValue(2,1));
  res.axesEllipse = axes;

  double r1 = res.mvalues[0] / (res.mvalues[0] + res.mvalues[1]);
  axes.setValue(0,0, res.center[0] + sx1 * sqrt((r1 * res.inertia) / (1. + e2)));
  axes.setValue(0,1, res.center[1] + sy1 * sqrt((r1 * res.inertia) / (1. + (1./e2))));
  axes.setValue(1,0, 2. * res.center[0] - axes.getValue(0,0));
  axes.setValue(1,1, 2. * res.center[1] - axes.getValue(0,1));
  axes.setValue(2,0, res.center[0] + sx2 * sqrt(((1.-r1) * res.inertia) / (1.+(1./e2))));
  axes.setValue(2,1, res.center[1] + sy2 * sqrt(((1.-r1) * res.inertia) / (1.+e2)));
  axes.setValue(3,0, 2. * res.center[0] - axes.getValue(2,0));
  axes.setValue(3,1, 2. * res.center[1] - axes.getValue(2,1));
  res.axesInertia = axes;

  double dx1 = res.axesInertia.getValue(1,0) - res.axesInertia.getValue(0,0);
  double dy1 = res.axesInertia.getValue(1,1) - res.axesInertia.getValue(0,1);
  res.theta = atan(dy1/dx1) * 180. / GV_PI;
  double dx2 = res.axesInertia.getValue(3,0) - res.axesInertia.getValue(2,0);
  double dy2 = res.axesInertia.getValue(3,1) - res.axesInertia.getValue(2,1);
  res.ra = 0.5 * sqrt(dx1*dx1 + dy1*dy1);
  res.rb = 0.5 * sqrt(dx2*dx2 + dy2*dy2);

  res.xl1 = {res.axesInertia.getValue(0,0), res.axesInertia.getValue(1,0)};
  res.yl1 = {res.axesInertia.getValue(0,1), res.axesInertia.getValue(1,1)};
  res.xl2 = {res.axesInertia.getValue(2,0), res.axesInertia.getValue(3,0)};
  res.yl2 = {res.axesInertia.getValue(2,1), res.axesInertia.getValue(3,1)};

  return res;
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
  double maille = 1.;
  if (db->isGrid())
  {
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(db);
    maille = dbgrid->getCellSize();
  }

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

void cgiPrint(const cgiOutput& cgioutput)
{
  VH::display("Gravity center", cgioutput.center);
  message("Number of active samples = %d\n", cgioutput.nvalid);
  message("Inertia = %lf\n", cgioutput.inertia);
  message("Isotropy = %lf\n", cgioutput.iso);
}

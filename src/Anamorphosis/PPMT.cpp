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
#include "geoslib_define.h"

#include "Anamorphosis/PPMT.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Db/Db.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Law.hpp"

PPMT::PPMT(int nvar)
    : AStringable(),
      _nvar(nvar),
      _anams(),
      _directions()
{
}

PPMT::PPMT(const PPMT &m)
    : AStringable(m),
      _nvar(m._nvar),
      _anams(m._anams),
      _directions(m._directions)
{
}

PPMT& PPMT::operator=(const PPMT &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _nvar = m._nvar;
    _anams = m._anams;
    _directions = m._directions;
  }
  return *this;
}

PPMT::~PPMT()
{
}

String PPMT::toString(const AStringFormat* strfmt) const
{
  SYMBOL_UNUSED(strfmt);

  std::stringstream sstr;

  int niter = (int) _anams.size();
  for (int iter = 0; iter < niter; iter++)
  {
    sstr << _anams[iter].toString(strfmt);
    sstr << "Direction = " << _directions[iter] << std::endl;
  }

  return sstr.str();
}

MatrixRectangular PPMT::fillLegendre(const VectorDouble& r, int n) const
{
  int nrow = (int) r.size();
  int ncol = n + 1;
  MatrixRectangular lp(nrow, ncol);

  // Initialization

  for (int i = 0; i < nrow; i++)
  {
    lp.setValue(i, 0, 1.);
    lp.setValue(i, 1, r[i]);
  }

  // Recursion

  for (int j = 1; j < n; j++)
    for (int i = 0; i < nrow; i++)
    {
      lp.setValue(i, j+1,
                  ((2*j+1) * r[i] * lp.getValue(i,j) -
                  (j) * lp.getValue(i,j-1))/(j+1));
    }
  return lp;
}

MatrixRectangular PPMT::sphering(const MatrixRectangular& X)
{
  if (X.isEmpty()) return MatrixRectangular();
  int nech = X.getNRows();
  int nvar = X.getNCols();

  AMatrix* TX = X.transpose();
  AMatrix* prod = prodMatrix(TX, &X);
  prod->prodScalar(1. / (double) nech);

  VectorDouble eigen_values(nvar);
  VectorDouble eigen_vectors(nvar * nvar);
  if (matrix_eigen(prod->getValues().data(), nvar,
                   eigen_values.data(), eigen_vectors.data()))
    return MatrixRectangular();

  // Invert the sign of the second Eigen vector (for compatibility with R output)
  MatrixSquareGeneral S(nvar);
  _S.reset(nvar,  nvar);
  S.setValues(eigen_vectors.data(),true);
  for (int ivar = 0; ivar < nvar ; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
    {
      double signe = (jvar < nvar-1) ? 1 : -1;
      S.setValue(ivar, jvar,
                 signe * S.getValue(ivar, jvar) / sqrt(eigen_values[jvar]));
    }

  AMatrix* YY = prodMatrix(&X, &S);
  MatrixRectangular Y(nech, nvar);
  Y.setValues(YY->getValues());

  delete TX;
  delete prod;
  delete YY;

  return Y;
}

VectorDouble PPMT::generateDirection(double angle) const
{
  VectorDouble direction(_nvar);
  direction[0] = cos(angle);
  direction[1] = sin(angle);
  return direction;
}

double PPMT::getIndex(const MatrixRectangular &X,
                      const VectorDouble &direction,
                      int j) const
{
  MatrixRectangular dirmat(_nvar, 1);
  dirmat.setValues(direction.data());
  AMatrix* XPP = prodMatrix(&X, &dirmat);

  VectorDouble r = XPP->getColumn(0);
  VH::divideConstant(r, 3.5);  // Normation empirique

  MatrixRectangular lp = fillLegendre(r, j);

  double idx = 0.;
  for (int l = 0; l < j; l++)
  {
    double mean = lp.getMeanByColumn(1+l);
    idx += (2*l+3)/2. * mean * mean;
  }
  delete XPP;
  return idx;
}

VectorDouble PPMT::optimize(const MatrixRectangular &X, int j, int N) const
{
  double idx_max = -1.;
  double ang_max = -1.;

  for (int i = 0; i < N; i++)
  {
    // Draw a direction at random

    double angle = (double)(i+1) * GV_PI / (double) N; // DR to fit R code
    VectorDouble direction = generateDirection(angle);
    double idx = getIndex(X, direction, j);
    if (idx < idx_max) continue;
    idx_max = idx;
    ang_max = angle;
  }

  VectorDouble result(2);
  result[0] = idx_max;
  result[1] = ang_max;
  return result;
}

MatrixRectangular PPMT::rotate(const MatrixRectangular &X,
                               double alpha,
                               bool direct) const
{
  double cs = cos(alpha);
  double sn = sin(alpha);
  MatrixSquareGeneral rotation = MatrixSquareGeneral(2);
  if (direct)
  {
    rotation.setValue(0, 0, cs);
    rotation.setValue(1, 0, sn);
    rotation.setValue(0, 1, -sn);
    rotation.setValue(1, 1, cs);
  }
  else
  {
    rotation.setValue(0, 0, cs);
    rotation.setValue(1, 0, -sn);
    rotation.setValue(0, 1, sn);
    rotation.setValue(1, 1, cs);
  }

  AMatrix* Xtemp = prodMatrix(&X, &rotation);

  MatrixRectangular XP = X;
  XP.setValues(Xtemp->getValues());
  return XP;
}

void PPMT::fit(const MatrixRectangular &X, int niter, int j, int N, int nbpoly)
{
  _anams.clear();
  _directions.clear();

  MatrixRectangular XP = sphering(X);

  for (int iter = 0; iter < niter; iter++)
  {
    VectorDouble result = optimize(XP, j, N);
    message("Iteration %d: angle=%lf max=%lf\n",iter+1,result[1],result[0]);

    double angle = result[1];
    MatrixRectangular XR = rotate(XP, angle, true);

    _directions.push_back(result);

    AnamHermite anam = AnamHermite(nbpoly);
    VectorDouble Z = XR.getColumn(0);
    anam.fitFromArray(Z);
    VectorDouble Y = anam.RawToTransformVec(Z);
    XR.setColumn(0, Y);
    MatrixRectangular XB = rotate(XR, angle, false);
    _anams.push_back(anam);

    XP = XB;
  }
}

MatrixRectangular PPMT::eval(const MatrixRectangular& X)
{
  MatrixRectangular XP = X;
  int nech = X.getNRows();
  int niter = getNiter();

  for (int iter = niter-1; iter >= 0; iter--)
  {
    double angle = _directions[iter][1];
    MatrixRectangular XR = rotate(XP, angle, true);
    VectorDouble Y = XR.getColumn(0);
    VectorDouble Z = _anams[iter].TransformToRawVec(Y);
    XR.setColumn(0, Z);
    MatrixRectangular XB = rotate(XR, angle, false);

    XP = XB;
  }

  // Inverse sphering

  MatrixSquareGeneral invS = _S;
  invS.invert();

  AMatrix* YY = prodMatrix(&XP, &invS);
  MatrixRectangular Y(nech, _nvar);
  Y.setValues(YY->getValues());

  return Y;
}

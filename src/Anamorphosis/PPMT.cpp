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
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Db/Db.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Law.hpp"

PPMT::PPMT(int nbpoly, int ndir, int legendre_order)
    : AStringable(),
      _nbpoly(nbpoly),
      _ndir(ndir),
      _legendreOrder(legendre_order),
      _anams(),
      _directions(),
      _nvar(0)
{
}

PPMT::PPMT(const PPMT &m)
    : AStringable(m),
      _nbpoly(m._nbpoly),
      _ndir(m._ndir),
      _legendreOrder(m._legendreOrder),
      _anams(m._anams),
      _directions(m._directions),
      _nvar(m._nvar)
{
}

PPMT& PPMT::operator=(const PPMT &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _nbpoly = m._nbpoly;
    _ndir = m._ndir;
    _legendreOrder = m._legendreOrder;
    _anams = m._anams;
    _directions = m._directions;
    _nvar = m._nvar;
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

MatrixRectangular PPMT::fillLegendre(const VectorDouble& r) const
{
  int nrow = (int) r.size();
  int ncol = _legendreOrder + 1;
  MatrixRectangular lp(nrow, ncol);

  // Initialization

  for (int i = 0; i < nrow; i++)
  {
    lp.setValue(i, 0, 1.);
    lp.setValue(i, 1, r[i]);
  }

  // Recursion

  for (int j = 1; j < _legendreOrder; j++)
    for (int i = 0; i < nrow; i++)
    {
      lp.setValue(i, j+1,
                  ((2*j+1) * r[i] * lp.getValue(i,j) -
                  (j) * lp.getValue(i,j-1))/(j+1));
    }
  return lp;
}

AMatrix* PPMT::sphering(const AMatrix* X)
{
  if (X->isEmpty()) return nullptr;
  int nech = X->getNRows();
  int nvar = X->getNCols();

  AMatrix* TX = X->transpose();
  AMatrix* prod = prodMatrix(TX, X);
  prod->prodScalar(1. / (double) nech);

  VectorDouble eigen_values(nvar);
  VectorDouble eigen_vectors(nvar * nvar);
  if (matrix_eigen(prod->getValues().data(), nvar,
                   eigen_values.data(), eigen_vectors.data()))
    return nullptr;

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

  AMatrix* Y = prodMatrix(X, &S);

  delete TX;
  delete prod;

  return Y;
}

VectorDouble PPMT::generateDirection(double angle) const
{
  VectorDouble direction(_nvar);
  direction[0] = cos(angle);
  direction[1] = sin(angle);
  return direction;
}

double PPMT::getIndex(const AMatrix *X, const VectorDouble &direction) const
{
  MatrixRectangular dirmat(_nvar, 1);
  dirmat.setValues(direction.data());
  AMatrix* XPP = prodMatrix(X, &dirmat);

  VectorDouble r = XPP->getColumn(0);
  delete XPP;

  VH::normalizeFromGaussianDistribution(r,-1.,1.);

  MatrixRectangular lp = fillLegendre(r);

  double idx = 0.;
  for (int l = 0; l < _legendreOrder; l++)
  {
    double mean = lp.getMeanByColumn(1+l);
    idx += (2*l+3)/2. * mean * mean;
  }
  return idx;
}

VectorDouble PPMT::optimize(const AMatrix* X) const
{
  double idx_max = -1.;
  double ang_max = -1.;

  for (int i = 0; i < _ndir; i++)
  {
    // Draw a direction at random

    double angle = (double)(i+1) * GV_PI / (double) _ndir;
    VectorDouble direction = generateDirection(angle);
    double idx = getIndex(X, direction);
    if (idx > idx_max)
    {
      idx_max = idx;
      ang_max = angle;
    }
  }

  VectorDouble result(2);
  result[0] = idx_max;
  result[1] = ang_max;
  return result;
}

AMatrix* PPMT::rotate(const AMatrix *X, double alpha, bool direct) const
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

  AMatrix* XR = prodMatrix(X, &rotation);

  return XR;
}

int PPMT::fit(const AMatrix* X, int niter)
{
  // Preliminary check
  if (_nbpoly <= 0)
  {
    messerr("The number of Hermite polynomials is not defined");
    return 1;
  }
  if (_ndir <= 0)
  {
    messerr("The number of Directions is not defined");
    return 1;
  }
  if (_legendreOrder <= 2)
  {
    messerr("The order of Legendre Polynomials is not defined");
    return 1;
  }

  // Cleaning
  _anams.clear();
  _directions.clear();
  _nvar = X->getNCols();

  // Processing
  AMatrix* XP = sphering(X);
  for (int iter = 0; iter < niter; iter++)
  {
    VectorDouble result = optimize(XP);
    message("Iteration %d: angle=%lf max=%lf\n",iter+1,result[1],result[0]);

    double angle = result[1];
    AMatrix* XR = rotate(XP, angle, true);
    delete XP;

    _directions.push_back(result);

    AnamHermite anam = AnamHermite(_nbpoly);
    VectorDouble Z = XR->getColumn(0);
    anam.fitFromArray(Z);
    VectorDouble Y = anam.RawToTransformVec(Z);
    XR->setColumn(0, Y);
    AMatrix* XB = rotate(XR, angle, false);
    delete XR;

    _anams.push_back(anam);

    XP = XB;
  }
  return 0;
}

AMatrix* PPMT::RawToTransform(const AMatrix* X)
{
  if (X->getNCols() != _nvar)
  {
    messerr("The input array has %d columns",X->getNCols());
    messerr("The number of variables in the PPMT is %d",_nvar);
    messerr("This is not correct");
    return nullptr;
  }
  int niter = getNiter();

  AMatrix* XP = dynamic_cast<AMatrix*>(X->clone());
  for (int iter = niter-1; iter >= 0; iter--)
  {
    double angle = _directions[iter][1];
    AMatrix* XR = rotate(XP, angle, true);
    delete XP;

    VectorDouble Y = XR->getColumn(0);
    VectorDouble Z = _anams[iter].TransformToRawVec(Y);
    XR->setColumn(0, Z);
    AMatrix* XB = rotate(XR, angle, false);
    delete XR;

    XP = XB;
  }

  // Inverse sphering

  MatrixSquareGeneral invS = _S;
  invS.invert();

  AMatrix* Y = prodMatrix(XP, &invS);

  return Y;
}

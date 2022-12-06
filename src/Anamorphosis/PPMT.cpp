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
#include "Matrix/AMatrix.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Db/Db.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Law.hpp"
#include "Basic/MathFunc.hpp"

#include <math.h>

PPMT::PPMT(int ndir, int niter, double alpha, const String &method)
    : AStringable(),
      _niter(niter),
      _ndir(ndir),
      _alpha(alpha),
      _method(method)
{
}

PPMT::PPMT(const PPMT &m)
    : AStringable(m),
      _niter(m._niter),
      _ndir(m._ndir),
      _alpha(m._alpha),
      _method(m._method)
{
}

PPMT& PPMT::operator=(const PPMT &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _niter = m._niter;
    _ndir = m._ndir;
    _ndim = m._ndim;
    _alpha = m._alpha;
    _method = m._method;
  }
  return *this;
}

PPMT::~PPMT()
{
}

PPMT* PPMT::create(int ndir, int niter, double alpha, const String &method)
{
  return new PPMT(ndir, niter, alpha, method);
}

String PPMT::toString(const AStringFormat* strfmt) const
{
  SYMBOL_UNUSED(strfmt);

  std::stringstream sstr;

  mestitle(1, "PPMT Method");
  if (getMethod() == "vdc")
    sstr << "- Using Van der Corput method" << std::endl;
  else
    sstr << "- Using Uniform method" << std::endl;
  sstr << "- Number of Directions =" << getNdir() << std::endl;
  sstr << "- Number of iterations =" << getNiter() << std::endl;
  sstr << "- Exponent value =" << getAlpha() << std::endl;

  return sstr.str();
}

MatrixRectangular PPMT::_fillLegendre(const VectorDouble& r, int legendreOrder) const
{
  int nrow = (int) r.size();
  int ncol = legendreOrder + 1;
  MatrixRectangular lp(nrow, ncol);

  // Initialization

  for (int i = 0; i < nrow; i++)
  {
    lp.setValue(i, 0, 1.);
    lp.setValue(i, 1, r[i]);
  }

  // Recursion

  for (int j = 1; j < legendreOrder; j++)
    for (int i = 0; i < nrow; i++)
    {
      lp.setValue(i, j+1,
                  ((2*j+1) * r[i] * lp.getValue(i,j) -
                  (j) * lp.getValue(i,j-1))/(j+1));
    }
  return lp;
}

AMatrix* PPMT::_sphering(const AMatrix* X)
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
  MatrixRectangular S(nvar, nvar);
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

void PPMT::_iteration(AMatrix *Y, const AMatrix* dir, int iter)
{
  int np   = Y->getNRows();
  int ndim = getNdim();
  int nd   = getNdir();
  double alpha = getAlpha();

  // Initialization
  VectorDouble sequence = VH::sequence(1., np, 1., 1. + np);
  VectorDouble N0 = VH::qnormVec(sequence);

  VectorDouble Y0(np, TEST);
  VectorInt R0(np, ITEST);
  int    idmax = -1;
  double ddmax = -1.e30;

  // Loop on directions
  VectorDouble Yi(np);
  for (int id = 0; id < nd; id++)
  {
    for (int ip = 0; ip < np; ip++)
    {
      double value = 0.;
      for (int idim = 0; idim < ndim; idim++)
        value += Y->getValue(ip, idim) * dir->getValue(id, idim);
      Yi[ip] = value;
    }

    VectorInt Ri = VH::sortRanks(Yi);
    double di = 0.;
    for (int ip = 0; ip < np; ip++)
    {
      double value = ABS(Yi[ip] - N0[Ri[ip]]);
      di += pow(value, alpha);
    }
    di /= (double) np;

    if (ddmax < di)
    {
      idmax = id;
      ddmax = di;
      Y0 = Yi;
      R0 = Ri;
    }
  }

  for (int ip = 0; ip < np; ip++)
  {
    double value = 0.;
    double scale = N0[R0[ip]] - Y0[ip];
    for (int idim = 0; idim < ndim; idim++)
    {
      value = Y->getValue(ip, idim) + scale * dir->getValue(idmax, idim);
      Y->setValue(ip, idim, value);
    }
  }

  // Returning arguments
  _serieAngle.push_back(idmax);
  _serieScore.push_back(ddmax);
  _directions.push_back(dir->getRow(idmax));
}

void PPMT::fit(AMatrix *Y, bool verbose)
{
  if (Y == nullptr)
  {
    messerr("Input Argument 'Y' (matrix) should be provided. Nothing is done");
  }
  _ndim  = Y->getNCols();
  int ndir = getNdir();
  int niter = getNiter();

  // Creating the directions

  MatrixRectangular* Umat;
  if (getMethod() == "vdc")
  {
    Umat = vanDerCorput(ndir, _ndim);
  }
  else
  {
    VectorDouble X = VH::simulateUniform(ndir * _ndim);
    Umat = MatrixRectangular::createFromVD(X, ndir, _ndim);
  }
  MatrixRectangular* dirmat = GeometryHelper::getDirectionsInRn(Umat);
  delete Umat;

  // Clearing the storage
  _serieAngle.clear();
  _serieScore.clear();
  _directions.clear();

  // Loop on the iterations
  for (int iter = 0; iter < niter; iter++)
    _iteration(Y, dirmat, iter);
  delete dirmat;

  // Optional printout
  if (verbose)
  {
    mestitle(1, "PPMT Method");
    if (getMethod() == "vdc")
      message("- Using Van der Corput method\n");
    else
      message("- Using Uniform method\n");
    message("- Space dimension = %d\n",getNdim());
    message("- Number of Directions = %d\n",getNdir());
    message("- Number of Iterations = %d\n",getNiter());
    message("- Exponent value = %lf\n", getAlpha());
  }
}

int PPMT::rawToGaussian(Db *db,
                        const VectorString &names,
                        bool useSel,
                        bool verbose,
                        const NamingConvention &namconv)
{
  // Extract the relevant information

  MatrixRectangular Y = db->getColumnsAsMatrix(names, useSel);
  if (Y.isEmpty())
  {
    messerr("This method requires several variables to be defined");
    return 1;
  }

  // Perform the PPMT procedure and transform the argument in place

  fit(&Y, verbose);

  // Add the newly created information in the Db

  int iptr = db->addColumns(Y.getValues());

  namconv.setNamesAndLocators(names, db, iptr);

  return 0;
}

VectorDouble PPMT::getSerieScore(bool flagLog) const
{
  VectorDouble vec;
  for (int iter = 0; iter < getNiter(); iter++)
  {
    double value = _serieScore[iter];
    if (flagLog) value = log(value);
    vec.push_back(value);
  }
  return vec;
}

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

#include "Anamorphosis/PPMT.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Db/Db.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorHelper.hpp"

PPMT::PPMT()
    : AStringable(),
      _anams(),
      _directions()
{
}

PPMT::PPMT(const PPMT &m)
    : AStringable(m),
      _anams(m._anams),
      _directions(m._directions)
{
}

PPMT& PPMT::operator=(const PPMT &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
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

void PPMT::addIteration(const AnamHermite& anam, const VectorDouble& dir)
{
  _anams.push_back(anam);
  _directions.push_back(dir);
}

MatrixRectangular PPMT::fillLegendre(const VectorDouble& r, int n)
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

MatrixSquareSymmetric PPMT::sphering(const MatrixRectangular& X)
{
  if (X.isEmpty()) return MatrixSquareSymmetric();
  int nech = X.getNRows();
  int nvar = X.getNCols();

  AMatrix* TX = X.transpose();
  AMatrix* prod = prodMatrix(TX, &X);
  prod->prodScalar(1. / (double) nech);
  prod->display();

//  if (matrix_eigen(prod.data(), ncol, _eigen.data(), _Z2F.data())) return (1);
  return MatrixSquareSymmetric();
}

/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include <LinearOp/ProdMatVect.hpp>

ProdMatVect::ProdMatVect(int nx, int ny, double *A)
    :
    ALinearOp(),
    _nx(nx),
    _ny(ny),
    _A(A)
{
}

ProdMatVect::~ProdMatVect()
{
}

void ProdMatVect::_evalDirect(const VectorDouble &inv, VectorDouble &outv) const
{
  double s = 0;
  for (int j = 0; j < _ny; j++)
    outv[j] = 0;

  for (int i = 0, nx = _nx; i < nx; i++)
  {
    s = 0;
    for (int j = 0; j < _ny; j++)
      s += _A[i + _ny * j] * inv[j];
    outv[i] = s;
  }
}


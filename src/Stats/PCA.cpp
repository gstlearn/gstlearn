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
#include "Stats/PCA.hpp"
#include "geoslib_f.h"
#include "geoslib_old_f.h"

PCA::PCA()
  : _nVar(0),
    _mean(),
    _sigma(),
    _eigen(),
    _Z2F(),
    _F2Z()
{
}

PCA::PCA(const PCA &m)
    : _nVar(m._nVar),
      _mean(m._mean),
      _sigma(m._sigma),
      _eigen(m._eigen),
      _Z2F(m._Z2F),
      _F2Z(m._F2Z)
{

}

PCA& PCA::operator=(const PCA &m)
{
  if (this != &m)
  {
    _nVar = m._nVar;
    _mean = m._mean;
    _sigma = m._sigma;
    _eigen = m._eigen;
    _Z2F = m._Z2F;
    _F2Z = m._F2Z;
  }
  return *this;
}

PCA::~PCA()
{

}

void PCA::init(int nvar)
{
  _nVar   = nvar;
  _mean.resize(nvar);
  _sigma.resize(nvar,0);
  _eigen.resize(nvar,0);
  _Z2F.resize(nvar * nvar,0);
  _F2Z.resize(nvar * nvar,0);
}

void PCA::clean()
{
  for (int ivar=0; ivar<_nVar; ivar++)
  {
    _mean[ivar]  = 0.;
    _sigma[ivar] = 1.;
    _eigen[ivar] = 0.;
    for (int jvar=0; jvar<_nVar; jvar++)
    {
      _Z2F[ivar * _nVar + jvar] = 0.;
      _F2Z[ivar * _nVar + jvar] = 0.;
    }
  }
}

int PCA::calculateEigen(int nvar, VectorDouble& c0)
{
  _nVar = nvar;
  _eigen.resize (nvar * nvar,0);
  _Z2F.resize(nvar * nvar,0);
  _F2Z.resize(nvar * nvar,0);

  /* Eigen decomposition */

  if (matrix_eigen(c0.data(), nvar, _eigen.data(),
                   _Z2F.data())) return (1);
  if (matrix_invert_copy(_Z2F.data(), nvar,
                         _F2Z.data())) return (1);

  return(0);
}

void PCA::display(int flag_center, int flag_stats)

{
  mestitle(1, "PCA Transform");
  if (flag_center)
  {
    print_matrix("Means", 0, 1, 1, _nVar, NULL, _mean.data());
    print_matrix("Standard deviations", 0, 1, 1, _nVar, NULL, _sigma.data());
  }
  if (flag_stats)
    print_matrix("Eigen Values", 0, 1, 1, _nVar, NULL, _eigen.data());

  message("\n");
  print_matrix("Matrix M to transform standardized Variables Z into Factors Y",
               0, 1, _nVar, _nVar, NULL, _Z2F.data());
  message("Y = Z * M (columns  = eigen vectors)\n");
  message("\n");
  print_matrix(
      "Matrix t(M) to back-transform Factors Y into standardized Variables Z",
      0, 1, _nVar, _nVar, NULL, _F2Z.data());
  message("Z = Y * t(M) (rows  = eigen vectors)\n");
}


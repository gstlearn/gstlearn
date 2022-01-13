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

#include "Stats/PCAStringFormat.hpp"
#include "geoslib_f.h"
#include "geoslib_old_f.h"

PCA::PCA(int nvar)
  : AStringable(),
    _nVar(0),
    _mean(),
    _sigma(),
    _eigen(),
    _Z2F(),
    _F2Z()
{
  init(nvar);
}

PCA::PCA(const PCA &m)
    : AStringable(m),
      _nVar(m._nVar),
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
    AStringable::operator=(m);
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

String PCA::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  const PCAStringFormat* pcafmt = dynamic_cast<const PCAStringFormat*>(strfmt);
  PCAStringFormat dsf;
  if (pcafmt != nullptr) dsf = *pcafmt;

  sstr << toTitle(1, "PCA Transform");

  if (dsf.getflagCenter())
  {
    sstr << toMatrix("Means", VectorString(), VectorString(), true, 1, _nVar,
                    _mean);
    sstr << toMatrix("Standard deviations", VectorString(), VectorString(), true,
                    1, _nVar, _sigma);
  }
  if (dsf.getflagStats())
    sstr << toMatrix("Eigen Values", VectorString(), VectorString(), true, 1,
                    _nVar, _eigen);

  sstr << toMatrix("Matrix M to transform standardized Variables Z into Factors Y",
                   VectorString(), VectorString(), true,
                   _nVar, _nVar, _Z2F);
  sstr << "Y = Z * M (columns  = eigen vectors)" << std::endl;
  sstr << toMatrix("Matrix t(M) to back-transform Factors Y into standardized Variables Z",
                   VectorString(), VectorString(), true,
                   _nVar, _nVar, _F2Z);
  sstr << "Z = Y * t(M) (rows  = eigen vectors)" << std::endl;

  return sstr.str();
}

int PCA::compute(const Db *db, bool verbose)
{
  return pca_compute(db, verbose, this);
}
int PCA::dbZ2F(Db* db,
               bool flag_norm,
               bool verbose,
               const NamingConvention& namconv)
{
  return pca_z2f(db, this, flag_norm, verbose, namconv);
}

int PCA::dbF2Z(Db* db,
               bool flag_norm,
               bool verbose,
               const NamingConvention& namconv)
{
  return pca_f2z(db, this, flag_norm, verbose, namconv);
}

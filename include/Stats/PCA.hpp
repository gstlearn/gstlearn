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
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Vector.hpp"
#include "Basic/NamingConvention.hpp"

class Db;

class Db;

class GSTLEARN_EXPORT PCA: public AStringable
{
public:
  PCA(int nvar = 1);
  PCA(const PCA &m);
  PCA& operator= (const PCA &m);
  virtual ~PCA();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  void init(int nvar);
  void clean();
  int calculateEigen(int nvar, VectorDouble& c0);

  const VectorDouble& getEigen() const { return _eigen; }
  double getEigen(int ivar) const { return _eigen[ivar]; }
  const VectorDouble& getMean() const { return _mean; }
  double getMean(int ivar) const { return _mean[ivar]; }
  int getNVar() const { return _nVar; }
  const VectorDouble& getF2Z() const { return _F2Z; }
  double getF2Z(int ivar, int jvar) const { return _F2Z[_getAddress(ivar,jvar)]; }
  const VectorDouble& getZ2F() const { return _Z2F; }
  double getZ2F(int ivar, int jvar) const { return _Z2F[_getAddress(ivar, jvar)]; }
  const VectorDouble& getSigma() const { return _sigma; }
  double getSigma(int ivar) const { return _sigma[ivar]; }

  void setPcaZ2F(VectorDouble& pcaz2f) { _Z2F = pcaz2f; }
  void setPcaZ2F(int ivar, int jvar, double pcaz2f) { _Z2F[_getAddress(ivar,jvar)] = pcaz2f; }
  void setPcaF2Z(VectorDouble& pcaf2z) { _F2Z = pcaf2z; }
  void setEigen(VectorDouble& eigen) { _eigen = eigen; }
  void setEigen(int ivar, double eigen) { _eigen[ivar] = eigen; }
  void setMean(VectorDouble& mean) { _mean = mean; }
  void setSigma(VectorDouble& sigma) { _sigma = sigma; }

  int compute(const Db *db, bool verbose = false);
  int dbZ2F(Db* db,
            bool flag_norm = true,
            bool verbose = false,
            const NamingConvention& namconv = NamingConvention("Z2F"));
  int dbF2Z(Db* db,
            bool flag_norm = true,
            bool verbose = false,
            const NamingConvention& namconv = NamingConvention("F2Z"));

private:
  int _getAddress(int ivar, int jvar) const { return (ivar * _nVar + jvar); }

private:
  int          _nVar;
  VectorDouble _mean;
  VectorDouble _sigma;
  VectorDouble _eigen;
  VectorDouble _Z2F;
  VectorDouble _F2Z;
};

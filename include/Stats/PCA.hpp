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
#include "Variogram/DirParam.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Vector.hpp"
#include "Basic/NamingConvention.hpp"

class Db;

class GSTLEARN_EXPORT PCA: public AStringable
{
public:
  PCA(int nvar = 0);
  PCA(const Db *db, bool verbose = false);
  PCA(Db *db, double h0, double dh, const DirParam& dirparam = DirParam(), bool verbose = false);
  PCA(const PCA &m);
  PCA& operator= (const PCA &m);
  virtual ~PCA();

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  void init(int nvar);

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

  void setPcaZ2F(const VectorDouble& pcaz2f) { _Z2F = pcaz2f; }
  void setPcaZ2F(int ivar, int jvar, double pcaz2f) { _Z2F[_getAddress(ivar,jvar)] = pcaz2f; }
  void setPcaF2Z(VectorDouble& pcaf2z) { _F2Z = pcaf2z; }
  void setEigen(VectorDouble& eigen) { _eigen = eigen; }
  void setEigen(int ivar, double eigen) { _eigen[ivar] = eigen; }
  void setMean(const VectorDouble& mean) { _mean = mean; }
  void setSigma(const VectorDouble& sigma) { _sigma = sigma; }

  int pca_compute(const Db *db, bool verbose = false);
  int maf_compute(Db *db,
                  double h0,
                  double dh,
                  const DirParam& dirparam,
                  bool verbose = false);
  int dbZ2F(Db* db,
            bool verbose = false,
            const NamingConvention& namconv = NamingConvention("Z2F"));
  int dbF2Z(Db* db,
            bool verbose = false,
            const NamingConvention& namconv = NamingConvention("F2Z"));
  VectorDouble mafOfIndex() const;

private:
  int _getAddress(int ivar, int jvar) const { return (ivar * _nVar + jvar); }
  VectorBool _getVectorIsotopic(const Db* db);
  void _loadData(const Db* db, int iech, VectorDouble& data);
  int  _calculateEigen(VectorDouble& c0);
  int _pcaCalculate(const Db *db, const VectorBool& isoFlag, bool verbose);
  int _normalization(const Db *db,
                     const VectorBool& isoFlag,
                     VectorDouble& mean,
                     VectorDouble& sigma,
                     bool verbose);
  void _covariance0(const Db *db,
                    const VectorBool& isoFlag,
                    const VectorDouble& mean,
                    const VectorDouble& sigma,
                    VectorDouble& c0,
                    bool verbose);
  int _covarianceh(Db *db,
                   double h0,
                   double dh,
                   const DirParam& dirparam,
                   const VectorBool& isoFlag,
                   VectorDouble& ch,
                   bool verbose);
  void _center(VectorDouble& data,
               const VectorDouble &mean,
               const VectorDouble &sigma);
  void _uncenter(VectorDouble& data,
                 const VectorDouble &mean,
                 const VectorDouble &sigma);
  void _pcaZ2F(bool flag_norm,
               int iptr,
               Db *db,
               const VectorBool isoFlag,
               const VectorDouble& mean,
               const VectorDouble& sigma);
  void _pcaF2Z(int iptr, Db *db, const VectorBool& isoFlag);

private:
  int          _nVar;
  VectorDouble _mean;
  VectorDouble _sigma;
  VectorDouble _eigen;
  VectorDouble _Z2F;
  VectorDouble _F2Z;
};

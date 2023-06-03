/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Variogram/VarioParam.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/NamingConvention.hpp"

class Db;

class GSTLEARN_EXPORT PCA: public AStringable
{
public:
  PCA(int nvar = 0);
  PCA(const PCA &m);
  PCA& operator= (const PCA &m);
  virtual ~PCA();

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  void init(int nvar);

  const VectorDouble& getEigVals() const { return _eigval; }
  double getEigVal(int ivar) const { return _eigval[ivar]; }
  const VectorDouble& getEigVecs() const { return _eigvec; }
  double getEigVec(int ivar, int jvar) const { return _eigvec[_getAddress(ivar,jvar)]; }
  VectorDouble getVarianceRatio() const;
  const VectorDouble& getMeans() const { return _mean; }
  double getMean(int ivar) const { return _mean[ivar]; }
  int getNVar() const { return _nVar; }
  const VectorDouble& getF2Zs() const { return _F2Z; }
  double getF2Z(int ivar, int jvar) const { return _F2Z[_getAddress(ivar,jvar)]; }
  const VectorDouble& getZ2Fs() const { return _Z2F; }
  double getZ2F(int ivar, int jvar) const { return _Z2F[_getAddress(ivar, jvar)]; }
  const VectorDouble& getSigmas() const { return _sigma; }
  double getSigma(int ivar) const { return _sigma[ivar]; }

  void setMeans(const VectorDouble &mean) { _mean = mean; }
  void setSigmas(const VectorDouble &sigma) { _sigma = sigma; }
  void setZ2Fs(const VectorDouble& z2f) { _Z2F = z2f; }
  void setZ2F(int ivar, int jvar, double z2f) { _Z2F[_getAddress(ivar,jvar)] = z2f; }
  void setF2Zs(VectorDouble& f2z) { _F2Z = f2z; }
  void setF2Z(int ivar, int jvar, double f2z) { _F2Z[_getAddress(ivar,jvar)] = f2z; }
  void setEigVals(VectorDouble& eigval) { _eigval = eigval; }
  void setEigVal(int ivar, double eigval) { _eigval[ivar] = eigval; }
  void setEigVecs(const VectorDouble& eigvec) { _eigvec = eigvec; }
  void setEigVec(int ivar, int jvar, double eigvec) { _eigvec[_getAddress(ivar,jvar)] = eigvec; }

  int pca_compute(const Db *db, bool verbose = false);
  int maf_compute(Db *db,
                  const VarioParam &varioparam,
                  int ilag0 = 1,
                  int idir0 = 0,
                  bool verbose = false);
  int maf_compute_interval(Db *db, double hmin, double hmax, bool verbose = false);
  int dbZ2F(Db* db,
            bool verbose = false,
            const NamingConvention& namconv = NamingConvention("F", false));
  int dbF2Z(Db* db,
            bool verbose = false,
            const NamingConvention& namconv = NamingConvention("Z", false));
  VectorDouble mafOfIndex() const;

private:
  int _getAddress(int ivar, int jvar) const { return (ivar * _nVar + jvar); }
  VectorBool _getVectorIsotopic(const Db* db);
  void _loadData(const Db* db, int iech, VectorDouble& data);
  int  _calculateEigen(bool verbose = false);
  int  _calculateGEigen(bool verbose);
  void _calculateNormalization(const Db *db,
                               const VectorBool &isoFlag,
                               bool verbose = false,
                               bool flag_nm1 = false);
  void _covariance0(const Db *db,
                    const VectorBool &isoFlag,
                    bool verbose = false,
                    bool flag_nm1 = false);
  void _variogramh(Db *db,
                   const VarioParam &varioparam,
                   int ilag0,
                   int idir0,
                   double hmin,
                   double hmax,
                   const VectorBool &isoFlag,
                   bool verbose);
  void _center(VectorDouble& data,
               const VectorDouble &mean,
               const VectorDouble &sigma,
               bool flag_center = true,
               bool flag_scale = false);
  void _uncenter(VectorDouble& data,
                 const VectorDouble &mean,
                 const VectorDouble &sigma,
                 bool flag_center = true,
                 bool flag_scale = false);
  void _pcaZ2F(int iptr,
               Db *db,
               const VectorBool isoFlag,
               const VectorDouble &mean,
               const VectorDouble &sigma);
  void _pcaF2Z(int iptr,
               Db *db,
               const VectorBool &isoFlag,
               const VectorDouble &mean,
               const VectorDouble &sigma);
  void _pcaFunctions(bool verbose);
  void _mafFunctions(bool verbose);
  int _mafCompute(Db *db,
                  const VarioParam &varioparam,
                  int ilag0,
                  int idir0,
                  double hmin,
                  double hmax,
                  bool verbose);

private:
  int          _nVar;
  VectorDouble _mean;
  VectorDouble _sigma;
  VectorDouble _eigval;
  VectorDouble _eigvec;
  VectorDouble _c0;
  VectorDouble _gh;
  VectorDouble _Z2F;
  VectorDouble _F2Z;
};

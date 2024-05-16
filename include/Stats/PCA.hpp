/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Variogram/VarioParam.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/MatrixRectangular.hpp"
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
  const MatrixRectangular& getEigVecs() const { return _eigvec; }
  double getEigVec(int ivar, int jvar) const { return _eigvec.getValue(ivar,jvar); }
  const VectorDouble getVarianceRatio() const;
  const VectorDouble& getMeans() const { return _mean; }
  double getMean(int ivar) const { return _mean[ivar]; }
  const MatrixSquareSymmetric& getC0() const { return _c0; }
  int getNVar() const { return _nVar; }
  const MatrixSquareGeneral& getF2Zs() const { return _F2Z; }
  const MatrixSquareGeneral& getZ2Fs() const { return _Z2F; }
  const VectorDouble& getSigmas() const { return _sigma; }
  double getSigma(int ivar) const { return _sigma[ivar]; }

  void setMeans(const VectorDouble &mean) { _mean = mean; }
  void setSigmas(const VectorDouble &sigma) { _sigma = sigma; }
  void setZ2Fs(const MatrixSquareGeneral& z2f) { _Z2F = z2f; }
  void setF2Zs(MatrixSquareGeneral& f2z) { _F2Z = f2z; }

  void setEigVals(VectorDouble& eigval) { _eigval = eigval; }
  void setEigVal(int ivar, double eigval) { _eigval[ivar] = eigval; }
  void setEigVecs(const MatrixRectangular& eigvec) { _eigvec = eigvec; }
  void setEigVec(int ivar, int jvar, double eigvec) { _eigvec.setValue(ivar,jvar,eigvec); }

  int pca_compute(const Db *db, bool verbose = false, bool optionPositive = true);
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
  double _getF2Z(int ivar, int ifac) const { return _F2Z.getValue(ivar, ifac); }
  void   _setF2Z(int ivar, int ifac, double f2z) {  _F2Z.setValue(ivar, ifac, f2z); }
  double _getZ2F(int ifac, int ivar) const { return _Z2F.getValue(ifac, ivar); }
  void   _setZ2F(int ifac, int ivar, double z2f) {  _Z2F.setValue(ifac, ivar, z2f); }

  VectorBool _getVectorIsotopic(const Db* db);
  void _loadData(const Db* db, int iech, VectorDouble& data);
  int  _calculateEigen(bool verbose = false, bool optionPositive = true);
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
  VectorDouble          _eigval;
  MatrixRectangular     _eigvec;
  MatrixSquareSymmetric _c0;
  MatrixSquareSymmetric _gh;
  MatrixSquareGeneral   _Z2F;
  MatrixSquareGeneral   _F2Z;
};

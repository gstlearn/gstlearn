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
#include "geoslib_define.h"

#include "Space/ASpaceObject.hpp"

#include "Matrix/MatrixSquareGeneral.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Space/SpacePoint.hpp"

#include <vector>

class Db;
class DbGrid;
class MatrixRectangular;
class MatrixSparse;

/**
 * \brief
 * Class containing the Covariance part of the Model.
 *
 * It is the uppermost class of the Covariance Tree and is conceived as simple as possible on purpose
 * (in order to let the user defined its own version if necessary): it must simply be able to return its value
 * between two end-point (see eval method).
 *
 * It is mainly implemented in CovAniso.hpp or ACovAnisoList.hpp
 */
class GSTLEARN_EXPORT ACov : public ASpaceObject
{
public:
  ACov(const ASpace* space = nullptr);
  ACov(const ACov &r);
  ACov& operator=(const ACov &r);
  virtual ~ACov();

  /// ACov Interface
  virtual int getNVariables() const = 0;
  virtual bool isIndexable() const { return false; }
  virtual bool isNoStat() const { return false; }

 
  /// Calculate the covariance between two variables for 0-distance (stationary case)
  virtual double eval0(int ivar = 0,
                       int jvar = 0,
                       const CovCalcMode* mode = nullptr) const = 0;
  /// Calculate the matrix of covariances for 0-distance (stationary case)
  virtual void eval0MatInPlace(MatrixSquareGeneral &mat,
                               const CovCalcMode *mode = nullptr) const;
  /// Calculate the covariance between two variables and two points (general case)
  virtual double eval(const SpacePoint& p1,
                      const SpacePoint& p2,
                      int ivar = 0,
                      int jvar = 0,
                      const CovCalcMode* mode = nullptr) const = 0;
  /// Calculate the matrix of covariances between two points (general case)
  virtual void evalMatInPlace(const SpacePoint &p1,
                              const SpacePoint &p2,
                              MatrixSquareGeneral &mat,
                              const CovCalcMode *mode = nullptr) const;
  /// Calculate the matrix of covariances between two points given by indices (optim)
  virtual void evalMatOptimInPlace(int icas1,
                                   int iech1,
                                   int icas2,
                                   int iech2,
                                   MatrixSquareGeneral &mat,
                                   const CovCalcMode *mode = nullptr) const = 0;
  /// Tell if the use of Optimization is enabled or not
  virtual bool isOptimEnabled() const { return _isOptimEnabled; }

  virtual double evalCovOnSphere(double alpha,
                                 int degree = 50,
                                 bool flagScaleDistance = false,
                                 const CovCalcMode* mode = nullptr) const
  {
    DECLARE_UNUSED(alpha);
    DECLARE_UNUSED(degree);
    DECLARE_UNUSED(flagScaleDistance);
    DECLARE_UNUSED(mode);
    return TEST;
  }
  virtual VectorDouble evalSpectrumOnSphere(int n,
                                            bool flagNormDistance = false,
                                            bool flagCumul = false) const
  {
    DECLARE_UNUSED(n);
    DECLARE_UNUSED(flagNormDistance);
    DECLARE_UNUSED(flagCumul);
    return VectorDouble();
  }
  virtual double evalSpectrum(const VectorDouble &freq,
                              int ivar,
                              int jvar) const
  {
    DECLARE_UNUSED(freq);
    DECLARE_UNUSED(ivar);
    DECLARE_UNUSED(jvar);
    return TEST;
  }

  virtual void updateCovByPoints(int icas1, int iech1, int icas2, int iech2)
  {
    DECLARE_UNUSED(icas1);
    DECLARE_UNUSED(iech1);
    DECLARE_UNUSED(icas2);
    DECLARE_UNUSED(iech2);
  }

  /////////////////////////////////////////////////////////////////////////////////
  ///
  void setOptimEnabled(bool isOptimEnabled) { _isOptimEnabled = isOptimEnabled; }
  VectorDouble eval(const std::vector<SpacePoint>& vec_p1,
                    const std::vector<SpacePoint>& vec_p2,
                    int ivar = 0,
                    int jvar = 0,
                    const CovCalcMode* mode = nullptr) const;
  MatrixSquareGeneral eval0Mat(const CovCalcMode* mode = nullptr) const;
  MatrixSquareGeneral evalMat(const SpacePoint& p1,
                              const SpacePoint& p2,
                              const CovCalcMode* mode = nullptr) const;
  double evalIvarIpas(double step,
                      const VectorDouble& dir = VectorDouble(),
                      int ivar = 0,
                      int jvar = 0,
                      const CovCalcMode* mode = nullptr) const;
  double evalIvarIpasIncr(const VectorDouble &dincr,
                          int ivar = 0,
                          int jvar = 0,
                          const CovCalcMode* mode = nullptr) const;
  VectorDouble evalIvarNpas(const VectorDouble& vec_step,
                            const VectorDouble& dir = VectorDouble(),
                            int ivar = 0,
                            int jvar = 0,
                            const CovCalcMode* mode = nullptr) const;
  MatrixSquareGeneral evalNvarIpas(double step,
                                   const VectorDouble& dir = VectorDouble(),
                                   const CovCalcMode* mode = nullptr) const;
  MatrixSquareGeneral evalNvarIpasIncr(const VectorDouble &dincr,
                                       const CovCalcMode* mode = nullptr) const;
  double evalIsoIvarIpas(double step,
                         int ivar = 0,
                         int jvar = 0,
                         const CovCalcMode* mode = nullptr) const;
  VectorDouble evalIsoIvarNpas(const VectorDouble& vec_step,
                               int ivar = 0,
                               int jvar = 0,
                               const CovCalcMode* mode = nullptr) const;
  MatrixSquareGeneral evalIsoNvarIpas(double step,
                                      const CovCalcMode* mode = nullptr) const;

  double evalCvv(const VectorDouble& ext,
                 const VectorInt& ndisc,
                 const VectorDouble& angles = VectorDouble(),
                 int ivar = 0,
                 int jvar = 0,
                 const CovCalcMode* mode = nullptr) const;
  double evalCvvShift(const VectorDouble& ext,
                      const VectorInt& ndisc,
                      const VectorDouble& shift,
                      const VectorDouble& angles = VectorDouble(),
                      int ivar = 0,
                      int jvar = 0,
                      const CovCalcMode* mode = nullptr) const;
  MatrixSquareGeneral evalCvvM(const VectorDouble& ext,
                               const VectorInt& ndisc,
                               const VectorDouble& angles = VectorDouble(),
                               const CovCalcMode* mode = nullptr) const;
  double evalCxv(const SpacePoint& p1,
                 const VectorDouble& ext,
                 const VectorInt& ndisc,
                 const VectorDouble& angles = VectorDouble(),
                 const VectorDouble& x0 = VectorDouble(),
                 int ivar = 0,
                 int jvar = 0,
                 const CovCalcMode* mode = nullptr) const;
  double evalCxv(const Db* db,
                 const VectorDouble& ext,
                 const VectorInt& ndisc,
                 const VectorDouble& angles = VectorDouble(),
                 const VectorDouble& x0 = VectorDouble(),
                 int ivar = 0,
                 int jvar = 0,
                 const CovCalcMode* mode = nullptr) const;
  MatrixSquareGeneral evalCxvM(const SpacePoint& p1,
                               const VectorDouble& ext,
                               const VectorInt& ndisc,
                               const VectorDouble& angles = VectorDouble(),
                               const VectorDouble& x0 = VectorDouble(),
                               const CovCalcMode* mode = nullptr) const;
  VectorDouble evalPointToDb(const SpacePoint& p1,
                             const Db* db2,
                             int ivar = 0,
                             int jvar = 0,
                             bool useSel = true,
                             const VectorInt& nbgh2 = VectorInt(),
                             const CovCalcMode* mode = nullptr) const;
  VectorDouble evalPointToDbAsSP(const std::vector<SpacePoint>& p1s,
                                 const SpacePoint& p2,
                                 int ivar = 0,
                                 int jvar = 0,
                                 const CovCalcMode* mode = nullptr) const;
  double evalAverageDbToDb(const Db* db1,
                           const Db* db2,
                           int ivar = 0,
                           int jvar = 0,
                           double eps = 0.,
                           int seed = 434132,
                           const CovCalcMode* mode = nullptr) const;
  double evalAverageIncrToIncr(const VectorVectorDouble& d1,
                               const VectorVectorDouble& d2,
                               int ivar = 0,
                               int jvar = 0,
                               const CovCalcMode* mode = nullptr) const;
  double evalAveragePointToDb(const SpacePoint& p1,
                              const Db* db2,
                              int ivar = 0,
                              int jvar = 0,
                              const CovCalcMode* mode = nullptr) const;
  MatrixRectangular evalCovMatrix(const Db* db1_arg,
                                  const Db* db2_arg = nullptr,
                                  int ivar0 = -1,
                                  int jvar0 = -1,
                                  const VectorInt& nbgh1 = VectorInt(),
                                  const VectorInt& nbgh2 = VectorInt(),
                                  const CovCalcMode* mode = nullptr);
  MatrixSquareSymmetric evalCovMatrixSymmetric(const Db *db1,
                                               int ivar0,
                                               const VectorInt &nbgh1,
                                               const CovCalcMode *mode);
  MatrixSparse* evalCovMatrixSparse(const Db *db1_arg,
                                    const Db *db2_arg = nullptr,
                                    int ivar0 = -1,
                                    int jvar0 = -1,
                                    const VectorInt &nbgh1 = VectorInt(),
                                    const VectorInt &nbgh2 = VectorInt(),
                                    const CovCalcMode *mode = nullptr,
                                    double eps = EPSILON3);
  double extensionVariance(const Db* db,
                           const VectorDouble& ext,
                           const VectorInt& ndisc,
                           const VectorDouble& angles = VectorDouble(),
                           const VectorDouble& x0 = VectorDouble(),
                           int ivar = 0,
                           int jvar = 0) const;
  double samplingDensityVariance(const Db* db,
                                 const VectorDouble& ext,
                                 const VectorInt& ndisc,
                                 const VectorDouble& angles = VectorDouble(),
                                 const VectorDouble& x0 = VectorDouble(),
                                 int ivar = 0,
                                 int jvar = 0) const;
  double specificVolume(const Db *db,
                        double mean,
                        const VectorDouble &ext,
                        const VectorInt &ndisc,
                        const VectorDouble &angles = VectorDouble(),
                        const VectorDouble &x0 = VectorDouble(),
                        int ivar = 0,
                        int jvar = 0) const;
  double coefficientOfVariation(const Db *db,
                                double volume,
                                double mean,
                                const VectorDouble &ext,
                                const VectorInt &ndisc,
                                const VectorDouble &angles = VectorDouble(),
                                const VectorDouble &x0 = VectorDouble(),
                                int ivar = 0,
                                int jvar = 0) const;
  double specificVolumeFromCoV(Db *db,
                               double cov,
                               double mean,
                               const VectorDouble &ext,
                               const VectorInt &ndisc,
                               const VectorDouble &angles = VectorDouble(),
                               const VectorDouble &x0 = VectorDouble(),
                               int ivar = 0,
                               int jvar = 0) const;


  void manage(const Db* db1,const Db* db2) const
  {
      _manage(db1, db2);
  }
protected:

  VectorInt _getActiveVariables(int ivar0) const;
  static void _updateCovMatrixSymmetricVerr(const Db* db1,
                                            AMatrix* mat,
                                            const VectorInt& ivars,
                                            const VectorVectorInt& index1);

private:
  virtual void _manage(const Db* db1,const Db* db2) const 
  {
    DECLARE_UNUSED(db1)
    DECLARE_UNUSED(db2)
  }

  DbGrid* _discretizeBlock(const VectorDouble& ext,
                           const VectorInt& ndisc,
                           const VectorDouble& angles = VectorDouble(),
                           const VectorDouble& x0 = VectorDouble()) const;
  Db* _discretizeBlockRandom(const DbGrid* dbgrid, int seed = 34131) const;
  double _getVolume(const VectorDouble& ext) const;

protected:
  bool _isOptimEnabled;

  // These temporary information is used to speed up processing (optimization functions)
  // They are in a protected section as they may be modified by class hierarchy
  mutable std::vector<SpacePoint> _p1As;
  mutable SpacePoint _p2A;
};

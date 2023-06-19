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
#include "geoslib_define.h"

#include "Space/ASpaceObject.hpp"

#include "Matrix/MatrixSquareGeneral.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Space/SpacePoint.hpp"

#include <vector>

class Db;
class DbGrid;
class MatrixRectangular;

class GSTLEARN_EXPORT ACov : public ASpaceObject
{
public:
  ACov(const ASpace* space = nullptr);
  ACov(const ACov &r);
  ACov& operator=(const ACov &r);
  virtual ~ACov();

  /// ACov Interface
  virtual int getNVariables() const = 0;
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

  virtual void evalOptim(const SpacePoint &p1,
                         VectorDouble &res,
                         VectorDouble &temp,
                         SpacePoint &pttr,
                         int ivar = 0,
                         int jvar = 0,
                         const CovCalcMode* mode = nullptr) const { };

  virtual double evalCovOnSphere(double /*alpha*/,
                                 int /*degree*/,
                                 bool /*normalize*/) const { return TEST; }
  virtual double evalSpectrum(const VectorDouble& /*freq*/,
                              int /*ivar*/, int /*jvar*/) const { return TEST; }

  virtual void 	preProcess(const std::vector<SpacePoint>& vec) const {};
  virtual void  cleanPreProcessInfo() const {}

  /////////////////////////////////////////////////////////////////////////////////

  VectorDouble eval(const std::vector<SpacePoint>& vec_p1,
                    const std::vector<SpacePoint>& vec_p2,
                    int ivar = 0,
                    int jvar = 0,
                    const CovCalcMode* mode = nullptr) const;
  void evalVect(VectorDouble &res,
                const SpacePoint &p1,
                const std::vector<SpacePoint> &vec_p2,
                int ivar = 0,
                int jvar = 0,
                const CovCalcMode* mode = nullptr) const;
  MatrixSquareGeneral eval0Mat(const CovCalcMode* mode = nullptr) const;
  MatrixSquareGeneral evalMat(const SpacePoint& p1,
                              const SpacePoint& p2,
                              const CovCalcMode* mode = nullptr) const;
  double evalIvarIpas(double step,
                      const VectorDouble& dir,
                      int ivar = 0,
                      int jvar = 0,
                      const VectorDouble& center = VectorDouble(),
                      const CovCalcMode* mode = nullptr) const;
  double evalIvarIpasIncr(const VectorDouble &dincr,
                          int ivar = 0,
                          int jvar = 0,
                          const CovCalcMode* mode = nullptr) const;
  VectorDouble evalIvarNpas(const VectorDouble& vec_step,
                            const VectorDouble& dir = VectorDouble(),
                            int ivar = 0,
                            int jvar = 0,
                            const VectorDouble& center = VectorDouble(),
                            const CovCalcMode* mode = nullptr) const;
  MatrixSquareGeneral evalNvarIpas(double step,
                                   const VectorDouble& dir,
                                   const VectorDouble& center = VectorDouble(),
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
  VectorDouble evalPointToDbAsSP(const SpacePoint& p1,
                                 const std::vector<SpacePoint>& p2s,
                                 int ivar = 0,
                                 int jvar = 0,
                                 const CovCalcMode* mode = nullptr) const;
  double evalAverageDbToDb(const Db* db1,
                           const Db* db2,
                           int ivar = 0,
                           int jvar = 0,
                           const CovCalcMode* mode = nullptr) const;
  double evalAveragePointToDb(const SpacePoint& p1,
                              const Db* db2,
                              int ivar = 0,
                              int jvar = 0,
                              const CovCalcMode* mode = nullptr) const;
  double evalAveragePointToDbAsSP(const SpacePoint &p1,
                                  const std::vector<SpacePoint> &p2s,
                                  int ivar = 0,
                                  int jvar = 0,
                                  const CovCalcMode *mode = nullptr) const;
  MatrixRectangular evalCovMatrix(const Db* db1,
                                  const Db* db2 = nullptr,
                                  int ivar = 0,
                                  int jvar = 0,
                                  const VectorInt& nbgh1 = VectorInt(),
                                  const VectorInt& nbgh2 = VectorInt(),
                                  const CovCalcMode* mode = nullptr) const;
  VectorVectorDouble evalCovMatrixOptim(const Db *db1,
                                        const Db *db2 = nullptr,
                                        int ivar = 0,
                                        int jvar = 0,
                                        const CovCalcMode* mode = nullptr) const;
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

private:
  DbGrid* _discretizeBlock(const VectorDouble& ext,
                           const VectorInt& ndisc,
                           const VectorDouble& angles = VectorDouble(),
                           const VectorDouble& x0 = VectorDouble()) const;
  Db* _discretizeBlockRandom(const DbGrid* dbgrid, int seed = 34131) const;
  double _getVolume(const VectorDouble& ext) const;
};

/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Space/ASpaceObject.hpp"

#include "Matrix/MatrixSquareGeneral.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Space/SpacePoint.hpp"

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
  virtual double eval0(int ivar,
                       int jvar,
                       const CovCalcMode& mode = CovCalcMode()) const = 0;
  virtual double eval(int ivar,
                      int jvar,
                      const SpacePoint& p1,
                      const SpacePoint& p2,
                      const CovCalcMode& mode = CovCalcMode()) const = 0;
  virtual double evalCovOnSphere(double /*alpha*/,
                                 int /*degree*/,
                                 bool /*normalize*/) const { return TEST; }
  virtual double evalSpectrum(const VectorDouble& /*freq*/,
                              int /*ivar*/, int /*jvar*/) const { return TEST; }
  /////////////////////////////////////////////////////////////////////////////////

  MatrixSquareGeneral eval0Nvar(const CovCalcMode& mode = CovCalcMode()) const;
  VectorDouble eval(int ivar,
                    int jvar,
                    const std::vector<SpacePoint>& vec_p1,
                    const std::vector<SpacePoint>& vec_p2,
                    const CovCalcMode& mode = CovCalcMode()) const;
  MatrixSquareGeneral eval(const SpacePoint& p1,
                           const SpacePoint& p2,
                           const CovCalcMode& mode = CovCalcMode()) const;
  double evalIvarIpas(int ivar,
                      int jvar,
                      double step,
                      const VectorDouble& dir,
                      const VectorDouble& center = VectorDouble(),
                      const CovCalcMode& mode = CovCalcMode()) const;
  double evalIvarIpas(int ivar,
                      int jvar,
                      const VectorDouble &dincr,
                      const CovCalcMode &mode) const;
  VectorDouble evalIvarNpas(int ivar,
                            int jvar,
                            const VectorDouble& vec_step,
                            const VectorDouble& dir = VectorDouble(),
                            const VectorDouble& center = VectorDouble(),
                            const CovCalcMode& mode = CovCalcMode()) const;
  MatrixSquareGeneral evalNvarIpas(double step,
                                   const VectorDouble& dir,
                                   const VectorDouble& center = VectorDouble(),
                                   const CovCalcMode& mode = CovCalcMode()) const;
  MatrixSquareGeneral evalNvarIpas(const VectorDouble &dincr,
                                   const CovCalcMode &mode) const;
  double evalIsoIvarIpas(int ivar,
                         int jvar,
                         double step,
                         const CovCalcMode& mode = CovCalcMode()) const;
  VectorDouble evalIsoIvarNpas(int ivar,
                               int jvar,
                               const VectorDouble& vec_step,
                               const CovCalcMode& mode = CovCalcMode()) const;
  MatrixSquareGeneral evalIsoNvarIpas(double step,
                                      const CovCalcMode& mode = CovCalcMode()) const;

  double evalCvv(const VectorDouble& ext,
                 const VectorInt& ndisc,
                 const VectorDouble& angles = VectorDouble(),
                 int ivar = 0,
                 int jvar = 0,
                 const CovCalcMode& mode = CovCalcMode()) const;
  double evalCvvShift(const VectorDouble& ext,
                      const VectorInt& ndisc,
                      const VectorDouble& shift,
                      const VectorDouble& angles = VectorDouble(),
                      int ivar = 0,
                      int jvar = 0,
                      const CovCalcMode& mode = CovCalcMode()) const;
  MatrixSquareGeneral evalCvvM(const VectorDouble& ext,
                               const VectorInt& ndisc,
                               const VectorDouble& angles = VectorDouble(),
                               const CovCalcMode& mode = CovCalcMode()) const;
  double evalCxv(const SpacePoint& p1,
                 const VectorDouble& ext,
                 const VectorInt& ndisc,
                 const VectorDouble& angles = VectorDouble(),
                 const VectorDouble& x0 = VectorDouble(),
                 int ivar = 0,
                 int jvar = 0,
                 const CovCalcMode& mode = CovCalcMode()) const;
  double evalCxv(const Db* db,
                 const VectorDouble& ext,
                 const VectorInt& ndisc,
                 const VectorDouble& angles = VectorDouble(),
                 const VectorDouble& x0 = VectorDouble(),
                 int ivar = 0,
                 int jvar = 0,
                 const CovCalcMode& mode = CovCalcMode()) const;
  MatrixSquareGeneral evalCxvM(const SpacePoint& p1,
                               const VectorDouble& ext,
                               const VectorInt& ndisc,
                               const VectorDouble& angles = VectorDouble(),
                               const VectorDouble& x0 = VectorDouble(),
                               const CovCalcMode& mode = CovCalcMode()) const;
  VectorDouble evalPointToDb(const SpacePoint& p1,
                             const Db* db2,
                             int ivar = 0,
                             int jvar = 0,
                             bool useSel = true,
                             const CovCalcMode& mode = CovCalcMode()) const;
  double evalAverageDbToDb(const Db* db1,
                           const Db* db2,
                           int ivar = 0,
                           int jvar = 0,
                           const CovCalcMode& mode = CovCalcMode()) const;
  double evalAveragePointToDb(const SpacePoint& p1,
                              const Db* db2,
                              int ivar = 0,
                              int jvar = 0,
                              const CovCalcMode& mode = CovCalcMode()) const;
  MatrixRectangular evalCovMatrix(const Db* db1,
                                  const Db* db2 = nullptr,
                                  int ivar = 0,
                                  int jvar = 0,
                                  const CovCalcMode& mode = CovCalcMode()) const;

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

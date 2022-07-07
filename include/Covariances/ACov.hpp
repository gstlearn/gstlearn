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
#include "geoslib_define.h"

#include "Space/ASpaceObject.hpp"

#include "Matrix/MatrixSquareGeneral.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Space/SpacePoint.hpp"

class Db;
class DbGrid;

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
                 const CovCalcMode& mode = CovCalcMode());
  double evalCxv(const SpacePoint& p1,
                 const VectorDouble& est,
                 const VectorInt& ndisc,
                 const VectorDouble& angle = VectorDouble(),
                 int ivar = 0,
                 int jvar = 0,
                 const CovCalcMode& mode = CovCalcMode());
  VectorDouble evalPointToDb(const SpacePoint& p1,
                             const Db* db2,
                             int ivar = 0,
                             int jvar = 0,
                             bool useSel = false,
                             const CovCalcMode& mode = CovCalcMode());
  double evalAverageDbToDb(const Db* db1,
                           const Db* db2,
                           int ivar = 0,
                           int jvar = 0,
                           const CovCalcMode& mode = CovCalcMode());
  double evalAveragePointToDb(const SpacePoint& p1,
                              const Db* db2,
                              int ivar = 0,
                              int jvar = 0,
                              const CovCalcMode& mode = CovCalcMode());

private:
  DbGrid* _discretizeBlock(const VectorDouble& ext,
                           const VectorInt& ndisc,
                           const VectorDouble& angles = VectorDouble());
  Db* _discretizeBlockRandom(const DbGrid* dbgrid);
};

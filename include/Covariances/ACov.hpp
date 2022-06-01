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
  /////////////////////////////////////////////////////////////////////////////////

  virtual MatrixSquareGeneral eval0(const CovCalcMode& mode = CovCalcMode()) const;
  virtual VectorDouble eval(int ivar,
                            int jvar,
                            const std::vector<SpacePoint>& vec_p1,
                            const std::vector<SpacePoint>& vec_p2,
                            const CovCalcMode& mode = CovCalcMode()) const;
  virtual MatrixSquareGeneral eval(const SpacePoint& p1,
                                   const SpacePoint& p2,
                                   const CovCalcMode& mode = CovCalcMode()) const;

  /// Covariance from a given point (center) in a given direction (dir * step)
  /// for a pair of variables and a single step
  virtual double evalIvarIpas(int ivar,
                              int jvar,
                              double step,
                              const VectorDouble& dir,
                              const VectorDouble& center = VectorDouble(),
                              const CovCalcMode& mode = CovCalcMode()) const;
  /// Covariance vector from a given point (center) in a given direction (dir * steps)
  /// for a pair of variables and a set of steps
  virtual VectorDouble evalIvarNpas(int ivar,
                                    int jvar,
                                    const VectorDouble& vec_step,
                                    const VectorDouble& dir = VectorDouble(),
                                    const VectorDouble& center = VectorDouble(),
                                    const CovCalcMode& mode = CovCalcMode()) const;
  /// Covariance Matrix from a given point (center) in a given direction (dir * step)
  /// for a set of variables and a given step
  virtual MatrixSquareGeneral evalNvarIpas(double step,
                                           const VectorDouble& dir,
                                           const VectorDouble& center = VectorDouble(),
                                           const CovCalcMode& mode = CovCalcMode()) const;

  /// Covariance for a given unit global distance (without anisotropy)
  /// for a pair of variables and a single step
  virtual double evalIsoIvarIpas(int ivar,
                                 int jvar,
                                 double step,
                                 const CovCalcMode& mode = CovCalcMode()) const;
  /// Covariance for a given unit global distance (without anisotropy)
  /// for a pair of variables and a set of steps
  virtual VectorDouble evalIsoIvarNpas(int ivar,
                                       int jvar,
                                       const VectorDouble& vec_step,
                                       const CovCalcMode& mode = CovCalcMode()) const;
  /// Covariance for a given unit global distance (without anisotropy)
  /// for a set of variables and a single step
  virtual MatrixSquareGeneral evalIsoNvarIpas(double step,
                                              const CovCalcMode& mode = CovCalcMode()) const;

};

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

  /// Number of variables (used for covariance matrix size)
  virtual int getNVariables() const = 0;
  /// Covariance at the origin
  virtual double eval0(int ivar,
                       int jvar,
                       const CovCalcMode& mode = CovCalcMode()) const = 0;
  /// Covariance between two points for two given variables
  virtual double eval(int ivar,
                      int jvar,
                      const SpacePoint& p1,
                      const SpacePoint& p2,
                      const CovCalcMode& mode = CovCalcMode()) const = 0;
  /// Compatibility with gradient calculations
  virtual bool isGradientCompatible() const;

  /////////////////////////////////////////////////////////////////////////
  /// Convenient shortcut methods

  /// Covariance at the origin
  virtual MatrixSquareGeneral eval0(const CovCalcMode& mode = CovCalcMode()) const;

  /// Covariance between two points
  virtual VectorDouble eval(int ivar,
                            int jvar,
                            const std::vector<SpacePoint>& vec_p1,
                            const std::vector<SpacePoint>& vec_p2,
                            const CovCalcMode& mode = CovCalcMode()) const;
  virtual MatrixSquareGeneral eval(const SpacePoint& p1,
                                   const SpacePoint& p2,
                                   const CovCalcMode& mode = CovCalcMode()) const;

  /// Covariance from a given point (center) in a given direction (dir * step)
  virtual double eval(int ivar,
                      int jvar,
                      double step,
                      const VectorDouble& dir,
                      const VectorDouble& center = VectorDouble(),
                      const CovCalcMode& mode = CovCalcMode()) const;
  virtual VectorDouble eval(int ivar,
                            int jvar,
                            const VectorDouble& vec_step,
                            const VectorDouble& dir,
                            const VectorDouble& center = VectorDouble(),
                            const CovCalcMode& mode = CovCalcMode()) const;
  virtual MatrixSquareGeneral eval(double step,
                                   const VectorDouble& dir,
                                   const VectorDouble& center = VectorDouble(),
                                   const CovCalcMode& mode = CovCalcMode()) const;

  /// Covariance for a given unit global distance (without anisotropy)
  virtual double eval(int ivar,
                      int jvar,
                      double step,
                      const CovCalcMode& mode = CovCalcMode()) const;
  virtual VectorDouble eval(int ivar,
                            int jvar,
                            const VectorDouble& vec_step,
                            const CovCalcMode& mode = CovCalcMode()) const;
  virtual MatrixSquareGeneral eval(double step,
                                   const CovCalcMode& mode = CovCalcMode()) const;

};

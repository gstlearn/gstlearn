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

#include "Covariances/CovAniso.hpp"
#include "Enum/EConsElem.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"
#include "Basic/ICloneable.hpp"
#include "Covariances/CovContext.hpp"
#include "Space/SpacePoint.hpp"
#include <vector>


class ACov;
/**
 * \brief
 * This class describes an **elementary covariance**.
 *
 * This covariance is described through the following list of parameters:
 * - the covariance **type**: the list of these types is provided in ECov.hpp
 * - the largest set of parameters for any covariance: **range(s)**, **anisotropy angle(s)**, **third parameter**. Some of these parameters
 * do not make sense, depending on the covariance type: e.g. the range for nugget effect, the third parameter for a spherical
 * structure, ...
 * All these parameters are processed and stored as a **tensor** in order to avoid repetitive calculations.
 * - the **sill**. This comes as a square symmetric matrix whose dimension is equal to the number of variables.
 */
class GSTLEARN_EXPORT CovGneiting: public ACov, public ICloneable//, public ICloneable
{
public:
  CovGneiting(const CovAniso* covS, const CovAniso* covTemp, double separability = 1.0);
  CovGneiting(const CovGneiting& r);
  CovGneiting& operator=(const CovGneiting& r);
  virtual ~CovGneiting();
  IMPLEMENT_CLONING(CovGneiting)

  bool isConsistent(const ASpace* space) const override 
  {
    DECLARE_UNUSED(space)
    return true; 
  }
  /// ACov Interface
 
  virtual double eval(const SpacePoint& p1,
                      const SpacePoint& p2,
                      int ivar = 0,
                      int jvar = 0,
                      const CovCalcMode* mode = nullptr) const override;

  virtual int getNVariables() const override { return 1; }
  void optimizationSetTargetByIndex(int iech) const override;
protected:
    void _optimizationSetTarget(const SpacePoint &pt) const override;

private:
  void _optimizationPreProcess(const std::vector<SpacePoint>& p) const override;
  void _optimizationPostProcess() const override;

private:
  CovContext _ctxt;                    /// Context (space, number of variables, ...) // TODO : Really store a copy ?
  const CovAniso* _covS;
  const CovAniso* _covTemp;
  double _separability;
  mutable CovAniso _covSCopy;


};


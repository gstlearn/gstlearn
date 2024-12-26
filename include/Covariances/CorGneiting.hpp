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
 * This class describes the Gneiting correlation function.
 *
 */
class GSTLEARN_EXPORT CorGneiting: public ACor, public ICloneable//, public ICloneable
{
public:
  CorGneiting(const CovAniso* covS, const CovAniso* covTemp, double separability = 1.0);
  CorGneiting(const CorGneiting& r);
  CorGneiting& operator=(const CorGneiting& r);
  virtual ~CorGneiting();
  IMPLEMENT_CLONING(CorGneiting)

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
  void optimizationSetTargetByIndex(int iech) const;
protected:
    void _optimizationSetTarget(const SpacePoint &pt) const;

private:
  void _optimizationPreProcess(const std::vector<SpacePoint>& p) const;
  void optimizationPostProcess() const override;

private:
  CovContext _ctxt;                    /// Context (space, number of variables, ...) // TODO : Really store a copy ?
  const CovAniso* _covS;
  const CovAniso* _covTemp;
  double _separability;
  mutable CovAniso _covSCopy;


};


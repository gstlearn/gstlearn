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

#include "Covariances/CorAniso.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"
#include "Basic/ICloneable.hpp"
#include "Space/SpacePoint.hpp"
#include <vector>

class ACov;
class CorAniso;
/**
 * \brief
 * This class describes the Gneiting correlation function.
 *
 */
class GSTLEARN_EXPORT CorGneiting: public ACov, public ICloneable//, public ICloneable
{
public:
  CorGneiting(const CorAniso* covS, const CorAniso* covTemp, double separability = 1.0);
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

  virtual int getNVar() const override { return 1; }
  void optimizationSetTargetByIndex(int iech) const override;
protected:
    void _optimizationSetTarget(const SpacePoint &pt) const override;

private:
  void _optimizationPreProcess(int mode, const std::vector<SpacePoint>& p) const override;
  void _optimizationPostProcess() const override;

private:
  CovContext _ctxt;                    /// Context (space, number of variables, ...) // TODO : Really store a copy ?
  const CorAniso* _covS;
  const CorAniso* _covTemp;
  double _separability;
  mutable CorAniso _covSCopy;
};


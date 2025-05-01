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

#include "Basic/AStringable.hpp"
#include "Space/SpaceTarget.hpp"

class Db;
/**
 * This class performs the test between two Space Points (usually the target point first and the tentative point second)
 * The function isOK returns TRUE if the two space points are "comparable" with respect to the current criterion
 * Examples:
 * - if the criterion is based on FAULTS, the method isOK returns TRUE if the two space points
 * are not separated by any fault segment
 * - if the criterion is based on DISTANCE, the method isOK returns TRUE if the distance between the two space points
 * is smaller than a threshold distance (registered as a member of the BiTargetCheckDistance class)
 */
class GSTLEARN_EXPORT ABiTargetCheck: public AStringable, public ICloneable
{
public:
  ABiTargetCheck();
  ABiTargetCheck(const ABiTargetCheck& r);
  ABiTargetCheck& operator=(const ABiTargetCheck& r);
  virtual ~ABiTargetCheck();

  virtual bool isOK(const SpaceTarget &T1,
                    const SpaceTarget &T2) const = 0;

  // TODO : isValid is not const (should be renamed)
  virtual bool isValid(const Db* dbin, const Db* dbout) { DECLARE_UNUSED(dbin, dbout); return true; }
/*
#ifndef SWIG
  // Implement clone for permitting director feature
  virtual ICloneable* clone() const override { return nullptr; } ;
#endif
*/
};

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

#include "Basic/ICloneable.hpp"
#include "Basic/AStringable.hpp"
#include "Space/SpaceTarget.hpp"

class Db;

class GSTLEARN_EXPORT ABiTargetCheck: public AStringable//, public ICloneable
{
public:
  ABiTargetCheck();
  ABiTargetCheck(const ABiTargetCheck& r);
  ABiTargetCheck& operator=(const ABiTargetCheck& r);
  virtual ~ABiTargetCheck();

  virtual bool isOK(const SpaceTarget &T1,
                    const SpaceTarget &T2) const = 0;

  // TODO : isValid is not const (should be renamed)
  virtual bool isValid(const Db* dbin, const Db* dbout) { return true; }
/*
#ifndef SWIG
  // Implement clone for permitting director feature
  virtual ICloneable* clone() const override { return nullptr; } ;
#endif
*/
};

/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include <Geometry/ABiTargetCheck.hpp>
#include "gstlearn_export.hpp"

#include "Faults/Faults.hpp"

class GSTLEARN_EXPORT BiTargetCheckCode: public ABiTargetCheck
{
public:
  BiTargetCheckCode(int optcode=1, double tolcode=EPSILON6);
  BiTargetCheckCode(const BiTargetCheckCode& r);
  BiTargetCheckCode& operator=(const BiTargetCheckCode& r);
  virtual ~BiTargetCheckCode();

  /// ICloneable Interface
  IMPLEMENT_CLONING(BiTargetCheckCode)

  virtual bool isOK(const SpaceTarget &T1,
                    const SpaceTarget &T2) const override;

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int getOptCode() const { return _optCode; }
  void setOptCode(int optCode) { _optCode = optCode; }
  double getTolCode() const { return _tolCode; }
  void setTolCode(double tolCode) { _tolCode = tolCode; }

private:
  int _optCode;
  double _tolCode;
};

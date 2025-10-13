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
#include "LinearOp/IProj.hpp"

namespace gstlrn
{
class GSTLEARN_EXPORT TestInheritance: public AStringable
{
public:
  TestInheritance();
  TestInheritance(const TestInheritance& r)            = delete;
  TestInheritance& operator=(const TestInheritance& r) = delete;

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  void setIproj(IProj* ipr) { _iproj = ipr; }
  virtual ~TestInheritance();

private:
  IProj* _iproj;
};
} // namespace gstlrn
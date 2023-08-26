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

#include "LinearOp/IProjMatrix.hpp"
#include "Basic/AStringable.hpp"

class GSTLEARN_EXPORT TestInheritance : public AStringable
{
public:
  TestInheritance();
  TestInheritance(const TestInheritance& r) = delete;
  TestInheritance& operator=(const TestInheritance& r) = delete;

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  void setIproj(IProjMatrix* ipr){ _iproj = ipr;}
  virtual ~TestInheritance();

private:
 IProjMatrix* _iproj;
};

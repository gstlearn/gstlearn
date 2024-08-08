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
#include <assert.h>

/**
 * Inherits from this interface to make your class cloneable.
 * You must use IMPLEMENT_CLONING macro in concrete classes only.
 */
class GSTLEARN_EXPORT ICloneable
{
public:
  ICloneable() {};
  virtual ~ICloneable() {};

  virtual ICloneable* clone() const = 0;
};

// Thanks to here (macro way):
// https://alfps.wordpress.com/2010/06/12/cppx-3-ways-to-mix-in-a-generic-cloning-implementation/
#define IMPLEMENT_CLONING(Class)                   \
public:                                            \
  inline virtual Class* clone() const override     \
  {                                                \
    static_assert(                                 \
      ! std::is_abstract<Class>::value,            \
      "Class cannot be cloned as it is abstract"   \
    );                                             \
    assert(typeid(*this) == typeid(Class));        \
    return (new Class(*this));                     \
  }


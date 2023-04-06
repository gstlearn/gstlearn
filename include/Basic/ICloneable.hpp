/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include <typeinfo>
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
    assert(typeid(*this) == typeid(Class));        \
    return (new Class(*this));                     \
  }


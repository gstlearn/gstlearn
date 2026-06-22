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
#include <cassert>
#include <memory>

namespace gstlrn
{

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
#ifndef SWIG
    std::shared_ptr<ICloneable> cloneShared() const
    {
      return std::shared_ptr<ICloneable>(clone());
    }

    std::unique_ptr<ICloneable> cloneUnique() const
    {
      return std::unique_ptr<ICloneable>(clone());
    }
#endif
  };

// std::remove_cvref_t not defined in C++17, however std::decay_t provides a
// close-enough equivalent since this is not used here with functions or arrays
// https://devblogs.microsoft.com/cppblog/cpp17-20-features-and-fixes-in-vs-2019/
#ifdef USE_BOOST_SPAN
#define REMOVE_CVREF_T std::decay_t
#else
#define REMOVE_CVREF_T std::remove_cvref_t
#endif

// Thanks to here (macro way):
// https://alfps.wordpress.com/2010/06/12/cppx-3-ways-to-mix-in-a-generic-cloning-implementation/
#define IMPLEMENT_CLONING(Class)                                               \
                                                                               \
public:                                                                        \
  inline virtual Class* clone() const override                                 \
  {                                                                            \
    static_assert(                                                             \
      !std::is_abstract<Class>::value,                                         \
      "Class cannot be cloned as it is abstract");                             \
    static_assert(std::is_base_of_v<Class, REMOVE_CVREF_T<decltype(*this)>>);  \
    return (new Class(*this));                                                 \
  }
} // namespace gstlrn

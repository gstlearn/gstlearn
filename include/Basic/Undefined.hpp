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

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include <cmath>

namespace gstlrn
{
  GSTLEARN_EXPORT bool FFFF(double value); // TODO isNA<double>
  GSTLEARN_EXPORT bool IFFFF(Id value); // TODO isNA<Id>
  GSTLEARN_EXPORT double getTEST(); // TODO getNA<double>
  GSTLEARN_EXPORT Id getITEST(); // TODO getNA<Id>

// No need this stuff through SWIG (because we use target language NAs)
#ifndef SWIG

#define DOUBLE_NA TEST
#define INT_NA ITEST
#define STRING_NA "NA"
#define FLOAT_NA                                                               \
  static_cast<float>(                                                          \
    TEST) // 1.234e30 is ok for 4 bytes but needs a cast for Windows
#define UCHAR_NA 255 // Warning: not to be used

  template<typename T>
  inline T getNA();

  template<>
  inline bool getNA()
  {
    return false;
  }

  template<>
  inline double getNA()
  {
    return DOUBLE_NA;
  }

  template<>
  inline Id getNA()
  {
    return INT_NA;
  }

  template<>
  inline int getNA()
  {
    return INT_NA;
  }

  template<>
  inline long getNA()
  {
    return INT_NA;
  }

  template<>
  inline String getNA()
  {
    return STRING_NA;
  }

  template<>
  inline float getNA()
  {
    return FLOAT_NA;
  }

  template<>
  inline unsigned char getNA()
  {
    return UCHAR_NA;
  }

  template<typename T>
  inline bool isNA(const T& v)
  {
    if constexpr (std::is_arithmetic_v<T>)
    {
      // fallback pour les scalaires non spécialisés (rare)
      return v == getNA<T>();
    }
    else
    {
      // Pour tout type non scalaire : jamais NA
      return false;
    }
  }

  template<>
  inline bool isNA(const double& v)
  {
    return (v == getNA<double>() || std::isnan(v) || std::isinf(v));
  }

  template<>
  inline bool isNA(const Id& v)
  {
    return (v == getNA<Id>());
  }

  template<>
  inline bool isNA(const String& v)
  {
    return (v == getNA<String>());
  }

  template<>
  inline bool isNA(const float& v)
  {
    return (v == getNA<float>() || std::isnan(v) || std::isinf(v));
  }

  template<>
  inline bool isNA(const unsigned char& v)
  {
    return (v == getNA<unsigned char>());
  }

#endif // SWIG
}; // namespace gstlrn

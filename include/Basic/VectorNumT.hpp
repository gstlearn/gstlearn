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

#include "Basic/Undefined.hpp"
#include "Basic/VectorT.hpp"
#include "geoslib_define.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace gstlrn
{
/***************************************************************************
 **
 ** Vector of T values (numerical type).
 ** T type must define copy constructor and assignment operator
 ** T type must override numerical operators (+, -, *, /)
 **
 ***************************************************************************/
template<typename T>
class VectorNumT: public VectorT<T>
{
public:
  typedef VectorT<T> Parent;
  typedef std::vector<T> Vector;
  typedef typename Vector::value_type value_type;
  typedef typename Vector::size_type size_type;
  typedef typename Vector::iterator iterator;
  typedef typename Vector::const_iterator const_iterator;
  typedef typename Vector::reverse_iterator reverse_iterator;
  typedef typename Vector::const_reverse_iterator const_reverse_iterator;

public:
  inline VectorNumT()
    : Parent()
  {
  }
  inline VectorNumT(const Vector& vec)
    : Parent(vec)
  {
  }
  inline VectorNumT(size_type count, const T& value = {})
    : Parent(count, value)
  {
  }
  template<class InputIt>
  inline VectorNumT(InputIt first, InputIt last)
    : Parent(first, last)
  {
  }
  inline VectorNumT(const VectorNumT& other) = default;
#ifndef SWIG
  inline VectorNumT(std::initializer_list<T> init)
    : Parent(init)
  {
  }
#endif
  inline ~VectorNumT() = default;

public:
  inline bool isEqual(const VectorNumT& other, double eps = 1.e-10) const;
  inline bool isConstant();

  inline T sum() const;
  inline T minimum(bool flagAbs = false) const;
  inline T maximum(bool flagAbs = false) const;
  inline double mean() const;
  inline double norm(Id size = 0) const;
  inline double norm2(Id size = 0) const;
  inline double normL1() const;
  inline double normInf() const;
  inline T prod() const;

  inline void identify() const;
  inline void dump(const String& title) const;

  inline double innerProduct(const VectorNumT<T>& v, Id size = 0) const;

  /**
   * \defgroup VectorNumT: Operation performed on Vector of Numerical values
   *
   **/

  /** @addtogroup Func_V
   * \ingroup VectorNumT
   *
   * Syntax: const func(vec)
   *
   * Action: 'this' = 'this' %func% 'vec'
   *  @{
   */
  inline const VectorNumT<T>& add(const VectorNumT<T>& v);
  inline const VectorNumT<T>& subtract(const VectorNumT<T>& v);
  inline const VectorNumT<T>& multiply(const VectorNumT<T>& v);
  inline const VectorNumT<T>& divide(const VectorNumT<T>& v);
  /**@}*/

  /** @addtogroup Func_Cst
   * \ingroup VectorNumT
   *
   * Syntax: const func(val)
   *
   * Action: 'this' = 'this' %func% 'val'
   *  @{
   */
  inline const VectorNumT<T>& addCst(const T& v);
  inline const VectorNumT<T>& subtractCst(const T& v);
  inline const VectorNumT<T>& multiplyCst(const T& v);
  inline const VectorNumT<T>& divideCst(const T& v);
  /**@}*/

  inline void normalizeInPlace(Id normType = 2);
};

template<typename T>
bool VectorNumT<T>::isEqual(const VectorNumT& other, double eps) const
{
  if (other.size() != VectorNumT::size()) return false;
  for (size_type i = 0, n = VectorNumT::size(); i < n; i++)
    if (abs(VectorNumT::at(i) - other.at(i)) > eps) return false;
  return true;
}

template<typename T>
bool VectorNumT<T>::isConstant()
{
  if (VectorNumT::size() <= 0) return false;
  T refval = VectorNumT::_v.at(0);
  for (auto v: VectorNumT::_v)
  {
    if (v != refval) return false;
  }
  return true;
}

template<typename T>
T VectorNumT<T>::sum() const
{
  if (VectorNumT::size() <= 0) return T();
  T sum = 0;
  for (auto v: VectorNumT::_v)
  {
    if (isNA(v)) continue;
    sum += v;
  }
  return (sum);
}

template<typename T>
T VectorNumT<T>::maximum(bool flagAbs) const // Prevent using max and min keywords (Visual)
{
  if (VectorNumT::size() <= 0) return 0;
  T mymax = (std::numeric_limits<T>::min)(); // https://stackoverflow.com/a/27443191/3952924
  for (auto v: VectorNumT::_v)
  {
    if (isNA(v)) continue;
    if (flagAbs) v = ABS(v);
    if (v > mymax) mymax = v;
  }
  return (mymax);
}

template<typename T>
T VectorNumT<T>::minimum(bool flagAbs) const // Prevent using max and min keywords (Visual)
{
  if (VectorNumT::size() <= 0) return 0;
  T mymin = (std::numeric_limits<T>::max)(); // https://stackoverflow.com/a/27443191/3952924
  for (auto v: VectorNumT::_v)
  {
    if (isNA(v)) continue;
    if (flagAbs) v = ABS(v);
    if (v < mymin) mymin = v;
  }
  return (mymin);
}

template<typename T>
double VectorNumT<T>::mean() const
{
  if (VectorNumT::size() <= 0) return static_cast<T>(NAN);
  double mean   = 0.;
  double number = 0.;
  for (auto v: VectorNumT::_v)
  {
    if (isNA(v)) continue;
    mean += v;
    number += 1.;
  }
  if (number > 0.)
    mean /= number;
  else
    mean = TEST;
  return (mean);
}

template<typename T>
double VectorNumT<T>::norm(Id size) const
{
  double ip = innerProduct(*this, size);
  return sqrt(ip);
}

template<typename T>
double VectorNumT<T>::norm2(Id size) const
{
  double ip = innerProduct(*this, size);
  return ip;
}

template<typename T>
double VectorNumT<T>::normL1() const
{
  double normL1 = 0.;
  for (auto v: VectorNumT::_v)
  {
    if (isNA(v)) continue;
    T value = ABS(v);
    normL1 += value;
  }
  return normL1;
}

template<typename T>
double VectorNumT<T>::normInf() const
{
  double norminf = 0.;
  for (auto v: VectorNumT::_v)
  {
    if (isNA(v)) continue;
    T value = ABS(v);
    if (value > norminf) norminf = value;
  }
  return norminf;
}

template<typename T>
T VectorNumT<T>::prod() const
{
  if (VectorNumT::size() <= 0) return 0;
  T prod = 1;
  for (auto v: VectorNumT::_v)
  {
    if (isNA(v)) continue;
    prod *= v;
  }
  return prod;
}

template<typename T>
double VectorNumT<T>::innerProduct(const VectorNumT<T>& v, Id size) const
{
  if (size <= 0) size = v.size();
  if (size != static_cast<Id>(VectorNumT::size())) return 0.;
  double prod = 0.;
  for (size_type i = 0, n = VectorNumT::size(); i < n; i++)
  {
    auto v1 = static_cast<double>(VectorNumT::at(i));
    auto v2 = static_cast<double>(v[i]);
    if (isNA(v1) || isNA(v2)) continue;
    prod += v1 * v2;
  }
  return prod;
}

template<typename T>
const VectorNumT<T>& VectorNumT<T>::add(const VectorNumT<T>& v)
{
  if (v.size() != VectorNumT::size()) return *this;
  for (size_type i = 0, n = VectorNumT::size(); i < n; i++)
  {
    auto v1 = VectorNumT::at(i);
    if (isNA(v1)) continue;
    if (isNA(v[i])) continue;
    VectorNumT::operator[](i) = v1 + v[i];
  }
  return *this;
}

template<typename T>
const VectorNumT<T>& VectorNumT<T>::subtract(const VectorNumT<T>& v)
{
  if (v.size() != VectorNumT::size()) return *this;
  for (size_type i = 0, n = VectorNumT::size(); i < n; i++)
  {
    auto v1 = VectorNumT::at(i);
    if (isNA(v1)) continue;
    if (isNA(v[i])) continue;
    VectorNumT::operator[](i) = v1 - v[i];
  }
  return *this;
}

template<typename T>
const VectorNumT<T>& VectorNumT<T>::multiply(const VectorNumT<T>& v)
{
  if (v.size() != VectorNumT::size()) return *this;
  for (size_type i = 0, n = VectorNumT::size(); i < n; i++)
  {
    auto v1 = VectorNumT::at(i);
    if (isNA(v1)) continue;
    if (isNA(v[i])) continue;
    VectorNumT::operator[](i) = v1 * v[i];
  }
  return *this;
}

template<typename T>
const VectorNumT<T>& VectorNumT<T>::divide(const VectorNumT<T>& v)
{
  if (v.size() != VectorNumT::size()) return *this;
  for (size_type i = 0, n = VectorNumT::size(); i < n; i++)
  {
    auto v1 = VectorNumT::at(i);
    if (isNA(v1)) continue;
    if (isNA(v[i])) continue;
    if (abs(v[i]) < 1.e-10) continue;
    VectorNumT::operator[](i) = v1 / v[i];
  }
  return *this;
}

template<typename T>
void VectorNumT<T>::normalizeInPlace(Id normType)
{
  double normValue;
  if (normType == 1)
    normValue = normL1();
  else
    normValue = norm();
  for (size_type i = 0, n = VectorNumT::size(); i < n; i++)
  {
    auto v1 = VectorNumT::at(i);
    if (isNA(v1)) continue;
    VectorNumT::operator[](i) = v1 / normValue;
  }
}

template<typename T>
const VectorNumT<T>& VectorNumT<T>::addCst(const T& v)
{
  std::for_each(VectorNumT::begin(), VectorNumT::end(), [v](T& d)
                { if (!isNA(d)) d += v; });
  return *this;
}

template<typename T>
const VectorNumT<T>& VectorNumT<T>::subtractCst(const T& v)
{
  std::for_each(VectorNumT::begin(), VectorNumT::end(), [v](T& d)
                { if (!isNA(d)) d -= v; });
  return *this;
}

template<typename T>
const VectorNumT<T>& VectorNumT<T>::multiplyCst(const T& v)
{
  std::for_each(VectorNumT::begin(), VectorNumT::end(), [v](T& d)
                { if (!isNA(d)) d *= v; });
  return *this;
}

template<typename T>
const VectorNumT<T>& VectorNumT<T>::divideCst(const T& v)
{
  if (abs(v) < 1.e-10)
    throw("VectorNumT<T>::divide: division by 0");
  std::for_each(VectorNumT::begin(), VectorNumT::end(), [v](T& d)
                { if (!isNA(d)) d /= v; });
  return *this;
}

/**
 * @brief Identify the VectorNumT
 *
 * @tparam T Input argument (template)
 */
template<typename T>
void VectorNumT<T>::identify() const
{
  if (VectorNumT::size() <= 0) return;

  if constexpr (std::is_same_v<T, Id>)
    std::cout << "--> This is a Vector of <Id>\n";
  else if constexpr (std::is_same_v<T, int>)
    std::cout << "--> This is a Vector of <int>\n";
  else if constexpr (std::is_same_v<T, long>)
    std::cout << "--> This is a Vector of <long>\n";
  else if constexpr (std::is_same_v<T, double>)
    std::cout << "--> This is a Vector of <double>\n";
  else
    std::cout << "??? This is a vector of unexpected type\n";
}

template<typename T>
void VectorNumT<T>::dump(const String& title) const
{
  if (VectorNumT::size() <= 0) return;

  if constexpr (std::is_same_v<T, double>)
  {
    std::cout << "Array of Double: " << title << std::endl;
    VectorT<T>::display();
  }
  else if constexpr (std::is_same_v<T, Id>)
  {
    std::cout << "Array of Id: " << title << std::endl;
    VectorT<T>::display();
  }
  else
  {
    std::cout << "Unknown type" << title << std::endl;
  }
}

#ifndef SWIG
template<typename T>
std::ostream& operator<<(std::ostream& os, const VectorT<VectorNumT<T>>& vec)
{
  os << "[";
  for (Id i = 0, n = static_cast<Id>(vec.size()); i < n; i++)
  {
    os << vec.at(i).toString();
    if (i != n - 1) os << " ";
  }
  os << "]";
  return os;
}
#endif

typedef VectorNumT<Id> VectorInt;
typedef VectorNumT<double> VectorDouble;
typedef VectorNumT<float> VectorFloat;
typedef VectorNumT<UChar> VectorUChar; // Use typedef because swig doesn't like 'unsigned char' in two words
typedef VectorT<VectorInt> VectorVectorInt;
typedef VectorT<VectorDouble> VectorVectorDouble;
typedef VectorT<VectorFloat> VectorVectorFloat;
} // namespace gstlrn

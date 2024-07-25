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

#include "VectorT.hpp"

#include <vector>
#include <limits>
#include <algorithm>
#include <cmath>

/***************************************************************************
 **
 ** Vector of T values (numerical type).
 ** T type must define copy constructor and assignment operator
 ** T type must override numerical operators (+, -, *, /)
 **
 ***************************************************************************/
template <typename T>
class VectorNumT : public VectorT<T>
{
public:
  typedef VectorT<T> Parent;
  typedef std::vector<T> Vector;
  typedef typename Vector::value_type               value_type;
  typedef typename Vector::size_type                size_type;
  typedef typename Vector::iterator                 iterator;
  typedef typename Vector::const_iterator           const_iterator;
  typedef typename Vector::reverse_iterator         reverse_iterator;
  typedef typename Vector::const_reverse_iterator   const_reverse_iterator;

public:
  inline VectorNumT()                                              : Parent() { }
  inline VectorNumT(const Vector& vec)                             : Parent(vec) { }
  inline VectorNumT(size_type count, const T& value = T())         : Parent(count, value) { }
  template< class InputIt >
  inline VectorNumT(InputIt first, InputIt last)                   : Parent(first, last) { }
  inline VectorNumT(const VectorNumT& other) = default;
#ifndef SWIG
  inline VectorNumT(std::initializer_list<T> init)                 : Parent(init) { }
#endif
  inline ~VectorNumT() = default;

// Only for C++ users
// These functions are not available in target language
// because numerical vectors are converted in target language vectors
public:
  inline bool isSame(const VectorNumT& other, double eps = 1.e-10) const;

  inline T sum() const;
  inline T minimum() const;
  inline T maximum() const;
  inline double mean() const;
  inline double norm() const;

  inline double innerProduct(const VectorNumT<T>& v) const;

  inline const VectorNumT<T>& add(const VectorNumT<T>& v);
  inline const VectorNumT<T>& subtract(const VectorNumT<T>& v);
  inline const VectorNumT<T>& multiply(const VectorNumT<T>& v);
  inline const VectorNumT<T>& divide(const VectorNumT<T>& v);

  inline const VectorNumT<T>& add(const T& v);
  inline const VectorNumT<T>& subtract(const T& v);
  inline const VectorNumT<T>& multiply(const T& v);
  inline const VectorNumT<T>& divide(const T& v);
};

template <typename T>
bool VectorNumT<T>::isSame(const VectorNumT& other,
                           double eps) const
{
  if (other.size() != VectorNumT::size()) return false;
  for (size_type i = 0, n = VectorNumT::size(); i < n; i++)
    if (abs(VectorNumT::at(i) - other.at(i)) > eps) return false;
  return true;
}

template <typename T>
T VectorNumT<T>::sum() const
{
  if (VectorNumT::size() <= 0) return T();
  T sum = 0;
  for (size_type i = 0, n = VectorNumT::size(); i < n; i++)
    sum += VectorNumT::_v->at(i);
  return (sum);
}

template <typename T>
T VectorNumT<T>::maximum() const  // Prevent using max and min keywords (Visual)
{
  if (VectorNumT::size() <= 0) return 0;
  T mymax = (std::numeric_limits<T>::min)(); // https://stackoverflow.com/a/27443191/3952924
  for (auto v : *VectorNumT::_v)
    if (v > mymax) mymax = v;
  return (mymax);
}

template <typename T>
T VectorNumT<T>::minimum() const  // Prevent using max and min keywords (Visual)
{
  if (VectorNumT::size() <= 0) return 0;
  T mymin = (std::numeric_limits<T>::max)(); // https://stackoverflow.com/a/27443191/3952924
  for (auto v : *VectorNumT::_v)
    if (v < mymin) mymin = v;
  return (mymin);
}

template <typename T>
double VectorNumT<T>::mean() const
{
  if (VectorNumT::size() <= 0) return static_cast<T>(NAN);
  double s = static_cast<double>(sum());
  return (s / static_cast<double>(VectorNumT::_v->size()));
}

template <typename T>
double VectorNumT<T>::norm() const
{
  double ip = innerProduct(*this);
  return sqrt(ip);
}

template <typename T>
double VectorNumT<T>::innerProduct(const VectorNumT<T>& v) const
{
  if (v.size() != VectorNumT::size())
    throw("VectorNumT<T>::innerProduct: Wrong size");
  double prod = 0.;
  for (size_type i = 0, n = VectorNumT::size(); i < n; i++)
    prod += static_cast<double>(VectorNumT::at(i)) * static_cast<double>(v[i]);
  return prod;
}

template <typename T>
const VectorNumT<T>& VectorNumT<T>::add(const VectorNumT<T>& v)
{
  if (v.size() != VectorNumT::size())
    throw("VectorNumT<T>::add: Wrong size");
  for (size_type i = 0, n = VectorNumT::size(); i < n; i++)
    VectorNumT::operator[](i) = VectorNumT::at(i)+v[i];
  return *this;
}

template <typename T>
const VectorNumT<T>& VectorNumT<T>::subtract(const VectorNumT<T>& v)
{
  if (v.size() != VectorNumT::size())
    throw("VectorNumT<T>::subtract: Wrong size");
  for (size_type i = 0, n = VectorNumT::size(); i < n; i++)
    VectorNumT::operator[](i) = VectorNumT::at(i)-v[i];
  return *this;
}

template <typename T>
const VectorNumT<T>& VectorNumT<T>::multiply(const VectorNumT<T>& v)
{
  if (v.size() != VectorNumT::size())
    throw("VectorNumT<T>::multiply: Wrong size");
  for (size_type i = 0, n = VectorNumT::size(); i < n; i++)
    VectorNumT::operator[](i) = VectorNumT::at(i)*v[i];
  return *this;
}

template <typename T>
const VectorNumT<T>& VectorNumT<T>::divide(const VectorNumT<T>& v)
{
  if (v.size() != VectorNumT::size())
    throw("VectorNumT<T>::divide: Wrong size");
  for (size_type i = 0, n = VectorNumT::size(); i < n; i++)
  {
    if (abs(v[i]) < 1.e-10)
      throw("VectorNumT<T>::divide: division by 0");
    VectorNumT::operator[](i) = VectorNumT::at(i)/v[i];
  }
  return *this;
}

template <typename T>
const VectorNumT<T>& VectorNumT<T>::add(const T& v)
{
  std::for_each(VectorNumT::begin(), VectorNumT::end(), [v](T& d) { d += v;});
  return *this;
}

template <typename T>
const VectorNumT<T>& VectorNumT<T>::subtract(const T& v)
{
  std::for_each(VectorNumT::begin(), VectorNumT::end(), [v](T& d) { d -= v;});
  return *this;
}

template <typename T>
const VectorNumT<T>& VectorNumT<T>::multiply(const T& v)
{
  std::for_each(VectorNumT::begin(), VectorNumT::end(), [v](T& d) { d *= v;});
  return *this;
}

template <typename T>
const VectorNumT<T>& VectorNumT<T>::divide(const T& v)
{
  if (abs(v) < 1.e-10)
    throw("VectorNumT<T>::divide: division by 0");
  std::for_each(VectorNumT::begin(), VectorNumT::end(), [v](T& d) { d /= v;});
  return *this;
}

template <typename T>
std::ostream& operator<<(std::ostream& os,
                         const VectorT<VectorNumT<T>>& vec)
{
  os << "[";
  for (int i = 0, n = (int)vec.size(); i < n; i++)
  {
    os << vec.at(i).toString();
    if (i != n - 1) os << " ";
  }
  os << "]";
  return os;
}

typedef VectorNumT<int>       VectorInt;
typedef VectorNumT<double>    VectorDouble;
typedef VectorNumT<float>     VectorFloat;
typedef VectorNumT<UChar>     VectorUChar; // Use typedef because swig doesn't like 'unsigned char' in two words
typedef VectorT<VectorInt>    VectorVectorInt;
typedef VectorT<VectorDouble> VectorVectorDouble;
typedef VectorT<VectorFloat>  VectorVectorFloat;

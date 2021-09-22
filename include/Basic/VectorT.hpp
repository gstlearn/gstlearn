/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "Basic/AStringable.hpp"
#include "Basic/AException.hpp"

#include "Basic/Vector.hpp"

#include <vector>
#include <sstream>
#include <memory>
#include <limits>
#include <algorithm>
#include <cmath>

/***************************************************************************
 **
 ** Numerical vector of T values (int, double, float...).
 ** T type must override numerical operators +*-/
 **
 ***************************************************************************/
template <typename T>
class VectorT : public AStringable
{
public:
  typedef std::vector<T> Vector;

  typedef typename Vector::size_type                size_type;
  typedef typename Vector::iterator                 iterator;
  typedef typename Vector::const_iterator           const_iterator;
  typedef typename Vector::reverse_iterator         reverse_iterator;
  typedef typename Vector::const_reverse_iterator   const_reverse_iterator;

public:
  inline VectorT()                                                 : _v(std::make_shared<Vector>()) { }
  inline VectorT(const Vector& vec)                                : _v(std::make_shared<Vector>(vec)) { }
  inline VectorT(size_type count, const T& value = T())            : _v(std::make_shared<Vector>(count, value)) { }
  inline VectorT(const VectorT& other)                             : _v(other._v) { }
#ifndef SWIG
  inline VectorT(std::initializer_list<T> init)                    : _v(std::make_shared<Vector>(init)) { }
  inline VectorT(VectorT&& other)                                  { _v.swap(other._v); }
#endif
  inline ~VectorT()                                                { }

#ifndef SWIG
  inline operator const Vector&() const                            { return *_v; } /// TODO : while we mix both
#endif

  inline Vector& getVector() const                                 { return *_v; } /// TODO : while we mix both

  inline VectorT& operator=(const Vector& vec)                     { _detach(); *_v = vec; return (*this); }
  inline VectorT& operator=(const VectorT& other)                  { _detach(); _v = other._v; return (*this); }
  inline VectorT& operator=(VectorT&& other)                       { _v.swap(other._v); return (*this); }
  inline VectorT& operator=(std::initializer_list<T> init)         { _detach(); (*_v) = init; return (*this); }

  inline bool operator==(const VectorT& other) const               { return *_v == *other._v; }
  inline bool operator!=(const VectorT& other) const               { return *_v != *other._v; }
  inline bool operator <(const VectorT& other) const               { return *_v  < *other._v; }
  inline bool operator<=(const VectorT& other) const               { return *_v <= *other._v; }
  inline bool operator >(const VectorT& other) const               { return *_v  > *other._v; }
  inline bool operator>=(const VectorT& other) const               { return *_v >= *other._v; }

  inline const T& at(size_type pos) const                          { if (pos >= _v->size()) my_throw("VectorT<T>::at: index out of range");                    return _v->operator[](pos); }
#ifndef SWIG
  inline T& operator[](size_type pos)                              { if (pos >= _v->size()) my_throw("VectorT<T>::operator[]: index out of range"); _detach(); return _v->operator[](pos); }
  inline const T& operator[](size_type pos) const                  { if (pos >= _v->size()) my_throw("VectorT<T>::operator[]: index out of range");            return _v->operator[](pos); }
#endif
  inline T& front()                                                { _detach(); return _v->front(); }
  inline const T& front() const                                    { return _v->front(); }
  inline T& back()                                                 { _detach(); return _v->back(); }
  inline const T& back() const                                     { return _v->back(); }
  inline T* data()                                                 { _detach(); return _v->data(); }
  inline const T* data() const                                     { return _v->data(); }
  inline const T* constData() const                                { return _v->data(); }

  inline bool empty() const                                        { return _v->empty(); }
  inline size_type size() const                                    { return _v->size(); }
  inline void reserve(size_type new_cap)                           { _v->reserve(new_cap); }
  inline size_type capacity() const                                { return _v->capacity(); }
  inline void clear()                                              { _detach(); _v->clear(); }

  inline void insert(size_type i, const T& value)                  { _detach(); _v->insert(begin() + i, value); }
  inline void insert(size_type i, size_type count, const T& value) { _detach(); _v->insert(begin() + i, count, value); }
  inline void remove(size_type i)                                  { _detach(); _v->erase(begin() + i); }
  inline void remove(size_type i, size_type count)                 { _detach(); _v->erase(begin() + i, begin() + i + count); }

  inline void push_back(const T& value)                            { _detach(); _v->push_back(value); }
  inline void push_back(const T&& value)                           { _detach(); _v->push_back(value); }
  inline void push_front(const T& value)                           { _detach(); _v->insert(begin(), value); }
  inline void push_front(const T&& value)                          { _detach(); _v->insert(begin(), value); }
  inline void resize(size_type count)                              { if (count == size()) return; _detach(); _v->resize(count); }
  inline void resize(size_type count, const T& value)              { if (count == size()) return; _detach(); _v->resize(count, value); }

  inline iterator begin()                                          { _detach(); return _v->begin(); }
  inline const_iterator begin() const                              { return _v->begin(); }
  inline const_iterator cbegin() const                             { return _v->cbegin(); }
  inline iterator end()                                            { _detach(); return _v->end(); }
  inline const_iterator end() const                                { return _v->end(); }
  inline const_iterator cend() const                               { return _v->cend(); }

  inline reverse_iterator rbegin()                                 { _detach(); return _v->rbegin(); }
  inline const_reverse_iterator crbegin() const                    { return _v->crbegin(); }
  inline reverse_iterator rend()                                   { _detach(); return _v->rend(); }
  inline const_reverse_iterator crend() const                      { return _v->crend(); }

#ifndef SWIG
  inline VectorT& operator<<(const T& value)                       { _detach(); _v->push_back(value); return (*this); }
  inline VectorT& operator<<(const VectorT<T>& v)                  { _detach(); reserve(size() + v.size()); std::for_each(v.cbegin(), v.cend(), [&](const T& value) { _v->push_back(value); }); return (*this); }
#endif

  inline void swap(VectorT& other);
  inline bool contains(const T& value) const;
  inline void fill(const T& value, size_type size = -1);
  inline void assign(const T* tab, size_type size);
  inline void set(size_type i, const T& value);
  inline bool isSame(const VectorT& v, double eps = EPSILON10) const;

  inline T sum() const;
  inline T min() const;
  inline T max() const;
  inline double mean() const;
  inline double norm() const;

  inline T innerProduct(const VectorT<T>& v) const;

  inline const VectorT<T>& add(const VectorT<T>& v);
  inline const VectorT<T>& subtract(const VectorT<T>& v);
  inline const VectorT<T>& multiply(const VectorT<T>& v);
  inline const VectorT<T>& divide(const VectorT<T>& v);

  inline const VectorT<T>& add(const T& v);
  inline const VectorT<T>& subtract(const T& v);
  inline const VectorT<T>& multiply(const T& v);
  inline const VectorT<T>& divide(const T& v);

  /*! Conversion to a string */
  inline std::string toString() const override;

private:
  inline void _detach();

  std::shared_ptr<Vector> _v;
};

template <typename T>
void VectorT<T>::swap(VectorT& other)
{
  std::swap(_v, other._v);
}

template <typename T>
bool VectorT<T>::contains(const T& value) const
{
  return (std::find(cbegin(), cend(), value) != cend());
}

template <typename T>
void VectorT<T>::fill(const T& value, size_type size)
{
  _detach();
  _v->resize(size);
  std::fill(begin(), end(), value);
}

template <typename T>
void VectorT<T>::assign(const T* tab, size_type size)
{
  _detach();
  _v->assign(tab, tab+size);
}

template <typename T>
void VectorT<T>::set(size_type i, const T& value)
{
  _detach();
  operator[](i) = value;
}

template <typename T>
bool VectorT<T>::isSame(const VectorT& other,
                        double eps) const
{
  if (other.size() != size()) return false;
  for (size_type i = 0, n = size(); i < n; i++)
    if (ABS(at(i) - other.at(i)) > eps) return false;
  return true;
}

template <typename T>
T VectorT<T>::sum() const
{
  if (size() <= 0) return T();
  T sum = 0.; /// TODO : always possible ?
  for (size_type i = 0, n = size(); i < n; i++)
    sum += _v->at(i);
  return (sum);
}

template <typename T>
T VectorT<T>::max() const
{
  if (size() <= 0) return 0.;
  T max = std::numeric_limits<T>::min();
  for (auto v : *_v)
    if (v > max) max = v;
  return (max);
}

template <typename T>
T VectorT<T>::min() const
{
  if (size() <= 0) return 0.;
  T min = std::numeric_limits<T>::max();
  for (auto v : *_v)
    if (v < min) min = v;
  return (min);
}

template <typename T>
double VectorT<T>::mean() const
{
  if (size() <= 0) return static_cast<T>(TEST);
  double s = static_cast<double>(sum());
  return (s / static_cast<double>(_v->size()));
}

template <typename T>
double VectorT<T>::norm() const
{
  double ip = static_cast<double>(innerProduct(*this));
  return sqrt(ip);
}

template <typename T>
T VectorT<T>::innerProduct(const VectorT<T>& v) const
{
  if (v.size() != size())
    my_throw("VectorT<T>::innerProduct: Wrong size");
  T prod = 0.; /// TODO : always possible ?
  for (size_type i = 0, n = size(); i < n; i++)
    prod += at(i) * v[i];
  return prod;
}

template <typename T>
const VectorT<T>& VectorT<T>::add(const VectorT<T>& v)
{
  if (v.size() != size())
    my_throw("VectorT<T>::add: Wrong size");
  for (size_type i = 0, n = size(); i < n; i++)
    set(i, at(i)+v[i]);
  return *this;
}

template <typename T>
const VectorT<T>& VectorT<T>::subtract(const VectorT<T>& v)
{
  if (v.size() != size())
    my_throw("VectorT<T>::subtract: Wrong size");
  for (size_type i = 0, n = size(); i < n; i++)
    set(i, at(i)-v[i]);
  return *this;
}

template <typename T>
const VectorT<T>& VectorT<T>::multiply(const VectorT<T>& v)
{
  if (v.size() != size())
    my_throw("VectorT<T>::multiply: Wrong size");
  for (size_type i = 0, n = size(); i < n; i++)
    set(i, at(i)*v[i]);
  return *this;
}

template <typename T>
const VectorT<T>& VectorT<T>::divide(const VectorT<T>& v)
{
  if (v.size() != size())
    my_throw("VectorT<T>::divide: Wrong size");
  for (size_type i = 0, n = size(); i < n; i++)
  {
    if (abs(v[i]) < EPSILON10)
      my_throw("VectorT<T>::divide: division by 0");
    set(i, at(i)/v[i]);
  }
  return *this;
}

template <typename T>
const VectorT<T>& VectorT<T>::add(const T& v)
{
  std::for_each(begin(), end(), [v](T& d) { d += v;});
  return *this;
}

template <typename T>
const VectorT<T>& VectorT<T>::subtract(const T& v)
{
  std::for_each(begin(), end(), [v](T& d) { d -= v;});
  return *this;
}

template <typename T>
const VectorT<T>& VectorT<T>::multiply(const T& v)
{
  std::for_each(begin(), end(), [v](T& d) { d *= v;});
  return *this;
}

template <typename T>
const VectorT<T>& VectorT<T>::divide(const T& v)
{
  if (abs(v) < EPSILON10)
    my_throw("VectorT<T>::divide: division by 0");
  std::for_each(begin(), end(), [v](T& d) { d /= v;});
  return *this;
}

template <typename T>
std::string VectorT<T>::toString() const
{
  std::stringstream sstr;
  sstr << "[";
  for (size_type i = 0, n = size(); i < n; i++)
  {
    sstr << _v->at(i);
    if (i < n-1) sstr << ", ";
  }
  sstr << "]";
  return sstr.str();
}

template <typename T>
void VectorT<T>::_detach()
{
  if (_v.unique())
    return;
  _v = std::make_shared<Vector>(*_v);
}


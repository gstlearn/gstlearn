#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"
#include "Basic/AException.hpp"

#include <vector>
#include <sstream>
#include <memory>
#include <limits>
#include <algorithm>
#include <cmath>

/***************************************************************************
 **
 ** Vector of T values (any type).
 ** T type must define copy constructor and assignment operator
 **
 ***************************************************************************/
template <typename T>
class GSTLEARN_EXPORT VectorT
{
public:
  typedef std::vector<T> Vector;
  typedef typename Vector::value_type               value_type;
  typedef typename Vector::size_type                size_type;
  typedef typename Vector::iterator                 iterator;
  typedef typename Vector::const_iterator           const_iterator;
  typedef typename Vector::reverse_iterator         reverse_iterator;
  typedef typename Vector::const_reverse_iterator   const_reverse_iterator;

public:
  inline VectorT()                                                    : _v(std::make_shared<Vector>()) { }
  inline VectorT(const Vector& vec)                                   : _v(std::make_shared<Vector>(vec)) { }
  inline VectorT(size_type count, const T& value = T())               : _v(std::make_shared<Vector>(count, value)) { }
  inline VectorT(const T* first, const T* last)                       : _v(std::make_shared<Vector>()) { _v->assign(first, last); }
  inline VectorT(const VectorT& other)                                : _v(other._v) { }
#ifndef SWIG
  inline VectorT(std::initializer_list<T> init)                       : _v(std::make_shared<Vector>(init)) { }
  inline VectorT(VectorT&& other)                                     { _v.swap(other._v); }
#endif
  inline ~VectorT()                                                   { }

#ifndef SWIG
  inline operator const Vector&() const                               { return *_v; }
#endif

  inline Vector& getVector() const                                    { return *_v; }
  inline Vector* getVectorPtr() const                                 { return _v.get(); }

#ifndef SWIG
  inline VectorT& operator=(const Vector& vec)                        { _detach(); *_v = vec; return (*this); }
  inline VectorT& operator=(const VectorT& other)                     { _detach(); _v = other._v; return (*this); }
  inline VectorT& operator=(VectorT&& other)                          { _v.swap(other._v); return (*this); }
  inline VectorT& operator=(std::initializer_list<T> init)            { _detach(); (*_v) = init; return (*this); }
#endif

  inline bool operator==(const VectorT& other) const                  { return *_v == *other._v; }
  inline bool operator!=(const VectorT& other) const                  { return *_v != *other._v; }
  inline bool operator <(const VectorT& other) const                  { return *_v  < *other._v; }
  inline bool operator<=(const VectorT& other) const                  { return *_v <= *other._v; }
  inline bool operator >(const VectorT& other) const                  { return *_v  > *other._v; }
  inline bool operator>=(const VectorT& other) const                  { return *_v >= *other._v; }

  inline const T& at(size_type pos) const                             { if (pos >= _v->size()) my_throw("VectorT<T>::at: index out of range");                    return _v->operator[](pos); }
#ifndef SWIG
  inline T& operator[](size_type pos)                                 { if (pos >= _v->size()) my_throw("VectorT<T>::operator[]: index out of range"); _detach(); return _v->operator[](pos); }
  inline const T& operator[](size_type pos) const                     { if (pos >= _v->size()) my_throw("VectorT<T>::operator[]: index out of range");            return _v->operator[](pos); }
#endif
  inline T& front()                                                   { _detach(); return _v->front(); }
  inline const T& front() const                                       { return _v->front(); }
  inline T& back()                                                    { _detach(); return _v->back(); }
  inline const T& back() const                                        { return _v->back(); }
  inline T* data()                                                    { _detach(); return _v->data(); }
  inline const T* data() const                                        { return _v->data(); }
  inline const T* constData() const                                   { return _v->data(); }

  inline bool empty() const                                           { return _v->empty(); }
  inline size_type size() const                                       { return _v->size(); }
  inline void reserve(size_type new_cap)                              { _v->reserve(new_cap); }
  inline size_type capacity() const                                   { return _v->capacity(); }
  inline void clear()                                                 { _detach(); _v->clear(); }

  inline void insert(size_type i, const T& value)                     { _detach(); _v->insert(begin() + i, value); }
  inline void insert(size_type i, size_type count, const T& value)    { _detach(); _v->insert(begin() + i, count, value); }
  inline iterator insert(const_iterator pos,
                         const_iterator first,
                         const_iterator last )                        { _detach(); return _v->insert(pos, first, last); }
  inline void remove(size_type i)                                     { _detach(); _v->erase(begin() + i); }
  inline void remove(size_type i, size_type count)                    { _detach(); _v->erase(begin() + i, begin() + i + count); }
  inline iterator erase( const_iterator pos )                         { _detach(); return _v->erase(pos); }
  inline iterator erase( const_iterator first, const_iterator last)   { _detach(); return _v->erase(first, last); }

  inline void push_back(const T& value)                               { _detach(); _v->push_back(value); }
  inline void push_front(const T& value)                              { _detach(); _v->insert(begin(), value); }
#ifndef SWIG
  inline void push_back(const T&& value)                              { _detach(); _v->push_back(value); }
  inline void push_front(const T&& value)                             { _detach(); _v->insert(begin(), value); }
#endif
  inline void resize(size_type count)                                 { if (count == size()) return; _detach(); _v->resize(count); }
  inline void resize(size_type count, const T& value)                 { if (count == size()) return; _detach(); _v->resize(count, value); }

  inline iterator begin()                                             { _detach(); return _v->begin(); }
  inline const_iterator begin() const                                 { return _v->begin(); }
  inline const_iterator cbegin() const                                { return _v->cbegin(); }
  inline iterator end()                                               { _detach(); return _v->end(); }
  inline const_iterator end() const                                   { return _v->end(); }
  inline const_iterator cend() const                                  { return _v->cend(); }

  inline reverse_iterator rbegin()                                    { _detach(); return _v->rbegin(); }
  inline const_reverse_iterator crbegin() const                       { return _v->crbegin(); }
  inline reverse_iterator rend()                                      { _detach(); return _v->rend(); }
  inline const_reverse_iterator crend() const                         { return _v->crend(); }

#ifndef SWIG
  inline VectorT& operator<<(const T& value)                          { _detach(); _v->push_back(value); return (*this); }
  inline VectorT& operator<<(const VectorT<T>& v)                     { _detach(); reserve(size() + v.size()); std::for_each(v.cbegin(), v.cend(), [&](const T& value) { _v->push_back(value); }); return (*this); }
#endif

  inline void swap(VectorT& other);
  inline bool contains(const T& value) const;
  inline void fill(const T& value, size_type size = -1);
  inline void assign(const T* first, const T* last);
  inline void set(size_type i, const T& value);

  inline String toString() const;

protected:
  std::shared_ptr<Vector> _v;

private:
  inline void _detach();
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
void VectorT<T>::assign(const T* first, const T* last)
{
  _detach();
  _v->assign(first, last);
}

template <typename T>
void VectorT<T>::set(size_type i, const T& value)
{
  _detach();
  operator[](i) = value;
}

template <typename T>
String VectorT<T>::toString() const
{
  std::stringstream sstr;
  sstr << "[";
  for (size_type i = 0, n = size(); i < n; i++)
  {
    sstr << at(i);
    if (i != n-1)
      sstr << " ";
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


// Force instantiation for VectorT (for windows export)
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
  // Do not export VectorXXX to SWIG (no more instantiation needed)
  #ifndef SWIG
    GSTLEARN_TEMPLATE_EXPORT template class VectorT<int>;
    GSTLEARN_TEMPLATE_EXPORT template class VectorT<double>;
    GSTLEARN_TEMPLATE_EXPORT template class VectorT<float>;
    GSTLEARN_TEMPLATE_EXPORT template class VectorT<unsigned char>;
//    GSTLEARN_TEMPLATE_EXPORT template class VectorT<bool>;
    GSTLEARN_TEMPLATE_EXPORT template class VectorT<String>;
  #endif
#endif

//typedef VectorT<bool> VectorBool;
// https://stackoverflow.com/a/61158013/3952924
typedef VectorT<unsigned char> VectorBool;
typedef VectorT<String> VectorString;

template <typename T>
std::ostream& operator<<(std::ostream& os, const VectorT<T>& vec);


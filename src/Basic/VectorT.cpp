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
#include "Basic/VectorT.hpp"

template <>
std::ostream& operator<<(std::ostream& os, const VectorT<int>& vec)
{
  os << vec.toString();
  return os;
}
template <>
std::ostream& operator<<(std::ostream& os, const VectorT<double>& vec)
{
  os << vec.toString();
  return os;
}
template <>
std::ostream& operator<<(std::ostream& os, const VectorT<String>& vec)
{
  os << vec.toString();
  return os;
}
template <>
std::ostream& operator<<(std::ostream& os, const VectorT<float>& vec)
{
  os << vec.toString();
  return os;
}
template <>
std::ostream& operator<<(std::ostream& os, const VectorT<UChar>& vec)
{
  os << vec.toString();
  return os;
}

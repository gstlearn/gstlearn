/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Basic/VectorNumT.hpp"

template <>
std::ostream& operator<<(std::ostream& os, const VectorT<VectorNumT<int> >& vec)
{
  os << "[";
  for (int i = 0, n = (int)vec.size(); i < n; i++)
  {
    os << vec.at(i).toString();
    if (i != n-1)
      os << " ";
  }
  os << "]";
  return os;
}

template <>
std::ostream& operator<<(std::ostream& os, const VectorT<VectorNumT<double> >& vec)
{
  os << "[";
  for (int i = 0, n = (int)vec.size(); i < n; i++)
  {
    os << vec.at(i).toString();
    if (i != n-1)
      os << " ";
  }
  os << "]";
  return os;
}

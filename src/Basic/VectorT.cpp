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

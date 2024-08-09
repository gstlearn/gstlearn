/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2024) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Matrix/VectorEigen.hpp"
#include "Basic/AException.hpp"
#include <Eigen/src/Core/Matrix.h>

VectorEigen::VectorEigen(int size) : _eigenVector(size) {}

VectorEigen::VectorEigen(const VectorEigen &v) : _eigenVector(v._eigenVector) {}

VectorEigen::VectorEigen(const VectorDouble &v)
  : _eigenVector(Eigen::VectorXd::Map(v.data(), v.size()))
{
}

#ifndef SWIG
VectorEigen::VectorEigen(const Eigen::VectorXd& v)
  : _eigenVector(v)
{
}
#endif

VectorEigen& VectorEigen::operator= (const VectorEigen &r)
{
  if (this != &r)
  {
    _eigenVector = r._eigenVector;
  }
  return *this;
}

VectorEigen::~VectorEigen() {}

/**
 * @brief Set the value at a given position in the vector 
 * 
 * @param i index position
 * @param value new value
 * @param flagCheck true to check index position consistency
 */
void VectorEigen::setValue(int i, double value, bool flagCheck)
{
  if (flagCheck && (i < 0 || i >= _eigenVector.size()))
    my_throw("Wrong vector index");
  _eigenVector[i] = value;
}

/**
 * @brief Get the value at a given position
 *
 * @param i index position
 * @param flagCheck true to check index position consistency
 * @return the value
 */
double VectorEigen::getValue(int i, bool flagCheck) const
{
  if (flagCheck && (i < 0 || i >= _eigenVector.size()))
    my_throw("Wrong vector index");
  return _eigenVector[i];
}

/**
 * @brief Get all values in a VectorDouble
 *
 * @return VectorDouble
 */
VectorDouble VectorEigen::getValues() const
{
  VectorDouble vec(_eigenVector.size());
  Eigen::Map<Eigen::VectorXd>(vec.data(), vec.size()) = _eigenVector;
  return vec;
}

/**
 * @brief Set all the values of the Vector at once
 *
 * @param value value to be filled
 */
void VectorEigen::fill(double value)
{
  _eigenVector.fill(value);
}

std::ostream& operator<<(std::ostream& os, const VectorEigen& vec)
{
  os << vec.getValues().toString();
  return os;
}
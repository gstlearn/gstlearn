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
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"

#include <iostream>

#ifndef SWIG
#include <Eigen/Core>
#endif

/**
 * Eigen vector wrapper class
 */

class GSTLEARN_EXPORT VectorEigen {

public:
  VectorEigen(int size = 0);
  VectorEigen(const VectorEigen &v);
  VectorEigen(const VectorDouble& v);
#ifndef SWIG
  VectorEigen(const Eigen::VectorXd& v);
#endif
  VectorEigen& operator= (const VectorEigen &r);
	virtual ~VectorEigen();

  /*! Set the value at a given position in the vector */
  void setValue(int i, double value, bool flagCheck = false);
  /*! Get the value at a given position */
  double getValue(int i, bool flagCheck = false) const;
  /*! Get all values in a VectorDouble */
  VectorDouble getValues() const;

  /*! Set all the values of the Vector at once */
  void fill(double value);

#ifndef SWIG
  /*! Get underlying Eigen vector */
  const Eigen::VectorXd& getVector() const { return _eigenVector; }
  /*! Get underlying Eigen vector */
  Eigen::VectorXd& getVector() { return _eigenVector; }

private:
  Eigen::VectorXd _eigenVector; /// Eigen storage for vector in Eigen Library
#endif
  
};

GSTLEARN_EXPORT std::ostream& operator<<(std::ostream& os,
                                         const VectorEigen& vec);
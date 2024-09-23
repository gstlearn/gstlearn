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

#ifndef SWIG
  #include <Eigen/Core>
  #include <Eigen/Dense>
  #include <Eigen/src/Core/Matrix.h>
#endif

#include <iostream>



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
  int size() const { return _eigenVector.size();}
  void resize(int n ) { _eigenVector.resize(n);}
  bool empty() const { return _eigenVector.size()==0; }
#ifndef SWIG
  static double maximum(const std::vector<Eigen::VectorXd>& vect);
  static void simulateGaussianInPlace(Eigen::VectorXd& vect, double mean = 0, double var = 1);
  static void fill(Eigen::VectorXd &vect, double val = 0.);
  static void addMultiplyConstantInPlace(double val ,const Eigen::VectorXd& in, Eigen::VectorXd& res, int iad);
  static void addMultiplyVectVectInPlace(const Eigen::VectorXd&in1 ,const Eigen::VectorXd& in2, Eigen::VectorXd& res, int iad);
  static void addInPlace(const Eigen::VectorXd& in, Eigen::VectorXd& res, int iad);
  static VectorDouble copyIntoVD(const Eigen::VectorXd& in);
  static double innerProduct(const std::vector<Eigen::VectorXd> &in1,const std::vector<Eigen::VectorXd> &in2 );
  static void fill(std::vector<Eigen::VectorXd> &vect, double val = 0.);
  static void copy(const Eigen::VectorXd& in, Eigen::VectorXd& dest);
  static void copy(const Eigen::VectorXd& in, Eigen::Map<Eigen::VectorXd>& dest); //TODO this function shouldn't be used
                                                                                  //use template to avoid it
  static void copy(const Eigen::VectorXd&in, VectorDouble& dest);
  static void  linearCombinationVVDInPlace(double coeff1, const std::vector<Eigen::VectorXd> &in1,
                                           double coeff2, const std::vector<Eigen::VectorXd> &in2,
                                           std::vector<Eigen::VectorXd> &res);
  static void  substractInPlace(const std::vector<Eigen::VectorXd> &in1,
                                const std::vector<Eigen::VectorXd> &in2,
                                std::vector<Eigen::VectorXd> &res);

  static void copy(const std::vector<Eigen::VectorXd> &in,std::vector<Eigen::VectorXd>& dest);

  static void addInPlace(const Eigen::VectorXd& in, Eigen::VectorXd& out);
  static void divideInPlace(const Eigen::VectorXd& in, Eigen::VectorXd& out);

  /*! Get underlying Eigen vector */
  static void addInPlace(const Eigen::VectorXd& t1, const Eigen::VectorXd& t2,Eigen::VectorXd& res);
  static void addInPlace(const std::vector<Eigen::VectorXd>& t1, 
                         const std::vector<Eigen::VectorXd>& t2,
                         std::vector<Eigen::VectorXd>& res);
  static void multiplyConstant(Eigen::VectorXd vect, double val = 1.);
  static Eigen::VectorXd flatten(const std::vector<Eigen::VectorXd>& vvd);
  static void flattenInPlace(const std::vector<Eigen::VectorXd>& vvd, Eigen::VectorXd& vd);
  static std::vector<Eigen::VectorXd> unflatten(const Eigen::VectorXd& vd, const VectorInt& sizes);
  static void unflattenInPlace(const Eigen::VectorXd& vd, std::vector<Eigen::VectorXd>& vvd);
 
  const Eigen::VectorXd& getVector() const { return _eigenVector; }
  /*! Get underlying Eigen vector */
  Eigen::VectorXd& getVector() { return _eigenVector; }
  /*! Get map to underlying Eigen vector */
  Eigen::Map<const Eigen::VectorXd> getMap() const
  {
    return {_eigenVector.data(), _eigenVector.size()};
  }
  /*! Get map to underlying Eigen vector */
  Eigen::Map<Eigen::VectorXd> getMap()
  {
    return {_eigenVector.data(), _eigenVector.size()};
  }

private:
  Eigen::VectorXd _eigenVector; /// Eigen storage for vector in Eigen Library
#endif
  
};

#ifndef SWIG
GSTLEARN_EXPORT std::ostream& operator<<(std::ostream& os,
                                         const VectorEigen& vec);
#endif
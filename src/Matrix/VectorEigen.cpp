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
#include "Basic/Law.hpp"
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

void VectorEigen::fill(Eigen::VectorXd &vect, double val)
  { 
    for (int i = 0; i < (int)vect.size(); i++)
    {
      vect[i] = val;
    }
  }

void VectorEigen::addMultiplyConstantInPlace(double val1,
                                              const Eigen::VectorXd &in,
                                              Eigen::VectorXd &out,
                                              int iad)
{
    double * outp = out.data() + iad;
    const double* inp = in.data();
    for (int i = 0; i < (int)in.size();i++)
    {
      *(outp++) += val1 * *(inp++);
    }
}

double VectorEigen::maximum(const std::vector<Eigen::VectorXd>& vect)
{
  double max = vect[0].maxCoeff();
  for (int i = 1; i < (int)vect.size(); i++)
  {
    double v = vect[i].maxCoeff();
    max = (v>max)?v:max;
  }
  return max;
}

void VectorEigen::addInPlace(const Eigen::VectorXd& in, Eigen::VectorXd& out)
{
  out += in;
}

void VectorEigen::addInPlace(const std::vector<Eigen::VectorXd>& t1,
                             const std::vector<Eigen::VectorXd>& t2,
                             std::vector<Eigen::VectorXd>& res)
{
   for (int i = 0; i < (int)t1.size(); i++)
   {
    res[i] = t1[i] + t2[i];
   }
}

void VectorEigen::divideInPlace(const Eigen::VectorXd& in, Eigen::VectorXd& out)
{
  out.array() /= in.array();
}

double VectorEigen::innerProduct(const std::vector<Eigen::VectorXd> &in1,const std::vector<Eigen::VectorXd> &in2)
{
  double res = 0.;
  for (int i = 0; i < (int)in1.size(); i++)
  {
    res+= in1[i].adjoint()*in2[i];
  }
  return res;
}

void VectorEigen::fill(std::vector<Eigen::VectorXd> &vect, double val)
{
  for(auto &e: vect)
  {
    e.fill(val);
  }
}
void VectorEigen::linearCombinationVVDInPlace(double coeff1, const std::vector<Eigen::VectorXd> &in1,
                                              double coeff2, const std::vector<Eigen::VectorXd> &in2,
                                              std::vector<Eigen::VectorXd> &res)
{
  for(int i = 0; i < (int)in1.size(); i++)
  {
    res[i] = coeff1 * in1[i] + coeff2 * in2[i];
  }
}

void VectorEigen::substractInPlace(const std::vector<Eigen::VectorXd> &in1,
                                   const std::vector<Eigen::VectorXd> &in2,
                                   std::vector<Eigen::VectorXd> &res)
{
  for(int i = 0; i < (int)in1.size(); i++)
  {
    res[i] = in2[i] - in2[i];
  }
}

void VectorEigen::copy(const Eigen::VectorXd& in, Eigen::VectorXd& dest)
{
  for (int i = 0; i < (int)in.size(); i++)
  {
    dest[i] = in[i];
  }
}

void VectorEigen::copy(const std::vector<Eigen::VectorXd> &in,std::vector<Eigen::VectorXd>& dest)
{
  for(int i = 0; i < (int)in.size(); i++)
  {
    VectorEigen::copy(in[i],dest[i]);
  }
}  
void VectorEigen::simulateGaussianInPlace(Eigen::VectorXd& vect, double mean, double var)
{
  for(int i = 0; i < vect.size(); i++)
    vect[i] = law_gaussian(mean,var);
}

void VectorEigen::addInPlace(const Eigen::VectorXd& t1, const Eigen::VectorXd& t2,Eigen::VectorXd& res)
{
  res = t1 + t2;
}

/**
 * Method which flattens a std::vector<Eigen::VectorXd> into an Eigen::VectorXd
 * @param vvd Input std::vector<Eigen::VectorXd>
 * @return Returned Eigen::VectorXd
 */
Eigen::VectorXd VectorEigen::flatten(const std::vector<Eigen::VectorXd>& vvd)
{
  int sizetot = 0;
  for(auto &e : vvd)
  {
    sizetot += e.size(); 
  }

  Eigen::VectorXd vd(sizetot);
  
  flattenInPlace(vvd,vd);

  return vd;
}

void VectorEigen::multiplyConstant(Eigen::VectorXd vect, double val)
{
  vect*= val;
}


void VectorEigen::flattenInPlace(const std::vector<Eigen::VectorXd>& vvd, Eigen::VectorXd& vd)
{
  int s = 0;
  for (int i = 0; i < (int) vvd.size(); i++)
    for (int j = 0; j < (int) vvd[i].size(); j++)
      vd[s++] = vvd[i][j];
}

std::vector<Eigen::VectorXd> VectorEigen::unflatten(const Eigen::VectorXd& vd, const VectorInt& sizes)
{
  std::vector<Eigen::VectorXd> vvd;

  int lec = 0;
  for (int i = 0, n = (int) sizes.size(); i < n; i++)
  {
    int lng = sizes[i];
    Eigen::VectorXd local(lng);
    vvd.push_back(local);
  }
  unflattenInPlace(vd,vvd);
  return vvd;
}

void VectorEigen::unflattenInPlace(const Eigen::VectorXd& vd, std::vector<Eigen::VectorXd>& vvd)
{
  int lec = 0;
  for (int i = 0, n = (int) vvd.size(); i < n; i++)
    for (int j = 0; j < (int) vvd[i].size(); j++)
      vvd[i][j] = vd[lec++];
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
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
#include "LinearOp/ASimulable.hpp"
#include "Basic/AStringable.hpp"

VectorDouble ASimulable::evalSimulate(const VectorDouble& whitenoise) const
{
  VectorDouble res(whitenoise.size());
  evalSimulate(whitenoise, res);
  return res;
}

int ASimulable::addSimulateToDest(const VectorDouble& whitenoise,
                                  VectorDouble& outv) const
{
  try
  {
    Eigen::Map<const Eigen::VectorXd> myInv(whitenoise.data(),
                                            whitenoise.size());
    Eigen::VectorXd myOut(outv.size());
    VectorEigen::fill(myOut, 0.);
    // Assume outv has the good size
    if (_addSimulateToDest(myInv, myOut)) return 1;

    VectorEigen::copy(myOut, outv);
  }
  catch (const std::string& str)
  {
    // TODO : Check if std::exception can be used
    messerr("%s", str.c_str());
  }
  return 0;
}

int ASimulable::addSimulateToDest(const VectorEigen& whitenoise,
                                  VectorEigen& outv) const
{
  return _addSimulateToDest(whitenoise.getVector(), outv.getVector());
}

int ASimulable::evalSimulate(const Eigen::VectorXd& whitenoise,
                             Eigen::VectorXd& outv) const
{
  int n = (int)outv.size();
  for (int i = 0; i < n; i++)
  {
    outv[i] = 0.;
  }
  return _addSimulateToDest(whitenoise, outv);
}

int ASimulable::addSimulateToDest(const Eigen::VectorXd& whitenoise,
                                  Eigen::VectorXd& outv) const
{
  return _addSimulateToDest(whitenoise, outv);
}

/*****************************************************************************/
/*!
**  Evaluate the product: 'outv' = Q * 'inv'
**
** \param[in]  whitenoise     Array of input values
**
** \param[out] outv           Array of output values
**
*****************************************************************************/
int ASimulable::evalSimulate(const VectorDouble& whitenoise,
                             VectorDouble& outv) const
{
  return addSimulateToDest(whitenoise, outv);
}

/*****************************************************************************/
/*!
**  Evaluate the product: 'outv' = Q * 'inv'
**
** \param[in]  whitenoise     Array of input values
**
** \param[out] outv    Array of output values
**
*****************************************************************************/
int ASimulable::evalSimulate(const VectorEigen& whitenoise,
                             VectorEigen& outv) const
{
  return evalSimulate(whitenoise.getVector(), outv.getVector());
}

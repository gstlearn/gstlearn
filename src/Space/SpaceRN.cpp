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
#include "Space/SpaceRN.hpp"
#include "Space/SpacePoint.hpp"
#include "Basic/Tensor.hpp"
#include "Basic/VectorHelper.hpp"

#include <math.h>

SpaceRN::SpaceRN(unsigned int ndim, bool addtime)
    : ASpace(ndim, addtime)
{
}

SpaceRN::SpaceRN(const SpaceRN &r)
    : ASpace(r)
{
}

SpaceRN& SpaceRN::operator=(const SpaceRN &r)
{
  if (this != &r)
  {
    ASpace::operator=(r);
  }
  return *this;
}

SpaceRN::~SpaceRN()
{
}

SpaceRN* SpaceRN::create(unsigned int ndim, bool addtime)
{
  return new SpaceRN(ndim, addtime);
}

void SpaceRN::_move(SpacePoint &p1, const VectorDouble &vec) const
{
  for (unsigned int i = _iDimOffset; i < _nDim + _iDimOffset; i++)
  {
    p1.setCoord(i, p1.getCoord(i) + vec[i]);
  }
}

/**
 * Return the distance between two space points in RN Space
 The distance between \f$p1=(x_1,y_1)\f$ and \f$p2=(x_2,y_2)\f$ is \f$\sqrt{(x_2-x_1)^2+(y_2-y_1)^2}\f$.
 \param[in] p1 First point
 \param[in] p2 Second point
 \return The distance between p1 and p2

 \note The code has been optimized in order to avoid using '_work1' for storing
 \note temporary results
 */
double SpaceRN::_getDistance(const SpacePoint &p1, const SpacePoint &p2) const
{
  double dist = 0.;
  double delta = 0.;
  for (unsigned int i = _iDimOffset; i < _nDim + _iDimOffset; i++)
  {
    delta = p2.getCoord(i) - p1.getCoord(i);
    dist += delta * delta;
  }
  return sqrt(dist);
}

double SpaceRN::_getDistance(const SpacePoint &p1,
                             const SpacePoint &p2,
                             const Tensor &tensor) const
{
  _getIncrementInPlace(p1, p2, _work1);

  if (!tensor.isFlagDefinedByInverse2())
  {
    tensor.applyInverseInPlace(_work1, _work2);
    return VH::norm(_work2);
  }
  tensor.applyInverse2InPlace(_work1, _work2);
  return sqrt(VH::innerProduct(_work1, _work2));
}

double SpaceRN::_getDistance1D(double c1, double c2) const
{
  return c2 - c1;
}

double SpaceRN::_getFrequentialDistance(const SpacePoint &p1,
                                        const SpacePoint &p2,
                                        const Tensor &tensor) const
{
  _getIncrementInPlace(p1, p2, _work1);
  tensor.applyDirectSwapInPlace(_work1, _work2);
  return VH::norm(_work2);
}

VectorDouble SpaceRN::_getIncrement(const SpacePoint &p1,
                                    const SpacePoint &p2) const
{
  _getIncrementInPlace(p1, p2, _work1);
  return _work1;
}

void SpaceRN::_getIncrementInPlace(const SpacePoint &p1,
                                   const SpacePoint &p2,
                                   VectorDouble &ptemp) const
{
  int j = 0;
  for (unsigned int i = _iDimOffset; i < _nDim + _iDimOffset; i++)
    ptemp[j++] = p2.getCoord(i) - p1.getCoord(i);
}

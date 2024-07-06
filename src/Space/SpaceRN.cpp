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

SpaceRN::SpaceRN(unsigned int ndim)
    : ASpace(ndim)
{
  if (ndim == 0)
  {
    _nDim = 2;
  }
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

SpaceRN* SpaceRN::create(unsigned int ndim)
{
  return new SpaceRN(ndim);
}

void SpaceRN::move(SpacePoint &p1, const VectorDouble &vec) const
{
  p1.setCoord(VH::add(p1.getCoord(), vec));
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
double SpaceRN::getDistance(const SpacePoint &p1, const SpacePoint &p2) const
{
//  _getIncrementInPlace(p1, p2, _work1);
//  return VH::norm(_work1);
  double dist = 0.;
  double delta = 0.;
  for (unsigned int i = 0; i < _nDim; i++)
  {
    delta = p2.getCoord(i) - p1.getCoord(i);
    dist += delta * delta;
  }
  return sqrt(dist);
}

double SpaceRN::getDistance(const SpacePoint &p1,
                            const SpacePoint &p2,
                            const Tensor &tensor) const
{
  _getIncrementInPlace(p1, p2, _work1);

  if (! tensor.isFlagDefinedByInverse2())
  {
    tensor.applyInverseInPlace(_work1, _work2);
    return VH::norm(_work2);
  }
  else
  {
    tensor.applyInverse2InPlace(_work1, _work2);
    return sqrt(VH::innerProduct(_work1, _work2));
  }
}

double SpaceRN::getDistance1D(const SpacePoint &p1,
                              const SpacePoint &p2,
                              int idim) const
{
  _getIncrementInPlace(p1, p2, _work1);
  if (idim > (int) getNDim())
    return TEST;
  else
    return _work1[idim];
}

void SpaceRN::_getIncrementInPlace(const SpacePoint &p1,
                                   const SpacePoint &p2,
                                   VectorDouble &ptemp) const
{
  for (unsigned int i = 0; i < _nDim; i++)
    ptemp[i] = p2.getCoord(i) - p1.getCoord(i);
}

double SpaceRN::getFrequentialDistance(const SpacePoint &p1,
                                       const SpacePoint &p2,
                                       const Tensor &tensor) const
{
  _getIncrementInPlace(p1, p2, _work1);
  tensor.applyDirectSwapInPlace(_work1, _work2);
  return VH::norm(_work2);
}

VectorDouble SpaceRN::getIncrement(const SpacePoint &p1,
                                   const SpacePoint &p2) const
{
  return VH::subtract(p1.getCoord(), p2.getCoord());
}

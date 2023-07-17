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
    messerr("Wrong dimension = %d when creating SpaceRN (ndim set to 2)", ndim);
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
 */
double SpaceRN::getDistance(const SpacePoint &p1, const SpacePoint &p2) const
{
  _getIncrementInPlace(p1, p2, _work1);
  return VH::norm(_work1);
}

double SpaceRN::getDistance(const SpacePoint &p1,
                            const SpacePoint &p2,
                            const Tensor &tensor) const
{
  _getIncrementInPlace(p1, p2, _work1);
  tensor.applyInverseInPlace(_work1, _work2);
  return VH::norm(_work2);
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
  tensor.applyDirectInPlace(_work1, _work2);
  return VH::norm(_work2);
}

VectorDouble SpaceRN::getIncrement(const SpacePoint &p1,
                                   const SpacePoint &p2) const
{
  return VH::subtract(p1.getCoord(), p2.getCoord());
}

String SpaceRN::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  sstr << "Space Type      = " << getType().getKey() << std::endl;
  sstr << "Space Dimension = " << getNDim() << std::endl;
  return sstr.str();
}

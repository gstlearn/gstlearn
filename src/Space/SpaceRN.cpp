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
#include "Space/ASpace.hpp"
#include "Space/SpacePoint.hpp"
#include "Basic/Tensor.hpp"
#include "Basic/VectorHelper.hpp"

#include <math.h>
#include <memory>

SpaceRN::SpaceRN(unsigned int ndim)
  : ASpace(ndim)
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
ASpaceSharedPtr SpaceRN::create(int ndim)
{
  return std::shared_ptr<const SpaceRN>(new SpaceRN(ndim));
}

void SpaceRN::_move(SpacePoint &p1, const VectorDouble &vec) const
{
  unsigned int offset = getOffset();
  unsigned int ndim = getNDim();
  for (unsigned int i = offset; i < ndim + offset; i++)
  {
    p1.setCoord(i, p1.getCoord(i) + vec[i]);
  }
}

/**
 * Return the distance between two space points in RN Space
 The distance between \f$p1=(x_1,y_1)\f$ and \f$p2=(x_2,y_2)\f$ is \f$\sqrt{(x_2-x_1)^2+(y_2-y_1)^2}\f$.
 \param[in] p1 First point
 \param[in] p2 Second point
 \param[in] ispace Rank of the sub-space (default: -1)
 \return The distance between p1 and p2

 \note The code has been optimized in order to avoid using '_work1' for storing
 \note temporary results
 */
double SpaceRN::_getDistance(const SpacePoint& p1,
                             const SpacePoint& p2,
                             int ispace) const
{
  DECLARE_UNUSED(ispace);
  double dist = 0.;
  double delta = 0.;
  unsigned int offset = getOffset();
  unsigned int ndim   = getNDim();
  for (unsigned int i = offset; i < ndim + offset; i++)
  {
    delta = p2.getCoord(i) - p1.getCoord(i);
    dist += delta * delta;
  }
  return sqrt(dist);
}

double SpaceRN::_getDistance(const SpacePoint& p1,
                             const SpacePoint& p2,
                             const Tensor& tensor,
                             int ispace) const
{
  DECLARE_UNUSED(ispace);
  _getIncrementInPlace(p1, p2, _work1);

  if (!tensor.isFlagDefinedByInverse2())
  {
    tensor.applyInverseInPlace(_work1, _work2);
    return VH::norm(_work2);
  }
  tensor.applyInverse2InPlace(_work1, _work2);
  return sqrt(VH::innerProduct(_work1, _work2));
}

double SpaceRN::_getFrequentialDistance(const SpacePoint& p1,
                                        const SpacePoint& p2,
                                        const Tensor& tensor,
                                        int ispace) const
{
  DECLARE_UNUSED(ispace);
  _getIncrementInPlace(p1, p2, _work1);
  tensor.applyDirectSwapInPlace(_work1, _work2);
  return VH::norm(_work2);
}

VectorDouble SpaceRN::_getIncrement(const SpacePoint& p1,
                                    const SpacePoint& p2,
                                    int ispace) const
{
  DECLARE_UNUSED(ispace);
  _getIncrementInPlace(p1, p2, _work1);
  return _work1;
}

void SpaceRN::_getIncrementInPlace(const SpacePoint& p1,
                                   const SpacePoint& p2,
                                   VectorDouble& ptemp,
                                   int ispace) const
{
  DECLARE_UNUSED(ispace);
  int j = 0;
  unsigned int offset = getOffset();
  unsigned int ndim   = getNDim();
  unsigned int maxlength = ndim + offset;
  for (unsigned int i = offset; i < maxlength; i++)
    ptemp[j++] = p2.getCoord(i) - p1.getCoord(i);
}


void SpaceRN::getDistancePointVectInPlace(const SpacePoint& p1,
										  const std::vector<SpacePoint>& p2,
	                                      VectorDouble& res) const
{
	double ti;
	double s;
	int nbp = res.size();
	for(int i = 0; i<nbp;i++)
	{
		s = 0.;

		for(unsigned int idim = 0;idim<_nDim;idim++)
		{
			ti = p1.getCoord(idim) - p2[i].getCoord(idim);
			s+= ti * ti;
		}

		res[i] = sqrt(s);
	}
}
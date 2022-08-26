#include "Space/SpaceRN.hpp"
#include "Space/SpacePoint.hpp"
#include "Basic/Tensor.hpp"
#include "Basic/Vector.hpp"

SpaceRN::SpaceRN(unsigned int ndim)
 : ASpace(ndim)
{
  if (ndim == 0)
  {
    messerr("Wrong dimension = %d when creating SpaceRN (ndim set to 2)", ndim);
    _nDim = 2;
  }
}

SpaceRN::SpaceRN(const SpaceRN& r)
: ASpace(r)
{
}

SpaceRN& SpaceRN::operator=(const SpaceRN& r)
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

void SpaceRN::move(SpacePoint& p1,
                   const VectorDouble& vec) const
{
  p1.setCoord(ut_vector_add(p1.getCoord(), vec));
}

double SpaceRN::getDistance(const SpacePoint& p1,
                            const SpacePoint& p2) const
{
  return ut_vector_norm(getIncrement(p1, p2));
}

double SpaceRN::getDistance(const SpacePoint& p1,
                            const SpacePoint& p2,
                            const Tensor& tensor) const
{
  return ut_vector_norm(tensor.applyInverse(getIncrement(p1, p2)));
}

double SpaceRN::getFrequentialDistance(const SpacePoint& p1,
                                       const SpacePoint& p2,
                                       const Tensor& tensor) const
{
  return ut_vector_norm(tensor.applyDirect(getIncrement(p1, p2),0));
}

VectorDouble SpaceRN::getIncrement(const SpacePoint& p1,
                                   const SpacePoint& p2) const
{
  return ut_vector_subtract(p1.getCoord(), p2.getCoord());
}

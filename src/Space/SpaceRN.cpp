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

SpaceRN* SpaceRN::create(unsigned int ndim)
{
  return new SpaceRN(ndim);
}

void SpaceRN::move(SpacePoint& p1,
                   const VectorDouble& vec) const
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
double SpaceRN::getDistance(const SpacePoint& p1,
                            const SpacePoint& p2) const
{
  return VH::norm(getIncrement(p1, p2));
}

double SpaceRN::getDistance(const SpacePoint& p1,
                            const SpacePoint& p2,
                            const Tensor& tensor) const
{
  return VH::norm(tensor.applyInverse(getIncrement(p1, p2)));
}


double SpaceRN::_getDistance(const SpacePoint& p1,
							 const SpacePoint & p2,
							 VectorDouble& ptemp,
       const Tensor& tensor,
	   VectorDouble& temp) const
{
	_getIncrementInPlace(p1, p2,ptemp);
	tensor.applyInverseInPlace(ptemp,temp);
	double s=0;
	double ti;
	for(int i = 0;i<(int)temp.size();i++)
	{
		ti = temp[i];
		s+= ti * ti;
	}
	return sqrt(s);

}

void SpaceRN::_getIncrementInPlace(const SpacePoint& p1,
							 const SpacePoint & p2,
							 VectorDouble& ptemp)const
{

	for(int i = 0; i<(int)ptemp.size();i++)
	{

		ptemp[i] = p2.getCoord(i)-p1.getCoord(i);
	}

}

void SpaceRN::norm(const SpacePoint& p1, const std::vector<SpacePoint>& p2v,const Tensor& tens, VectorDouble& res,VectorDouble & temp,VectorDouble& out) const
{
	double t;
	for(int i = 0; i < (int)res.size(); i++)
	{
		for (int idim = 0; idim < (int)temp.size();idim++)
			temp[i] = p1.getCoord(idim) - p2v[i].getCoord(idim);
		tens.applyInverseInPlace(temp, out);
		res[i] = 0.;
		for (int idim = 0;idim <(int)temp.size();idim++)
		{
			t = out[idim];
			res[i]+= t * t;
		}
			res[i] = sqrt(res[i]);
	}
}
double SpaceRN::getFrequentialDistance(const SpacePoint& p1,
                                       const SpacePoint& p2,
                                       const Tensor& tensor) const
{
  return VH::norm(tensor.applyDirect(getIncrement(p1, p2),0));
}

VectorDouble SpaceRN::getIncrement(const SpacePoint& p1,
                                   const SpacePoint& p2) const
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

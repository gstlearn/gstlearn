/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "Boolean/TokenHalfEllipsoid.hpp"
#include "Boolean/Object.hpp"

TokenHalfEllipsoid::TokenHalfEllipsoid(double proportion,
                                       double xext,
                                       double yext,
                                       double zext,
                                       double theta)
    : AToken()
{
  initParams(getNParams());
  setParamDefault(0, "X-Extension", xext);
  setParamDefault(1, "Y-Extension", yext);
  setParamDefault(2, "Z-Extension", zext);
  setParamDefault(3, "Orientation Angle", theta);
  setProportion(proportion);
}

TokenHalfEllipsoid::TokenHalfEllipsoid(const TokenHalfEllipsoid &r)
    : AToken(r)
{
}

TokenHalfEllipsoid& TokenHalfEllipsoid::operator=(const TokenHalfEllipsoid &r)
{
  if (this != &r)
  {
    AToken::operator =(r);
  }
  return *this;
}

TokenHalfEllipsoid::~TokenHalfEllipsoid()
{
}

/*****************************************************************************/
/*!
 **  Generate the geometry of the object
 **
 ** \param[in]  ndim    Space dimension
 **
 *****************************************************************************/
Object* TokenHalfEllipsoid::generateObject(int ndim)

{
  Object* object = new Object(this);
  if (ndim >= 1) object->setExtension(0, generateParam(0));
  if (ndim >= 2) object->setExtension(1, generateParam(1));
  if (ndim >= 3) object->setExtension(2, generateParam(2));
  object->setOrientation(generateParam(3));
  return object;
}

/****************************************************************************/
/*!
 **  Check if the pixel (x,y,z) belongs to the object
 **
 ** \return  1 if the pixel is in the grain, 0 if it is in the pore
 **
 *****************************************************************************/
bool TokenHalfEllipsoid::belongObject(const VectorDouble& coor,
                                      const Object* object) const
{
  int ndim = (int) coor.size();
  double dx = (ndim >= 1) ? coor[0] / (object->getExtension(0) / 2.) : 0.;
  double dy = (ndim >= 2) ? coor[1] / (object->getExtension(1) / 2.) : 0.;
  double dz = (ndim >= 3) ? coor[2] / (object->getExtension(2))      : 0.;
  if (dx * dx + dy * dy + dz * dz > 1) return false;
  return true;
}


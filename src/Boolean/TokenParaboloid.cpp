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
#include "Boolean/TokenParaboloid.hpp"
#include "Boolean/Object.hpp"

TokenParaboloid::TokenParaboloid()
    : AToken()
{
  initParams(4);
  setParamName(0, "X-Extension");
  setParamName(1, "Y-Extension");
  setParamName(2, "Z-Extension");
  setParamName(3, "Orientation Angle");
}

TokenParaboloid::TokenParaboloid(const TokenParaboloid &r)
    : AToken(r)
{
}

TokenParaboloid& TokenParaboloid::operator=(const TokenParaboloid &r)
{
  if (this != &r)
  {
    AToken::operator =(r);
  }
  return *this;
}

TokenParaboloid::~TokenParaboloid()
{
}

/*****************************************************************************/
/*!
 **  Generate the geometry of the object
 **
 ** \param[in]  ndim    Space dimension
 **
 *****************************************************************************/
Object* TokenParaboloid::generateObject(int ndim)

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
bool TokenParaboloid::belongObject(const VectorDouble& coor,
                                   const Object* object) const
{
  int ndim = (int) coor.size();
  double dx = (ndim >= 1) ? coor[0] / (object->getExtension(0) / 2.) : 0.;
  double dy = (ndim >= 2) ? coor[1] / (object->getExtension(1) / 2.) : 0.;
  double dz = (ndim >= 3) ? coor[2] / (object->getExtension(2))      : 0.;
  if (dx * dx + dy * dy - dz > 1) return false;
  if (dx * dx + dy * dy + dz > 1) return false;
  return true;
}

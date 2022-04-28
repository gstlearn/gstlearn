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
#include "Boolean/TokenHalfSinusoid.hpp"
#include "Boolean/Object.hpp"

#include <math.h>

TokenHalfSinusoid::TokenHalfSinusoid()
    : AToken()
{
  initParams(4);
  setParamName(0, "X-Extension");
  setParamName(1, "Y-Extension");
  setParamName(2, "Z-Extension");
  setParamName(3, "Orientation Angle");
}

TokenHalfSinusoid::TokenHalfSinusoid(const TokenHalfSinusoid &r)
    : AToken(r)
{
}

TokenHalfSinusoid& TokenHalfSinusoid::operator=(const TokenHalfSinusoid &r)
{
  if (this != &r)
  {
    AToken::operator =(r);
  }
  return *this;
}

TokenHalfSinusoid::~TokenHalfSinusoid()
{
}

/*****************************************************************************/
/*!
 **  Generate the geometry of the object
 **
 ** \param[in]  ndim    Space dimension
 **
 *****************************************************************************/
Object* TokenHalfSinusoid::generateObject(int ndim)

{
  Object* object = new Object(this);
  if (ndim >= 1) object->setValue(0, generateParam(0));
  if (ndim >= 2) object->setValue(1, generateParam(1));
  if (ndim >= 3) object->setValue(2, generateParam(2));
  if (ndim >= 1) object->setExtension(0, generateParam(3));
  if (ndim >= 2) object->setExtension(1, object->getValue(1) + object->getValue(2));
  if (ndim >= 3) object->setExtension(2, generateParam(4));
  object->setOrientation(generateParam(5));
  return object;
}

/****************************************************************************/
/*!
 **  Check if the pixel (x,y,z) belongs to the object
 **
 ** \return  1 if the pixel is in the grain, 0 if it is in the pore
 **
 *****************************************************************************/
bool TokenHalfSinusoid::belongObject(const VectorDouble& coor,
                                     const Object* object) const
{
  int ndim = (int) coor.size();
  double dx = (ndim >= 1) ? coor[0] / object->getValue(0) : 0.;
  double dz = (ndim >= 3) ? coor[2] / object->getExtension(2) : 0.;
  double yloc = object->getValue(1) * cos(2. * GV_PI * dx) / 2.;
  double dy = (ndim >= 2) ? (coor[1] - yloc) / (object->getValue(2) / 2.) :  0.;

  if (dx * dx + dy * dy + dz * dz > 1) return false;
  return true;
}

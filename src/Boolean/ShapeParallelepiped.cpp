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
#include "Boolean/ShapeParallelepiped.hpp"
#include "Simulation/BooleanObject.hpp"

ShapeParallelepiped::ShapeParallelepiped(double proportion,
                                         double xext,
                                         double yext,
                                         double zext,
                                         double theta)
    : AShape()
{
  initParams(getNParams());
  setParamDefault(0, "X-Extension", xext);
  setParamDefault(1, "Y-Extension", yext);
  setParamDefault(2, "Z-Extension", zext);
  setParamDefault(3, "Orientation Angle", theta);
  setProportion(proportion);
}

ShapeParallelepiped::ShapeParallelepiped(const ShapeParallelepiped &r)
    : AShape(r)
{
}

ShapeParallelepiped& ShapeParallelepiped::operator=(const ShapeParallelepiped &r)
{
  if (this != &r)
  {
    AShape::operator =(r);
  }
  return *this;
}

ShapeParallelepiped::~ShapeParallelepiped()
{
}

/*****************************************************************************/
/*!
 **  Generate the geometry of the object
 **
 ** \param[in]  ndim    Space dimension
 **
 *****************************************************************************/
BooleanObject* ShapeParallelepiped::generateObject(int ndim)

{
  BooleanObject* object = new BooleanObject(this);
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
bool ShapeParallelepiped::belongObject(const VectorDouble& /*coor*/,
                                       const BooleanObject* /*object*/) const
{
  return true;
}


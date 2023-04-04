/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Boolean/ShapeParaboloid.hpp"
#include "Simulation/BooleanObject.hpp"

ShapeParaboloid::ShapeParaboloid(double proportion,
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

ShapeParaboloid::ShapeParaboloid(const ShapeParaboloid &r)
    : AShape(r)
{
}

ShapeParaboloid& ShapeParaboloid::operator=(const ShapeParaboloid &r)
{
  if (this != &r)
  {
    AShape::operator =(r);
  }
  return *this;
}

ShapeParaboloid::~ShapeParaboloid()
{
}

/*****************************************************************************/
/*!
 **  Generate the geometry of the object
 **
 ** \param[in]  ndim    Space dimension
 **
 *****************************************************************************/
BooleanObject* ShapeParaboloid::generateObject(int ndim)

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
bool ShapeParaboloid::belongObject(const VectorDouble& coor,
                                   const BooleanObject* object) const
{
  int ndim = (int) coor.size();
  double dx = (ndim >= 1) ? coor[0] / (object->getExtension(0) / 2.) : 0.;
  double dy = (ndim >= 2) ? coor[1] / (object->getExtension(1) / 2.) : 0.;
  double dz = (ndim >= 3) ? coor[2] / (object->getExtension(2))      : 0.;
  if (dx * dx + dy * dy - dz > 1) return false;
  if (dx * dx + dy * dy + dz > 1) return false;
  return true;
}

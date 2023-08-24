/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "Boolean/ShapeHalfSinusoid.hpp"
#include "Simulation/BooleanObject.hpp"

#include <math.h>

ShapeHalfSinusoid::ShapeHalfSinusoid(double proportion,
                                     double period,
                                     double amplitude,
                                     double thickness,
                                     double xext,
                                     double zext,
                                     double theta)
    : AShape()
{
  initParams(getNParams());
  setParamDefault(0, "Period", period);
  setParamDefault(1, "Amplitude", amplitude);
  setParamDefault(2, "Thickness", thickness);
  setParamDefault(3, "X-Extension", xext);
  setParamDefault(4, "Z-Extension", zext);
  setParamDefault(5, "Orientation Angle", theta);
  setProportion(proportion);
}

ShapeHalfSinusoid::ShapeHalfSinusoid(const ShapeHalfSinusoid &r)
    : AShape(r)
{
}

ShapeHalfSinusoid& ShapeHalfSinusoid::operator=(const ShapeHalfSinusoid &r)
{
  if (this != &r)
  {
    AShape::operator =(r);
  }
  return *this;
}

ShapeHalfSinusoid::~ShapeHalfSinusoid()
{
}

/*****************************************************************************/
/*!
 **  Generate the geometry of the object
 **
 ** \param[in]  ndim    Space dimension
 **
 *****************************************************************************/
BooleanObject* ShapeHalfSinusoid::generateObject(int ndim)

{
  BooleanObject* object = new BooleanObject(this);
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
bool ShapeHalfSinusoid::belongObject(const VectorDouble& coor,
                                     const BooleanObject* object) const
{
  int ndim = (int) coor.size();
  double dx = (ndim >= 1) ? coor[0] / object->getValue(0) : 0.;
  double dz = (ndim >= 3) ? coor[2] / object->getExtension(2) : 0.;
  double yloc = object->getValue(1) * cos(2. * GV_PI * dx) / 2.;
  double dy = (ndim >= 2) ? (coor[1] - yloc) / (object->getValue(2) / 2.) :  0.;

  if (dx * dx + dy * dy + dz * dz > 1) return false;
  return true;
}

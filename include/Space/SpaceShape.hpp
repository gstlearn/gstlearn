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
#pragma once

#include "gstlearn_export.hpp"

/// TODO : to be kept ?

/****************************************************************************
**
****************************************************************************/
class GSTLEARN_EXPORT ASpaceShape
{
  public:
  ASpaceShape(){}
  virtual ~ASpaceShape(){}
  
  //bool is_inside(SpacePoint pt) const = 0;
};

/****************************************************************************
** Description of a Cone
****************************************************************************/
class GSTLEARN_EXPORT Cone : public ASpaceShape
{
  public:
    Cone() : ASpaceShape(),angle(90){}
    ~Cone(){}
    //bool is_inside(SpacePoint pt) const override;
    double angle;
};

/****************************************************************************
**
****************************************************************************/
class GSTLEARN_EXPORT Cylinder : public ASpaceShape
{
  public:
    Cylinder() : ASpaceShape(),radius(TEST){}
    ~Cylinder(){}
    //bool is_inside(SpacePoint pt) const override;
    double radius;
};

/****************************************************************************
**
****************************************************************************/
class GSTLEARN_EXPORT Pencil : public ASpaceShape
{
  public:

    Pencil() : ASpaceShape(), angle(90), radius(TEST){}
    ~Pencil(){}
    Pencil(const Pencil& ref) : ASpaceShape(ref), angle(ref.angle),radius(ref.radius){}

    Pencil& operator=(const Pencil &ref)
    {
      if (this != &ref)
      {
        angle = ref.angle;
        radius= ref.radius;
      }
      return(*this);
    }
    //bool is_inside(SpacePoint pt) const override{};
    void setAngle(double ang){angle = ang;}
    double angle;
    double radius;
};


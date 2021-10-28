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
#pragma once

/// TODO : to be kept ?

/****************************************************************************
**
****************************************************************************/
class ASpaceShape
{
  public:
  ASpaceShape(){}
  virtual ~ASpaceShape(){}
  
  //bool is_inside(SpacePoint pt) const = 0;
};

/****************************************************************************
** Description of a Cone
****************************************************************************/
class Cone : public ASpaceShape
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
class Cylinder : public ASpaceShape
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
class Pencil : public ASpaceShape
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


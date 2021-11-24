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
#include "Model/ConsItem.hpp"
#include "Basic/Utilities.hpp"
#include "geoslib_f.h"
#include "geoslib_old_f.h"


ConsItem::ConsItem(bool authAssign,
                   const EConsType& icase,
                   int igrf,
                   int icov,
                   const EConsElem& type,
                   int iv1,
                   int iv2,
                   double value)
    : _icase(icase),
      _igrf(igrf),
      _icov(icov),
      _elemType(type),
      _iv1(iv1),
      _iv2(iv2),
      _value(value),
      _authAssign(authAssign)
{
}

ConsItem::ConsItem(const ConsItem &m)
    : _icase(m._icase),
      _igrf(m._igrf),
      _icov(m._icov),
      _elemType(m._elemType),
      _iv1(m._iv1),
      _iv2(m._iv2),
      _value(m._value),
      _authAssign(m._authAssign)
{

}

ConsItem& ConsItem::operator=(const ConsItem &m)
{
  if (this != &m)
  {
    _icase = m._icase;
    _igrf = m._igrf;
    _icov = m._icov;
    _elemType = m._elemType;
    _iv1   = m._iv1;
    _iv2   = m._iv2;
    _value = m._value;
    _authAssign = m._authAssign;
  }
  return *this;
}

ConsItem::~ConsItem()
{

}

int ConsItem::init(const EConsType& icase,
                   int igrf,
                   int icov,
                   const EConsElem& type,
                   int v1,
                   int v2,
                   double value)
{
  _icase = icase;
  _igrf  = igrf;
  _icov  = icov;
  _elemType  = type;
  _iv1   = v1;
  _iv2   = v2;
  _value = value;

  // Check to avoid rotation of a Model defined on the sphere
  int flag_sphere;
  variety_query(&flag_sphere);
  if (flag_sphere && type == EConsElem::ANGLE)
  {
    messerr("When working on the Sphere Geometry");
    messerr("Rotation must be specified using 'I' constraints (not 'A')");
    return 1;
  }
  return 0;
}

String ConsItem::toString(int /*level*/) const
{
  std::stringstream sstr;

  if (_authAssign)
  {
    switch (_icase.toEnum())
    {
      case EConsType::E_LOWER:
        sstr << "Constraint Type = Lower Bound" << std::endl;
        break;

      case EConsType::E_DEFAULT:
        sstr << "Constraint Type = Default Parameter" << std::endl;
        break;

      case EConsType::E_UPPER:
        sstr << "Constraint Type = Upper Bound" << std::endl;
        break;

      case EConsType::E_EQUAL:
        sstr << "Constraint Type = Equality" << std::endl;
        break;

      default:
        sstr << "Constraint Type = UNKNOWN!!" << std::endl;
        break;
    }
  }

  switch (_elemType.toEnum())
  {
    case EConsElem::E_RANGE:
      sstr << "Range     ";
      break;

    case EConsElem::E_ANGLE:
      sstr << "Angle     ";
      break;

    case EConsElem::E_PARAM:
      sstr << "Param     ";
      break;

    case EConsElem::E_SILL:
      sstr << "Sill      ";
      break;

    case EConsElem::E_SCALE:
      sstr << "Scale     ";
      break;

    case EConsElem::E_T_RANGE:
      sstr << "Tapering  ";
      break;

    case EConsElem::E_VELOCITY:
      sstr << "Velocity  ";
      break;

    case EConsElem::E_SPHEROT:
      sstr << "S-Rotation";
      break;

    case EConsElem::E_TENSOR:
      sstr << "Anis-Matrix";
      break;

    default:
      break;
  }
  sstr << " :";
  sstr << " GRF=" << _igrf + 1;
  sstr << " Str=" << _icov + 1;
  sstr << " V#1=" << _iv1  + 1;
  sstr << " V#2=" << _iv2  + 1;

  if (_authAssign)
  {
    if (FFFF(_value))
      sstr << " Value=NA" << std::endl;
    else
      sstr << " Value=" << _value << std::endl;
  }
  sstr << std::endl;
  return sstr.str();
}

IClonable* ConsItem::clone() const
{
  return new ConsItem(*this);
}

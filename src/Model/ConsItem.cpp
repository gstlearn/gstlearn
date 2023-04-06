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
#include "geoslib_f.h"
#include "geoslib_old_f.h"

#include "Enum/EConsElem.hpp"

#include "Space/ASpaceObject.hpp"
#include "Model/ConsItem.hpp"
#include "Basic/Utilities.hpp"

ConsItem::ConsItem(const CovParamId& paramid,
                   const EConsType& type,
                   double value)
    : AStringable(),
      _paramId(),
      _type(),
      _value()
{
  (void) _init(paramid, type, value);
}

ConsItem::ConsItem(const ConsItem &m)
    : AStringable(m),
      _paramId(m._paramId),
      _type(m._type),
      _value(m._value)
{
}

ConsItem& ConsItem::operator=(const ConsItem &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _paramId = m._paramId;
    _type = m._type;
    _value = m._value;
  }
  return *this;
}

ConsItem::~ConsItem()
{

}

ConsItem* ConsItem::create(const CovParamId &paramid,
                           const EConsType &type,
                           double value)
{
  return new ConsItem(paramid, type, value);
}

ConsItem* ConsItem::createFromParamId(int icov,
                                      const EConsElem &elem,
                                      const EConsType &type,
                                      double value,
                                      int igrf,
                                      int iv1,
                                      int iv2)
{
  CovParamId paramid(igrf, icov, elem, iv1, iv2);
  return new ConsItem(paramid, type, value);
}

int ConsItem::_init(const CovParamId& paramid,
                   const EConsType& type,
                   double value)
{
  _paramId = paramid;
  _type  = type;
  _value = value;

  // Check to avoid rotation of a Model defined on the sphere
  bool flag_sphere = (getDefaultSpaceType() == ESpaceType::SN);
  if (flag_sphere && type == EConsElem::ANGLE)
  {
    messerr("When working on the Sphere Geometry");
    messerr("Rotation must be specified using 'I' constraints (not 'A')");
    return 1;
  }
  return 0;
}

String ConsItem::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  switch (_type.toEnum())
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

  sstr << _paramId.toString(strfmt);

  if (FFFF(_value))
    sstr << " Value=NA" << std::endl;
  else
    sstr << " Value=" << _value << std::endl;
  sstr << std::endl;
  return sstr.str();
}

ConsItem ConsItem::define(const EConsElem& elem,
                          int icov,
                          int iv1,
                          int iv2,
                          const EConsType& type,
                          double value)
{
  CovParamId parid(0, icov, elem, iv1, iv2);
  return ConsItem(parid, type, value);
}

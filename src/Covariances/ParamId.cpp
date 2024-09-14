/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Space/ASpaceObject.hpp"
#include "Covariances/ParamId.hpp"

ParamId::ParamId(const EConsElem& elem,
                 int iv1,
                 int iv2)
    : AStringable(),
      _elemType(elem),
      _iv1(iv1),
      _iv2(iv2)
{
}

ParamId::ParamId(const ParamId &m)
    : AStringable(m),
      _elemType(m._elemType),
      _iv1(m._iv1),
      _iv2(m._iv2)
{

}

ParamId& ParamId::operator=(const ParamId &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _elemType = m._elemType;
    _iv1   = m._iv1;
    _iv2   = m._iv2;
  }
  return *this;
}

ParamId::~ParamId()
{

}

ParamId* ParamId::create(const EConsElem &elem,
                         int iv1,
                         int iv2)
{
  return new ParamId(elem, iv1, iv2);
}

int ParamId::init(const EConsElem& type,
                  int v1,
                  int v2)
{
  _elemType  = type;
  _iv1   = v1;
  _iv2   = v2;

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

String ParamId::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  switch (_elemType.toEnum())
  {
    case EConsElem::E_RANGE:
      sstr << "Range      :"
           << " IDir=" << _iv1 + 1;
      break;

    case EConsElem::E_ANGLE:
      sstr << "Angle      :"
           << " IdAngle=" << _iv1 + 1;
      break;

    case EConsElem::E_PARAM:
      sstr << "Param      ";
      break;

    case EConsElem::E_SILL:
      sstr << "Sill       :"
           << " Ivar=" << _iv1 << " Jvar=" << _iv2;
      break;

    case EConsElem::E_SCALE:
      sstr << "Scale      :"
           << " IDir=" << _iv1 + 1;
      break;

    case EConsElem::E_T_RANGE:
      sstr << "Tapering   :"
           << " IDir=" << _iv1 + 1;
      break;

    case EConsElem::E_VELOCITY:
      sstr << "Velocity   :"
           << " Ivar=" << _iv1 + 1
           << " Jvar=" << _iv2 + 1;
      break;

    case EConsElem::E_SPHEROT:
      sstr << "S-Rotation :"
           << " Ivar=" << _iv1 + 1
           << " Jvar=" << _iv2 + 1;
      break;

    case EConsElem::E_TENSOR:
      sstr << "Anis-Matrix :"
           << " Ivar=" << _iv1 + 1
           << " Jvar=" << _iv2 + 1;
      break;

    default:
      return sstr.str();
  }
  sstr << std::endl;

  return sstr.str();
}

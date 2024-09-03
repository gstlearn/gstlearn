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
#include "Model/CovParamId.hpp"

CovParamId::CovParamId(int igrf,
                       int icov,
                       const EConsElem& elem,
                       int iv1,
                       int iv2)
    : AStringable(),
      _igrf(igrf),
      _icov(icov),
      _elemType(elem),
      _iv1(iv1),
      _iv2(iv2)
{
}

CovParamId::CovParamId(const CovParamId &m)
    : AStringable(m),
      _igrf(m._igrf),
      _icov(m._icov),
      _elemType(m._elemType),
      _iv1(m._iv1),
      _iv2(m._iv2)
{

}

CovParamId& CovParamId::operator=(const CovParamId &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _igrf = m._igrf;
    _icov = m._icov;
    _elemType = m._elemType;
    _iv1   = m._iv1;
    _iv2   = m._iv2;
  }
  return *this;
}

CovParamId::~CovParamId()
{

}

CovParamId* CovParamId::create(int igrf,
                               int icov,
                               const EConsElem &elem,
                               int iv1,
                               int iv2)
{
  return new CovParamId(igrf, icov, elem, iv1, iv2);
}

int CovParamId::init(int igrf,
                     int icov,
                     const EConsElem& type,
                     int v1,
                     int v2)
{
  _igrf  = igrf;
  _icov  = icov;
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

String CovParamId::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  switch (_elemType.toEnum())
  {
    case EConsElem::E_RANGE:
      sstr << "Range      :" << "IdCov=" << _icov+1 << "IDir=" << _iv1+1;
      break;

    case EConsElem::E_ANGLE:
      sstr << "Angle      :" << "IdCov=" << _icov+1 << "IdAngle=" << _iv1;
      break;

    case EConsElem::E_PARAM:
      sstr << "Param      :" << "IdCov" << _icov+1;
      break;

    case EConsElem::E_SILL:
      sstr << "Sill       :" << "IdCov=" << _icov+1 << "Ivar=" << _iv1 << "Jvar=" << _iv2;
      break;

    case EConsElem::E_SCALE:
      sstr << "Scale      :" << "IdCov=" << _icov + 1 << "IDir=" << _iv1 + 1;
      break;

    case EConsElem::E_T_RANGE:
      sstr << "Tapering   :" << "IdCov=" << _icov + 1 << "IDir=" << _iv1 + 1;
      break;

    case EConsElem::E_VELOCITY:
      sstr << "Velocity   :" << "IdCov=" << _icov + 1 << "Ivar=" << _iv1 << "Jvar=" << _iv2;
      break;

    case EConsElem::E_SPHEROT:
      sstr << "S-Rotation :" << "IdCov=" << _icov + 1 << "Ivar=" << _iv1 << "Jvar=" << _iv2;
      break;

    case EConsElem::E_TENSOR:
      sstr << "Anis-Matrix :" << "IdCov=" << _icov + 1 << "Ivar=" << _iv1 << "Jvar=" << _iv2;
      break;

    default:
      return sstr.str();
  }
  if (_igrf > 0)
    sstr << " (GRF=" << _igrf + 1 << ")";
  sstr << std::endl;

  return sstr.str();
}


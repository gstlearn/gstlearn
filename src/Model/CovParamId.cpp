/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "geoslib_f.h"
#include "geoslib_old_f.h"

#include "Space/ASpaceObject.hpp"
#include "Model/CovParamId.hpp"
#include "Basic/Utilities.hpp"

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
  bool flag_V2 = true;

  switch (_elemType.toEnum())
  {
    case EConsElem::E_RANGE:
      flag_V2 = false;
      sstr << "Range     ";
      break;

    case EConsElem::E_ANGLE:
      flag_V2 = false;
      sstr << "Angle     ";
      break;

    case EConsElem::E_PARAM:
      flag_V2 = false;
      sstr << "Param     ";
      break;

    case EConsElem::E_SILL:
      flag_V2 = true;
      sstr << "Sill      ";
      break;

    case EConsElem::E_SCALE:
      flag_V2 = false;
      sstr << "Scale     ";
      break;

    case EConsElem::E_T_RANGE:
      flag_V2 = false;
      sstr << "Tapering  ";
      break;

    case EConsElem::E_VELOCITY:
      flag_V2 = true;  // TODO: to be validated
      sstr << "Velocity  ";
      break;

    case EConsElem::E_SPHEROT:
      flag_V2 = true;  // TODO: to be validated
      sstr << "S-Rotation";
      break;

    case EConsElem::E_TENSOR:
      flag_V2 = true;
      sstr << "Anis-Matrix";
      break;

    default:
      return sstr.str();
  }
  sstr << " :";
  sstr << " GRF=" << _igrf + 1;
  sstr << " Str=" << _icov + 1;
  sstr << " V#1=" << _iv1  + 1;
  if (flag_V2) sstr << " V#2=" << _iv2  + 1;
  sstr << std::endl;

  return sstr.str();
}


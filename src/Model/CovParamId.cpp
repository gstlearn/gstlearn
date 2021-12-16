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
#include "Model/CovParamId.hpp"
#include "Basic/Utilities.hpp"
#include "geoslib_f.h"
#include "geoslib_old_f.h"

CovParamId::CovParamId(int igrf,
                       int icov,
                       const EConsElem& elem,
                       int iv1,
                       int iv2)
    : _igrf(igrf),
      _icov(icov),
      _elemType(elem),
      _iv1(iv1),
      _iv2(iv2)
{
}

CovParamId::CovParamId(const CovParamId &m)
    : _igrf(m._igrf),
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

String CovParamId::toString(int /*level*/) const
{
  std::stringstream sstr;

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

  return sstr.str();
}

IClonable* CovParamId::clone() const
{
  return new CovParamId(*this);
}

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
#include "Model/ElemNostat.hpp"
#include "Basic/Utilities.hpp"
#include "geoslib_enum.h"
#include "geoslib_f.h"

ElemNostat::ElemNostat()
    : _locType(EConsElem::UNKNOWN),
      _rankGRF(0),
      _rankStr(0),
      _rankV1(0),
      _rankV2(0),
      _val1(0),
      _val2(0)
{

}

ElemNostat::ElemNostat(const ElemNostat &m)
    : _locType(m._locType),
      _rankGRF(m._rankGRF),
      _rankStr(m._rankStr),
      _rankV1(m._rankV1),
      _rankV2(m._rankV2),
      _val1(m._val1),
      _val2(m._val2)
{

}

ElemNostat& ElemNostat::operator=(const ElemNostat &m)
{
  if (this != &m)
  {
    _locType = m._locType;
    _rankGRF = m._rankGRF;
    _rankStr = m._rankStr;
    _rankV1 = m._rankV1;
    _rankV2 = m._rankV2;
    _val1 = m._val1;
    _val2 = m._val2;
  }
  return *this;
}

ElemNostat::~ElemNostat()
{

}

void ElemNostat::init(const EConsElem& loctype,
                      int rank_grf,
                      int rank_str,
                      int rank_v1,
                      int rank_v2)
{
  _locType = loctype;
  _rankGRF = rank_grf;
  _rankStr = rank_str;
  _rankV1 = rank_v1;
  _rankV2 = rank_v2;
  _val1 = TEST;
  _val2 = TEST;
}

String ElemNostat::toString(int level) const
{
  std::stringstream sstr;
  switch (getLocType().toEnum())
  {
    case EConsElem::E_RANGE:
      sstr << "Type = Range";
      break;

    case EConsElem::E_SCALE:
      sstr << "Type = Scale";
      break;

    case EConsElem::E_ANGLE:
      sstr << "Type = Angle";
      break;

    case EConsElem::E_PARAM:
      sstr << "Type = Third";
      break;

    case EConsElem::E_SILL:
      sstr << "Type = Sill";
      break;

    case EConsElem::E_SPHEROT:
      sstr << "Type = Spherical Rotation";
      break;

    case EConsElem::E_T_RANGE:
      sstr << "Type = T Range";
      break;

    case EConsElem::E_VELOCITY:
      sstr << "Type = Velocity";
      break;

    case EConsElem::E_ROTMAT:
       sstr << "Anis-Matrix";
       break;

    default:
      sstr << "Type = UNKNOWN";
      break;
  }
  sstr << " - Structure = " << getRankStr() + 1 << " - Variable = "
       << getRankV1() + 1;

  if (!IFFFF(getRankV2())) sstr << " - " << getRankV2() + 1;

  return sstr.str();
}

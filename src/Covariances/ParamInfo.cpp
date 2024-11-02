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


#include "Covariances/ParamInfo.hpp"
#include "Basic/AStringable.hpp"
#include "geoslib_define.h"

ParamInfo::ParamInfo(double valuser, double valmin, double valmax, bool mute)
  : _valUser(valuser),
    _valMin(valmin),
    _valMax(valmax),
    _valMinUser(valmin),
    _valMaxUser(valmax),
    _mutable(mute)
{
}

ParamInfo::ParamInfo(const ParamInfo& m)
  : AStringable(m)
{
  this->_valUser = m._valUser;
  this->_valMin  = m._valMin;
  this->_valMax  = m._valMax;
  this->_valMinUser = m._valMinUser;
  this->_valMaxUser = m._valMaxUser;  
  this->_mutable = m._mutable;
}

ParamInfo& ParamInfo::operator=(const ParamInfo& m)
{
  if (this != &m)
  {
    _valUser = m._valUser;
    _valMin = m._valMin;
    _valMax = m._valMax;
    _valMinUser = m._valMinUser;
    _valMaxUser = m._valMaxUser;
    _mutable = m._mutable;
  }
  return *this;
}




String ParamInfo::toString(const AStringFormat* strfmt) const
{
    DECLARE_UNUSED(strfmt)
    std::stringstream sstr;
  
    sstr << toTitle(1, "Information on covariance parameter");
    
    sstr << "User value = " << std::to_string(_valUser) << " - ";
    sstr << "Min value = " << std::to_string(_valMin) << " - "; 
    sstr << "Max value = " << std::to_string(_valMax) << " - ";
    sstr << "Parameter is ";
    if (!_mutable)
        sstr << "not ";
    sstr << "mutable";
    return sstr.str();
}

ParamInfo::~ParamInfo()
{

}

void ParamInfo::setValMax(double valmax)
{
  if (valmax > _valMax)
  {
    messerr("The maximum value of the parameter is greater than ");
    messerr("the maximum value of the parameter in the model");
    messerr("which is fixed to %f", _valMax);
    return;
  }
    _valMaxUser = valmax;
}

void ParamInfo::setValMin(double valmin)
{
  if (valmin < _valMin)
  {
    messerr("The minimum value of the parameter is lower than ");
    messerr("the minimum value of the parameter in the model");
    messerr("which is fixed to %f", _valMin);
    return;
  }
    _valMinUser = valmin;
}
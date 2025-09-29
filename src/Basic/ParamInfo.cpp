#include "Basic/ParamInfo.hpp"
#include "Basic/AStringable.hpp"
#include "geoslib_define.h"

#include <sstream>
namespace gstlrn
{ 
ParamInfo::ParamInfo(const String& name,
                     double value,
                     const std::array<double, 2>& absoluteBounds,
                     const String& description,
                     bool isfixed)
  : AStringable()
  , _name(name)
  , _value(value)
  , _currentValue(value)
  , _absoluteBounds(absoluteBounds)
  , _userBounds(absoluteBounds)
  , _isFixed(isfixed)
  , _description(description)
  , _address(ITEST)
{
}

ParamInfo::ParamInfo(const ParamInfo& other)
  : AStringable(other)
  , _name(other._name)
  , _value(other._value)
  , _currentValue(other._currentValue)
  , _absoluteBounds(other._absoluteBounds)
  , _userBounds(other._userBounds)
  , _isFixed(other._isFixed)
  , _description(other._description)
  , _address(ITEST)
{
}

ParamInfo& ParamInfo::operator=(const ParamInfo& other)
{
  if (this != &other)
  {
    AStringable::operator=(other);
    _name           = other._name;
    _value          = other._value;
    _currentValue   = other._currentValue;
    _absoluteBounds = other._absoluteBounds;
    _userBounds     = other._userBounds;
    _isFixed        = other._isFixed;
    _description    = other._description;
    _address        = ITEST;
  }
  return *this;
}

ParamInfo::~ParamInfo()
{
}

void ParamInfo::increaseMin(double value)
{
  _userBounds[0] = std::max(value, _userBounds[0]);
  _value         = std::max(_userBounds[0], _value);
  _currentValue  = _value;
}


void ParamInfo::decreaseMax(double value)
{
  _userBounds[1] = std::min(value, _userBounds[1]);
  _value         = std::min(_userBounds[1], _value);
  _currentValue  = _value;
}

String ParamInfo::toString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);
  std::stringstream sstr;
  sstr << _description << std::endl;

  sstr << "    Value: " << std::to_string(_value);
  if (_isFixed) sstr << " (fixed)";
  sstr << std::endl;

  sstr << "    Absolute Bounds: ";
  for (const auto& bound: _absoluteBounds)
  {
    sstr << bound << " ";
  }
  sstr << std::endl;

  sstr << "    User Bounds: ";
  for (const auto& bound: _userBounds)
  {
    sstr << bound << " ";
  }
  return sstr.str();
}

void ParamInfo::setMinValue(double value)
{
  if (value < _absoluteBounds[0])
  {
    _userBounds[0] = value;
    _value         = std::max(_userBounds[0], _value);
    _currentValue  = _value;
  }
  else
  {
    messerr("Warning for parameter: %s", _name.c_str());
    messerr("Minimum value %f is less than the minimum authorized value %f", value, _absoluteBounds[0]);
    messerr("The minimum user value is set to the minimum authorized value");
    _userBounds[0] = _absoluteBounds[0];
  }
}

void ParamInfo::setValue(double value)
{

  if (value < _userBounds[0] || value > _userBounds[1])
  {
    String parenthesisleft  = "[";
    String parenthesisright = "]";
    if (_userBounds[0] == -INF)
      parenthesisleft = "(";
    if (_userBounds[1] == INF)
      parenthesisright = ")";
    messerr("Warning for parameter: %s", _name.c_str());
    messerr("Input value: %f", value);
    messerr(" is out of user bounds: %s %f, %f %s", parenthesisleft.c_str(), _userBounds[0], _userBounds[1], parenthesisright.c_str());
    messerr("The value is set to the closest bound");
    _value = std::clamp(value, _userBounds[0], _userBounds[1]);
  }
  else
  {
    _value = value;
  }
  _currentValue = _value;
}

void ParamInfo::setMaxValue(double value)
{
  if (value < _absoluteBounds[1])
  {
    _userBounds[1] = value;
    _value         = std::min(_userBounds[1], _value);
    _currentValue  = _value;
  }
  else
  {
    messerr("Warning for parameter: %s", _name.c_str());
    messerr("Maximum value %f is greater than the maximum authorized value %f", value, _absoluteBounds[1]);
    messerr("The maximum user value is set to the maximum authorized value");
    _userBounds[1] = value;
  }
}
}
#include "Basic/ParamInfo.hpp"
#include "Basic/AStringable.hpp"
#include "geoslib_define.h"

#include <sstream>


ParamInfo::ParamInfo(const String& name, 
                    double value,
                    const std::array<double,2>& absoluteBounds,
                    const String& description)
: AStringable()
, _name(name)
, _value(value)
, _currentValue(value)
, _absoluteBounds(absoluteBounds)
, _userBounds(absoluteBounds)
, _isFixed(false)
, _description(description)
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
{
}

ParamInfo::~ParamInfo()
{

}

String ParamInfo::toString(const AStringFormat* strfmt) const {
    DECLARE_UNUSED(strfmt);
    std::stringstream sstr;
    sstr << " Description of parameter " << _name << std::endl;
    sstr <<  _description << std::endl;
    sstr << "  Value: " <<  std::to_string(_value)  << std::endl;
    sstr << "  Absolute Bounds: ";
    for (const auto& bound : _absoluteBounds) {
        sstr << bound << " ";
    }
    sstr << std::endl;
    sstr << "  User Bounds: ";
    for (const auto& bound : _userBounds) {
        sstr << bound << " ";
    }
    sstr << std::endl;
    sstr << "  Is Fixed: " << (_isFixed ? "true" : "false") << std::endl;
    return sstr.str();
}

void ParamInfo::setMinValue(double value)
{
    if (value < _absoluteBounds[0]) 
    {
        _userBounds[0] = value;
    }
    else 
    {
        messerr("Value is less than the minimum authorized value");
        messerr("Setting the minimum user value to the minimum authorized value");
        _userBounds[0] = _absoluteBounds[0];
    } 

}
void ParamInfo::setMaxValue(double value)
{
    if (value > _absoluteBounds[1]) 
    {
        _userBounds[1] = value;
    }
    else 
    {
        messerr("Value is greater than the maximum authorized value");
        messerr("Setting the maximum user value to the maximum authorized value");
        _userBounds[1] = value;
    } 
}
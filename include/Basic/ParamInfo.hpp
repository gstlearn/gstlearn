#pragma once

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/AStringFormat.hpp"

#include <array>
#include <limits>
#include <string>


#define INF std::numeric_limits<double>::infinity()

/**
 * @brief Definition of the generic parameter
 * 
 * @param _name Name assigned to the parameter
 * @param _value Value (initial) assigned to the parameter
 * @param _currentValue Value currently assigned to the parameter
 * @param _absoluteBounds Absolute bounds within which the parameter may vary
 * @param _userBounds Current bounds for the parameter
 * @param _isFixed True if the parameter value may not vary
 * @param _description String describing the role of this parameter
 */
class GSTLEARN_EXPORT ParamInfo : public AStringable {
public:
    ParamInfo(const String& name = "", 
              double value = TEST,
              const std::array<double,2>& absoluteBounds = {-INF, INF},
              const std::string& description = "",
              bool  isfixed = false);
    ParamInfo(const ParamInfo& other);
    ParamInfo& operator=(const ParamInfo& other);
    virtual ~ParamInfo();

    double getValue() const {return _value;};
    double getAbsoluteMinValue() const {return _absoluteBounds[0];};
    double getAbsoluteMaxValue() const {return _absoluteBounds[1];};

    void setMinValue(double value);
    void setMaxValue(double value);
    void decreaseMax(double value);
    void increaseMin(double value);
    void setValue(double value)
    {
      _value = value;
      _currentValue = value;
    }
    double getUserMin() const {return _userBounds[0];};
    double getUserMax() const {return _userBounds[1];};
    void setFixed(bool isFixed) {_isFixed = isFixed;};
    bool isFixed() const {return _isFixed;};

    String toString(const AStringFormat* strfmt = nullptr) const override;

private:
  String _name;
  double _value;
  mutable double _currentValue;
  std::array<double, 2> _absoluteBounds;
  std::array<double, 2> _userBounds;
  bool _isFixed;
  String _description;
};

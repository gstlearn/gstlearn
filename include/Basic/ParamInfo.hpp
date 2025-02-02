#pragma once

#include <string>
#include "Basic/AStringable.hpp"
#include "Basic/AStringFormat.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"


// Définition d'un type générique pour les paramètres

class GSTLEARN_EXPORT ParamInfo : public AStringable {
public:
    ParamInfo(const std::string& name ="", 
              double value = TEST,
              const std::array<double,2>& absoluteBounds = {0.0, 0.0},
              const std::string& description = "");
    ParamInfo(const ParamInfo& other);
    ParamInfo& operator=(const ParamInfo& other) = delete;
    virtual ~ParamInfo();

    double getAbsoluteMinValue() const {return _absoluteBounds[0];};
    double getAbsoluteMaxValue() const {return _absoluteBounds[1];};
    void setMinValue(double value);
    void setMaxValue(double value);

    String toString(const AStringFormat* strfmt = nullptr) const override;

private:
    String _name;
    double _value;
    mutable double _currentValue;
    const std::array<double,2> _absoluteBounds;
          std::array<double,2> _userBounds;
    bool  _isFixed;
    const String _description;
};


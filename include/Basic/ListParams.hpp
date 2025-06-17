#pragma once

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/AStringFormat.hpp"
#include "Basic/ParamInfo.hpp"
#include <vector>
#include <functional>

class GSTLEARN_EXPORT ListParams: public AStringable
{
public:
  ListParams();
  ListParams(const ListParams& other)            = delete;
  ListParams& operator=(const ListParams& other) = delete;
  virtual ~ListParams()                          = default;

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  void addParam(ParamInfo& param);
  void addParams(std::vector<ParamInfo>& params)
  {
    for (auto& param: params)
    {
      _params.push_back(param);
    }
  }
  double getValue(int index) const;
  void setValue(int index, double value);
  std::vector<double> getValues() const;
  void setValues(const std::vector<double>& values);
  std::vector<double> getMinValues() const;
  std::vector<double> getMaxValues() const;

private:
  std::vector<std::reference_wrapper<ParamInfo>> _params; // List of parameters
};

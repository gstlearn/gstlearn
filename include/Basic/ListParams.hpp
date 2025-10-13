#pragma once

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/AStringFormat.hpp"
#include "Basic/ParamInfo.hpp"
#include <cstddef>
#include <vector>
#include <functional>

namespace gstlrn
{
class GSTLEARN_EXPORT ListParams: public AStringable
{
public:
  ListParams();
  ListParams(const ListParams& other)            = delete;
  ListParams& operator=(const ListParams& other) = delete;
  virtual ~ListParams()                          = default;

  String toString(const AStringFormat* strfmt = nullptr) const override;

  void addParam(ParamInfo& param);
  void addParams(std::vector<ParamInfo>& params);
  void clear();
  double getValue(Id index) const;
  void setValue(Id index, double value);
  std::vector<double> getOptimizableValues() const;
  void setValues(const std::vector<double>& values);
  std::vector<double> getMinValues(double epsilon = 0.) const;
  std::vector<double> getMaxValues(double epsilon = 0.) const;
  void makeDispatchIndexFromDispatch();
  double getOptimizableValue(size_t index) const;

  size_t getNOptimizableParams() const { return _dispatchIndex.size(); }
  size_t getNParams() const { return _params.size(); }
  void updateDispatch();
  std::vector<size_t> getDispatch() const { return _dispatch; }
  std::vector<size_t> getDispatchIndex() const { return _dispatchIndex; }
private:
  
  // Example: The model has 6 parameters, but only 4 are used in the optimization.
  // [0, 1, 2, 3, 4, 5] with 1 = 3 and 4 = 5
  // _params will contain references to all 6 parameters (0, 1, 2, 3, 4, 5),
  // but the optimizable parameters will be [0, 1, 2, 4]
  // _dispatch = {0, 1, 2, 1, 3, 3} says how to get the value inside the vector of optimizable parameters
  // _dispatchIndex = {0, 1, 2, 4} 
  std::vector<std::reference_wrapper<ParamInfo>> _params; // List of parameters
  std::vector<size_t> _dispatch;
  std::vector<size_t> _dispatchIndex; // Indexes for dispatching parameters
};
}
#include "Basic/ListParams.hpp"
#include "Basic/AStringable.hpp"
#include "geoslib_define.h"
#include <cstddef>
#include <sstream>

ListParams::ListParams()
  : AStringable()
{
}

void ListParams::updateDispatch()
{
  _dispatch.clear();
  _dispatchIndex.clear();
  size_t nmax = 0;
  for (size_t i = 0; i < _params.size(); ++i)
  {
    size_t index = _params[i].get().getAddress();
    if (index > nmax)
    {
      nmax = nmax + 1;
      _dispatch.push_back(nmax);
    }
    else
    {
      _dispatch.push_back(index);
    }
    
  }
  makeDispatchIndexFromDispatch();
}

void ListParams::addParam(ParamInfo& param)
{
  if (param.isFixed()) return;
  _params.push_back(param);
  _dispatch.push_back(_params.size() - 1);
  _dispatchIndex.push_back(_params.size() - 1);
  param.setAddress(_params.size() - 1);
}

void ListParams::addParams(std::vector<ParamInfo>& params)
{
  for (auto& param: params)
  {
    addParam(param);
  }
}

void ListParams::clear()
{
  _params.clear();
  _dispatch.clear();
  _dispatchIndex.clear();
}

double ListParams::getValue(int index) const
{
  if (index < 0 || index >= static_cast<int>(_params.size()))
  {
    messerr("Index out of range in ListParams::getValue");
    return TEST;
  }
  return _params[index].get().getValue();
}

double ListParams::getOptimizableValue(size_t index) const
{
  if (index >= getNOptimizableParams())
  {
    messerr("Index out of range in ListParams::getOptimizableValue");
    return TEST;
  }

  return _params[_dispatchIndex[index]].get().getValue();
}

void ListParams::setValue(int index, double value)
{
  if (index < 0 || index >= static_cast<int>(_params.size()))
  {
    messerr("Index out of range in ListParams::setValue");
    return;
  }
  _params[index].get().setValue(value);
}


String ListParams::toString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);
  std::stringstream result;
  result << toTitle(1,"List of Parameters:");
  for (int ipar = 0, jpar = 0, npar = (int) _dispatchIndex.size(); ipar < npar; ipar++)
  {
    jpar++;
    result << jpar << " - " << _params[_dispatchIndex[ipar]].get().toString() << std::endl;
  }
  return result.str();
}


void ListParams::makeDispatchIndexFromDispatch()
{
  _dispatchIndex.clear();
  for (size_t i = 0; i < getNParams(); ++i)
  {
    if (std::find(_dispatchIndex.begin(), _dispatchIndex.end(), _dispatch[i]) == _dispatchIndex.end())
    {
      _dispatchIndex.push_back(_dispatch[i]);
    }
  }
}
std::vector<double> ListParams::getOptimizableValues() const
{
  size_t nparam = getNOptimizableParams();
  std::vector<double> values(nparam);
  for (size_t i = 0; i < nparam; ++i)
  {
    values[i] = getOptimizableValue(i);
  }
  return values;
}

std::vector<double> ListParams::getMinValues() const
{
  size_t nparam = _params.size();
  std::vector<double> values(nparam);
  for (size_t i = 0; i < nparam; ++i)
  {
    values[i] = _params[i].get().getUserMin();
  }
  return values;
}

std::vector<double> ListParams::getMaxValues() const
{
  size_t nparam = _params.size();
  std::vector<double> values(nparam);
  for (size_t i = 0; i < nparam; ++i)
  {
    values[i] = _params[i].get().getUserMax();
  }
  return values;
}

void ListParams::setValues(const std::vector<double>& values)
{
  size_t size = _dispatch.size();
  for (size_t i = 0; i < size; i++)
  {
    _params[i].get().setValue(values[_dispatch[i]]);
  }
}
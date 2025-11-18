#include "Basic/ListParams.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"
#include "geoslib_define.h"
#include <algorithm>
#include <cstddef>
#include <sstream>
#include <unordered_map>

namespace gstlrn
{

struct DSU
{
  std::vector<Id> parent;
  DSU(Id n)
    : parent(n)
  {
    for (Id i = 0; i < n; i++) parent[i] = i;
  }
  Id find(Id x)
  {
    return parent[x] == x ? x : parent[x] = find(parent[x]);
  }
  void unite(Id x, Id y)
  {
    Id rx = find(x), ry = find(y);
    if (rx != ry) parent[ry] = rx;
  }
};

void reindex(const VectorInt& v, std::vector<size_t>& result)
{
  if (v.empty())
  {
    result.clear();
    return;
  }

  auto maxVal = *std::max_element(v.begin(), v.end());
  DSU dsu(maxVal + 1);

  for (size_t i = 0; i < v.size(); i++)
  {
    dsu.unite(v[i], v[i]);
  }

  std::unordered_map<Id, Id> newIndex;
  Id next = 0;
  result.resize(v.size());

  for (size_t i = 0; i < v.size(); i++)
  {
    Id root = dsu.find(v[i]);
    if (!newIndex.count(root))
    {
      newIndex[root] = next++;
    }
    result[i] = newIndex[root];
  }
}

ListParams::ListParams()
  : AStringable()
{
}

void ListParams::updateDispatch()
{
  _dispatch.clear();
  _dispatchIndex.clear();
  VectorInt adresses(_params.size());
  for (size_t i = 0; i < _params.size(); ++i)
    adresses[i] = static_cast<Id>(_params[i].get().getAddress());

  reindex(adresses, _dispatch);
  makeDispatchIndexFromDispatch();
}

void ListParams::addParam(ParamInfo& param)
{
  if (param.isFixed()) return;
  _params.push_back(param);

  _dispatch.push_back(_params.size() - 1);
  _dispatchIndex.push_back(_params.size() - 1);
  param.setAddress(static_cast<Id>(_params.size()) - 1);
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

double ListParams::getValue(Id index) const
{
  if (index < 0 || index >= static_cast<Id>(_params.size()))
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

void ListParams::setValue(Id index, double value)
{
  if (index < 0 || index >= static_cast<Id>(_params.size()))
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
  result << toStrTitle(1, "List of Parameters:");
  for (Id ipar = 0, jpar = 0, npar = static_cast<Id>(_dispatchIndex.size()); ipar < npar; ipar++)
  {
    jpar++;
    result << jpar << " - " << _params[_dispatchIndex[ipar]].get().toString() << std::endl;
  }
  return result.str();
}

void ListParams::makeDispatchIndexFromDispatch()
{
  _dispatchIndex.clear();
  size_t nmax = 0;
  for (size_t i = 0; i < _dispatch.size(); ++i)
  {
    if (_dispatch[i] >= nmax)
    {
      _dispatchIndex.push_back(i);
      nmax = _dispatch[i] + 1;
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

std::vector<double> ListParams::getMinValues(double epsilon) const
{
  size_t nparam = _params.size();
  std::vector<double> values(nparam);
  for (size_t i = 0; i < nparam; ++i)
  {
    values[i] = _params[i].get().getUserMin() + epsilon;
  }
  return values;
}

std::vector<double> ListParams::getMaxValues(double epsilon) const
{
  size_t nparam = _params.size();
  std::vector<double> values(nparam);
  for (size_t i = 0; i < nparam; ++i)
  {
    values[i] = _params[i].get().getUserMax() - epsilon;
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
} // namespace gstlrn
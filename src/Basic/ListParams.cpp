#include "Basic/ListParams.hpp"
#include "geoslib_define.h"
#include <cstddef>
#include <sstream>

ListParams::ListParams()
: AStringable()
{
}

void ListParams::addParam(ParamInfo& param)
{
    _params.push_back(param);
}

String ListParams::toString(const AStringFormat* strfmt) const
{
    DECLARE_UNUSED(strfmt);
    std::stringstream result;
    result << "List of Parameters:\n";
    result << "---------------------\n";
    for (const auto& param : _params) {
        result << param.get().toString() + "\n";
    }
    return result.str();
}

std::vector<double> ListParams::getValues() const
{
    size_t nparam = _params.size();
    std::vector<double> values(nparam);
    for (size_t i = 0; i < nparam; ++i) {
        values[i] = _params[i].get().getValue();
    }
    return values;
}

void ListParams::setValues(const std::vector<double>& values)
{
    size_t size = values.size();
    for (size_t i = 0; i < size; i ++)
    {
        _params[i].get().setValue(values[i]);
    }
}
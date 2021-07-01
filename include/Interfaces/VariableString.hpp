#pragma once

#include <string>
#include <vector>
#include <iostream>

#include "Interfaces/AVariableTemplate.hpp"
#include "Interfaces/interface_d.hpp"

class VariableString: public AVariableTemplate<String>
{
  public:
    VariableString();
    VariableString(const String &name);
    VariableString(const VectorString& values);
    VariableString(const VariableString& ref);
    virtual ~VariableString();
    VariableString& operator= (const VariableString& ref);

    virtual VariableString* clone() const override;

    VectorDouble            getValues() const override;

#ifdef _USE_NETCDF
    netCDF::NcVar serialize(netCDF::NcFile& file, std::vector<netCDF::NcDim>& dims) const override;
    void deserialize(const netCDF::NcFile& file, const netCDF::NcVar& var) override;
#endif

    //bool                  isUndefined(int i) const override;

  private:
};

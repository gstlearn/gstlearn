#pragma once

#include "gstlearn_export.hpp"
#include "Interfaces/AVariableTemplate.hpp"
#include "Interfaces/interface_d.hpp"

#include <string>
#include <vector>
#include <iostream>

class GSTLEARN_EXPORT VariableString: public AVariableTemplate<String>
{
  public:
    VariableString();
    VariableString(const String &name);
    VariableString(const VectorString& values);
    VariableString(const VariableString& ref);
    virtual ~VariableString();
    VariableString& operator= (const VariableString& ref);

    /// Cloneable interface
    IMPLEMENT_CLONING(VariableString)

    VectorDouble            getValues() const override;

#ifdef _USE_NETCDF
    netCDF::NcVar serialize(netCDF::NcFile& file, std::vector<netCDF::NcDim>& dims) const override;
    void deserialize(const netCDF::NcFile& file, const netCDF::NcVar& var) override;
#endif

    //bool                  isUndefined(int i) const override;

  private:
};

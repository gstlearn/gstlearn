#pragma once

#include "DataBase/AVariableTemplate.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
  class GSTLEARN_EXPORT VariableString: public AVariableTemplate<String>
  {
  public:
    VariableString();
    VariableString(
      const String& name,
      size_t ncols = 1,
      const String& unit = "",
      bool naallowed = true,
      const VariableFormat& format = VariableFormat());
    VariableString(const VectorString& values);
    VariableString(const VariableString& ref);
    virtual ~VariableString();
    VariableString& operator=(const VariableString& ref);

    /// Cloneable interface
    IMPLEMENT_CLONING(VariableString)

    String getType() const override { return ("string"); }

    VectorDouble getValues(Id icol = 0) const override;

#ifdef _USE_NETCDF
    netCDF::NcVar serialize(
      netCDF::NcFile& file,
      std::vector<netCDF::NcDim>& dims) const override;
    void deserialize(const netCDF::NcFile& file, const netCDF::NcVar& var)
      override;
#endif

    bool isUndefined(size_t i) const override;

  private:
  };

} // namespace gstlrn

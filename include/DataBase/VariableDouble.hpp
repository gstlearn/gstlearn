/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "DataBase/AVariableTemplate.hpp"
#include "DataBase/VariableFormat.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
  class GSTLEARN_EXPORT VariableDouble: public AVariableTemplate<double>
  {
  public:
    VariableDouble();
    VariableDouble(
      const String& name,
      size_t ncols = 1,
      const String& unit = "",
      bool naallowed = true,
      const VariableFormat& format = VariableFormat());
    VariableDouble(const VectorDouble& values);
    VariableDouble(const VariableDouble& ref);
    virtual ~VariableDouble();
    VariableDouble& operator=(const VariableDouble& ref);

    /// Cloneable interface
    IMPLEMENT_CLONING(VariableDouble)

    String getType() const override { return ("double"); }

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

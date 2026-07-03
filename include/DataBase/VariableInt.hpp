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
#include "gstlearn_export.hpp"

namespace gstlrn
{

  class GSTLEARN_EXPORT VariableInt: public AVariableTemplate<Id>
  {
  public:
    VariableInt();
    VariableInt(
      const String& name,
      size_t ncols = 1,
      const String& unit = "",
      bool naallowed = true,
      const VariableFormat& format = VariableFormat());
    VariableInt(const VectorInt& values);
    VariableInt(const VariableInt& ref);
    ~VariableInt() override;
    VariableInt& operator=(const VariableInt& ref);

    /// Cloneable interface
    IMPLEMENT_CLONING(VariableInt)

    String getType() const override { return ("int"); }

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

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

  class GSTLEARN_EXPORT VariableBool: public AVariableTemplate<UChar>
  {
  public:
    VariableBool();
    VariableBool(const String& name, size_t ncols = 1);
    VariableBool(const VectorBool& values);
    VariableBool(const VariableBool& ref);
    virtual ~VariableBool();
    VariableBool& operator=(const VariableBool& ref);

    /// Cloneable interface
    IMPLEMENT_CLONING(VariableBool)

    String getType() const override { return ("bool"); }

    VectorDouble getValues(Id icol = 0) const override;

    bool getBool(size_t i) const;

#ifdef _USE_NETCDF
    netCDF::NcVar serialize(
      netCDF::NcFile& file,
      std::vector<netCDF::NcDim>& dims) const override;
    void deserialize(const netCDF::NcFile& file, const netCDF::NcVar& var)
      override;
#endif

    bool isUndefined(size_t i) const override;

    Id getNTrue(Id icol = 0) const;

    Id getNFalse(Id icol = 0) const { return (getNSamples() - getNTrue(icol)); }

  private:
  };

} // namespace gstlrn

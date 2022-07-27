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

#include "gstlearn_export.hpp"
#include "Interfaces/AVariableTemplate.hpp"

#include <string>
#include <vector>
#include <iostream>

class GSTLEARN_EXPORT VariableBool : public AVariableTemplate<bool>
{
  public:
    VariableBool();
    VariableBool(const String &name);
    VariableBool(const VectorBool& values);
    VariableBool(const VariableBool& ref);
    virtual ~VariableBool();
    VariableBool& operator= (const VariableBool& ref);

    /// Cloneable interface
    IMPLEMENT_CLONING(VariableBool)

    VectorDouble getValues() const override;

#ifdef _USE_NETCDF
    netCDF::NcVar serialize(netCDF::NcFile& file, std::vector<netCDF::NcDim>& dims) const override;
    void deserialize(const netCDF::NcFile& file,const netCDF::NcVar& var) override;
#endif

    //bool                    isUndefined(int i) const override;
  private:
};

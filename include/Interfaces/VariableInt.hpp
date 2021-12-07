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
#include "Interfaces/interface_d.hpp"

#include <string>
#include <vector>
#include <iostream>

class GSTLEARN_EXPORT VariableInt: public AVariableTemplate<int>
{
  public:
    VariableInt();
    VariableInt(const String &name);
    VariableInt(const VectorInt& values);
    VariableInt(const VariableInt& ref);
    virtual ~VariableInt();
    VariableInt& operator=(const VariableInt& ref);
    
    virtual VariableInt* clone() const override;
    VectorDouble            getValues() const override;
    //bool                  isUndefined(int i) const override;

#ifdef _USE_NETCDF
    netCDF::NcVar serialize(netCDF::NcFile& file, std::vector<netCDF::NcDim>& dims) const override;
    void deserialize(const netCDF::NcFile& file, const netCDF::NcVar &var) override;
#endif

  private:
};

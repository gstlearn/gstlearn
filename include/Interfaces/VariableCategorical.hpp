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
#include "Interfaces/Dictionary.hpp"

#include <string>

/// TODO : to be terminated
class GSTLEARN_EXPORT VariableCategorical: public AVariableTemplate<Category>
{
public:
  VariableCategorical(const String& name, const Dictionary& dico);
  VariableCategorical(const VariableCategorical &ref);
  virtual ~VariableCategorical();

  /// Cloneable interface
  IMPLEMENT_CLONING(VariableCategorical)

  virtual void resize(int n, const Category& val) override;
  VectorDouble getValues() const override;
  virtual void setValue(int i, const Category& cat) override;
  void setValue(int i, int value);
  void setValue(int i, const String& label);

#ifdef _USE_NETCDF
  netCDF::NcVar serialize(netCDF::NcFile& file, std::vector<netCDF::NcDim>& dims) const override;
  void deserialize(const netCDF::NcFile& file,const netCDF::NcVar& var) override;
#endif

private:
  Dictionary _dico;
};

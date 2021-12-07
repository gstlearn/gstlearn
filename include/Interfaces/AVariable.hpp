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
#include "Interfaces/interface_d.hpp"
#include "Basic/AStringable.hpp"

#include <string>
#include <vector>

#ifdef _USE_NETCDF
#include <netcdf>
#endif

class VariableDouble;
class VariableBool;
class VariableInt;
class VariableString;

class GSTLEARN_EXPORT AVariable : public AStringable
{
public:
  AVariable();
  AVariable(const String& name);
  AVariable(const AVariable& ref);
  virtual ~AVariable();
  AVariable& operator=(const AVariable& ref);

  virtual AVariable* clone() const = 0; 
  virtual unsigned int getNValues() const = 0;
  virtual VectorDouble getValues() const = 0;

#ifdef _USE_NETCDF
  virtual netCDF::NcVar serialize(netCDF::NcFile& file, std::vector<netCDF::NcDim>& dims) const = 0;
  virtual void deserialize(const netCDF::NcFile& file, const netCDF::NcVar &var) = 0;
#endif

  static AVariable* createVariable(const String& type, const String& name);
  void setName(const String& name);

  const String& getName() const;
protected:
  String _name;
};

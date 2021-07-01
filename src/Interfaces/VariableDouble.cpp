#include "Interfaces/VariableDouble.hpp"

VariableDouble::VariableDouble()
: AVariableTemplate()
{
}

VariableDouble::VariableDouble(const String &name)
: AVariableTemplate(name)
{
}

VariableDouble::VariableDouble(const VectorDouble& values)
:  AVariableTemplate("", values)
{
}

VariableDouble::VariableDouble(const VariableDouble& ref)
:  AVariableTemplate(ref)
{
}

VariableDouble& VariableDouble::operator=(const VariableDouble& ref)
{
  if (this != &ref)
  {
    _name = ref._name;
    _values = ref._values;
  }
  return(*this);
}

VariableDouble::~VariableDouble()
{
}

/**
 * Return a copy of the internal values vector
 */
VectorDouble  VariableDouble::getValues() const
{
  return (_values);
}

/**
 * Clone
 */
VariableDouble* VariableDouble::clone() const
{
  return (new VariableDouble(*this));
}

#ifdef _USE_NETCDF
/**
 *
 *  Create a NcVar object, add it to a netCDF file, create an attribute
 *  with the correct type(double) and return the NcVar object
 *
 * \param[in]  file   netCDF file to write in
 * \param[in]  dims   Vector of dimensions associate to Variable
 *
 */
netCDF::NcVar VariableDouble::serialize(netCDF::NcFile& file, std::vector<netCDF::NcDim>& dims) const
{
  netCDF::NcVar var_ncdf;
  var_ncdf = file.addVar(_name, netCDF::ncDouble, dims);
  var_ncdf.putVar(getValuesAsType().data());
  var_ncdf.putAtt("type", "double");
  return (var_ncdf);
}

/**
 *  TODO: First Param: netCDF::NcFile is unused, delete it in All Var.cpp et.hpp
 *  Read the content of a NcVar and set values accordingly
 *
 * @param[in]  file   netCDF file we read from
 *
 */
void VariableDouble::deserialize(const netCDF::NcFile& file, const netCDF::NcVar &var)
{
  std::vector<netCDF::NcDim> dims;
  int size = 1;

  dims= var.getDims();
  for (const auto& dim:dims)
    size *= dim.getSize();
  
  VectorDouble val(size);
  var.getVar(val.data());
  setValues(val);
}
#endif

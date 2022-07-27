#include "Interfaces/VariableInt.hpp"

VariableInt::VariableInt()
    : AVariableTemplate()
{
}

VariableInt::VariableInt(const String &name)
    : AVariableTemplate(name)
{
}

VariableInt::VariableInt(const VectorInt& values)
    : AVariableTemplate("", values)
{
}

VariableInt::VariableInt(const VariableInt& ref)
    : AVariableTemplate(ref)
{
}

VariableInt& VariableInt::operator=(const VariableInt& ref)
{
  if (this != &ref)
  {
    _name = ref._name;
    _values = ref._values;
  }
  return (*this);
}

VariableInt::~VariableInt()
{
}

/**
 * Return a copy of the values converted into double
 */
VectorDouble VariableInt::getValues() const
{
  VectorDouble res;
  for (VectorDouble::size_type i = 0; i < getNValues(); i++)
    res.push_back(static_cast<double>(_values.at(i)));
  return (res);
}

#ifdef _USE_NETCDF
/**
 *  Create a NcVar object, add it to a netCDF file, create an attribute
 *  with the correct type(int) and return the NcVar object
 *
 * @param[in]  file   netCDF file to write in
 * @param[in]  dims   Vector of dimensions associate to Variable
 */
netCDF::NcVar VariableInt::serialize(netCDF::NcFile& file, std::vector<netCDF::NcDim>& dims) const
{
  netCDF::NcVar var_ncdf;
  var_ncdf = file.addVar(_name, netCDF::ncInt, dims);
  var_ncdf.putVar(getValuesAsType().data());
  var_ncdf.putAtt("type", "int");
  return (var_ncdf);
}

/**
 * Read the content of a NcVar and set values accordingly
 *
 * @param[in]  file   netCDF file we read from
 * @param[in]  dims   Vector of dimensions associate to Variable
 */
void VariableInt::deserialize(const netCDF::NcFile& file,const netCDF::NcVar &var)
{
  std::vector<netCDF::NcDim> dims;
  int size = 1;

  dims= var.getDims();
  for (const auto& dim:dims)
    size *= dim.getSize();

  VectorInt val(size);
  var.getVar(val.data());
  setValues(val);
}
#endif


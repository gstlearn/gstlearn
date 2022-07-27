#include "Interfaces/VariableBool.hpp"

VariableBool::VariableBool()
    : AVariableTemplate()
{
}

VariableBool::VariableBool(const String &name)
    : AVariableTemplate(name)
{
}

VariableBool::VariableBool(const VectorBool& values)
    : AVariableTemplate("", values)
{
}

VariableBool::VariableBool(const VariableBool& ref)
    : AVariableTemplate(ref)
{
}

VariableBool& VariableBool::operator=(const VariableBool& ref)
{
  if (this != &ref)
  {
    _name = ref._name;
    _values = ref._values;
  }
  return (*this);
}

VariableBool::~VariableBool()
{
}

/**
 * Return a copy of the values converted into double
 */
VectorDouble VariableBool::getValues() const
{
  VectorDouble res;
  for (VectorDouble::size_type i = 0; i < getNValues(); i++)
    res.push_back(static_cast<double>(_values.at(i)));
  return (res);
}

#ifdef _USE_NETCDF
/**
 *
 *  Create a NcVar object, add it to a netCDF file, create an attribute
 *  with the correct type(bool) and return the NcVar object
 *
 * @param[in]  file   netCDF file to write in
 * @param[in]  dims   Vector of dimensions associate to Variable
 *
 * @remark    Netcdf file doesn't have a boolean type.
 *            Values are store as Byte (ncByte).
 */
netCDF::NcVar VariableBool::serialize(netCDF::NcFile& file, std::vector<netCDF::NcDim>& dims) const
{
  netCDF::NcVar var_ncdf;

  var_ncdf = file.addVar(_name, netCDF::ncByte, dims);
  var_ncdf.putVar(getValues().data());
  var_ncdf.putAtt("type", "bool");
  return (var_ncdf);
}

/*
 *  Read the content of a NcVar and set values accordingly
 *
 * @param[in]  file   netCDF file we read from
 * @param[in]  dims   Vector of dimensions associate to Variable
 *
 * @remark  Cannot read values directly from NcVar to VectorBool:
 *          Values are read as VectorDouble then  convert to VectorBool
 */
void VariableBool::deserialize(const netCDF::NcFile& file, const netCDF::NcVar& var)
{
  std::vector<netCDF::NcDim> dims;
  int size = 1;

  dims= var.getDims();
  for (const auto& dim:dims)
    size *= dim.getSize();

  VectorDouble val(size);
  var.getVar(val.data());
  VectorBool values(val.begin(), val.end());
  setValues(values);
}
#endif

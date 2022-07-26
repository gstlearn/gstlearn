#include "Interfaces/VariableString.hpp"

#include <sstream>

VariableString::VariableString()
    : AVariableTemplate()
{
}

VariableString::VariableString(const String &name)
    : AVariableTemplate(name)
{
}

VariableString::VariableString(const VectorString& values)
    : AVariableTemplate("", values)
{
}

VariableString::VariableString(const VariableString& ref)
    : AVariableTemplate(ref)
{
}

VariableString& VariableString::operator=(const VariableString& ref)
{
  if (this != &ref)
  {
    _name = ref._name;
    _values = ref._values;
  }
  return (*this);
}

VariableString::~VariableString()
{
}

/**
 * Return the values into a vector of double (provided that each string
 * can be converted to a double
 */
VectorDouble VariableString::getValues() const
{
  VectorDouble vals;
  for (const auto& val : _values)
  {
    double dval = UNDEF_DOUBLE;
    std::stringstream sstr(val);
    sstr >> dval;
    if (!sstr.good())
    {
      // No error (just ignore the string value)
      dval = UNDEF_DOUBLE;
    }
    vals.push_back(dval);
  }
  return vals;
}

#ifdef _USE_NETCDF
/**
 *  Create a NcVar object, add it to a netCDF file, create an attribute
 *  with the correct type(string) and return the NcVar object
 *
 * @param[in]  file   netCDF file to write in
 * @param[in]  dims   Vector of dimensions associate to Variable
 *
 * @remark     Values as to be convert to Vector<char*> before getting add to
 *             the NcVar
 */
netCDF::NcVar VariableString::serialize(netCDF::NcFile &file, std::vector<netCDF::NcDim>& dims) const
{
  // Convert VectorString to Vector<char*>
  std::vector<char*> cstring;
  VectorString v_string = getValuesAsType();
  for (auto& str : v_string)
  {
    cstring.push_back(&str.front());
  }

  netCDF::NcVar var_ncdf;
  var_ncdf = file.addVar(_name, netCDF::ncString, dims);
  var_ncdf.putVar(cstring.data());
  var_ncdf.putAtt("type", "string");
  return (var_ncdf);
}

/**
 * Read the content of a NcVar and set values accordingly
 *
 * @param[in]  file   netCDF file we read from
 * @param[in]  dims   Vector of dimensions associate to Variable
 *
 * @remark    Cannot read values from Netcdf as VectorString:
 *            Values are read as vector<char*> then converted to VectorString
 */
void VariableString::deserialize(const netCDF::NcFile& file, const netCDF::NcVar& var)
{
  std::vector<netCDF::NcDim> dims;
  int size = 1;

  dims = var.getDims();
  for (const auto& dim:dims)
    size *= dim.getSize();

  std::vector<char*> cstring(size);
  var.getVar(cstring.data());

  VectorString value(cstring.begin(), cstring.end());
  setValues(value);
}
#endif


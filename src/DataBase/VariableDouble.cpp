/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "DataBase/VariableDouble.hpp"
#include "Basic/Message.hpp"

using namespace gstlrn;

VariableDouble::VariableDouble()
  : AVariableTemplate()
{
}

VariableDouble::VariableDouble(
  const String& name,
  size_t ncols,
  const String& unit,
  bool naallowed,
  const VariableFormat& format)
  : AVariableTemplate(name, ncols, unit, naallowed, format)
{
}

VariableDouble::VariableDouble(const VectorDouble& values)
  : AVariableTemplate("", values)
{
}

VariableDouble::VariableDouble(const VariableDouble& ref)
  : AVariableTemplate(ref)
{
}

VariableDouble& VariableDouble::operator=(const VariableDouble& ref)
{
  if (this != &ref)
  {
    _name = ref._name;
    _values = ref._values;
  }
  return (*this);
}

VariableDouble::~VariableDouble() {}

VectorDouble VariableDouble::getValues(Id icol) const
{
  if (_getFlagVariableCheck()
      && !checkArg("VariableDouble::getValues", icol, getNCols()))
    return VectorDouble();
  const size_t nsamples = getNSamples();
  const size_t offset = icol * nsamples;

  VectorDouble res(nsamples);
  for (size_t i = 0; i < nsamples; i++)
    res[i] = static_cast<double>(_values.at(offset + i));

  return res;
}

bool VariableDouble::isUndefined(size_t i) const
{
  return (_values.at(i) == getNA<double>());
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
netCDF::NcVar VariableDouble::serialize(
  netCDF::NcFile& file,
  std::vector<netCDF::NcDim>& dims) const
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
void VariableDouble::deserialize(
  const netCDF::NcFile& file,
  const netCDF::NcVar& var)
{
  std::vector<netCDF::NcDim> dims;
  int size = 1;

  dims = var.getDims();
  for (const auto& dim: dims) size *= dim.getSize();

  VectorDouble val(size);
  var.getVar(val.data());
  setValues(val);
}
#endif

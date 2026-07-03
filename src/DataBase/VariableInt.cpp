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
#include "DataBase/VariableInt.hpp"
#include "Basic/Message.hpp"
#include "DataBase/AVariable.hpp"

using namespace gstlrn;

VariableInt::VariableInt()
  : AVariableTemplate()
{
}

VariableInt::VariableInt(
  const String& name,
  size_t ncols,
  const String& unit,
  bool naallowed,
  const VariableFormat& format)
  : AVariableTemplate(name, ncols, unit, naallowed, format)
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

VariableInt::~VariableInt() {}

VectorDouble VariableInt::getValues(Id icol) const
{
  if (_getFlagVariableCheck()
      && !checkArg("VariableInt::getValues", icol, getNCols()))
    return VectorDouble();
  const size_t nsamples = getNSamples();
  const size_t offset = icol * nsamples;

  VectorDouble res(nsamples);
  for (size_t i = 0; i < nsamples; i++)
    res[i] = static_cast<double>(_values.at(offset + i));

  return res;
}

bool VariableInt::isUndefined(size_t i) const
{
  return (_values.at(i) == getNA<Id>());
}

#ifdef _USE_NETCDF
/**
 *  Create a NcVar object, add it to a netCDF file, create an attribute
 *  with the correct type(int) and return the NcVar object
 *
 * @param[in]  file   netCDF file to write in
 * @param[in]  dims   Vector of dimensions associate to Variable
 */
netCDF::NcVar VariableInt::serialize(
  netCDF::NcFile& file,
  std::vector<netCDF::NcDim>& dims) const
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
void
  VariableInt::deserialize(const netCDF::NcFile& file, const netCDF::NcVar& var)
{
  std::vector<netCDF::NcDim> dims;
  int size = 1;

  dims = var.getDims();
  for (const auto& dim: dims) size *= dim.getSize();

  VectorInt val(size);
  var.getVar(val.data());
  setValues(val);
}
#endif

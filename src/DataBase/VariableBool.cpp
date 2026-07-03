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
#include "DataBase/VariableBool.hpp"
#include "Basic/Message.hpp"

namespace gstlrn
{
  VariableBool::VariableBool()
    : AVariableTemplate()
  {
  }

  VariableBool::VariableBool(const String& name, size_t ncols)
    : AVariableTemplate(name, ncols)
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

  VariableBool::~VariableBool() {}

  VectorDouble VariableBool::getValues(Id icol) const
  {
    if (_getFlagVariableCheck()
        && !checkArg("VariableBool::getValues", icol, getNCols()))
      return VectorDouble();
    const size_t nsamples = getNSamples();
    const size_t offset = icol * nsamples;

    VectorDouble res(nsamples);
    for (size_t i = 0; i < nsamples; i++)
      res[i] = static_cast<double>(_values.at(offset + i));

    return res;
  }

  bool VariableBool::isUndefined(size_t i) const
  {
    DECLARE_UNUSED(i);
    return false;
  }

  Id VariableBool::getNTrue(Id icol) const
  {
    if (_getFlagVariableCheck()
        && !checkArg("VariableBool::getNTrue", icol, getNCols()))
      return 0;
    Id ntrue = 0;
    const size_t nsamples = getNSamples();
    const size_t offset = icol * nsamples;
    for (size_t i = 0; i < nsamples; i++) ntrue += _values.at(offset + i);
    return (ntrue);
  }

  bool VariableBool::getBool(size_t i) const
  {
    return _values.at(i) != 0;
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
  netCDF::NcVar VariableBool::serialize(
    netCDF::NcFile& file,
    std::vector<netCDF::NcDim>& dims) const
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
  void VariableBool::deserialize(
    const netCDF::NcFile& file,
    const netCDF::NcVar& var)
  {
    std::vector<netCDF::NcDim> dims;
    int size = 1;

    dims = var.getDims();
    for (const auto& dim: dims) size *= dim.getSize();

    VectorDouble val(size);
    var.getVar(val.data());
    VectorBool values(val.begin(), val.end());
    setValues(values);
  }
#endif
} // namespace gstlrn

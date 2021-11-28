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
#include "Basic/HDF5format.hpp"

#include "geoslib_old_f.h"
#include <malloc.h>
#include <stdio.h>
#include <typeinfo>

#define DEBUG 0

HDF5format::HDF5format(const String& filename,
                       const String& varname)
  : _filename(filename)
  , _varname(varname)
{
}

HDF5format::HDF5format(const HDF5format &r)
  : _filename(r._filename)
  , _varname(r._varname)
{
}

HDF5format& HDF5format::operator= (const HDF5format &r)
{
  if (this != &r)
  {
    _filename = r._filename;
    _varname  = r._varname;
  }
  return *this;
}

HDF5format::~HDF5format()
{
}

/****************************************************************************/
/*!
**  Create a HDF5 file and wirte an integer array into it
**
** \param[in]  type      Data type (H5T_NATIVE_XXX)
** \param[in]  ndim      Space dimension
** \param[in]  dims      Array giving grid dimension in all space directions
**                       (Dimension: ndim)
** \param[in]  wdata     Array of values to be written
**                       (Dimension: product of dims)
**
** \remarks Any message has been suppressed as HDF5 provides error information
**
*****************************************************************************/
int HDF5format::createRegular(hid_t type,
                              int ndim,
                              hsize_t *dims,
                              void *wdata)
{
  try
  {
    hsize_t start0[ndim];
    for (int idim = 0; idim < ndim; idim++)
      start0[idim] = 0;

    H5File datafile(_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    DataSpace dataspace(ndim, dims);
    DataType datatype(H5::PredType::NATIVE_DOUBLE);
    DataSet dataset = datafile.createDataSet(_varname.c_str(),datatype, dataspace);
    dataspace.selectHyperslab(H5S_SELECT_SET, dims, start0);
    DataSpace memspace(ndim, dims, NULL);
    memspace.selectHyperslab(H5S_SELECT_SET, dims, start0);

    dataset.write((double *) wdata, PredType::IEEE_F32LE, memspace, dataspace);

    dataspace.close();
    datatype.close();
    dataset.close();
    datafile.close();
    memspace.close();

    return 0;
  }
  catch (FileIException& error)
  {
    return 1;
  }
  catch (GroupIException& error)
  {
    return 1;
  }
}

/****************************************************************************/
/*!
**  Read an array from an HDF5 file
**
** \return The returned array or NULL
**
** \param[in]  flag_compress 1 if the returned array must be compressed
** \param[in]  type      Data type (H5T_NATIVE_XXX)
** \param[in]  ndim      Space dimension
** \param[in]  start     Array of starting position (from 0)
**                       (Dimension: ndim)
** \param[in]  stride    Array of number of elements between blocks (optional)
**                       (Dimension: ndim)
** \param[in]  count     Array of number of elements per block
**                       (Dimension: ndim)
** \param[in]  block     Array of size of block (optional)
**                       (Dimension: ndim)
**
** \param[out] dimout    Array of final dimensions
**
** \remarks The returned array is dimensioned to the product of dimout
** \remarks Any message has been suppressed as HDF5 provides error information
**
*****************************************************************************/
void* HDF5format::readRegular(int flag_compress,
                              hid_t type,
                              int ndim,
                              hsize_t *start,
                              hsize_t *stride,
                              hsize_t *count,
                              hsize_t *block,
                              hsize_t *dimout)
{
  try
  {
    hsize_t start0[ndim];
    for (int idim = 0; idim < ndim; idim++)
      start0[idim] = 0;
    hsize_t dims[ndim];

    H5File datafile(_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT, H5P_DEFAULT);
    DataSet dataset = datafile.openDataSet(_varname.c_str());
    DataType datatype = dataset.getDataType();
    DataSpace dataspace = dataset.getSpace();

    dataspace.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    // Core allocation for returned array

    if (!flag_compress)
    {
      (void) dataspace.getSimpleExtentDims( dims, NULL);
      for (int idim = 0; idim < ndim; idim++)
        dimout[idim] = dims[idim];
    }
    else
    {
      for (int idim = 0; idim < ndim; idim++)
      {
        dimout[idim] = count[idim];
        if (block != NULL) dimout[idim] *= block[idim];
      }
    }

    // Core allocation

    void* rdata = allocArray(type, ndim, dimout);
    if (rdata == NULL) return (rdata);

    DataSpace memspace(ndim, dimout, NULL);
    if (flag_compress)
      memspace.selectHyperslab(H5S_SELECT_SET, dimout, start0);
    else
      memspace.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    dataset.read((double*) rdata, PredType::IEEE_F32LE, memspace, dataspace);

    dataspace.close();
    datatype.close();
    dataset.close();
    datafile.close();
    memspace.close();

    return rdata;
  }
  catch (FileIException& error)
  {
    return nullptr;
  }
  catch (GroupIException& error)
  {
    return nullptr;
  }
}

/****************************************************************************/
/*!
**  Write an array in an HDF5 file
**
** \return Error returned argument
**
** \param[in]  type      Data type (H5T_NATIVE_XXX)
** \param[in]  ndim      Space dimension
** \param[in]  dims      Array giving grid dimension in all space directions
**                       (Dimension: ndim)
** \param[in]  start     Array of starting position (from 0)
**                       (Dimension: ndim)
** \param[in]  stride    Array of number of elements between blocks
**                       (Dimension: ndim)
** \param[in]  count     Array of number of elements per block
**                       (Dimension: ndim)
** \param[in]  block     Array of size of block
**                       (Dimension: ndim)
** \param[in]  wdata     Array to be written
**                       (Dimension: product of dims)
**
** \remarks Any message has been suppressed as HDF5 provides error information
**
*****************************************************************************/
int HDF5format::writeRegular(hid_t type,
                             int ndim,
                             hsize_t *dims,
                             hsize_t *start,
                             hsize_t *stride,
                             hsize_t *count,
                             hsize_t *block,
                             void *wdata)
{
  try
  {
    Exception::dontPrint();
    hsize_t start0[ndim];
    for (int idim = 0; idim < ndim; idim++)
      start0[idim] = 0;

    H5File datafile(_filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    DataSet dataset = datafile.openDataSet(_varname.c_str());
    DataType datatype = dataset.getDataType();
    DataSpace dataspace = dataset.getSpace();
    dataspace.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
    DataSpace memspace(ndim, dims, NULL);
    memspace.selectHyperslab(H5S_SELECT_SET, count, start0, block, block);

    dataset.write((double *) wdata, PredType::IEEE_F32LE, memspace, dataspace);

    dataspace.close();
    datatype.close();
    dataset.close();
    datafile.close();
    memspace.close();

    return 0;
  }

  catch (FileIException& error)
  {
    return 1;
  }
  catch (GroupIException& error)
  {
    return 1;
  }
}

/****************************************************************************/
/*!
**  Allocate an array
**
** \return Pointer to the allocated array or NULL
**
** \param[in]  type      Data type (H5T_NATIVE_XXX)
** \param[in]  ndim      Space dimension
** \param[in]  dims      Array giving grid dimension in all space directions
**
*****************************************************************************/
void* HDF5format::allocArray(hid_t type, int ndim, hsize_t *dims)
{
  hsize_t size = H5Tget_size(type);
  int ntot = 1;
  for (int idim=0; idim<ndim; idim++) ntot *= dims[idim];
  void* data = (void *) calloc(ntot,size);
  return data;
}

int HDF5format::deleteFile()
{
  if (_filename.empty()) return 0;
  if( remove(_filename.c_str() ) != 0 ) return 1;
  return 0;
}

int HDF5format::getSize() const
{
  try
  {
    Exception::dontPrint();
    H5File file(H5std_string (_filename.c_str()), H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(_varname.c_str());
    DataSpace dataspace = dataset.getSpace();
    const int npts = dataspace.getSimpleExtentNpoints();
    return npts;
  }
  catch (FileIException& error)
  {
    int err = -1;
    return err;
  }
  catch (GroupIException& error)
  {
    int err = -1;
    return err;
  }
}

int HDF5format::getDataInt() const
{
  try
  {
    Exception::dontPrint();
    H5File file(H5std_string (_filename.c_str()), H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(_varname.c_str());
    DataType datatype = dataset.getDataType();
    DataSpace dataspace = dataset.getSpace();
    hsize_t* dims = _getDims(dataspace);
    H5T_class_t classt = datatype.getClass();
    if (classt != 0)
    {
      messerr("%s is not an int... you can't save this as an int.",
              _varname.c_str());
      return -1;
    }
    int *data = new int;
    H5T_order_t order;
    size_t size;
    bool big_endian;
    _getOrderSize(dataset, &order, &size, &big_endian);

    if (order == 0 && size == 1)
      dataset.read(data, PredType::STD_I8LE); // Our standard integer
    else if (order == 0 && size == 2)
      dataset.read(data, PredType::STD_I16LE); // Our standard integer
    else if (order == 0 && size == 4)
      dataset.read(data, PredType::STD_I32LE); // Our standard integer
    else if (order == 0 && size == 8)
      dataset.read(data, PredType::STD_I64LE);
    else if (order == 1 && size == 1)
      dataset.read(data, PredType::STD_I8BE); // Our standard integer
    else if (order == 1 && size == 2)
      dataset.read(data, PredType::STD_I16BE); // Our standard integer
    else if (order == 1 && size == 4)
      dataset.read(data, PredType::STD_I32BE);
    else if (order == 1 && size == 8)
      dataset.read(data, PredType::STD_I64BE);
    else
      messageAbort("Did not find data type");

    // Manage our memory properly
    int v = *data;

    delete[] dims;
    delete data;
    dataspace.close();
    datatype.close();
    dataset.close();
    file.close();
    return v;
  }
  catch (FileIException& error)
  {
    int err = -1;
    return err;
  }
  catch (GroupIException& error)
  {
    int err = -1;
    return err;
  }
}

float HDF5format::getDataFloat() const
{
  try
  {
    Exception::dontPrint();
    H5File file(H5std_string (_filename.c_str()), H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(_varname.c_str());
    DataType datatype = dataset.getDataType();
    DataSpace dataspace = dataset.getSpace();
    hsize_t* dims = _getDims(dataspace);
    H5T_class_t classt = datatype.getClass();
    if (classt != 1)
    {
      messerr("%s is not a float... you can't save this as a float.", _varname.c_str());
      return -1;
    }
    H5T_order_t order;
    size_t size;
    bool big_endian;
    _getOrderSize(dataset, &order, &size, &big_endian);
    float *data = new float;

    if (order == 0 && size == 4)
      dataset.read(data, PredType::IEEE_F32LE); // Our standard integer
    else if (order == 0 && size == 8)
    {
      dataset.read((float*) data, PredType::IEEE_F64LE);
      message("NOTE: This is actually double data. We are casting to float\n");
    }
    else if (order == 1 && size == 4)
      dataset.read(data, PredType::IEEE_F32BE);
    else if ((big_endian || order == 1) && size == 8)
    {
      message("NOTE: This is actually double data. We are casting to float\n");
      dataset.read((float*) data, PredType::IEEE_F64BE);
    }
    else
      messageAbort("Did not find data type");

    float v = *data;

    delete[] dims;
    delete data;
    dataspace.close();
    datatype.close();
    dataset.close();
    file.close();
    return v;
  }
  catch (FileIException& error)
  {
    float err = -1.;
    return err;
  }
  catch (GroupIException& error)
  {
    float err = -1.;
    return err;
  }
}

double HDF5format::getDataDouble() const
{
  try
  {
    Exception::dontPrint();
    H5File file(H5std_string (_filename.c_str()), H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(_varname.c_str());
    DataType datatype = dataset.getDataType();
    DataSpace dataspace = dataset.getSpace();
    hsize_t* dims = _getDims(dataspace);
    H5T_class_t classt = datatype.getClass();
    if (classt != 1)
    {
      messerr("%s is not a float... you can't save this as a float.",_varname.c_str());
      return -1.;
    }
    H5T_order_t order;
    size_t size;
    bool big_endian;
    _getOrderSize(dataset, &order, &size, &big_endian);
    double *data = new double;

    if (order == 0 && size == 4)
    {
      message("NOTE: This is actually float data. We are casting to double\n");
      dataset.read((double*) data, PredType::IEEE_F32LE); // Our standard integer
    }
    else if (order == 0 && size == 8)
      dataset.read(data, PredType::IEEE_F64LE);
    else if (order == 1 && size == 4)
    {
      message("NOTE: This is actually float data. We are casting to double\n");
      dataset.read((double*) data, PredType::IEEE_F32BE);
    }
    else if (order == 1 && size == 8)
      dataset.read((double*) data, PredType::IEEE_F64BE);
    else
      messageAbort("Did not find data type");

    float v = *data;

    delete[] dims;
    delete data;
    dataspace.close();
    datatype.close();
    dataset.close();
    file.close();
    return v;
  }
  catch (FileIException& error)
  {
    double err = -1.;
    return err;
  }
  catch (GroupIException& error)
  {
    double err = -1.;
    return err;
  }
}

VectorInt HDF5format::getDataVInt() const
{
  try
  {
    Exception::dontPrint();
    H5File file(H5std_string (_filename.c_str()), H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(_varname.c_str());
    DataType datatype = dataset.getDataType();
    DataSpace dataspace = dataset.getSpace();
    hsize_t* dims = _getDims(dataspace);
    const int npts = dataspace.getSimpleExtentNpoints(); // Gets length of data
    H5T_class_t classt = datatype.getClass(); // Gets the data type of the data
    // Let's make a quick error check
    if (classt != 0)
    {
      messerr(" %s is not an int... you can't save this as an int.", _varname.c_str());
      VectorInt err { 1, -1 };
      return err;
    }
    int *data = new int[npts]; // allocate at run time what the size will be
    IntType itype = dataset.getIntType();
    H5T_order_t order;
    size_t size;
    bool big_endian;
    _getOrderSize(dataset, &order, &size, &big_endian);

    if (order == 0 && size == 1)
      dataset.read(data, PredType::STD_I8LE); // Our standard integer
    else if (order == 0 && size == 2)
      dataset.read(data, PredType::STD_I16LE); // Our standard integer
    else if (order == 0 && size == 4)
      dataset.read(data, PredType::STD_I32LE); // Our standard integer
    else if (order == 0 && size == 8)
      dataset.read(data, PredType::STD_I64LE);
    else if (order == 1 && size == 1)
      dataset.read(data, PredType::STD_I8BE); // Our standard integer
    else if (order == 1 && size == 2)
      dataset.read(data, PredType::STD_I16BE); // Our standard integer
    else if (order == 1 && size == 4)
      dataset.read(data, PredType::STD_I32BE);
    else if (order == 1 && size == 8)
      dataset.read(data, PredType::STD_I64BE);
    else
      messageAbort("Did not find data type");

    VectorInt v(data, data + npts); // Arrays are nice, but vectors are better

    delete[] dims;
    delete[] data;
    dataspace.close();
    datatype.close();
    dataset.close();
    file.close();
    return v;
  }
  catch (FileIException& error)
  {
    VectorInt err { 1, -1 };
    return err;
  }
  catch (GroupIException& error)
  {
    VectorInt err { 1, -1 };
    return err;
  }
}

VectorFloat HDF5format::getDataVFloat() const
{
  try
  {
    Exception::dontPrint();
    H5File file(H5std_string (_filename.c_str()), H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(_varname.c_str());
    DataType datatype = dataset.getDataType();
    DataSpace dataspace = dataset.getSpace();
    hsize_t* dims = _getDims(dataspace);

    const int npts = dataspace.getSimpleExtentNpoints();
    H5T_class_t classt = datatype.getClass();
    if (classt != 1)
    {
      messerr("%s is not a float... you can't save this as a float.",
              _varname.c_str());
      VectorFloat err { 1, -1. };
      return err;
    }
    H5T_order_t order;
    size_t size;
    bool big_endian;
    _getOrderSize(dataset, &order, &size, &big_endian);
    float *data = new float[npts];

    if (order == 0 && size == 4)
      dataset.read(data, PredType::IEEE_F32LE); // Our standard integer
    else if (order == 0 && size == 8)
    {
      message("NOTE: This is actually double data. We are casting to float\n");
      dataset.read((float*) data, PredType::IEEE_F64LE);
    }
    else if (order == 1 && size == 4)
      dataset.read(data, PredType::IEEE_F32BE);
    else if ((big_endian || order == 1) && size == 8)
    {
      message("NOTE: This is actually double data We are casting to float\n");
      dataset.read((float*) data, PredType::IEEE_F64BE);
    }
    else
      messageAbort("Did not find data type");

    VectorFloat v(data, data + npts);

    delete[] dims;
    delete[] data;
    dataspace.close();
    datatype.close();
    dataset.close();
    file.close();
    return v;
  }
  catch (FileIException& error)
  {
    VectorFloat err { 1, -1. };
    return err;
  }
  catch (GroupIException& error)
  {
    VectorFloat err { 1, -1. };
    return err;
  }
}

VectorDouble HDF5format::getDataVDouble() const
{
  try
  {
    Exception::dontPrint();
    H5File file(H5std_string(_filename.c_str()), H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(_varname.c_str());
    DataType datatype = dataset.getDataType();
    DataSpace dataspace = dataset.getSpace();
    hsize_t* dims = _getDims(dataspace);

    const int npts = dataspace.getSimpleExtentNpoints();
    H5T_class_t classt = datatype.getClass();
    if (classt != 1)
    {
      messerr(" is not a float... you can't save this as a float.",
              _varname.c_str());
      VectorDouble err { 1, -1. };
      return err;
    }
    H5T_order_t order;
    size_t size;
    bool big_endian;
    _getOrderSize(dataset, &order, &size, &big_endian);
    double *data = new double[npts];

    if (order == 0 && size == 4)
    {
      message("NOTE: This is actually float data. We are casting to double\n");
      dataset.read((double*) data, PredType::IEEE_F32LE); // Our standard integer
    }
    else if (order == 0 && size == 8)
      dataset.read(data, PredType::IEEE_F64LE);
    else if (order == 1 && size == 4)
    {
      message("NOTE: This is actually float data. We are casting to double\n");
      dataset.read((double*) data, PredType::IEEE_F32BE);
    }
    else if (order == 1 && size == 8)
      dataset.read((double*) data, PredType::IEEE_F64BE);
    else
      messageAbort("Did not find data type\n");

    VectorDouble v(data, data + npts);

    delete[] dims;
    delete[] data;
    dataspace.close();
    datatype.close();
    dataset.close();
    file.close();
    return v;
  }
  catch (FileIException& error)
  {
    VectorDouble err { 1, -1. };
    return err;
  }
  catch (GroupIException& error)
  {
    VectorDouble err { 1, -1. };
    return err;
  }
}

/**
 * Reading VectorVectorInt
 * @return
 */
VectorVectorInt HDF5format::getDataVVInt() const
{
  try
  {
    Exception::dontPrint();
    H5File file(H5std_string (_filename.c_str()), H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(_varname.c_str());
    DataType datatype = dataset.getDataType();
    DataSpace dataspace = dataset.getSpace();
    hsize_t* dims = _getDims(dataspace);
    H5T_order_t order;
    size_t size;
    bool big_endian;
    _getOrderSize(dataset, &order, &size, &big_endian);
    size_t dim1 = dims[0];
    size_t dim2 = dims[1];
    auto data = new int[dim1 * dim2];
    auto md = new int*[dim1];
    for (size_t i = 0; i < dim1; ++i)
      md[i] = data + i * dim2;

    if (order == 0 && size == 1)
      dataset.read(data, PredType::STD_I8LE); // Our standard integer
    else if (order == 1 && size == 1)
      dataset.read(data, PredType::STD_I8BE);
    else if (order == 0 && size == 2)
      dataset.read(data, PredType::STD_I16LE);
    else if (order == 1 && size == 2)
      dataset.read(data, PredType::STD_I16BE);
    else if (order == 0 && size == 4)
      dataset.read(data, PredType::STD_I32LE);
    else if (order == 1 && size == 4)
      dataset.read(data, PredType::STD_I32BE);
    else if (order == 0 && size == 8)
      dataset.read(data, PredType::STD_I64LE);
    else if (order == 1 && size == 8)
      dataset.read(data, PredType::STD_I64BE);
    else
      messageAbort("Did not find data type");

    VectorVectorInt v(dim1, VectorInt(dim2, 0)); //data, data + npts);
    for (size_t i = 0; i < dim1; ++i)
      for (size_t j = 0; j < dim2; ++j)
        v[i][j] = md[i][j];

    delete[] dims;
    delete[] md;
    delete data;
    dataspace.close();
    datatype.close();
    dataset.close();
    file.close();
    return v;
  }
  catch (FileIException& error)
  {
    VectorVectorInt err { 1, VectorInt(1, -1) };
    return err;
  }
  catch (GroupIException& error)
  {
    VectorVectorInt err { 1, VectorInt(1, -1) };
    return err;
  }
}

VectorVectorFloat HDF5format::getDataVVFloat() const
{
  try
  {
    Exception::dontPrint();
    H5File file(H5std_string (_filename.c_str()), H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(_varname.c_str());
    DataType datatype = dataset.getDataType();
    DataSpace dataspace = dataset.getSpace();
    hsize_t* dims = _getDims(dataspace);
    H5T_order_t order;
    size_t size;
    bool big_endian;
    _getOrderSize(dataset, &order, &size, &big_endian);
    size_t dim1 = dims[0];
    size_t dim2 = dims[1];
    //float data[dim1][dim2];
    auto data = new float[dim1 * dim2];
    auto md = new float*[dim1];
    for (size_t i = 0; i < dim1; ++i)
      md[i] = data + i * dim2;

    if (order == 0 && size == 4)
      dataset.read(data, PredType::IEEE_F32LE); // Our standard integer
    else if (order == 0 && size == 8)
      dataset.read(data, PredType::IEEE_F64LE);
    else if (order == 1 && size == 4)
      dataset.read(data, PredType::IEEE_F32BE);
    else if (order == 1 && size == 8)
      dataset.read(data, PredType::IEEE_F64BE);
    else
      messageAbort("Did not find data type");

    VectorVectorFloat v(dim1, VectorFloat(dim2, 0)); //data, data + npts);
    for (size_t i = 0; i < dim1; ++i)
      for (size_t j = 0; j < dim2; ++j)
        v[i][j] = md[i][j];

    delete[] dims;
    delete[] md;
    delete data;
    dataspace.close();
    datatype.close();
    dataset.close();
    file.close();
    return v;
  }
  catch (FileIException& error)
  {
    VectorVectorFloat err { 1, VectorFloat(1, -1.) };
    return err;
  }
  catch (GroupIException& error)
  {
    VectorVectorFloat err { 1, VectorFloat(1, -1.) };
    return err;
  }
}

VectorVectorDouble HDF5format::getDataVVDouble() const
{
  try
  {
    Exception::dontPrint();
    H5File file(H5std_string (_filename.c_str()), H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(_varname.c_str());
    DataType datatype = dataset.getDataType();
    DataSpace dataspace = dataset.getSpace();
    hsize_t* dims = _getDims(dataspace);
    H5T_order_t order;
    size_t size;
    bool big_endian;
    _getOrderSize(dataset, &order, &size, &big_endian);
    size_t dim1 = dims[0];
    size_t dim2 = dims[1];
    auto data = new double[dim1 * dim2];
    auto md = new double*[dim1];
    for (size_t i = 0; i < dim1; ++i)
      md[i] = data + i * dim2;

    if (order == 0 && size == 4)
      dataset.read(data, PredType::IEEE_F32LE); // Our standard integer
    else if (order == 0 && size == 8)
      dataset.read(data, PredType::IEEE_F64LE);
    else if (order == 1 && size == 4)
      dataset.read(data, PredType::IEEE_F32BE);
    else if (order == 1 && size == 8)
      dataset.read(data, PredType::IEEE_F64BE);
    else
      messageAbort("Did not find data type\n");

    VectorVectorDouble v(dim1, VectorDouble(dim2, 0)); //data, data + npts);
    for (size_t i = 0; i < dim1; ++i)
      for (size_t j = 0; j < dim2; ++j)
        v[i][j] = md[i][j];

    delete[] dims;
    delete[] md;
    delete data;
    dataspace.close();
    datatype.close();
    dataset.close();
    file.close();
    return v;
  }
  catch (FileIException& error)
  {
    VectorVectorDouble err { 1, VectorDouble(1, -1.) };
    return err;
  }
  catch (GroupIException& error)
  {
    VectorVectorDouble err { 1, VectorDouble(1, -1.) };
    return err;
  }
}

/**
 * This function extracts one VectorDouble from a data set constructed
 * with a set of VectorDouble (i.e. VectorVectorDouble)
 * @param myrank Rank of the extracted VectorDouble
 * @return The extracted VectorDouble
 */
VectorDouble HDF5format::getDataDoublePartial(int myrank) const
{
  try
  {
    H5File file(H5std_string (_filename.c_str()), H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(_varname.c_str());
    DataType datatype = dataset.getDataType();
    DataSpace dataspace = dataset.getSpace();
    hsize_t* dims = _getDims(dataspace);
    H5T_order_t order;
    size_t size;
    bool big_endian;
    _getOrderSize(dataset, &order, &size, &big_endian);
    size_t dim1 = dims[0];
    size_t dim2 = dims[1];
    if (myrank < 0 || myrank >= (int) dim1)
    {
      messerr("The argument 'myrank' should be smaller than %d",(int) dim1);
      return VectorDouble();
    }

    // Read a single element (VectorDouble) from the dataset
    // - First define memory dataspace, then define hyperslab and read it into output array.

    hsize_t mydims[1];
    mydims[0] = dim2;
    DataSpace myspace( 1, mydims );
    hsize_t  count[2] = { 1, dim2 };
    hsize_t offset[2] = { (hsize_t) myrank, 0};
    auto data = new double[dim2];
    dataspace.selectHyperslab( H5S_SELECT_SET, count, offset );

    if (order == 0 && size == 4)
    {
      message("NOTE: This is actually float data. We are casting to double\n");
      dataset.read((double*) data, PredType::IEEE_F32LE, myspace, dataspace);
    }
    else if (order == 0 && size == 8)
      dataset.read(data, PredType::IEEE_F64LE, myspace, dataspace);
    else if (order == 1 && size == 4)
    {
      message("NOTE: This is actually float data. We are casting to double\n");
      dataset.read((double*) data, PredType::IEEE_F32BE, myspace, dataspace);
    }
    else if (order == 1 && size == 8)
      dataset.read((double*) data, PredType::IEEE_F64BE, myspace, dataspace);
    else
      messageAbort("Did not find data type\n");

    VectorDouble v(data, data + dim2);

    delete[] dims;
    delete data;
    dataspace.close();
    datatype.close();
    dataset.close();
    myspace.close();
    file.close();
    return v;
  }
  catch (FileIException& error)
  {
    VectorDouble err { VectorDouble(1, -1.) };
    return err;
  }
  catch (GroupIException& error)
  {
    VectorDouble err { VectorDouble(1, -1.) };
    return err;
  }
}

int HDF5format::writeDataDoublePartial(int myrank, const VectorDouble& data)
{
  Exception::dontPrint();
  uint npts = data.size(); // size of our data
  auto *a = new double[npts]; // convert to an array
  // convert std::vector to array. H5 does not seem to like the pointer implementation
  for (size_t i = 0; i < npts; ++i)
    a[i] = data[i];

  try
  {
    H5File file(H5std_string (_filename.c_str()), H5F_ACC_RDWR);
    DataSet dataset = file.openDataSet(_varname.c_str());
    DataType datatype = dataset.getDataType();
    DataSpace dataspace = dataset.getSpace();
    hsize_t* dims = _getDims(dataspace);
    H5T_order_t order;
    size_t size;
    bool big_endian;
    _getOrderSize(dataset, &order, &size, &big_endian);
    size_t dim1 = dims[0];
    size_t dim2 = dims[1];
    if (myrank < 0 || myrank >= (int) dim1)
    {
      messerr("The argument 'myrank' should be smaller than %d",(int) dim1);
      return 1;
    }
    if (npts != dim2)
    {
      messerr("Wrong dimension of argument 'data'. It should be %d long",dim2);
      return 1;
    }

    // Read a single element (VectorDouble) from the dataset
    // - First define memory dataspace, then define hyperslab and read it into output array.

    hsize_t mydims[1];
    mydims[0] = dim2;
    DataSpace myspace( 1, mydims );
    hsize_t  count[2] = { 1, dim2 };
    hsize_t offset[2] = { (hsize_t) myrank, 0};
    dataspace.selectHyperslab( H5S_SELECT_SET, count, offset );

    if (order == 0 && size == 4)
    {
      message("NOTE: This is actually float data. We are casting to double\n");
      dataset.write((double*) a, PredType::IEEE_F32LE, myspace, dataspace);
    }
    else if (order == 0 && size == 8)
      dataset.write(a, PredType::IEEE_F64LE, myspace, dataspace);
    else if (order == 1 && size == 4)
    {
      message("NOTE: This is actually float data. We are casting to double\n");
      dataset.write((double*) a, PredType::IEEE_F32BE, myspace, dataspace);
    }
    else if (order == 1 && size == 8)
      dataset.write((double*) a, PredType::IEEE_F64BE, myspace, dataspace);
    else
      messageAbort("Did not find data type\n");

    delete[] dims;
    dataspace.close();
    datatype.close();
    dataset.close();
    file.close();
    return 0;
  }
  catch (FileIException& error)
  {
    return 1;
  }
  catch (GroupIException& error)
  {
    return 1;
  }
}

/**
 * Initialize an empty HDF5 file. The user must give the vector of dimension
 * and the type of storage
 * @param argdims Vector of dimensions
 * @param type 1 for Int, 2 for Uint, 3 for Float, 4 for Double
 */
void HDF5format::createData(const VectorInt& argdims,int type)
{
  Exception::dontPrint();
  uint itr = 0;
  int ndim = argdims.size();
  hsize_t dims[ndim];
  for (int i = 0; i < ndim; i++) dims[i] = argdims[i];

  while (true)
  {
    try
    {
      H5File file(H5std_string (_filename.c_str()), H5F_ACC_RDWR);
      DataSpace dsp = DataSpace(ndim, dims);
      // int
      if (type == 1)
      {

        DataSet dset = file.createDataSet(H5std_string (_varname.c_str()), PredType::STD_I32LE, dsp);
        dset.close();
      }
      // uint
      else if (type == 2)
      {
        DataSet dset = file.createDataSet(H5std_string (_varname.c_str()), PredType::STD_U32LE, dsp);
        dset.close();
      }
      // float
      else if (type == 3)
      {
        DataSet dset = file.createDataSet(H5std_string (_varname.c_str()), PredType::IEEE_F32LE, dsp);
        dset.close();
      }
      // double
      else if (type == 4)
      {
        DataSet dset = file.createDataSet(H5std_string (_varname.c_str()), PredType::IEEE_F64LE, dsp);
        dset.close();
      }
      else
      {
        messerr("Unknown data type! EXITING");
      }

      // remember to close everything and delete our arrays
      dsp.close();
      file.close();
      break;
    }
    // Here we are catching if the file does not exist. We will then create a new file and return
    // back to the try statement
    catch (FileIException& error)
    {
      H5File file(H5std_string (_filename.c_str()), H5F_ACC_TRUNC);
      file.close();
      itr++;
      if (itr > 3)
      {
        messerr("We've tried too many times in the Int writing sequence");
        break;
      }
    }
  }
}

herr_t
file_info(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata)
{
    hid_t group;
    group = H5Gopen2(loc_id, name, H5P_DEFAULT);
    message("Data Set Name : %s\n",name);
    H5Gclose(group);
    return 0;
}

int HDF5format::displayNames() const
{
  try
  {
    H5File file(H5std_string(_filename.c_str()), H5F_ACC_RDWR);
    DataSet dataset = file.openDataSet(_varname.c_str());
    DataType datatype = dataset.getDataType();
    DataSpace dataspace = dataset.getSpace();

    mestitle(2,"Liste of Data Sets contained in the File");
    (void) H5Literate(file.getId(), H5_INDEX_NAME, H5_ITER_INC, NULL, file_info, NULL);
    message("\n");

    return 0;
  }
  catch (FileIException& error)
  {
    return 1;
  }
  catch (GroupIException& error)
  {
    return 1;
  }
}

hsize_t* HDF5format::_getDims(const DataSpace& dataspace) const
{
  int rank = dataspace.getSimpleExtentNdims();
  hsize_t* dims = new hsize_t[rank];
  dataspace.getSimpleExtentDims(dims);

  if (DEBUG)
  {
    message("Dataset rank = %d. dimensions = %d", rank, dims[0]);
    for (int idim = 1; idim < rank; idim++)
      message(" x %d", dims[idim]);
    message("\n");
  }
  return dims;
}

void HDF5format::_getOrderSize(const DataSet& dataset,
                               H5T_order_t* order,
                               size_t* size,
                               bool* big_endian) const
{
  H5std_string order_string;
  FloatType ftype = dataset.getFloatType();
  *order = ftype.getOrder(order_string);
  *size = ftype.getSize();
  *big_endian = (order_string == "Big endian byte order_stringing (1)");

  if (DEBUG)
  {
    message("Order = %d - Size = %d",(int) *order,(int) *size);
    if (*big_endian) message(" - Big Endian");
    message("\n");
  }
}

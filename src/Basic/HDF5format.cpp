/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/******************************************************************************/
#include "geoslib_old_f.h"
#include "Basic/HDF5format.hpp"

#include <malloc.h>
#include <stdio.h>
#include <typeinfo>

#define DEBUG 0

HDF5format::HDF5format(const String& filename,
                       const String& varname)
  : _filename(filename)
  , _varname(varname)
  , _datafile()
  , _dataset()
  , _datatype()
  , _dataspace()
{
}

HDF5format::HDF5format(const HDF5format &r)
  : _filename(r._filename)
  , _varname(r._varname)
  , _datafile(r._datafile)
  , _dataset(r._dataset)
  , _datatype(r._datatype)
  , _dataspace(r._dataspace)
{
}

HDF5format& HDF5format::operator= (const HDF5format &r)
{
  if (this != &r)
  {
    _filename = r._filename;
    _varname  = r._varname;
    _datafile = r._datafile;
    _dataset = r._dataset;
    _datatype = r._datatype;
    _dataspace = r._dataspace;
  }
  return *this;
}

HDF5format::~HDF5format()
{
}

/****************************************************************************/
/*!
**  Read an array from an HDF5 file
**
** \return The returned array or NULL
**
** \param[in]  flag_compress 1 if the returned array must be compressed
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
                              hsize_t *start,
                              hsize_t *stride,
                              hsize_t *count,
                              hsize_t *block,
                              hsize_t *dimout)
{
  try
  {
    int ndim = _getNDim();
    std::vector<hsize_t> start0(ndim);
    std::vector<hsize_t> dims(ndim);

    H5::DataType datatype = _dataset.getDataType();
    H5::DataSpace dataspace = _dataset.getSpace();
    dataspace.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    // Core allocation for returned array

    if (!flag_compress)
    {
      (void) _dataspace.getSimpleExtentDims(dims.data(), NULL);
      for (int idim = 0; idim < ndim; idim++)
        dimout[idim] = dims[idim];
    }
    else
    {
      for (int idim = 0; idim < ndim; idim++)
      {
        start0[idim] = 0;
        dimout[idim] = count[idim];
        if (block != NULL) dimout[idim] *= block[idim];
      }
    }

    void* rdata = allocArray(datatype, ndim, dimout);
    if (rdata == NULL) return (rdata);

    H5::DataSpace memspace(ndim, dimout, NULL);
    if (flag_compress)
      memspace.selectHyperslab(H5S_SELECT_SET, dimout, start0.data());
    else
      memspace.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

    _dataset.read((double*) rdata, datatype, memspace, dataspace);

    dataspace.close();
    memspace.close();

    return rdata;
  }
  catch (H5::Exception& error)
  {
    messerr("---> Problem in readRegular. Operation aborted");
    EXCEPTION_PRINT_ERROR(error);
    return nullptr;
  }
}

/****************************************************************************/
/*!
**  Write an array in an HDF5 file
**
** \return Error returned argument
**
** \param[in]  start     Array of starting position (from 0)
**                       (Dimension: ndim)
** \param[in]  stride    Array of number of elements between blocks
**                       (Dimension: ndim)
** \param[in]  count     Array of number of elements per block
**                       (Dimension: ndim)
** \param[in]  block     Array of size of block
**                       (Dimension: ndim)
**
** \param[in]  wdata     Array of compressed values to be written
**                       (Dimension: product of dimout)
**
** \remarks Any message has been suppressed as HDF5 provides error information
**
*****************************************************************************/
int HDF5format::writeRegular(hsize_t *start,
                             hsize_t *stride,
                             hsize_t *count,
                             hsize_t *block,
                             void *wdata)
{
  try
  {
    H5::Exception::dontPrint();
    int ndim = _getNDim();
    std::vector<hsize_t> start0(ndim);
    std::vector<hsize_t> dimin(ndim);
    for (int idim = 0; idim < ndim; idim++)
    {
      start0[idim] = 0;
      dimin[idim] = count[idim];
      if (block != NULL) dimin[idim] *= block[idim];
    }

    H5::DataType datatype = _dataset.getDataType();
    H5::DataSpace dataspace = _dataset.getSpace();
    dataspace.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
    H5::DataSpace memspace(ndim, dimin.data(), NULL);
    memspace.selectHyperslab(H5S_SELECT_SET, dimin.data(), start0.data());

    _dataset.write((double *) wdata, datatype, memspace, dataspace);

    dataspace.close();
    memspace.close();

    return 0;
  }
  catch (H5::Exception& error)
  {
    messerr("---> Problem in writeRegular. Operation aborted");
    EXCEPTION_PRINT_ERROR(error);
    return 1;
  }
}

/****************************************************************************/
/*!
**  Allocate an array
**
** \return Pointer to the allocated array or NULL
**
** \param[in]  type      Diension of element
** \param[in]  ndim      Space dimension
** \param[in]  dims      Array giving grid dimension in all space directions
**
*****************************************************************************/
void* HDF5format::allocArray(H5::DataType type, int ndim, hsize_t *dims)
{
  hsize_t size = type.getSize();
  int ntot = 1;
  for (int idim=0; idim<ndim; idim++) ntot *= dims[idim];
  void* data = (void *) calloc(ntot,(size_t) size);
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
    H5::Exception::dontPrint();
    H5::DataSpace dataspace = _dataset.getSpace();
    return dataspace.getSimpleExtentNpoints();
  }
  catch (H5::Exception& error)
  {
    EXCEPTION_PRINT_ERROR(error);
    return 0;
  }
}

int HDF5format::getDataInt() const
{
  try
  {
    H5::Exception::dontPrint();
    if (_checkClass(0)) return 0;

    int *data = new int;
    _readInt(data);
    int v = *data;

    delete data;
    return v;
  }
  catch (H5::Exception& error)
  {
    messerr("---> Problem in getDataInt. Operation aborted");
    EXCEPTION_PRINT_ERROR(error);
    return -1;
  }
}

float HDF5format::getDataFloat() const
{
  try
  {
    H5::Exception::dontPrint();
    if (_checkClass(1)) return 0.;

    float *data = new float;
    _readFloat(data);
    float v = *data;

    delete data;
    return v;
  }
  catch (H5::Exception& error)
  {
    messerr("---> Problem in getDataFloat. Operation aborted");
    EXCEPTION_PRINT_ERROR(error);
    return -1.;
  }
}

double HDF5format::getDataDouble() const
{
  try
  {
    H5::Exception::dontPrint();
    if (_checkClass(1)) return 0.;

    double *data = new double;
    _readDouble(data);
    double v = *data;

    delete data;
    return v;
  }
  catch (H5::Exception& error)
  {
    messerr("---> Problem in getDataDouble. Operation aborted");
    EXCEPTION_PRINT_ERROR(error);
    return -1.;
  }
}

VectorInt HDF5format::getDataVInt() const
{
  try
  {
    H5::Exception::dontPrint();
    if (_checkClass(0)) return VectorInt{1,-1};

    const int npts = _dataset.getSpace().getSimpleExtentNpoints();
    int *data = new int[npts];
    _readInt(data);
    VectorInt v(data, data + npts);

    delete[] data;
    return v;
  }
  catch (H5::Exception& error)
  {
    messerr("---> Problem in getDataVInt. Operation aborted");
    EXCEPTION_PRINT_ERROR(error);
    return { 1, -1 };
  }
}

VectorFloat HDF5format::getDataVFloat() const
{
  try
  {
    H5::Exception::dontPrint();
    if (_checkClass(1)) return VectorFloat{1,-1.};

    const int npts = _dataset.getSpace().getSimpleExtentNpoints();
    float *data = new float[npts];
    _readFloat(data);
    VectorFloat v(data, data + npts);

    delete[] data;
    return v;
  }
  catch (H5::Exception& error)
  {
    messerr("---> Problem in getDataVFloat. Operation aborted");
    EXCEPTION_PRINT_ERROR(error);
    return { 1, -1. };
  }
}

VectorDouble HDF5format::getDataVDouble() const
{
  try
  {
    H5::Exception::dontPrint();
    if (_checkClass(1)) return VectorDouble{1, -1.};

    const int npts = _dataset.getSpace().getSimpleExtentNpoints();
    double *data = new double[npts];
    _readDouble(data);
    VectorDouble v(data, data + npts);

    delete[] data;
    return v;
  }
  catch (H5::Exception& error)
  {
    messerr("---> Problem in getDataVDouble. Operation aborted");
    EXCEPTION_PRINT_ERROR(error);
    return { 1, -1. };
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
    H5::Exception::dontPrint();
    if (_checkClass(0)) return { 1, VectorInt(1, -1) };

    hsize_t* dims = _getDims();
    size_t dim1 = dims[0];
    size_t dim2 = dims[1];

    auto data = new int[dim1 * dim2];
    auto md = new int*[dim1];
    for (size_t i = 0; i < dim1; ++i)
      md[i] = data + i * dim2;
    _readInt(data);
    VectorVectorInt v(dim1, VectorInt(dim2, 0));
    for (size_t i = 0; i < dim1; ++i)
      for (size_t j = 0; j < dim2; ++j)
        v[i][j] = md[i][j];

    delete[] dims;
    delete[] md;
    delete[] data;
    return v;
  }
  catch (H5::Exception& error)
  {
    messerr("---> Problem in getDataVVInt. Operation aborted");
    EXCEPTION_PRINT_ERROR(error);
    return { 1, VectorInt(1, -1) };
  }
}

VectorVectorFloat HDF5format::getDataVVFloat() const
{
  try
  {
    H5::Exception::dontPrint();
    if (_checkClass(1)) return { 1, VectorFloat(1, -1.) };

    hsize_t* dims = _getDims();
    size_t dim1 = dims[0];
    size_t dim2 = dims[1];

    auto data = new float[dim1 * dim2];
    auto md = new float*[dim1];
    for (size_t i = 0; i < dim1; ++i)
      md[i] = data + i * dim2;
    _readFloat(data);
    VectorVectorFloat v(dim1, VectorFloat(dim2, 0));
    for (size_t i = 0; i < dim1; ++i)
      for (size_t j = 0; j < dim2; ++j)
        v[i][j] = md[i][j];

    delete[] dims;
    delete[] md;
    delete data;
    return v;
  }
  catch (H5::Exception& error)
  {
    messerr("---> Problem in getDataVVFloat. Operation aborted");
    EXCEPTION_PRINT_ERROR(error);
    return { 1, VectorFloat(1, -1.) };
  }
}

VectorVectorDouble HDF5format::getDataVVDouble() const
{
  try
  {
    H5::Exception::dontPrint();
    if (_checkClass(1)) return { 1, VectorDouble(1, -1.) };

    hsize_t* dims = _getDims();
    size_t dim1 = dims[0];
    size_t dim2 = dims[1];

    auto data = new double[dim1 * dim2];
    auto md = new double*[dim1];
    for (size_t i = 0; i < dim1; ++i)
      md[i] = data + i * dim2;
    _readDouble(data);
    VectorVectorDouble v(dim1, VectorDouble(dim2, 0)); //data, data + npts);
    for (size_t i = 0; i < dim1; ++i)
      for (size_t j = 0; j < dim2; ++j)
        v[i][j] = md[i][j];

    delete[] dims;
    delete[] md;
    delete[] data;
    return v;
  }
  catch (H5::Exception& error)
  {
    messerr("---> Problem in getDataVVDouble. Operation aborted");
    EXCEPTION_PRINT_ERROR(error);
    return { 1, VectorDouble(1, -1.) };
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
    H5::DataSpace dataspace = _dataset.getSpace();

    hsize_t* dims = _getDims();
    H5T_order_t order;
    size_t size;
    bool big_endian;
    _getOrderSize(&order, &size, &big_endian);
    size_t dim1 = dims[0];
    size_t dim2 = dims[1];
    if (myrank < 0 || myrank >= (int) dim1)
    {
      messerr("The argument 'myrank' should be smaller than %d",(int) dim1);
      return VectorDouble();
    }

    hsize_t mydims[1];
    mydims[0] = dim2;
    H5::DataSpace memspace( 1, mydims );
    hsize_t  count[2] = { 1, dim2 };
    hsize_t offset[2] = { (hsize_t) myrank, 0};
    dataspace.selectHyperslab( H5S_SELECT_SET, count, offset );

    auto data = new double[dim2];
    _readDouble(data, memspace, dataspace);
    VectorDouble v(data, data + dim2);

    delete[] dims;
    delete[] data;
    dataspace.close();
    memspace.close();
    return v;
  }
  catch (H5::Exception& error)
  {
    messerr("---> Problem in getDataDoublePartial. Operation aborted");
    EXCEPTION_PRINT_ERROR(error);
    return { VectorDouble(1, -1.) };
  }
}

int HDF5format::writeDataDoublePartial(int myrank, const VectorDouble& data)
{
  H5::Exception::dontPrint();
  size_t npts = data.size(); // size of our data
  auto *a = new double[npts]; // convert to an array
  // convert std::vector to array. H5 does not seem to like the pointer implementation
  for (size_t i = 0; i < npts; ++i)
    a[i] = data[i];

  try
  {
    H5::DataSpace dataspace = _dataset.getSpace();
    hsize_t* dims = _getDims();
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
    H5::DataSpace memspace( 1, mydims );
    hsize_t  count[2] = { 1, dim2 };
    hsize_t offset[2] = { (hsize_t) myrank, 0};
    dataspace.selectHyperslab( H5S_SELECT_SET, count, offset );

    _writeDouble(a, memspace, dataspace);

    delete[] dims;
    delete[] a;
    memspace.close();
    dataspace.close();
    return 0;
  }
  catch (H5::Exception& error)
  {
    messerr("---> Problem in writeDataDoublePartial. Operation aborted");
    EXCEPTION_PRINT_ERROR(error);
    return 1;
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
    mestitle(2,"List of Data Sets contained in the File");
    (void) H5Literate(_datafile.getId(), H5_INDEX_NAME, H5_ITER_INC, NULL, file_info, NULL);
    message("\n");

    return 0;
  }
  catch (H5::Exception& error)
  {
    messerr("---> Problem in displayNames. Operation aborted");
    EXCEPTION_PRINT_ERROR(error);
    return 1;
  }
}

int HDF5format::_getNDim() const
{
  H5::DataSpace dataspace = _dataset.getSpace();
  int rank = dataspace.getSimpleExtentNdims();
  return rank;
}

hsize_t* HDF5format::_getDims() const
{
  H5::DataSpace dataspace = _dataset.getSpace();
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

void HDF5format::_getOrderSize(H5T_order_t* order,
                               size_t* size,
                               bool* big_endian) const
{
  H5std_string order_string;
  H5::FloatType ftype = _dataset.getFloatType();
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

void HDF5format::openFile(const String& filename) const
{
  // Define the File name
  if (! filename.empty()) _filename = filename;

  // Perform the Open operation
  try
  {
    _datafile = H5::H5File(H5std_string (_filename.c_str()), H5F_ACC_RDONLY);
  }
  catch (H5::Exception& error)
  {
    messerr("---> Problem in openFile. Operation aborted");
    EXCEPTION_PRINT_ERROR(error);
    return;
  }
}

void HDF5format::openNewFile(const String& filename) const
{
  // Define the File name
  if (! filename.empty()) _filename = filename;

  // Perform the Open operation
  try
  {
    _datafile = H5::H5File(_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT);

  }
  catch (H5::Exception& error)
  {
    messerr("---> Problem in openNewFile. Operation aborted");
    EXCEPTION_PRINT_ERROR(error);
    return;
  }
}

void HDF5format::openDataSet(const String& varname) const
{
  // Define the Data Set Name
  if (! varname.empty()) _varname = varname;

  // Perform the Open operation
  try
  {
    _dataset = _datafile.openDataSet(_varname.c_str());
    _datatype = _dataset.getDataType();
    _dataspace = _dataset.getSpace();
  }
  catch (H5::Exception& error)
  {
    messerr("---> Problem in openDataSet. Operation aborted");
    EXCEPTION_PRINT_ERROR(error);
    return;
  }
}

void HDF5format::openNewDataSet(const String& varname,
                                int ndim,
                                hsize_t *dims,
                                H5::DataType type) const
{
  // Define the Data Set Name
  if (! varname.empty()) _varname = varname;

  // Perform the Open operation
  try
  {
    _dataspace = H5::DataSpace(ndim, dims);
    _datatype = type;
    _dataset = _datafile.createDataSet(H5std_string(_varname.c_str()), _datatype, _dataspace);
  }
  catch (H5::Exception& error)
  {
    messerr("---> Problem in openNewDataSet. Operation aborted");
    EXCEPTION_PRINT_ERROR(error);
    return;
  }
}

void HDF5format::closeFile() const
{
  _datafile.close();
}

void HDF5format::closeDataSet() const
{
  _dataspace.close();
  _dataset.close();
  _datatype.close();
}

void HDF5format::_readInt(int *data,
                          const H5::DataSpace& memspace,
                          const H5::DataSpace& dataspace) const
{
  H5T_order_t order;
  size_t size;
  bool big_endian;
  _getOrderSize(&order, &size, &big_endian);

  if (order == 0 && size == 1)
    _dataset.read(data, H5::PredType::STD_I8LE, memspace, dataspace);
  else if (order == 0 && size == 2)
    _dataset.read(data, H5::PredType::STD_I16LE, memspace, dataspace);
  else if (order == 0 && size == 4)
    _dataset.read(data, H5::PredType::STD_I32LE, memspace, dataspace);
  else if (order == 0 && size == 8)
    _dataset.read(data, H5::PredType::STD_I64LE, memspace, dataspace);
  else if (order == 1 && size == 1)
    _dataset.read(data, H5::PredType::STD_I8BE, memspace, dataspace);
  else if (order == 1 && size == 2)
    _dataset.read(data, H5::PredType::STD_I16BE, memspace, dataspace);
  else if (order == 1 && size == 4)
    _dataset.read(data, H5::PredType::STD_I32BE, memspace, dataspace);
  else if (order == 1 && size == 8)
    _dataset.read(data, H5::PredType::STD_I64BE, memspace, dataspace);
  else
    messageAbort("Did not find data type");
}

void HDF5format::_readFloat(float *data,
                            const H5::DataSpace& memspace,
                            const H5::DataSpace& dataspace) const
{
  H5T_order_t order;
  size_t size;
  bool big_endian;
  _getOrderSize(&order, &size, &big_endian);

  if (order == 0 && size == 4)
    _dataset.read((float*) data, H5::PredType::IEEE_F32LE, memspace, dataspace);
  else if (order == 0 && size == 8)
    _dataset.read((float*) data, H5::PredType::IEEE_F64LE, memspace, dataspace);
  else if (order == 1 && size == 4)
    _dataset.read((float*) data, H5::PredType::IEEE_F32BE, memspace, dataspace);
  else if ((big_endian || order == 1) && size == 8)
    _dataset.read((float*) data, H5::PredType::IEEE_F64BE, memspace, dataspace);
  else
    messageAbort("Did not find data type");
}

void HDF5format::_readDouble(double *data,
                             const H5::DataSpace& memspace,
                             const H5::DataSpace& dataspace) const
{
  H5T_order_t order;
  size_t size;
  bool big_endian;
  _getOrderSize(&order, &size, &big_endian);

  if (order == 0 && size == 4)
    _dataset.read((double*) data, H5::PredType::IEEE_F32LE, memspace, dataspace);
  else if (order == 0 && size == 8)
    _dataset.read((double*) data, H5::PredType::IEEE_F64LE, memspace, dataspace);
  else if (order == 1 && size == 4)
    _dataset.read((double*) data, H5::PredType::IEEE_F32BE, memspace, dataspace);
  else if (order == 1 && size == 8)
    _dataset.read((double*) data, H5::PredType::IEEE_F64BE, memspace, dataspace);
  else
    messageAbort("Did not find data type");
}

int HDF5format::_checkClass(int value) const
{
  H5T_class_t classt = _datatype.getClass();
  if (classt != value)
  {
    messerr("%s is not a float... you can't save this in current format.",_varname.c_str());
    return 1;
  }
  return 0;
}

void HDF5format::_writeDouble(double *data,
                              const H5::DataSpace& memspace,
                              const H5::DataSpace& dataspace) const
{
  H5T_order_t order;
  size_t size;
  bool big_endian;
  _getOrderSize(&order, &size, &big_endian);

  if (order == 0 && size == 4)
    _dataset.write((double*) data, H5::PredType::IEEE_F32LE, memspace, dataspace);
  else if (order == 0 && size == 8)
    _dataset.write((double*) data, H5::PredType::IEEE_F64LE, memspace, dataspace);
  else if (order == 1 && size == 4)
    _dataset.write((double*) data, H5::PredType::IEEE_F32BE, memspace, dataspace);
  else if (order == 1 && size == 8)
    _dataset.write((double*) data, H5::PredType::IEEE_F64BE, memspace, dataspace);
  else
    messageAbort("Did not find data type\n");
}

void HDF5format::_writeAll(const char* type, void *a)
{
  try
  {
    if (type == (char*) typeid(int).name())
      _dataset.write(a, H5::PredType::STD_I32LE);
    else if (type == (char*) typeid(float).name())
      _dataset.write(a, H5::PredType::IEEE_F32LE);
    else if (type == (char*) typeid(double).name())
      _dataset.write(a, H5::PredType::IEEE_F64LE);
    else
      messerr("Unknown data type! EXITING");

    return;
  }
  catch (H5::Exception& error)
  {
    messerr("---> Problem in writeAll. Operation aborted");
    error.printError();
    return;
  }
}

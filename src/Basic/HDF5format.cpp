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

#include "H5Cpp.h"

using namespace H5;

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
  hid_t    datafile, dataspace, memspace, dataset, err;
  hsize_t *start0;
  int      error;

  // Initializations

  error    = 1;
  datafile = dataspace = memspace = dataset = -1;

  // Local core allocation

  start0 = (hsize_t *) mem_alloc(sizeof(hsize_t) * ndim,0);
  if (start0 == NULL) goto label_end;
  for (int idim=0; idim<ndim; idim++) start0[idim] = 0;

  // Open the Datafile

  datafile = H5Fcreate(_filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  if (datafile < 0) goto label_end;

  dataspace = H5Screate_simple(ndim, dims, NULL);
  if (dataspace < 0) goto label_end;

  err = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,start0,NULL,dims,NULL);
  if (err < 0) goto label_end;

  memspace = H5Screate_simple(ndim, dims, NULL);
  if (memspace < 0) goto label_end;

  err = H5Sselect_hyperslab(memspace,H5S_SELECT_SET,start0,NULL,dims,NULL);
  if (err < 0) goto label_end;

  dataset = H5Dcreate(datafile, _varname.c_str(), type,
                      dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset < 0) goto label_end;

  err = H5Dwrite(dataset, type, memspace, dataspace, H5P_DEFAULT, wdata);
  if (err < 0) goto label_end;

  // Set the error return code

  error = 0;

label_end:
  start0 = (hsize_t *) mem_free((char *) start0);
  if (memspace  >= 0) (void) H5Sclose(memspace);
  if (dataspace >= 0) (void) H5Sclose(dataspace);
  if (dataset   >= 0) (void) H5Dclose(dataset);
  if (datafile  >= 0) (void) H5Fclose(datafile);
  return(error);
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
  hid_t    datafile, dataset, dataspace, memspace, err;
  hsize_t *start0,*dims;
  int      error;
  void    *rdata;

  // Initializations

  error    = 1;
  datafile = dataset = dataspace = memspace = -1;
  rdata    = NULL;
  start0   = dims = (hsize_t *) NULL;

  // Local core allocation

  start0 = (hsize_t *) mem_alloc(sizeof(hsize_t) * ndim,0);
  if (start0 == NULL) goto label_end;
  dims   = (hsize_t *) mem_alloc(sizeof(hsize_t) * ndim,0);
  if (dims   == NULL) goto label_end;
  for (int idim=0; idim<ndim; idim++) start0[idim] = 0;

  // Open the Datafile

  datafile = H5Fopen(_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (datafile < 0) goto label_end;

  dataset = H5Dopen(datafile, _varname.c_str(), H5P_DEFAULT);
  if (dataset < 0) goto label_end;

  dataspace = H5Dget_space(dataset);
  if (dataspace < 0) goto label_end;

  // Core allocation for returned array

  if (! flag_compress)
  {
    err = H5Sget_simple_extent_dims(dataspace, dims, NULL);
    if (err < 0) goto label_end;
    for (int idim=0; idim<ndim; idim++) dimout[idim] = dims[idim];
  }
  else
  {
    for (int idim=0; idim<ndim; idim++)
    {
      dimout[idim] = count[idim];
      if (block != NULL) dimout[idim] *= block[idim];
    }
  }

  // Core allocation

  rdata = allocArray(type,ndim,dimout);
  if (rdata == NULL) return(rdata);

  err = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,start,stride,count,block);
  if (err < 0) goto label_end;

  memspace = H5Screate_simple(ndim,dimout,NULL);
  if (memspace < 0) goto label_end;

  // In the compression case, update the arguments

  if (flag_compress)
    err = H5Sselect_hyperslab(memspace,H5S_SELECT_SET,
                              start0,block,count,block);
  else
    err = H5Sselect_hyperslab(memspace,H5S_SELECT_SET,
                              start,stride,count,block);
  if (err < 0) goto label_end;

  err = H5Dread(dataset, type, memspace, dataspace,H5P_DEFAULT, rdata);
  if (err < 0) goto label_end;

  // Set the error return code

  error = 0;

label_end:
  start0 = (hsize_t *) mem_free((char *) start0);
  dims   = (hsize_t *) mem_free((char *) dims);
  if (dataset   >= 0) (void) H5Dclose(dataset);
  if (dataspace >= 0) (void) H5Sclose(dataspace);
  if (memspace  >= 0) (void) H5Sclose(memspace);
  if (datafile  >= 0) (void) H5Fclose(datafile);
  if (error)
  {
    rdata = (void *) mem_free((char *) rdata);
    rdata = NULL;
  }
  return rdata;
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
  hid_t    datafile, dataset, dataspace, memspace, err;
  hsize_t *start0;
  int      error;

  // Initializations

  error    = 1;
  datafile = dataset = dataspace = memspace = -1;

  // Local core allocation

  start0 = (hsize_t *) mem_alloc(sizeof(hsize_t) * ndim,0);
  if (start0 == NULL) goto label_end;
  for (int idim=0; idim<ndim; idim++) start0[idim] = 0;

  // Open the Datafile

  datafile = H5Fopen(_filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (datafile < 0) goto label_end;

  dataset = H5Dopen(datafile, _varname.c_str(), H5P_DEFAULT);
  if (dataset < 0) goto label_end;

  dataspace = H5Dget_space(dataset);
  if (dataspace < 0) goto label_end;

  err = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,start,stride,count,block);
  if (err < 0) goto label_end;

  memspace = H5Screate_simple(ndim,dims,NULL);
  if (memspace < 0) goto label_end;

  // Define the hyperslab in memory

  err = H5Sselect_hyperslab(memspace,H5S_SELECT_SET,start0,block,count,block);
  if (err < 0) goto label_end;

  err = H5Dwrite(dataset, type, memspace, dataspace ,H5P_DEFAULT, wdata);
  if (err < 0) goto label_end;

  // Set the error return code

  error = 0;

label_end:
  start0 = (hsize_t *) mem_free((char *) start0);
  if (dataset   >= 0) (void) H5Dclose(dataset);
  if (dataspace >= 0) (void) H5Sclose(dataspace);
  if (memspace  >= 0) (void) H5Sclose(memspace);
  if (datafile  >= 0) (void) H5Fclose(datafile);
  return error;
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

int HDF5format::deleteFile(const String& filename)
{
  if( remove( filename.c_str() ) != 0 ) return 1;
  return 0;
}

/**
 * Numeric implementation of our write data function
 * Only accepts numerical values. Integers, floats, or doubles
 * @param data
 */
template<typename T>
void HDF5format::writeData(const T &data)
{
  uint itr = 0;
  auto *a = new T { data };
  char* type = (char*) (typeid(T).name());
  int vrank = 0;
  hsize_t dims[1];
  dims[0] = 1;
  while (true)
  {
    try
    {
      H5File file(_filename.c_str(), H5F_ACC_RDWR);
      DataSpace dsp = DataSpace(vrank, dims);
      // int
      if (type == (char*) typeid(int).name())
      {
        DataSet dset = file.createDataSet(_varname.c_str(), PredType::STD_I32LE, dsp);
        dset.write(a, PredType::STD_I32LE);
        dset.close();
      }
      // float
      else if (type == (char*) typeid(float).name())
      {
        DataSet dset = file.createDataSet(_varname.c_str(), PredType::IEEE_F32LE, dsp);
        dset.write(a, PredType::IEEE_F32LE);
        dset.close();
      }
      else if (type == (char*) typeid(double).name())
      {
        DataSet dset = file.createDataSet(_varname.c_str(), PredType::IEEE_F64LE, dsp);
        dset.write(a, PredType::IEEE_F64LE);
        dset.close();
      }
      else
      {
        messerr("Unknown data type! EXITING");
      }
      dsp.close();
      file.close();
      delete a;
      break;
    }
    catch (FileIException error)
    {
      H5File file(_filename.c_str(), H5F_ACC_TRUNC);
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

template<typename T>
void HDF5format::writeData(const std::vector<T> &data)
{
  uint itr = 0; // Used to ensure we don't get stuck in an infinite loop
  uint npts = data.size(); // size of our data
  auto *a = new T[npts]; // convert to an array
  char* type = (char*) (typeid(a[0]).name());
  int vrank = 1; // since we are using std::vectors we are storing everything in one dimension

  // convert std::vector to array. H5 does not seem to like the pointer implementation
  for (size_t i = 0; i < npts; ++i)
    a[i] = data[i];
  // conventional syntax for H5 data writing
  hsize_t dims[1];
  dims[0] = npts;
  // Let's make sure we are doing what we want and output it to the std output

  // loop here will check if the file exists.
  while (true)
  {
    // This assumes that the file already exists and will then write to the file
    try
    {
      H5File file(_filename.c_str(), H5F_ACC_RDWR);
      DataSpace dsp = DataSpace(vrank, dims);
      // int
      if (type == (char*) typeid(int).name())
      {
        DataSet dset = file.createDataSet(_varname.c_str(), PredType::STD_I32LE, dsp);
        dset.write(a, PredType::STD_I32LE);
        dset.close();
      }
      // uint
      else if (type == (char*) typeid(uint).name())
      {
        DataSet dset = file.createDataSet(_varname.c_str(), PredType::STD_U32LE, dsp);
        dset.write(a, PredType::STD_U32LE);
        dset.close();
      }
      // float
      else if (type == (char*) typeid(float).name())
      {
        DataSet dset = file.createDataSet(_varname.c_str(), PredType::IEEE_F32LE, dsp);
        dset.write(a, PredType::IEEE_F32LE);
        dset.close();
      }
      // double
      else if (type == (char*) typeid(double).name())
      {
        DataSet dset = file.createDataSet(_varname.c_str(), PredType::IEEE_F64LE, dsp);
        dset.write(a, PredType::IEEE_F64LE);
        dset.close();
      }
      else
      {
        messerr("Unknown data type! EXITING");
      }

      // remember to close everything and delete our arrays
      dsp.close();
      file.close();
      delete[] a;
      break;
    }
    // Here we are catching if the file does not exist. We will then create a new file and return
    // back to the try statement
    catch (FileIException error)
    {
      H5File file(_filename.c_str(), H5F_ACC_TRUNC);
      file.close();
      // Just some warning that we have gone through this catch
      itr++;
      // This is to prevent us from getting caught in an infinite loop. While (true) loops
      // are useful, but they can be dangerous. Always ensure some escape sequence. Could
      // just use a for loop
      if (itr > 3)
      {
        messerr("We've tried too many times in the Int writing sequence");
        break;
      }
    }
  }
}

template<typename T>
void HDF5format::writeData(const std::vector<std::vector<T> > &data)
{
  uint itr = 0; // Used to ensure we don't get stuck in an infinite loop
  uint dim1 = data.size(); // size of our data
  uint dim2 = data[0].size();
  auto a = new T[dim1 * dim2]; // convert to an array
  auto md = new T*[dim1];
  for (size_t i = 0; i < dim1; ++i)
    md[i] = a + i * dim2;

  int vrank = 2; // since we are using std::vectors we are storing everything in one dimension

  // convert std::vector to array. H5 does not seem to like the pointer implementation
  for (size_t i = 0; i < dim1; ++i)
    for (size_t j = 0; j < dim2; ++j)
      md[i][j] = data[i][j];
  // conventional syntax for H5 data writing
  hsize_t dims[2];
  dims[0] = (int) dim1;
  dims[1] = (int) dim2;
  //hid_t memspace_id = H5Screate_simple(vrank, dims, NULL);
  // Let's make sure we are doing what we want and output it to the std output

  // loop here will check if the file exists.
  while (true)
  {
    // This assumes that the file already exists and will then write to the file
    try
    {
      H5File file(_filename.c_str(), H5F_ACC_RDWR);
      DataSpace dsp = DataSpace(vrank, dims);
      // int
      if (typeid(T).name() == typeid(int).name())
      {

        DataSet dset = file.createDataSet(_varname.c_str(), PredType::STD_I32LE, dsp);
        dset.write(a, PredType::STD_I32LE);
        dset.close();
      }
      // uint
      else if (typeid(T).name() == typeid(uint).name())
      {
        DataSet dset = file.createDataSet(_varname.c_str(), PredType::STD_U32LE, dsp);
        dset.write(a, PredType::STD_U32LE);
        dset.close();
      }
      // float
      else if (typeid(T).name() == typeid(float).name())
      {
        DataSet dset = file.createDataSet(_varname.c_str(), PredType::IEEE_F32LE, dsp);
        dset.write(a, PredType::IEEE_F32LE);
        dset.close();
      }
      // double
      else if (typeid(T).name() == typeid(double).name())
      {
        DataSet dset = file.createDataSet(_varname.c_str(), PredType::IEEE_F64LE, dsp);
        dset.write(a, PredType::IEEE_F64LE);
        dset.close();
      }
      else
      {
        messerr("Unknown data type! EXITING");
      }

      // remember to close everything and delete our arrays
      dsp.close();
      file.close();
      delete[] md;
      delete a;
      break;
    }
    // Here we are catching if the file does not exist. We will then create a new file and return
    // back to the try statement
    catch (FileIException error)
    {
      H5File file(_filename.c_str(), H5F_ACC_TRUNC);
      file.close();
      // Just some warning that we have gone through this catch
      itr++;
      // This is to prevent us from getting caught in an infinite loop. While (true) loops
      // are useful, but they can be dangerous. Always ensure some escape sequence. Could
      // just use a for loop
      if (itr > 3)
      {
        messerr("We've tried too many times in the Int writing sequence");
        break;
      }
    }
  }
}

int HDF5format::getSize() const
{
  try
  {
    H5File file(_filename.c_str(), H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(_varname);
    DataSpace dataspace = dataset.getSpace();
    const int npts = dataspace.getSimpleExtentNpoints();
    return npts;
  }
  catch (FileIException error)
  {
    int err = -1;
    return err;
  }
  catch (GroupIException error)
  {
    int err = -1;
    return err;
  }
}

int HDF5format::getDataint() const
{
  try
  {
    H5File file(_filename.c_str(), H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(_varname);
    DataType datatype = dataset.getDataType();
    DataSpace dataspace = dataset.getSpace();
    H5T_class_t classt = datatype.getClass();
    if (classt != 0)
    {
      messerr("%s is not an int... you can't save this as an int.",
              _varname.c_str());
      return -1;
    }
    int *data = new int;
    IntType itype = dataset.getIntType();

    H5std_string order_string;
    FloatType ftype = dataset.getFloatType();
    H5T_order_t order = ftype.getOrder(order_string);
    size_t size = ftype.getSize();
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
      messerr("Did not find data type");
    // Manage our memory properly
    int v = *data;
    delete data;
    dataspace.close();
    datatype.close();
    dataset.close();
    file.close();
    return v;
  }
  catch (FileIException error)
  {
    int err = -1;
    return err;
  }
  catch (GroupIException error)
  {
    int err = -1;
    return err;
  }
}

float HDF5format::getDatafloat() const
{
  try
  {
    H5File file(_filename.c_str(), H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(_varname);
    DataType datatype = dataset.getDataType();
    DataSpace dataspace = dataset.getSpace();
    H5T_class_t classt = datatype.getClass();
    if (classt != 1)
    {
      messerr("%s is not a float... you can't save this as a float.", _varname.c_str());
      return -1;
    }
    FloatType ftype = dataset.getFloatType();
    H5std_string order_string;
    H5T_order_t order = ftype.getOrder(order_string);
    size_t size = ftype.getSize();
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
    else if ((order_string == "Big endian byte order_stringing (1)"
        || order == 1)
             && size == 8)
    {
      message("NOTE: This is actually double data. We are casting to float\n");
      dataset.read((float*) data, PredType::IEEE_F64BE);
    }
    else
      messerr("Did not find data type");
    float v = *data;
    delete data;
    dataspace.close();
    datatype.close();
    dataset.close();
    file.close();
    return v;
  }
  catch (FileIException error)
  {
    float err = -1.;
    return err;
  }
  catch (GroupIException error)
  {
    float err = -1.;
    return err;
  }
}

double HDF5format::getDatadouble() const
{
  try
  {
    H5File file(_filename.c_str(), H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(_varname);
    DataType datatype = dataset.getDataType();
    DataSpace dataspace = dataset.getSpace();
    H5T_class_t classt = datatype.getClass();
    if (classt != 1)
    {
      messerr("%s is not a float... you can't save this as a float.",_varname.c_str());
      return -1.;
    }
    FloatType ftype = dataset.getFloatType();
    H5std_string order_string;
    H5T_order_t order = ftype.getOrder(order_string);
    size_t size = ftype.getSize();
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
      messerr("Did not find data type");
    float v = *data;
    delete data;
    dataspace.close();
    datatype.close();
    dataset.close();
    file.close();
    return v;
  }
  catch (FileIException error)
  {
    double err = -1.;
    return err;
  }
  catch (GroupIException error)
  {
    double err = -1.;
    return err;
  }
}

VectorInt HDF5format::getDataVint() const
{
  try
  {
    H5File file(_filename.c_str(), H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(_varname);
    DataType datatype = dataset.getDataType();
    DataSpace dataspace = dataset.getSpace();
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
    H5std_string order_string;
    H5T_order_t order = itype.getOrder(order_string);
    size_t size = itype.getSize();
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
      messerr("Did not find data type");
    VectorInt v(data, data + npts); // Arrays are nice, but vectors are better
    // Manage our memory properly
    delete[] data;
    dataspace.close();
    datatype.close();
    dataset.close();
    file.close();
    return v;
  }
  catch (FileIException error)
  {
    VectorInt err { 1, -1 };
    return err;
  }
  catch (GroupIException error)
  {
    VectorInt err { 1, -1 };
    return err;
  }
}

VectorFloat HDF5format::getDataVfloat() const
{
  try
  {
    H5File file(_filename.c_str(), H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(_varname);
    DataType datatype = dataset.getDataType();
    DataSpace dataspace = dataset.getSpace();
    const int npts = dataspace.getSimpleExtentNpoints();
    H5T_class_t classt = datatype.getClass();
    if (classt != 1)
    {
      messerr("%s is not a float... you can't save this as a float.",
              _varname.c_str());
      VectorFloat err { 1, -1. };
      return err;
    }
    FloatType ftype = dataset.getFloatType();
    H5std_string order_string;
    H5T_order_t order = ftype.getOrder(order_string);
    size_t size = ftype.getSize();
    float *data = new float[npts];
    if (order == 0 && size == 4)
      dataset.read(data, PredType::IEEE_F32LE); // Our standard integer
    else if (order == 0 && size == 8)
    {
      dataset.read((float*) data, PredType::IEEE_F64LE);
      message("NOTE: This is actually double data. We are casting to float\n");
    }
    else if (order == 1 && size == 4)
      dataset.read(data, PredType::IEEE_F32BE);
    else if ((order_string == "Big endian byte order_stringing (1)"
        || order == 1)
             && size == 8)
    {
      message("NOTE: This is actually double data We are casting to float\n");
      dataset.read((float*) data, PredType::IEEE_F64BE);
    }
    else
      messerr("Did not find data type");
    VectorFloat v(data, data + npts);
    delete[] data;
    dataspace.close();
    datatype.close();
    dataset.close();
    file.close();
    return v;
  }
  catch (FileIException error)
  {
    VectorFloat err { 1, -1. };
    return err;
  }
  catch (GroupIException error)
  {
    VectorFloat err { 1, -1. };
    return err;
  }
}

VectorDouble HDF5format::getDataVDouble() const
{
   try
   {
      H5File file(_filename.c_str(), H5F_ACC_RDONLY);
      DataSet dataset = file.openDataSet(_varname);
      DataType datatype = dataset.getDataType();
      DataSpace dataspace = dataset.getSpace();
      const int npts = dataspace.getSimpleExtentNpoints();
      H5T_class_t classt = datatype.getClass();
      if ( classt != 1 )
      {
         messerr(" is not a float... you can't save this as a float.",_varname.c_str());
         VectorDouble err{1,-1.};
         return err;
      }
      FloatType ftype = dataset.getFloatType();
      H5std_string order_string;
      H5T_order_t order = ftype.getOrder( order_string);
      size_t size = ftype.getSize();
      double *data = new double[npts];
      if ( order==0 && size == 4 )
      {
         message("NOTE: This is actually float data. We are casting to double\n");
         dataset.read((double*)data, PredType::IEEE_F32LE); // Our standard integer
      }
      else if ( order == 0 && size == 8 )
         dataset.read(data, PredType::IEEE_F64LE);
      else if ( order == 1 && size == 4 )
      {
         message("NOTE: This is actually float data. We are casting to double\n");
         dataset.read((double*)data, PredType::IEEE_F32BE);
      }
      else if ( order ==1 && size == 8 )
         dataset.read((double*)data, PredType::IEEE_F64BE);
      else
         message("Did not find data type\n");
      VectorDouble v(data, data + npts);
      delete[] data;
      dataspace.close();
      datatype.close();
      dataset.close();
      file.close();
      return v;
   }
   catch (FileIException error)
   {
      VectorDouble err{1,-1.};
      return err;
   }
   catch (GroupIException error)
   {
      VectorDouble err{1,-1.};
      return err;
   }
}

/**
 * Reading VectorVectorInt
 * @return
 */
VectorVectorInt HDF5format::getData2Dint() const
{
  try
  {
    H5File file(_filename.c_str(), H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(_varname);
    DataType datatype = dataset.getDataType();
    DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();
    hsize_t dims[rank];
    dataspace.getSimpleExtentDims(dims);
    FloatType ftype = dataset.getFloatType();
    H5std_string order_string;
    H5T_order_t order = ftype.getOrder(order_string);
    size_t size = ftype.getSize();
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
      messerr("Did not find data type");
    VectorVectorInt v(dim1, VectorInt(dim2, 0)); //data, data + npts);
    // Assign 2D vector
    for (size_t i = 0; i < dim1; ++i)
      for (size_t j = 0; j < dim2; ++j)
        v[i][j] = md[i][j];
    delete[] md;
    delete data;
    dataspace.close();
    datatype.close();
    dataset.close();
    file.close();
    return v;
  }
  catch (FileIException error)
  {
    VectorVectorInt err { 1, VectorInt(1, -1) };
    return err;
  }
  catch (GroupIException error)
  {
    VectorVectorInt err { 1, VectorInt(1, -1) };
    return err;
  }
}

VectorVectorFloat HDF5format::getData2Dfloat() const
{
  try
  {
    H5File file(_filename.c_str(), H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(_varname);
    DataType datatype = dataset.getDataType();
    DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();
    hsize_t dims[rank];
    dataspace.getSimpleExtentDims(dims);
    FloatType ftype = dataset.getFloatType();
    H5std_string order_string;
    H5T_order_t order = ftype.getOrder(order_string);
    size_t size = ftype.getSize();
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
      message("Did not find data typn");
    VectorVectorFloat v(dim1, VectorFloat(dim2, 0)); //data, data + npts);
    // Assign 2D vector
    for (size_t i = 0; i < dim1; ++i)
      for (size_t j = 0; j < dim2; ++j)
        v[i][j] = md[i][j];
    delete[] md;
    delete data;
    dataspace.close();
    datatype.close();
    dataset.close();
    file.close();
    return v;
  }
  catch (FileIException error)
  {
    VectorVectorFloat err { 1, VectorFloat(1, -1.) };
    return err;
  }
  catch (GroupIException error)
  {
    VectorVectorFloat err { 1, VectorFloat(1, -1.) };
    return err;
  }
}

VectorVectorDouble HDF5format::getData2Ddouble() const
{
  try
  {
    H5File file(_filename.c_str(), H5F_ACC_RDONLY);
    DataSet dataset = file.openDataSet(_varname);
    DataType datatype = dataset.getDataType();
    DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();
    hsize_t dims[rank];
    dataspace.getSimpleExtentDims(dims);
    FloatType ftype = dataset.getFloatType();
    H5std_string order_string;
    H5T_order_t order = ftype.getOrder(order_string);
    size_t size = ftype.getSize();
    size_t dim1 = dims[0];
    size_t dim2 = dims[1];
    //double data[dim1][dim2];
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
      messerr("Did not find data type\n");
    VectorVectorDouble v(dim1, VectorDouble(dim2, 0)); //data, data + npts);
    // Assign 2D vector
    for (size_t i = 0; i < dim1; ++i)
      for (size_t j = 0; j < dim2; ++j)
        v[i][j] = md[i][j];
    delete[] md;
    delete data;
    dataspace.close();
    datatype.close();
    dataset.close();
    file.close();
    return v;
  }
  catch (FileIException error)
  {
    VectorVectorDouble err { 1, VectorDouble(1, -1.) };
    return err;
  }
  catch (GroupIException error)
  {
    VectorVectorDouble err { 1, VectorDouble(1, -1.) };
    return err;
  }
}

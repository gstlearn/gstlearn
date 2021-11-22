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

HDF5format::HDF5format()
{
}

HDF5format::HDF5format(const HDF5format &r)
{
}

HDF5format& HDF5format::operator= (const HDF5format &r)
{
  if (this != &r)
  {
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
** \param[in]  filename  Name of the HDF5 file to be created
** \param[in]  dsname    Name of the Data Set
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
int HDF5format::create(const String& filename,
                       const String& dsname,
                       hid_t type,
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

  datafile = H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  if (datafile < 0) goto label_end;

  dataspace = H5Screate_simple(ndim, dims, NULL);
  if (dataspace < 0) goto label_end;

  err = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,start0,NULL,dims,NULL);
  if (err < 0) goto label_end;

  memspace = H5Screate_simple(ndim, dims, NULL);
  if (memspace < 0) goto label_end;

  err = H5Sselect_hyperslab(memspace,H5S_SELECT_SET,start0,NULL,dims,NULL);
  if (err < 0) goto label_end;

  dataset = H5Dcreate(datafile, dsname.c_str(), type,
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
** \param[in]  filename  Name of the HDF5 file to be created
** \param[in]  dsname    Name of the Data Set
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
void* HDF5format::read(const String& filename,
                       const String& dsname,
                       int flag_compress,
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

  datafile = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (datafile < 0) goto label_end;

  dataset = H5Dopen(datafile, dsname.c_str(), H5P_DEFAULT);
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
** \param[in]  filename  Name of the HDF5 file to be created
** \param[in]  dsname    Name of the Data Set
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
int HDF5format::write(const String& filename,
                      const String& dsname,
                      hid_t type,
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

  datafile = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (datafile < 0) goto label_end;

  dataset = H5Dopen(datafile, dsname.c_str(), H5P_DEFAULT);
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

int HDF5format::delfile(const String& filename)
{
  if( remove( filename.c_str() ) != 0 ) return 1;
  return 0;
}

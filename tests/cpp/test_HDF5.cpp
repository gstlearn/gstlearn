/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_old_f.h"
#include "Basic/HDF5format.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/File.hpp"

#include <stdlib.h>

#ifndef _USE_HDF5
int main (void)
{
  message("No HDF5 support: this test does nothing!");
  return 0;
}

#else
bool st_is_integer(H5::DataType type)
{
  return (type.getClass() == H5T_INTEGER);
}

bool st_is_float(H5::DataType type)
{
  return (type.getClass() == H5T_FLOAT && type.getSize() == 4);
}

bool st_is_double(H5::DataType type)
{
  return (type.getClass() == H5T_FLOAT && type.getSize() == 8);
}

/**
 * This test is meant to check the HDF5 read/write facility
 */
static void st_init(H5::DataType type,
                    int /*ndim*/,
                    hsize_t *dims,
                    void *data)
{
  int     jx,jy;
  int    *idata;
  float  *fdata;
  double *ddata;

  int base = 1000;
  int ecr = 0;
  if (st_is_integer(type))
  {
    idata = (int    *) data;
    for (int ix=0; ix<(int) dims[0]; ix++)
    {
      jx = base * base * (ix+1);
      for (int iy=0; iy<(int) dims[1]; iy++)
      {
        jy = base * (iy+1);
        for (int iz=0; iz<(int) dims[2]; iz++)
          idata[ecr++] = (iz+1) + jy + jx;
      }
    }
  }
  else if (st_is_float(type))
  {
    fdata = (float  *) data;
    for (int ix=0; ix<(int) dims[0]; ix++)
    {
      jx = base * base * (ix+1);
      for (int iy=0; iy<(int) dims[1]; iy++)
      {
        jy = base * (iy+1);
        for (int iz=0; iz<(int) dims[2]; iz++)
          fdata[ecr++] = (iz+1) + jy + jx;
      }
    }
  }
  else if (st_is_double(type))
  {
    ddata = (double *) data;
    for (int ix=0; ix<(int) dims[0]; ix++)
    {
      jx = base * base * (ix+1);
      for (int iy=0; iy<(int) dims[1]; iy++)
      {
        jy = base * (iy+1);
        for (int iz=0; iz<(int) dims[2]; iz++)
          ddata[ecr++] = (iz+1) + jy + jx;
      }
    }
  }
  else
  {
    messerr("Initialization has not been coded for this type of object");
  }
}

void* st_allocArray(H5::DataType type, int ndim, hsize_t *dims)
{
  hsize_t size = type.getSize();
  int ntot = 1;
  for (int idim=0; idim<ndim; idim++) ntot *= (int) dims[idim];
  void* data = (void *) calloc(ntot,(size_t) size);
  return data;
}


static void st_modify(H5::DataType type, int /*ndim*/, hsize_t *dims, void *data)
{
  int    *idata;
  float  *fdata;
  double *ddata;
  static int incr = 1;

  int ecr = 0;
  if (st_is_integer(type))
  {
    idata = (int    *) data;
    for (int ix=0; ix< (int) dims[0]; ix++)
      for (int iy=0; iy< (int) dims[1]; iy++)
        for (int iz=0; iz< (int) dims[2]; iz++)
          idata[ecr++] += incr;
  }
  else if (st_is_float(type))
  {
    fdata = (float  *) data;
    for (int ix=0; ix< (int) dims[0]; ix++)
      for (int iy=0; iy< (int) dims[1]; iy++)
        for (int iz=0; iz< (int) dims[2]; iz++)
          fdata[ecr++] += incr;
  }
  else if (st_is_double(type))
  {
    ddata = (double *) data;
    for (int ix=0; ix< (int) dims[0]; ix++)
      for (int iy=0; iy< (int) dims[1]; iy++)
        for (int iz=0; iz< (int) dims[2]; iz++)
          ddata[ecr++] += incr;
  }
  else
  {
    messerr("Modification has not been coded for this type of object");
  }
}

static void st_print(int verbose,
                     H5::DataType type,
                     int /*ndim*/,
                     hsize_t *dims,
                     void *data)
{
  if (data == NULL || ! verbose) return;
  message("\nNX(%d) panels of NY(%d) rows and NZ(%d) columns\n\n",
          (int) dims[0], (int) dims[1], (int) dims[2]);

  int lec = 0;

  if (st_is_integer(type))
  {
  int    *idata;
  idata = (int    *) data;
    for (int ix=0; ix< (int) dims[0]; ix++)
    {
      for (int iy=0; iy< (int) dims[1]; iy++)
      {
        for (int iz=0; iz< (int) dims[2]; iz++)
          message(" %8d",idata[lec++]);
        message ("\n");
      }
      message("\n");
    }
  }
  else if (st_is_float(type))
  {
    float  *fdata;
    fdata = (float  *) data;
    for (int ix=0; ix< (int) dims[0]; ix++)
    {
      for (int iy=0; iy< (int) dims[1]; iy++)
      {
        for (int iz=0; iz< (int) dims[2]; iz++)
          message(" %8.0f",fdata[lec++]);
        message ("\n");
      }
      message("\n");
    }
  }
  else if (st_is_double(type))
  {
    double *ddata;
    ddata = (double *) data;
    for (int ix=0; ix< (int) dims[0]; ix++)
    {
      for (int iy=0; iy< (int) dims[1]; iy++)
      {
        for (int iz=0; iz< (int) dims[2]; iz++)
          message(" %8.0lf",ddata[lec++]);
        message ("\n");
      }
      message("\n");
    }
  }
  else
  {
    messerr("Printout is not coded for this type of object");
  }
}

void st_print_condition_item(const String& title, int ndim, hsize_t* itab)
{
  message("%s : %d",title.c_str(),(int) itab[0]);
  for (int idim = 1; idim < ndim; idim++)
    message(" x %d",(int) itab[idim]);
  message("\n");
}

void st_print_condition(int verbose,
                        int ndim,
                        hsize_t* dims,
                        hsize_t* count,
                        hsize_t* start,
                        hsize_t* stride,
                        hsize_t* block)
{
  if (! verbose) return;
  message("Number of dimensions = %d\n",ndim);
  st_print_condition_item("Dims  ",ndim,dims);
  st_print_condition_item("Count ",ndim,count);
  st_print_condition_item("Start ",ndim,start);
  st_print_condition_item("Stride",ndim,stride);
  st_print_condition_item("Block ",ndim,block);
}

static void st_dimension(int icas,
                         int *ndim,
                         hsize_t *dims,
                         hsize_t *start,
                         hsize_t *count,
                         hsize_t *stride,
                         hsize_t *block,
                         hsize_t *start0,
                         int *flag_print)
{
  if (icas == 1)
  {
    *ndim       = 3;
    *flag_print = 1;

    dims[0]   = 3;
    dims[1]   = 7;
    dims[2]   = 8;

    start[0]  = 1;
    start[1]  = 2;
    start[2]  = 1;

    count[0]  = 2;
    count[1]  = 2;
    count[2]  = 3;

    stride[0] = 1;
    stride[1] = 2;
    stride[2] = 2;

    block[0]  = 1;
    block[1]  = 2;
    block[2]  = 1;
  }
  else
  {
    *ndim       = 3;
    *flag_print = 0;

    dims[0]   = 1000;
    dims[1]   = 1000;
    dims[2]   = 1000;

    start[0]  = 100;
    start[1]  = 110;
    start[2]  = 120;

    count[0]  = 3;
    count[1]  = 5;
    count[2]  = 7;

    stride[0] = 1;
    stride[1] = 1;
    stride[2] = 1;

    block[0]  = 1;
    block[1]  = 1;
    block[2]  = 1;
  }

  start0[0] = 0;
  start0[1] = 0;
  start0[2] = 0;
}

int main (void)
{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  hsize_t     dims[3],start[3],start0[3],stride[3],count0[3],count[3],block[3],dimout[3];
  int         ndim,flag_print;

  // Initializations

  int ipart = 0;

  // Main dispatch

  if (ipart == 0 || ipart == 1)
  {
    int icas  = 1;
    int nfois = 1;
    int niter = (icas == 1) ? 3 : 1000;
    double mult = 2.;

    // Do not use assignment operator here !
    H5::DataType type(H5::PredType::NATIVE_INT);
    bool verbose = (icas == 1);

    // Define the dimensions

    mestitle(1, "Read/Write for Regular File");
    st_dimension(icas, &ndim, dims, start, count0, stride, block, start0, &flag_print);

    // Define the HDF5 file and variable names

    HDF5format hdf5 = HDF5format();

    // Core allocation & filling the array
    void* wdata = st_allocArray(type, ndim, dims);
    if (wdata == NULL) return 1;
    st_init(type, ndim, dims, wdata);

    // Creating the HDF5 file

    if (verbose) message("Initial Array\n");
    hdf5.openNewFile("h5data1.h5");
    hdf5.openNewDataSetInt("DS1", ndim, dims);

    // Writing the Initial Information

    if (verbose) message("Creating the HDF5 file\n");
    hdf5.writeRegular(start0, NULL, dims, NULL, wdata);
    st_print(flag_print, type, ndim, dims, wdata);

    // Reading without compression

    if (verbose) message("Extraction without compression\n");
    st_print_condition(verbose,ndim,dims,count0,start,stride,block);
    void* rdata1 = hdf5.readRegular(0, start, stride, count0, block, dimout);
    st_print(flag_print, type, ndim, dimout, rdata1);
    rdata1 = (void *) mem_free((char * ) rdata1);

    // Reading with compression

    if (verbose) message("Extraction with compression\n");
    st_print_condition(verbose,ndim,dims,count0,start,stride,block);
    void* rdata2 = hdf5.readRegular(1, start, stride, count0, block, dimout);
    st_print(flag_print, type, ndim, dimout, rdata2);
    rdata2 = (void *) mem_free((char * ) rdata2);

    // Loop on multiple of chunk dimensions

    for (int ifois = 0; ifois < nfois; ifois++)
    {
      for (int idim = 0; idim < ndim; idim++)
        count[idim] = count0[idim] * (ifois+0.5) * mult;

      // Loop on iterations

      double total_read = 0.;
      double total_write = 0.;
      for (int iter = 0; iter < niter; iter++)
      {
        // Reading

        if (verbose) message("Modifying by adding a unit to all terms (%d times)\n",iter+1);
        void* rdata3 = hdf5.readRegular(1, start, stride, count, block, dimout);

        // Modifying the array

        st_modify(type, ndim, dimout, rdata3);

        // Writing

        hdf5.writeRegular(start, stride, count, block, rdata3);
        st_print(flag_print, type, ndim, dimout, rdata3);
        rdata3 = (void *) mem_free((char * ) rdata3);
      }
    }

    // Delete the file

    hdf5.deleteFile();
    free(wdata);
    wdata = NULL;
  }

  // Define the HDF5 file and variable names

  if (ipart == 0 || ipart == 2)
  {
#define dim0 10
#define dim1 7
#define dim2 5

    law_set_random_seed(32121.);
    VectorInt ival(dim0);
    for (size_t i = 0; i < dim0; i++)
      ival[i] = (i + 10);
    VectorFloat fval(dim0);
    for (size_t i = 0; i < dim0; i++)
      fval[i] = (float) law_uniform(0.,1.);
    VectorDouble dval(dim0);
    for (size_t i = 0; i < dim0; i++)
      dval[i] = law_uniform(0.,1.);
    VectorVectorInt vival(dim1, VectorInt(dim2));
    for (size_t i = 0; i < dim1; ++i)
      for (size_t j = 0; j < dim2; ++j)
        vival[i][j] = law_uniform(0.,1.);
    VectorVectorFloat vfval(dim1, VectorFloat(dim2));
    for (size_t i = 0; i < dim1; ++i)
      for (size_t j = 0; j < dim2; ++j)
        vfval[i][j] = law_uniform(0.,1.);
    VectorVectorDouble vdval(dim1, VectorDouble(dim2));
    for (size_t i = 0; i < dim1; ++i)
      for (size_t j = 0; j < dim2; ++j)
        vdval[i][j] = law_uniform(0.,1.);

    mestitle(1, "Read/Write series of same type");
    HDF5format hdf5b = HDF5format();

    // Write
    hdf5b.openNewFile("h5data2.h5");

    ndim = 0;
    dims[0] = 1;

    hdf5b.openNewDataSetInt("ValueInt", ndim, dims);
    hdf5b.writeData(ival[0]);
    hdf5b.closeDataSet();

    hdf5b.openNewDataSetFloat("ValueFloat", ndim, dims);
    hdf5b.writeData(fval[0]);
    hdf5b.closeDataSet();

    hdf5b.openNewDataSetDouble("ValueDouble", ndim, dims);
    hdf5b.writeData(dval[0]);
    hdf5b.closeDataSet();

    ndim = 1;
    dims[0] = dim0;

    hdf5b.openNewDataSetInt("VectorInt", ndim, dims);
    hdf5b.writeData(ival);
    hdf5b.closeDataSet();

    hdf5b.openNewDataSetFloat("VectorFloat", ndim, dims);
    hdf5b.writeData(fval);
    hdf5b.closeDataSet();

    hdf5b.openNewDataSetDouble("VectorDouble", ndim, dims);
    hdf5b.writeData(dval);
    hdf5b.closeDataSet();

    ndim = 2;
    dims[0] = dim1;
    dims[1] = dim2;

    hdf5b.openNewDataSetInt("VectorVectorInt", ndim, dims);
    hdf5b.writeData(vival);
    hdf5b.closeDataSet();

    hdf5b.openNewDataSetFloat("VectorVectorFloat", ndim, dims);
    hdf5b.writeData(vfval);
    hdf5b.closeDataSet();

    hdf5b.openNewDataSetDouble("VectorVectorDouble", ndim, dims);
    hdf5b.writeData(vdval);
    hdf5b.closeDataSet();

    // Get the list of data sets
    hdf5b.displayNames();

    // To Load Data
    hdf5b.openDataSet("ValueInt");
    int rival0 = hdf5b.getData();
    if (rival0 != ival[0]) messageAbort("Error when handling Int");
    hdf5b.closeDataSet();

    hdf5b.openDataSet("ValueFloat");
    float rfval0 = hdf5b.getData();
    if (rfval0 != fval[0]) messageAbort("Error when handling Float");
    hdf5b.closeDataSet();

    hdf5b.openDataSet("ValueDouble");
    double rdval0 = hdf5b.getData();
    if (rdval0 != dval[0]) messageAbort("Error when handling Double");
    hdf5b.closeDataSet();

    hdf5b.openDataSet("VectorInt");
    VectorInt rival = hdf5b.getData();
    if (ival != rival) messageAbort("Error when handling VectorInt");
    hdf5b.closeDataSet();

    hdf5b.openDataSet("VectorFloat");
    VectorFloat rfval = hdf5b.getData();
    if (fval != rfval) messageAbort("Error when handling VectorFloat");
    hdf5b.closeDataSet();

    hdf5b.openDataSet("VectorDouble");
    VectorDouble rdval = hdf5b.getData();
    if (dval != rdval) messageAbort("Error when handling VectorDouble");
    hdf5b.closeDataSet();

    hdf5b.openDataSet("VectorVectorInt");
    VectorVectorInt rvival = hdf5b.getData();
    if (vival != rvival) messageAbort("Error when handling VectorVectorInt");
    hdf5b.closeDataSet();

    hdf5b.openDataSet("VectorVectorFloat");
    VectorVectorFloat rvfval = hdf5b.getData();
    if (vfval != rvfval) messageAbort("Error when handling VectorVectorFloat");
    hdf5b.closeDataSet();

    hdf5b.openDataSet("VectorVectorDouble");
    VectorVectorDouble rvdval = hdf5b.getData();
    if (vdval != rvdval) messageAbort("Error when handling VectorVectorDouble");
    hdf5b.closeDataSet();

    // Extract a row (VectorDouble) of the VectorVectorDouble file

    hdf5b.openDataSet("VectorVectorDouble");
    rvdval = hdf5b.getData();
    VH::display("Ensemble of VectorDouble",rvdval);

    int myrank = 3;
    message("Extract Vector #%d (add 100) and replace\n\n",myrank);
    VectorDouble rpdval = hdf5b.getDataDoublePartial(myrank);
    VH::addConstant(rpdval, 100.);
    hdf5b.writeDataDoublePartial(myrank, rpdval);

    // Extract the whole file and print it
    rvdval = hdf5b.getData();
    VH::display("Modified VectorVectorDouble",rvdval);

    // Delete the file

    hdf5b.closeDataSet();
    hdf5b.deleteFile();

#undef dim0
#undef dim1
#undef dim2
  }

  // Defining a HDF5 file and file it incrementally

  if (ipart == 0 || ipart == 3)
  {
#define dim1 8
#define dim2 5

    VectorVectorDouble vdval(dim1, VectorDouble(dim2));
    for (size_t i = 0; i < dim1; ++i)
      for (size_t j = 0; j < dim2; ++j)
        vdval[i][j] = i + j;

    mestitle(1, "Read/Write VectorDouble in a file created incrementally");
    HDF5format hdf5c = HDF5format();

    // Create the empty file
    ndim = 2;
    dims[0] = dim1;
    dims[1] = dim2;
    hdf5c.openNewFile("h5data3.h5");
    hdf5c.openNewDataSetDouble("Set3", ndim, dims);

    // Store VectorDouble incrementally
    VectorDouble rowval(dim2);
    for (int irow = 0; irow < dim1; irow++)
    {
      // Define the values in the row
      for (int icol = 0; icol < dim2; icol++) rowval[icol] = 100 * (1+irow) + (1+icol);

      // Store the row (VectorDouble)
      hdf5c.writeDataDoublePartial(irow, rowval);
    }

    // Extract the whole file and print it
    VectorVectorDouble rpvdval = hdf5c.getData();
    VH::display("Initial VectorDouble",rpvdval);

    // Extracting one VectorDouble (at a given row number). Modify it by adding 1000

    int myrank = 5;
    message("\nExtracting the rank #%d (then add 1000) and replace\n\n",myrank);
    VectorDouble rpdval = hdf5c.getDataDoublePartial(myrank);
    VH::addConstant(rpdval, 1000.);

    hdf5c.writeDataDoublePartial(myrank, rpdval);

    // Extract the whole file and print it
    rpvdval = hdf5c.getData();
    VH::display("Modified VectorVectorDouble",rpvdval);

    // Delete the file

    hdf5c.closeDataSet();
    hdf5c.closeFile();
    hdf5c.deleteFile();

#undef dim1
#undef dim2
  }

  // Core deallocation

  return 0;
}
#endif

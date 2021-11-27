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
#include "Basic/AStringable.hpp"
#include "Basic/Timer.hpp"
#include "geoslib_old_f.h"
#include <malloc.h>

/**
 * This test is meant to check the HDF5 read/write facility
 */
static void st_init(int verbose,
                    hid_t type,
                    int ndim,
                    hsize_t *dims,
                    void *data)
{
  int     jx,jy;
  int    *idata;
  float  *fdata;
  double *ddata;

  int base = 1000;

  if (verbose)
  {
    message("  Rule: TAB(ix,iy,iz) = (iz+1) + (iy+1)*base + (ix+1)*base^2\n");
    message("  Base: %d\n",base);
  }

  int ecr = 0;
  if (type == H5T_NATIVE_INT)
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
  else if (type == H5T_NATIVE_FLOAT)
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
  else if (type == H5T_NATIVE_DOUBLE)
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

static void st_modify(hid_t type, int ndim, hsize_t *dims, void *data)
{
  int    *idata;
  float  *fdata;
  double *ddata;
  static int incr = 1;

  int ecr = 0;
  if (type == H5T_NATIVE_INT)
  {
    idata = (int    *) data;
    for (int ix=0; ix< (int) dims[0]; ix++)
      for (int iy=0; iy< (int) dims[1]; iy++)
        for (int iz=0; iz< (int) dims[2]; iz++)
          idata[ecr++] += incr;
  }
  else if (type == H5T_NATIVE_FLOAT)
  {
    fdata = (float  *) data;
    for (int ix=0; ix< (int) dims[0]; ix++)
      for (int iy=0; iy< (int) dims[1]; iy++)
        for (int iz=0; iz< (int) dims[2]; iz++)
          fdata[ecr++] += incr;
  }
  else if (type == H5T_NATIVE_DOUBLE)
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
                     hid_t type,
                     int ndim,
                     hsize_t *dims,
                     void *data)
{
  if (data == NULL || ! verbose) return;
  message("\nNX(%d) panels of NY(%d) rows and NZ(%d) columns\n\n",
          (int) dims[0], (int) dims[1], (int) dims[2]);

  int lec = 0;

  if (type == H5T_NATIVE_INT)
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
  else if (type == H5T_NATIVE_FLOAT)
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
  else if (type == H5T_NATIVE_DOUBLE)
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

static void st_dimension(int icas,
                         int *ndim,
                         hsize_t *dims,
                         hsize_t *start,
                         hsize_t *count,
                         hsize_t *stride,
                         hsize_t *block,
                         int *flag_print)
{
  if (icas == 1)
  {
    *ndim       = 3;
    *flag_print = 1;

    dims[0]   = 3;
    dims[1]   = 8;
    dims[2]   = 9;

    start[0]  = 1;
    start[1]  = 2;
    start[2]  = 1;

    count[0]  = 2;
    count[1]  = 2;
    count[2]  = 4;

    stride[0] = 1;
    stride[1] = 3;
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
}

int main (void)
{
  hsize_t     dims[3],start[3],stride[3],count0[3],count[3],block[3],dimout[3];
  int         ndim,flag_compress,flag_print;
  Timer       timer;

  // Initializations

  int verbose  = 0;
  int ipart    = 0;

  // Main dispatch

  if (ipart == 0 || ipart == 1)
  {
#define FILE    "h5data1.h5"
#define DATASET "DS1"

    int icas = 2;
    int nfois = 1;
    int niter = 1000;
    double mult = 10.;
    hid_t type = H5T_NATIVE_INT;

    // Define the dimensions

    mestitle(1, "Read/Write for Regular File");
    st_dimension(icas, &ndim, dims, start, count0, stride, block, &flag_print);

    // Define the HDF5 file and variable names

    HDF5format hdf5(FILE, DATASET);

    // Core allocation & filling the array

    void* wdata = hdf5.allocArray(type, ndim, dims);
    if (wdata == NULL) return 1;
    st_init(verbose, type, ndim, dims, wdata);

    // Creating the HDF5 file

    hdf5.createRegular(type, ndim, dims, wdata);
    st_print(flag_print, type, ndim, dims, wdata);
    timer.Interval("Creating the HDF5 file");

    // Reading without compression

    flag_compress = 0;
    void* rdata1 = hdf5.readRegular(flag_compress, type, ndim, start, stride,
                                    count0, block, dimout);
    st_print(flag_print, type, ndim, dimout, rdata1);
    rdata1 = (void *) mem_free((char * ) rdata1);
    timer.Interval("Reading HDF5 array (no compression)");

    // Reading with compression

    flag_compress = 1;
    void* rdata2 = hdf5.readRegular(flag_compress, type, ndim, start, stride,
                                    count0, block, dimout);
    rdata2 = (void *) mem_free((char * ) rdata2);
    timer.Interval("Reading HDF5 array (with compression)");

    // Loop on multiple of chunk dimensions

    for (int ifois = 0; ifois < nfois; ifois++)
    {
      for (int idim = 0; idim < ndim; idim++)
        count[idim] = count0[idim] * (ifois + 1) * mult;

      // Loop on iterations

      Timer timer_read;
      Timer timer_write;
      double total_read = 0.;
      double total_write = 0.;
      for (int iter = 0; iter < niter; iter++)
      {
        // Reading

        void* rdata3 = hdf5.readRegular(1, type, ndim, start, stride, count,
                                        block, dimout);
        st_print(flag_print, type, ndim, dimout, rdata3);
        total_read += timer_read.getInterval();

        // Modifying the array

        st_modify(type, ndim, dimout, rdata3);

        // Writing

        hdf5.writeRegular(type, ndim, dimout, start, stride, count, block,
                          rdata3);
        rdata3 = (void *) mem_free((char * ) rdata3);
        total_write += timer_write.getInterval();
      }
      timer_read.display("Reading HDF5 modified array", total_read);
      timer_write.display("Writing HDF5 modified array", total_write);

      // Print time scores

      double n1 = 1.;
      double n2 = 1.;
      for (int idim = 0; idim < ndim; idim++)
      {
        n1 *= dims[idim];
        n2 *= count[idim] * block[idim];
      }
      message("Complete Grid (uncompressed): %d %d %d -> %d\n", (int) dims[0],
              (int) dims[1], (int) dims[2], (int) n1);
      message("Local Grid (compressed): %d(x%d) %d(x%d) %d(x%d) -> %d\n",
              (int) count[0], (int) block[0], (int) count[1], (int) block[1],
              (int) count[2], (int) block[2], (int) n2);
    }

    // Delete the file

    hdf5.deleteFile();
    free(wdata);
    wdata = NULL;

#undef FILE
#undef DATASET
  }

  // Define the HDF5 file and variable names

  if (ipart == 0 || ipart == 2)
  {
#define FILE "h5data2.h5"
#define dim0 10
#define dim1 7
#define dim2 5

    VectorInt ival(dim0);
    for (size_t i = 0; i < dim0; i++)
      ival[i] = i;
    VectorFloat fval(dim0);
    for (size_t i = 0; i < dim0; i++)
      fval[i] = i * 0.1;
    VectorDouble dval(dim0);
    for (size_t i = 0; i < dim0; i++)
      dval[i] = i * 0.1;
    VectorVectorInt vival(dim1, VectorInt(dim2));
    for (size_t i = 0; i < dim1; ++i)
      for (size_t j = 0; j < dim2; ++j)
        vival[i][j] = i + j;
    VectorVectorDouble vdval(dim1, VectorDouble(dim2));
    for (size_t i = 0; i < dim1; ++i)
      for (size_t j = 0; j < dim2; ++j)
        vdval[i][j] = i + j;

    mestitle(1, "Read/Write series of same type");
    HDF5format hdf5b(FILE);

    // Write
    hdf5b.setVarName("ValueInt");
    hdf5b.writeData(ival[0]);
    hdf5b.setVarName("ValueFloat");
    hdf5b.writeData(fval[0]);
    hdf5b.setVarName("ValueDouble");
    hdf5b.writeData(dval[0]);
    hdf5b.setVarName("VectorInt");
    hdf5b.writeData(ival);
    hdf5b.setVarName("VectorFloat");
    hdf5b.writeData(fval);
    hdf5b.setVarName("VectorDouble");
    hdf5b.writeData(dval);
    hdf5b.setVarName("VectorVectorInt");
    hdf5b.writeData(vival);
    hdf5b.setVarName("VectorVectorDouble");
    hdf5b.writeData(vdval);

    // Get the list of data sets
    hdf5b.displayNames();

    // To Load Data
//    hdf5b.setVarName("ValueInt");
//    int rival0 = hdf5b.getData();
//    if (rival0 != ival[0]) messageAbort("Error when handling Int");
//    hdf5b.setVarName("ValueFloat");
//    float rfval0 = hdf5b.getData();
//    if (rfval0 != fval[0]) messageAbort("Error when handling Float");
    hdf5b.setVarName("ValueDouble");
    double rdval0 = hdf5b.getData();
    if (rdval0 != dval[0]) messageAbort("Error when handling Double");
    hdf5b.setVarName("VectorInt");
    VectorInt rival = hdf5b.getData();
    if (ival != rival) messageAbort("Error when handling VectorInt");
    hdf5b.setVarName("VectorFloat");
    VectorFloat rfval = hdf5b.getData();
    if (fval != rfval) messageAbort("Error when handling VectorFloat");
    hdf5b.setVarName("VectorDouble");
    VectorDouble rdval = hdf5b.getData();
    if (dval != rdval) messageAbort("Error when handling VectorDouble");
    hdf5b.setVarName("VectorVectorInt");
    VectorVectorInt rvival = hdf5b.getData();
    if (vival != rvival) messageAbort("Error when handling VectorVectorInt");
    hdf5b.setVarName("VectorVectorDouble");
    VectorVectorDouble rvdval = hdf5b.getData();
    if (vdval != rvdval) messageAbort("Error when handling VectorVectorDouble");

    // Extract a row (VectorDouble) of the VectorVectorDouble file
    hdf5b.setVarName("VectorVectorDouble");
    ut_vector_display("Ensemble of VectorDouble",rvdval);

    int myrank = 3;
    message("\nExtracting the rank #%d (then add 100)\n",myrank);
    VectorDouble rpdval = hdf5b.getDataDoublePartial(myrank);
    ut_vector_addval(rpdval, 100.);
    hdf5b.writeDataDoublePartial(myrank, rpdval);

    // Extract the whole file and print it
    VectorVectorDouble rpvdval = hdf5b.getData();
    ut_vector_display("Modified VectorVectorDouble",rpvdval);

    // Delete the file

    hdf5b.deleteFile();

#undef FILE
#undef dim0
#undef dim1
#undef dim2
  }

  // Defining a HDF5 file and file it incrementally

  if (ipart == 0 || ipart == 3)
  {
#define FILE    "h5data3.h5"
#define DATASET "Set3"
#define dim1 8
#define dim2 5

    VectorVectorDouble vdval(dim1, VectorDouble(dim2));
    for (size_t i = 0; i < dim1; ++i)
      for (size_t j = 0; j < dim2; ++j)
        vdval[i][j] = i + j;

    mestitle(1, "Read/Write VectorDouble in a file created incrementally");
    HDF5format hdf5c(FILE, DATASET);

    // Create the empty file
    VectorInt dims = { dim1, dim2 };
    hdf5c.createData(dims, 4);

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
    ut_vector_display("Initial VectorDouble",rpvdval);

    // Extracting one VectorDouble (at a given row number). Modify it by adding 1000

    int myrank = 5;
    message("\nExtracting the rank #%d (then add 1000)\n",myrank);
    VectorDouble rpdval = hdf5c.getDataDoublePartial(myrank);
    ut_vector_addval(rpdval, 1000.);

    hdf5c.writeDataDoublePartial(myrank, rpdval);

    // Extract the whole file and print it
    rpvdval = hdf5c.getData();
    ut_vector_display("Modified VectorVectorDouble",rpvdval);

    // Delete the file

    hdf5c.deleteFile();

#undef dim1
#undef dim2
  }

  // Core deallocation

  return 0;
}

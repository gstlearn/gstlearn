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

#define FILE            "h5data.h5"
#define DATASET         "DS1"

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
    *ndim     = 3;
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
    *ndim     = 3;
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
  hid_t       type;
  int         ndim,verbose,icas,flag_compress,flag_print,niter,nfois;
  void       *wdata,*rdata;
  double      n1,n2,mult;
  Timer       timer;
  HDF5format  hdf5;

  // Initializations

  verbose  = 0;
  icas     = 2;
  nfois    = 1;
  niter    = 1000;
  mult     = 10.;
  type     = H5T_NATIVE_INT;
  wdata    = rdata = (void *) NULL;

  // Define the dimensions

  st_dimension(icas,&ndim,dims,start,count0,stride,block,&flag_print);

  // Core allocation

  wdata = hdf5.allocArray(type,ndim,dims);
  if (wdata == NULL) goto label_end;
  timer.Interval("Allocating the array");

  // Filling the data

  st_init(verbose,type,ndim,dims,wdata);
  timer.Interval("Initializing the Data");

  // Creating the HDF5 file

  hdf5.create(FILE,DATASET,type,ndim,dims,wdata);
  st_print(flag_print,type,ndim,dims,wdata);
  timer.Interval("Creating the HDF5 file");

  // Reading without compression

  flag_compress = 0;
  rdata = hdf5.read(FILE,DATASET,flag_compress,type,
                    ndim,start,stride,count0,block,dimout);
  st_print(flag_print,type,ndim,dimout,rdata);
  timer.Interval("Reading HDF5 array (no compression)");

  // Reading with compression

  rdata = (void *) mem_free((char *) rdata);
  flag_compress = 1;
  rdata = hdf5.read(FILE,DATASET,flag_compress,type,
                    ndim,start,stride,count0,block,dimout);
  timer.Interval("Reading HDF5 array (with compression)");

  // Loop on multiple of chunk dimensions

  for (int ifois=0; ifois<nfois; ifois++)
  {
    for (int idim=0; idim<ndim; idim++)
      count[idim] = count0[idim] * (ifois+1) * mult;

    // Loop on iterations

    Timer timer_read;
    Timer timer_write;
    double total_read  = 0.;
    double total_write = 0.;
    for (int iter=0; iter<niter; iter++)
    {
      // Reading

      rdata = (void *) mem_free((char *) rdata);
      rdata = hdf5.read(FILE,DATASET,1,type,
                        ndim,start,stride,count,block,dimout);
      st_print(flag_print,type,ndim,dimout,rdata);
      total_read += timer_read.getInterval();

      // Modifying the array

      st_modify(type,ndim,dimout,rdata);

      // Writing

      if (hdf5.write(FILE,DATASET,type,ndim,dimout,start,stride,count,block,
                     rdata)) goto label_end;
      total_write += timer_write.getInterval();
    }
    timer_read.display("Reading HDF5 modified array",total_read);
    timer_write.display("Writing HDF5 modified array",total_write);

    // Print time scores

    n1 = n2 = 1.;
    for (int idim=0; idim<ndim; idim++)
    {
      n1 *= dims[idim];
      n2 *= count[idim] * block[idim];
    }
    message("Complete Grid (uncompressed): %d %d %d -> %d\n",
           (int) dims[0],(int) dims[1],(int) dims[2],(int) n1);
    message("Local Grid (compressed): %d(x%d) %d(x%d) %d(x%d) -> %d\n",
           (int) count[0],(int) block[0],
           (int) count[1],(int) block[1],
           (int) count[2],(int) block[2],(int) n2);
  }

  // Delete the file

  hdf5.delfile(FILE);

  // Core deallocation

label_end:
  free(wdata); wdata = NULL;
  return 0;
}

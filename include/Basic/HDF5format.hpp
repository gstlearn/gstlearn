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
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/String.hpp"
#include "hdf5.h"

class GSTLEARN_EXPORT HDF5format
{
public:
  HDF5format();
  HDF5format(const HDF5format &r);
  HDF5format& operator=(const HDF5format &r);
  virtual ~HDF5format();

public:
  int createRegular(const String& filename,
                    const String& dsname,
                    hid_t type,
                    int ndim,
                    hsize_t *dims,
                    void *wdata);
  void* readRegular(const String& filename,
                    const String& dsname,
                    int flag_compress,
                    hid_t type,
                    int ndim,
                    hsize_t *start,
                    hsize_t *stride,
                    hsize_t *count,
                    hsize_t *block,
                    hsize_t *dimout);
  int writeRegular(const String& filename,
                   const String& dsname,
                   hid_t type,
                   int ndim,
                   hsize_t *dims,
                   hsize_t *start,
                   hsize_t *stride,
                   hsize_t *count,
                   hsize_t *block,
                   void *wdata);
  int delfile(const String& filename);
  void* allocArray(hid_t type, int ndim, hsize_t *dims);
};

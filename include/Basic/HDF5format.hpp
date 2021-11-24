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

#include "Basic/String.hpp"
#include "hdf5.h"

class HDF5format
{
public:
  HDF5format(const String& filename = String(), const String& varname = String());
  HDF5format(const HDF5format &r);
  HDF5format& operator=(const HDF5format &r);
  virtual ~HDF5format();

public:
  int createRegular(hid_t type,
                    int ndim,
                    hsize_t *dims,
                    void *wdata);
  void* readRegular(int flag_compress,
                    hid_t type,
                    int ndim,
                    hsize_t *start,
                    hsize_t *stride,
                    hsize_t *count,
                    hsize_t *block,
                    hsize_t *dimout);
  int writeRegular(hid_t type,
                   int ndim,
                   hsize_t *dims,
                   hsize_t *start,
                   hsize_t *stride,
                   hsize_t *count,
                   hsize_t *block,
                   void *wdata);
  int deleteFile(const String& filename);
  void* allocArray(hid_t type, int ndim, hsize_t *dims);

  void setFilename(const String& filename) { _filename = filename; }
  void setVarname(const String& varname) { _varname = varname; }

  // Functions to be overloaded
  template<typename T>
  void writeData(const T&);
  template<typename T>
  void writeData(const std::vector<T>&);
  template<typename T>
  void writeData(const std::vector<std::vector<T> >&);

  int getDataint() const;
  float getDatafloat() const;
  double getDatadouble() const;
  VectorInt getDataVint() const;
  VectorFloat getDataVfloat() const;
  VectorDouble getDataVDouble() const;
  VectorVectorInt getData2Dint() const;
  VectorVectorFloat getData2Dfloat() const;
  VectorVectorDouble getData2Ddouble() const;
  // Return the size of the data
  // Note that for multi-dim arrays that it gets the total size and not the size of a single row.
  int getSize() const;

public:
  String _filename;
  String _varname;
};

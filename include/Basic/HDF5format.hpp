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
  int deleteFile();
  void* allocArray(hid_t type, int ndim, hsize_t *dims);

  void setFileName(const String& filename) { _filename = filename; }
  void setVarName(const String& varname) { _varname = varname; }

  // Functions to be overloaded
  template<typename T>
  void writeData(const T&);
  template<typename T>
  void writeData(const std::vector<T>&);
  template<typename T>
  void writeData(const std::vector<std::vector<T> >&);

  int getDataInt() const;
  float getDataFloat() const;
  double getDataDouble() const;
  VectorInt getDataVInt() const;
  VectorFloat getDataVFloat() const;
  VectorDouble getDataVDouble() const;
  VectorVectorInt getDataVVInt() const;
  VectorVectorFloat getDataVVFloat() const;
  VectorVectorDouble getDataVVDouble() const;
  // Return the size of the data
  // Note that for multi-dim arrays that it gets the total size and not the size of a single row.
  int getSize() const;

  // We now make a proxy class so that we can overload the return type and use a single
  // function to get data whether int or float. This could be made more advanced by
  // adding more data types (such as double).
  class Proxy
  {
  private:
    HDF5format const* myOwner;
  public:
    Proxy(const HDF5format* owner)
        : myOwner(owner)
    {
    }
    operator int() const
    {
      return myOwner->getDataInt();
    }
    operator float() const
    {
      return myOwner->getDataFloat();
    }
    operator double() const
    {
      return myOwner->getDataDouble();
    }
    operator std::vector<int>() const
    {
      return myOwner->getDataVInt();
    }
    operator std::vector<float>() const
    {
      return myOwner->getDataVFloat();
    }
    operator std::vector<double>() const
    {
      return myOwner->getDataVDouble();
    }
    operator std::vector<std::vector<int> >() const
    {
      return myOwner->getDataVVInt();
    }
    operator std::vector<std::vector<float> >() const
    {
      return myOwner->getDataVVFloat();
    }
    operator std::vector<std::vector<double> >() const
    {
      return myOwner->getDataVVDouble();
    }
  };
  // Here we use the Proxy class to have a single getData function
  Proxy getData() const
  {
    return Proxy(this);
  }

public:
  String _filename;
  String _varname;
};

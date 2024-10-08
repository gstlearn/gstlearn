/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Basic/VectorNumT.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#ifdef _USE_HDF5
#include <H5Cpp.h>

#if H5_VERSION_GE(1,8,20)
#define EXCEPTION_PRINT_ERROR(e) e.printErrorStack();
#else
#define EXCEPTION_PRINT_ERROR(e) e.printError();
#endif
#endif

class GSTLEARN_EXPORT HDF5format
{
public:
  HDF5format(const String& filename = "", const String& varname = "");
  HDF5format(const HDF5format &r);
  HDF5format& operator=(const HDF5format &r);
  virtual ~HDF5format();

public:
#ifdef _USE_HDF5
  void* readRegular(int flag_compress,
                    hsize_t *start,
                    hsize_t *stride,
                    hsize_t *count,
                    hsize_t *block,
                    hsize_t *dimout);
  int writeRegular(hsize_t *start,
                   hsize_t *stride,
                   hsize_t *count,
                   hsize_t *block,
                   void *wdata);
#endif

  int deleteFile() const;

  void setFileName(const String& filename) { _filename = filename; }
  void setVarName(const String& varname) { _varname = varname; }

  static int displayNames();

  void openFile(const String& filename = "");
  void openNewFile(const String& filename);
  void openDataSet(const String& varname = "");

#ifdef _USE_HDF5
  void openNewDataSetInt(const String& varname,
                         int ndim,
                         hsize_t *dims);
  void openNewDataSetFloat(const String& varname,
                           int ndim,
                           hsize_t *dims);
  void openNewDataSetDouble(const String& varname,
                            int ndim,
                            hsize_t *dims);
#endif

  void closeFile();
  void closeDataSet();

  // Functions to be overloaded

  template<typename T>
  void writeData(const T&);
  template<typename T>
  void writeData(const VectorT<T>&);
  template<typename T>
  void writeData(const VectorNumT<T>&);
  template<typename T>
  void writeData(const VectorT<VectorNumT<T> >&);

  static int getDataInt();
  static float getDataFloat();
  static double getDataDouble();
  static VectorInt getDataVInt();
  static VectorFloat getDataVFloat();
  static VectorDouble getDataVDouble();
  static VectorVectorInt getDataVVInt();
  static VectorVectorFloat getDataVVFloat();
  static VectorVectorDouble getDataVVDouble();

  static VectorDouble getDataDoublePartial(int myrank);
  static int writeDataDoublePartial(int myrank, const VectorDouble& data);

  // Return the size of the data
  // Note that for multi-dim arrays that it gets the total size and not the size of a single row.
  static int getSize();

  // We now make a proxy class so that we can overload the return type and use a single
  // function to get data whether int or float or double.
  class Proxy
  {
  public:
    Proxy(const HDF5format* owner)
    {
      DECLARE_UNUSED(owner);
    }
    operator int() const
    {
      return HDF5format::getDataInt();
    }
    operator float() const
    {
      return HDF5format::getDataFloat();
    }
    operator double() const
    {
      return HDF5format::getDataDouble();
    }
    operator VectorInt() const
    {
      return HDF5format::getDataVInt();
    }
    operator VectorFloat() const
    {
      return HDF5format::getDataVFloat();
    }
    operator VectorDouble() const
    {
      return HDF5format::getDataVDouble();
    }
    operator VectorVectorInt() const
    {
      return HDF5format::getDataVVInt();
    }
    operator VectorVectorFloat() const
    {
      return HDF5format::getDataVVFloat();
    }
    operator VectorVectorDouble() const
    {
      return HDF5format::getDataVVDouble();
    }
  };
  // Here we use the Proxy class to have a single getData function
  Proxy getData() const
  {
    return Proxy(this);
  }

private:
  static int _getNDim();

#ifdef _USE_HDF5
  hsize_t* _getDims() const;
  void _getOrderSize(H5T_order_t* order, size_t* size, bool* big_endian) const;
  void* _allocArray(const H5::DataType& myh5type, int ndim, hsize_t *dims);
  void _readInt(int *data,
                const H5::DataSpace& memspace = H5::DataSpace::ALL,
                const H5::DataSpace& dataspace = H5::DataSpace::ALL) const;
  void _readFloat(float *data,
                  const H5::DataSpace&  = H5::DataSpace::ALL,
                  const H5::DataSpace& dataspace = H5::DataSpace::ALL) const;
  void _readDouble(double *data,
                   const H5::DataSpace& memspace = H5::DataSpace::ALL,
                   const H5::DataSpace& dataspace = H5::DataSpace::ALL) const;
  void _writeDouble(double *data,
                    const H5::DataSpace& memspace = H5::DataSpace::ALL,
                    const H5::DataSpace& dataspace = H5::DataSpace::ALL) const;
#endif

  static int _checkClass(int value);
  static void _writeAll(const char* myh5type, void* a);

public:
  String        _filename;
  String        _varname;
#ifdef _USE_HDF5
  H5::H5File    _datafile;
  H5::DataSet   _dataset;
  H5::DataType  _datatype;
  H5::DataSpace _dataspace;
#endif
};

/**
 * Numeric implementation of our write data function
 * Only accepts numerical values. Integers, floats, or doubles
 * @param data
 */
template<typename T>
void HDF5format::writeData(const T &data)
{
#ifdef _USE_HDF5
  H5::Exception::dontPrint();
  char* myh5type = (char*) (typeid(T).name());
  auto *a = new T { data };
  _writeAll(myh5type, (void*) a);
  delete a;
#endif
}

template<typename T>
void HDF5format::writeData(const VectorT<T> &data)
{
#ifdef _USE_HDF5
  H5::Exception::dontPrint();
  size_t npts = data.size();
  auto *a = new T[npts];
  char* myh5type = (char*) (typeid(a[0]).name());
  for (size_t i = 0; i < npts; ++i)
    a[i] = data[i];
  _writeAll(myh5type, (void*) a);
  delete[] a;
#endif
 }

template<typename T>
void HDF5format::writeData(const VectorNumT<T> &data)
{
#ifdef _USE_HDF5
  H5::Exception::dontPrint();
  size_t npts = data.size();
  auto *a = new T[npts];
  char* myh5type = (char*) (typeid(a[0]).name());
  for (size_t i = 0; i < npts; ++i)
    a[i] = data[i];
  _writeAll(myh5type, (void*) a);
  delete[] a;
#endif
 }

template<typename T>
void HDF5format::writeData(const VectorT<VectorNumT<T> > &data)
{
#ifdef _USE_HDF5
  H5::Exception::dontPrint();
  size_t dim1 = data.size();
  size_t dim2 = data[0].size();
  auto a = new T[dim1 * dim2];
  auto md = new T*[dim1];
  for (size_t i = 0; i < dim1; ++i)
    md[i] = a + i * dim2;
  for (size_t i = 0; i < dim1; ++i)
    for (size_t j = 0; j < dim2; ++j)
      md[i][j] = data[i][j];
  char* myh5type = (char*) (typeid(a[0]).name());
  _writeAll(myh5type, (void* ) a);
  delete[] md;
  delete[] a;
#endif
}


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
#include "Basic/AStringable.hpp"
#include "hdf5.h"
#include "H5Cpp.h"
#include "typeinfo"

using namespace H5;

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

  int displayNames() const;

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

  void createData(const VectorInt& argdims,int type);
  VectorDouble getDataDoublePartial(int myrank) const;
  int writeDataDoublePartial(int myrank, const VectorDouble& data);

  // Return the size of the data
  // Note that for multi-dim arrays that it gets the total size and not the size of a single row.
  int getSize() const;

  // We now make a proxy class so that we can overload the return type and use a single
  // function to get data whether int or float or double.
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

private:
  hsize_t* _getDims(const DataSpace& dataspace) const;
  void _getOrderSize(const DataSet& dataset,
                     H5T_order_t* order,
                     size_t* size,
                     bool* big_endian) const;

public:
  String _filename;
  String _varname;
};

/**
 * Numeric implementation of our write data function
 * Only accepts numerical values. Integers, floats, or doubles
 * @param data
 */
template<typename T>
void HDF5format::writeData(const T &data)
{
  Exception::dontPrint();
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
      H5File file(H5std_string (_filename.c_str()), H5F_ACC_TRUNC);
      DataSpace dsp = DataSpace(vrank, dims);
      // int
      if (type == (char*) typeid(int).name())
      {
        DataSet dset = file.createDataSet(H5std_string (_varname.c_str()), PredType::STD_I32LE, dsp);
        dset.write(a, PredType::STD_I32LE);
        dset.close();
      }
      // float
      else if (type == (char*) typeid(float).name())
      {
        DataSet dset = file.createDataSet(H5std_string (_varname.c_str()), PredType::IEEE_F32LE, dsp);
        dset.write(a, PredType::IEEE_F32LE);
        dset.close();
      }
      else if (type == (char*) typeid(double).name())
      {
        DataSet dset = file.createDataSet(H5std_string (_varname.c_str()), PredType::IEEE_F64LE, dsp);
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
      H5File file(H5std_string (_filename.c_str()), H5F_ACC_TRUNC);
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
  Exception::dontPrint();
  uint itr = 0;
  uint npts = data.size();
  auto *a = new T[npts];
  char* type = (char*) (typeid(a[0]).name());
  int vrank = 1;
  for (size_t i = 0; i < npts; ++i)
    a[i] = data[i];
  hsize_t dims[1];
  dims[0] = npts;

  while (true)
  {
    // This assumes that the file already exists and will then write to the file
    try
    {
      H5File file(H5std_string (_filename.c_str()), H5F_ACC_RDWR);
      DataSpace dsp = DataSpace(vrank, dims);
      // int
      if (type == (char*) typeid(int).name())
      {
        DataSet dset = file.createDataSet(H5std_string (_varname.c_str()), PredType::STD_I32LE, dsp);
        dset.write(a, PredType::STD_I32LE);
        dset.close();
      }
      // uint
      else if (type == (char*) typeid(uint).name())
      {
        DataSet dset = file.createDataSet(H5std_string (_varname.c_str()), PredType::STD_U32LE, dsp);
        dset.write(a, PredType::STD_U32LE);
        dset.close();
      }
      // float
      else if (type == (char*) typeid(float).name())
      {
        DataSet dset = file.createDataSet(H5std_string (_varname.c_str()), PredType::IEEE_F32LE, dsp);
        dset.write(a, PredType::IEEE_F32LE);
        dset.close();
      }
      // double
      else if (type == (char*) typeid(double).name())
      {
        DataSet dset = file.createDataSet(H5std_string (_varname.c_str()), PredType::IEEE_F64LE, dsp);
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
      H5File file(H5std_string (_filename.c_str()), H5F_ACC_TRUNC);
      file.close();
      // Just some warning that we have gone through this catch
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
void HDF5format::writeData(const std::vector<std::vector<T> > &data)
{
  Exception::dontPrint();
  uint itr = 0;
  uint dim1 = data.size();
  uint dim2 = data[0].size();
  auto a = new T[dim1 * dim2];
  auto md = new T*[dim1];
  for (size_t i = 0; i < dim1; ++i)
    md[i] = a + i * dim2;
  int vrank = 2;
  for (size_t i = 0; i < dim1; ++i)
    for (size_t j = 0; j < dim2; ++j)
      md[i][j] = data[i][j];
  hsize_t dims[2];
  dims[0] = (int) dim1;
  dims[1] = (int) dim2;

  while (true)
  {
    // This assumes that the file already exists and will then write to the file
    try
    {
      H5File file(H5std_string (_filename.c_str()), H5F_ACC_RDWR);
      DataSpace dsp = DataSpace(vrank, dims);
      // int
      if (typeid(T).name() == typeid(int).name())
      {

        DataSet dset = file.createDataSet(H5std_string (_varname.c_str()), PredType::STD_I32LE, dsp);
        dset.write(a, PredType::STD_I32LE);
        dset.close();
      }
      // uint
      else if (typeid(T).name() == typeid(uint).name())
      {
        DataSet dset = file.createDataSet(H5std_string (_varname.c_str()), PredType::STD_U32LE, dsp);
        dset.write(a, PredType::STD_U32LE);
        dset.close();
      }
      // float
      else if (typeid(T).name() == typeid(float).name())
      {
        DataSet dset = file.createDataSet(H5std_string (_varname.c_str()), PredType::IEEE_F32LE, dsp);
        dset.write(a, PredType::IEEE_F32LE);
        dset.close();
      }
      // double
      else if (typeid(T).name() == typeid(double).name())
      {
        DataSet dset = file.createDataSet(H5std_string (_varname.c_str()), PredType::IEEE_F64LE, dsp);
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
      H5File file(H5std_string (_filename.c_str()), H5F_ACC_TRUNC);
      file.close();
      // Just some warning that we have gone through this catch
      itr++;
      if (itr > 3)
      {
        messerr("We've tried too many times in the Int writing sequence");
        break;
      }
    }
  }
}


/******************************************************************************/
/*                                                                            */
/*                          gstlearn Python package                           */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
// Keep sync with PYTHON_PACKAGE_NAME in CMakeLists.txt
%module(directors="1") gstlearn // TODO : configure this using CMake configure_file

%feature(director) AFunction; // TODO : director for AFunction

// Note : Keep order in this file!

// Swig for Python doesn't like our assignment operators 
%ignore *::operator=;

// https://stackoverflow.com/a/26035360/3952924
%import "doc/documentation.i"

//////////////////////////////////////////////////////////////
//    Specific typemaps and fragments for Python language   //
//////////////////////////////////////////////////////////////

%begin %{
// For converting NumPy integers to C++ integers
// https://github.com/swig/swig/issues/888
#define SWIG_PYTHON_CAST_MODE
// For isnan
#include <cmath>
// For isinf
#include <math.h>
// For numeric_limits
#include <limits>
// For cout
#include <iostream>

// Look below (python code) for integer NaN (inan definition)
#if defined(_WIN32) || defined(_WIN64)
  #define NPY_INT_TYPE       NPY_INT32
  #define NPY_INT_OUT_TYPE   int
  #define NPY_INT_NA         std::numeric_limits<NPY_INT_OUT_TYPE>::min()
#else // Linux or MacOS
  #define NPY_INT_TYPE       NPY_INT64
  #define NPY_INT_OUT_TYPE   int64_t
  #define NPY_INT_NA         std::numeric_limits<NPY_INT_OUT_TYPE>::min()
#endif

#define NPY_DOUBLE_NA        std::numeric_limits<double>::quiet_NaN()
#define NPY_FLOAT_NA         static_cast<float>(std::nan(""))
%}

%fragment("ToCpp", "header")
{
  int isNumericVector(PyObject* obj)
  {
    if (PySequence_Check(obj) || PyArray_CheckExact(obj))
    {
      int size = (int)PySequence_Length(obj);
      for (int i = 0; i < size; ++i)
      {
        PyObject* item = PySequence_GetItem(obj, i);
        if (!PyNumber_Check(item))
          return SWIG_TypeError;
      }
      return SWIG_OK;
    }
    return SWIG_TypeError;
  }
  int isStringVector(PyObject* obj)
  {
    if (PySequence_Check(obj) || PyArray_CheckExact(obj))
    {
      int size = (int)PySequence_Length(obj);
      for (int i = 0; i < size; ++i)
      {
        PyObject* item = PySequence_GetItem(obj, i);
        if (!PyUnicode_Check(item))
          return SWIG_TypeError;
      }
      return SWIG_OK;
    }
    return SWIG_TypeError;
  }

  template <typename Type> int convertToCpp(PyObject* obj, Type& value);
  
  template <> int convertToCpp(PyObject* obj, int& value)
  {
    // Test argument
    if (obj == NULL) return SWIG_TypeError;

    long long v = 0; // Biggest integer type whatever the platform
    int myres = SWIG_AsVal_long_SS_long(obj, &v);
    //std::cout << "convertToCpp(int): v=" << v << std::endl;
    if (SWIG_IsOK(myres) || myres == SWIG_OverflowError)
    {
      if (myres == SWIG_OverflowError || v == NPY_INT_NA) // NaN, Inf or out of bound value becomes NA
      {
        myres = SWIG_OK;
        value = getNA<int>();
      }
      else
        myres = SWIG_AsVal_int(obj, &value);
    }
    return myres;
  }
  template <> int convertToCpp(PyObject* obj, double& value)
  {
    // Test argument
    if (obj == NULL) return SWIG_TypeError;

    int myres = SWIG_AsVal_double(obj, &value);
    //std::cout << "convertToCpp(double): value=" << value << std::endl;
    if (SWIG_IsOK(myres))
    {
      if (std::isnan(value) || std::isinf(value))
        value = getNA<double>();
    }
    return myres; 
  }
  template <> int convertToCpp(PyObject* obj, String& value)
  {
    // Test argument
    if (obj == NULL) return SWIG_TypeError;
      
    int myres = SWIG_AsVal_std_string(obj, &value);
    //std::cout << "convertToCpp(String): value=" << value << std::endl;
    // No undefined
    return myres;
  }
  template <> int convertToCpp(PyObject* obj, float& value)
  {
    // Test argument
    if (obj == NULL) return SWIG_TypeError;
    
    int myres = SWIG_AsVal_float(obj, &value);
    //std::cout << "convertToCpp(float): value=" << value << std::endl;
    if (SWIG_IsOK(myres))
    {
      if (std::isnan(value) || std::isinf(value))
        value = getNA<float>();
    }
    return myres; 
  }
  template <> int convertToCpp(PyObject* obj, UChar& value)
  {
    // Test argument
    if (obj == NULL) return SWIG_TypeError;
    
    int v = 0;
    int myres = SWIG_AsVal_int(obj, &v);
    //std::cout << "convertToCpp(UChar): value=" << v << std::endl;
    if (myres == SWIG_OverflowError || 
        v < std::numeric_limits<UChar>::min() ||
        v > std::numeric_limits<UChar>::max()) // Out of bound value is error (no NA for UChar)
    {
      myres = SWIG_OverflowError;
    }
    else if (!SWIG_IsOK(myres) || v == NPY_INT_NA) // NaN or Inf is error (no NA for UChar)
    {
      myres = SWIG_TypeError;
    }
    else
    {
      value = static_cast<UChar>(v);
    }
    return myres;
  }
  template <> int convertToCpp(PyObject* obj, bool& value)
  {
    // Test argument
    if (obj == NULL) return SWIG_TypeError;
    
    int v = 0;
    int myres = SWIG_AsVal_int(obj, &v);
    //std::cout << "convertToCpp(bool): value=" << v << std::endl;
    if (v == 0)
      value = false;
    else
      value = true;
    return myres;
  }
  
  template <typename Vector>
  int vectorToCpp(PyObject* obj, Vector& vec)
  {
    // Type definitions
    using ValueType = typename Vector::value_type;
    vec.clear();

    // Test argument
    if (obj == NULL) return SWIG_TypeError;

    // Conversion
    int myres = SWIG_OK;
    int size = (int)PySequence_Length(obj);
    if (size < 0)
    {
      // Not a sequence (maybe a single value ?)
      ValueType value;
      // Clear Python error indicator
      PyErr_Restore(NULL, NULL, NULL);
      // Try to convert
      myres = convertToCpp(obj, value);
      if (SWIG_IsOK(myres))
        vec.push_back(value);
    }
    else if (size > 0)
    {
      // Real sequence 
      vec.reserve(size);
      for (int i = 0; i < size && SWIG_IsOK(myres); i++)
      {
        PyObject* item = PySequence_GetItem(obj, i);
        ValueType value;
        myres = convertToCpp(item, value);
        if (SWIG_IsOK(myres))
          vec.push_back(value);
      }
    }
    // else size is zero (empty vector)
    return myres;
  }

  template <typename VectorVector>
  int vectorVectorToCpp(PyObject* obj, VectorVector& vvec)
  {
    // Type definitions
    using InputVector = typename VectorVector::value_type;
    vvec.clear();

    // Test argument
    if (obj == NULL) return SWIG_TypeError;

    // Conversion
    int myres = SWIG_OK;
    int size = (int)PySequence_Length(obj);
    if (size < 0)
    {
      // Not a sequence
      InputVector vec;
      // Clear Python error indicator
      PyErr_Restore(NULL, NULL, NULL);
      // Try to convert
      myres = vectorToCpp(obj, vec);
      if (SWIG_IsOK(myres))
        vvec.push_back(vec);
    }
    else if (size > 0)
    {
      for (int i = 0; i < size && SWIG_IsOK(myres); i++)
      {
        PyObject* item = PySequence_GetItem(obj, i);
        InputVector vec;
        myres = vectorToCpp(item, vec);
        if (SWIG_IsOK(myres))
          vvec.push_back(vec);
      }
    }
    // else size is zero (empty vector)
    return myres;
  }
}

%typecheck(SWIG_TYPECHECK_UINT8) UChar {}

// Add numerical vector typecheck typemaps for dispatching functions
%typemap(typecheck, noblock=1, fragment="ToCpp", precedence=SWIG_TYPECHECK_DOUBLE_ARRAY) const VectorInt&,    VectorInt,
                                                                                         const VectorDouble&, VectorDouble,
                                                                                         const VectorUChar&,  VectorUChar,
                                                                                         const VectorBool&,   VectorBool
{
  $1 = SWIG_CheckState(isNumericVector($input));
}

// Add generic vector typecheck typemaps for dispatching functions
%typemap(typecheck, noblock=1, fragment="ToCpp", precedence=SWIG_TYPECHECK_STRING_ARRAY) const VectorString&, VectorString
{
  $1 = SWIG_CheckState(isStringVector($input));
}

// Include numpy interface for creating arrays

%{
  #define SWIG_FILE_WITH_INIT
%}
%include numpy.i
%init %{
  import_array(); // Mandatory for using PyArray_* functions
%}

%fragment("FromCpp", "header")
{
  template <typename Type> NPY_TYPES numpyType();
  template <> NPY_TYPES numpyType<int>()    { return NPY_INT_TYPE; }
  template <> NPY_TYPES numpyType<double>() { return NPY_DOUBLE; }
  template <> NPY_TYPES numpyType<String>() { return NPY_STRING; }
  template <> NPY_TYPES numpyType<float>()  { return NPY_FLOAT; }
  template <> NPY_TYPES numpyType<UChar>()  { return NPY_UBYTE; }
  template <> NPY_TYPES numpyType<bool>()   { return NPY_BOOL; }
  
  template<typename Type> struct TypeHelper;
  template <> struct TypeHelper<int>    { static bool hasFixedSize() { return true; } };
  template <> struct TypeHelper<double> { static bool hasFixedSize() { return true; } };
  template <> struct TypeHelper<String> { static bool hasFixedSize() { return false; } };
  template <> struct TypeHelper<float>  { static bool hasFixedSize() { return true; } };
  template <> struct TypeHelper<UChar>  { static bool hasFixedSize() { return true; } };
  template <> struct TypeHelper<bool>   { static bool hasFixedSize() { return true; } };
  template <typename Type> bool hasFixedSize() { return TypeHelper<Type>::hasFixedSize(); }
  
  template <typename InputType> struct OutTraits;
  template <> struct OutTraits<int>     { using OutputType = NPY_INT_OUT_TYPE; };
  template <> struct OutTraits<double>  { using OutputType = double; };
  template <> struct OutTraits<String>  { using OutputType = const char*; };
  template <> struct OutTraits<float>   { using OutputType = float; };
  template <> struct OutTraits<UChar>   { using OutputType = UChar; };
  template <> struct OutTraits<bool>    { using OutputType = bool; };
  
  template <typename Type> typename OutTraits<Type>::OutputType convertFromCpp(const Type& value);
  template <> NPY_INT_OUT_TYPE convertFromCpp(const int& value)
  {
    //std::cout << "convertFromCpp(int): value=" << value << std::endl;
    NPY_INT_OUT_TYPE vres = static_cast<NPY_INT_OUT_TYPE>(value);
    if (isNA<int>(value))
      vres = NPY_INT_NA;
    return vres;
  }
  template <> double convertFromCpp(const double& value)
  {
    //std::cout << "convertFromCpp(double): value=" << value << std::endl;
    if (isNA<double>(value))
      return NPY_DOUBLE_NA;
    return value;
  }
  template <> const char* convertFromCpp(const String& value)
  {
    //std::cout << "convertFromCpp(String): value=" << value << std::endl;
    return value.c_str(); // No special conversion provided
  }
  template <> float convertFromCpp(const float& value)
  {
    //std::cout << "convertFromCpp(float): value=" << value << std::endl;
    if (isNA<float>(value))
      return NPY_FLOAT_NA;
    return value;
  }
  template <> UChar convertFromCpp(const UChar& value)
  {
    //std::cout << "convertFromCpp(UChar): value=" << value << std::endl;
    return value; // No special conversion provided
  }
  template <> bool convertFromCpp(const bool& value)
  {
    //std::cout << "convertFromCpp(bool): value=" << value << std::endl;
    return value; // No special conversion provided
  }
  
  template <typename Type> PyObject* objectFromCpp(const Type& value);
  template <> PyObject* objectFromCpp(const int& value)
  {
    return PyLong_FromLongLong(convertFromCpp(value));
  }
  template <> PyObject* objectFromCpp(const double& value)
  {
    return PyFloat_FromDouble(convertFromCpp(value));
  }
  template <> PyObject* objectFromCpp(const String& value)
  {
    return PyUnicode_FromString(convertFromCpp(value));
  }
  template <> PyObject* objectFromCpp(const float& value)
  {
    return PyFloat_FromDouble(static_cast<double>(convertFromCpp(value)));
  }
  template <> PyObject* objectFromCpp(const UChar& value)
  {
    return PyLong_FromLong(static_cast<long>(convertFromCpp(value)));
  }
  template <> PyObject* objectFromCpp(const bool& value)
  {
    return PyBool_FromLong(static_cast<long>(convertFromCpp(value)));
  }
  
  template <typename Vector>
  int vectorFromCpp(PyObject** obj, const Vector& vec)
  {
    // Type definitions
    int myres = SWIG_TypeError;
    using SizeType = typename Vector::size_type;
    using InputType = typename Vector::value_type;
    
    // Conversion
    if (hasFixedSize<InputType>()) // Convert to 1D NumPy array
    {
      using OutputType = typename OutTraits<InputType>::OutputType;
      SizeType size = vec.size();
      npy_intp dims[1] = { (npy_intp)size };
      //*obj = PyArray_SimpleNew(1, dims, numpyType<InputType>());
      PyArray_Descr* descr = PyArray_DescrFromType(numpyType<InputType>());
      *obj = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, NULL, NULL, 0, NULL);
      if (*obj != NULL)
      {
        OutputType* array_ptr = (OutputType*) PyArray_DATA((PyArrayObject*)(*obj));
        std::transform(vec.cbegin(), vec.cend(), array_ptr, convertFromCpp<InputType>);
        myres = SWIG_OK;
      }
    }
    else // Convert to a tuple using standard std_vector
    {
      // Test NA values
      auto vec2 = vec.getVector();
      SizeType size = vec2.size();
      for(SizeType i = 0; i < size; i++)
      {
        vec2[i] = convertFromCpp(vec2[i]);
      }
      // Convert to tuple
      *obj = swig::from(vec.getVector());
      myres = (*obj) == NULL ? SWIG_TypeError : SWIG_OK;
    }
    return myres;
  }

  template <typename VectorVector>
  int vectorVectorFromCpp(PyObject** obj, const VectorVector& vec)
  {
    // Type definitions
    int myres = SWIG_TypeError;
    using SizeType = typename VectorVector::size_type;
    using VectorType = typename VectorVector::value_type;
    using InputType = typename VectorType::value_type;
    using OutputType = typename OutTraits<InputType>::OutputType;

    // Check if empty
    if (vec.empty())
    {
      VectorType v; // Create an empty 1D NumPy array
      return vectorFromCpp(obj, v);
    }
    
    // Check the size of sub-vectors
    bool same_size = true;
    SizeType size = vec.at(0).size();
    for(auto v : vec)
    {
      if (same_size)
        same_size = (v.size() == size);
    }
    
    // Conversion
    if (same_size && hasFixedSize<InputType>()) // Convert to a 2D NumPy array
    {
      SizeType size2 = vec.size();
      npy_intp dims[2] = { (npy_intp)size2, (npy_intp)size };
      *obj = PyArray_SimpleNew(2, dims, numpyType<InputType>());
      if (*obj != NULL)
      {
        OutputType* array_ptr = (OutputType*) PyArray_DATA((PyArrayObject*)(*obj));
        for (auto v : vec)
        {
          std::transform(v.cbegin(), v.cend(), array_ptr, convertFromCpp<InputType>);
          array_ptr += size;
        }
        myres = SWIG_OK;
      }
    }
    else // Convert to a list of 1D NumPy array (fixed item size) or tuple (no item size)
    {
      // https://stackoverflow.com/questions/36050713/using-py-buildvalue-to-create-a-list-of-tuples-in-c
      *obj = PyList_New(0);
      if(*obj != NULL)
      {
        myres = SWIG_OK;
        SizeType size2 = vec.size();
        for(SizeType i = 0; i < size2 && SWIG_IsOK(myres); i++)
        {
          PyObject* tuple;
          myres = vectorFromCpp(&tuple, vec.at(i));
          if (SWIG_IsOK(myres))
            myres = PyList_Append(*obj, tuple);
        }
      }
    }
    
    return myres;
  }
}

//////////////////////////////////////////////////////////////
//         C++ library SWIG includes and typemaps           //
//////////////////////////////////////////////////////////////

%include ../swig/swig_inc.i

//////////////////////////////////////////////////////////////
//                C++ library SWIG interface                //
//////////////////////////////////////////////////////////////

%include ../swig/swig_exp.i

//////////////////////////////////////////////////////////////
//                    Add C++ extension                     //
//////////////////////////////////////////////////////////////

%{
#include <algorithm>
#include <string>
#include <sstream>

// Use pybind11 for printing text directly to jupyter notebook output cells
#include <pybind11/pybind11.h>
namespace py = pybind11;
//https://github.com/pybind/pybind11/pull/372/files
using namespace py::literals; // To get access to _a literal below 

void write_f(const char * str)
{
  py::print(str, "end"_a=" ");  // pybind11 print wrapper with no new line
}

void warn_f(const char * str)
{
  py::print(str, "end"_a=" ");  // pybind11 print wrapper with no new line
}

void read_f(const char * str, char * out)
{
  message("%s :",str);
  char line[10000];
  while (fgets(line,10000,stdin) == NULL);
  (void) strcpy(out,line);
  out[strlen(out)-1] = '\0';
}

void exit_f(void)
{
  exit(0);
}
%}

%init %{
  redefine_message(write_f);
  redefine_error(warn_f);
  redefine_read(read_f);
  redefine_exit(exit_f);
%}

// Do not use VectorInt here
%extend VectorNumT<int> {
  std::string __repr__() {  return $self->toString(); }
}
%extend VectorNumT<double> {
  std::string __repr__() {  return $self->toString(); }
}
%extend VectorNumT<float> {
  std::string __repr__() {  return $self->toString(); }
}
%extend VectorNumT<UChar> {
  std::string __repr__() {  return $self->toString(); }
}
%extend VectorNumT<bool> {
  std::string __repr__() {  return $self->toString(); }
}
%extend VectorT<int> {
  std::string __repr__() {  return $self->toString(); }
}
%extend AMatrix {
  std::string __repr__() {  return $self->toString(); }
}
%extend AMatrixSquare {
  std::string __repr__() {  return $self->toString(); }
}
%extend MatrixInt {
  std::string __repr__() {  return $self->toString(); }
}
%extend MatrixRectangular {
  std::string __repr__() {  return $self->toString(); }
}
%extend MatrixSquareDiagonal {
  std::string __repr__() {  return $self->toString(); }
}
%extend MatrixSquareDiagonalCst {
  std::string __repr__() {  return $self->toString(); }
}
%extend MatrixSquareGeneral {
  std::string __repr__() {  return $self->toString(); }
}
%extend MatrixSquareSymmetric {
  std::string __repr__() {  return $self->toString(); }
}
%extend ASpace {
  std::string __repr__() {  return $self->toString(); }
}
%extend SpacePoint {
  std::string __repr__() {  return $self->toString(); }
}
%extend SpaceTarget {
  std::string __repr__() {  return $self->toString(); }
}
%extend Db {
  std::string __repr__() {  return $self->toString(); }
}
%extend DbGrid {
  std::string __repr__() {  return $self->toString(); }
}
%extend Vario {
  std::string __repr__() {  return $self->toString(); }
}
%extend DirParam {
  std::string __repr__() {  return $self->toString(); }
}
%extend VarioParam {
  std::string __repr__() {  return $self->toString(); }
}
%extend Option_AutoFit {
  std::string __repr__() {  return $self->toString();  }
}
%extend Option_VarioFit {
  std::string __repr__() {  return $self->toString();  }
}
%extend Constraints {
  std::string __repr__() {  return $self->toString();  }
}
%extend Polygons {
  std::string __repr__() {  return $self->toString(); }
}
%extend Model {
  std::string __repr__() {  return $self->toString(); }
}
%extend CovAniso {
  std::string __repr__() {  return $self->toString(); }
}
%extend SimuBooleanParam {
  std::string __repr__() {  return $self->toString(); }
}
%extend SimuSphericalParam {
  std::string __repr__() {  return $self->toString(); }
}
%extend SimuPartitionParam {
  std::string __repr__() {  return $self->toString(); }
}
%extend SimuSubstitutionParam {
  std::string __repr__() {  return $self->toString(); }
}
%extend SimuFFTParam {
  std::string __repr__() {  return $self->toString(); }
}
%extend SimuRefineParam {
  std::string __repr__() {  return $self->toString(); }
}
%extend Array {
  std::string __repr__() {  return $self->toString(); }
}
%extend BImage {
  std::string __repr__() {  return $self->toString(); }
}
%extend ProjMatrix {
  std::string __repr__() {  return $self->toString(); }
}
%extend FracEnviron {
  std::string __repr__() {  return $self->toString(); }
}
%extend PolyLine2D {
  std::string __repr__() {  return $self->toString(); }
}
%extend Table {
  std::string __repr__() {  return $self->toString(); }
}
%extend Selectivity {
  std::string __repr__() {  return $self->toString(); }
}
%extend Limits {
  std::string __repr__() {  return $self->toString(); }
}
%extend NeighMoving {
  std::string __repr__() {  return $self->toString(); }
}
%extend NeighUnique {
  std::string __repr__() {  return $self->toString(); }
}
%extend NeighImage {
  std::string __repr__() {  return $self->toString(); }
}
%extend NeighBench {
  std::string __repr__() {  return $self->toString(); }
}
%extend NeighCell {
  std::string __repr__() {  return $self->toString(); }
}
%extend MeshEStandard {
  std::string __repr__() {  return $self->toString(); }
}
%extend MeshETurbo {
  std::string __repr__() {  return $self->toString(); }
}
%extend MeshSpherical {
  std::string __repr__() {  return $self->toString(); }
}
%extend CSVformat {
  std::string __repr__() {  return $self->toString(); }
}
%extend PCA {
  std::string __repr__() {  return $self->toString(); }
}
%extend Rule {
  std::string __repr__() {  return $self->toString(); }
}
%extend PPMT {
  std::string __repr__() {  return $self->toString(); }
}
%extend AnamHermite {
  std::string __repr__() {  return $self->toString(); }
}

//////////////////////////////////////////////////////////////
//       Add target language additional features below      //
//////////////////////////////////////////////////////////////

%pythoncode %{

import gstlearn as gl
import numpy as np
import scipy.sparse as sc
import pandas as pd
import os
import sys

## Version and authors
from gstlearn.version import __version__
from gstlearn.version import __author__

## Numpy default integer type could be different according the OS (int64 or int32)
## https://stackoverflow.com/questions/36278590/numpy-array-dtype-is-coming-as-int32-by-default-in-a-windows-10-64-bit-machine
## Remind that there is no standard NaN for integers in Python
## We cannot use std::numeric_limits<int>::quiet_NaN() in C++
## But we can use the minimum signed integer value as follow

import os
## Integer NaN custom value
if os.name == 'nt':         # Windows
  inan = -2147483648
else:                       # Others
  inan = -sys.maxsize - 1

## isNaN custom function
def isNaN(value):
  if (type(value).__module__ == np.__name__): # Numpy type
    if (np.dtype(value) == 'intc' or np.dtype(value) == 'int64' or np.dtype(value) == 'int32'):
      return value == gl.inan
  else:
    if (type(value).__name__ == 'int'):
      return value == gl.inan
  return np.isnan(value)


## Add operator [] to VectorXXX R class [1-based index] ##
## ---------------------------------------------------- ##

def setitem(self, idx, item):
  if idx < 0 or idx >= self.length():
    raise IndexError("Index out or range")
  self.setAt(idx,item)
  
def getitem(self, idx):
  if idx < 0 or idx >= self.length():
    raise IndexError("Index out or range")
  return self.getAt(idx)

setattr(gl.VectorDouble,       "__getitem__", getitem)
setattr(gl.VectorDouble,       "__setitem__", setitem)
setattr(gl.VectorInt,          "__getitem__", getitem)
setattr(gl.VectorInt,          "__setitem__", setitem)
setattr(gl.VectorString,       "__getitem__", getitem)
setattr(gl.VectorString,       "__setitem__", setitem)
setattr(gl.VectorFloat,        "__getitem__", getitem)
setattr(gl.VectorFloat,        "__setitem__", setitem)
setattr(gl.VectorUChar,        "__getitem__", getitem)
setattr(gl.VectorUChar,        "__setitem__", setitem)
setattr(gl.VectorBool,         "__getitem__", getitem)
setattr(gl.VectorBool,         "__setitem__", setitem)
setattr(gl.VectorVectorDouble, "__getitem__", getitem)
setattr(gl.VectorVectorDouble, "__setitem__", setitem)
setattr(gl.VectorVectorInt,    "__getitem__", getitem)
setattr(gl.VectorVectorInt,    "__setitem__", setitem)
setattr(gl.VectorVectorFloat,  "__getitem__", getitem)
setattr(gl.VectorVectorFloat,  "__setitem__", setitem)


## Override operator [] for the Db class ##
## ------------------------------------- ##
# Thanks to Nicolas Desassis:

def is_list_type(mylist, types):
    """Check if an input is an iterable (tuple, list or numpy array) containing
       elements of only a given type"""
    all_type = True
    if not(isinstance(mylist,(tuple, list, np.ndarray))):
        all_type = False
    i = 0
    while all_type and i<len(mylist):
        if not(isinstance(mylist[i], types)):
               all_type = False
        i += 1
    return all_type

def check_nrows(db, nrows):
    """Check if a number of rows matches with the number of samples of a Db, and returns the flag
    for useSel (whether it matches the number of active samples or the total number of samples)"""
    if nrows == db.getActiveSampleNumber() :
        useSel = True
    elif nrows == db.getSampleNumber() or db.getSampleNumber()==0:
        useSel = False
    else:
        if db.getActiveSampleNumber() != db.getSampleNumber():
            raise ValueError("Error of dimension. Your number of lines ("+str(nrows)+") has to be equal to " +
                str(db.getActiveSampleNumber()) + " or " + str(db.getSampleNumber()))
        else :
            raise ValueError("Error of dimension. Your number of lines ("+str(nrows)+") has to be equal to " +
                  str(db.getActiveSampleNumber()))
    return useSel

def findColumnNames(self, columns):
    """Extract names of columns from Db, given different possible types of arguments: 
        names, indices, or locator"""
    if isinstance(columns, str) or is_list_type(columns, (str, np.str_)): #get variable(s) by name
        names = self.identifyNames(np.atleast_1d(columns))
    
    elif isinstance(columns, (int, np.int_)):
        names = self.getNameByColIdx(columns)
    
    elif isinstance(columns, slice):
        Nmax = self.getColumnNumber()
        names = []
        for i in range(Nmax)[columns]:
            names.append(self.getNameByColIdx(i))

    elif is_list_type(columns, (int, np.int_)):
        names = []
        Nfields = self.getColumnNumber()
        for i in columns:
            if i >= Nfields:
                print("Warning: the index {} is out of bounds with {}, this index is ignored".format(i,Nfields))
            else:
                names.append(self.getNameByColIdx(int(i)))
        
    else:
        raise ValueError("Argument for columns of wrong type: {}".format(type(columns)))
        
    return np.atleast_1d(names)

def has_row_selection(self, arg):
    """Check if the argument given contains a rows selection [rows,columns], 
    or only column selection [columns].
    If the argument is a tuple of length 2 and its first element is a valid argument
    for indexing rows, then the function returns True."""
    valid_row_indexing = False
    if isinstance(arg, tuple) and len(arg)==2:
        array_test = np.zeros(getNrows(self))
        try: # test if first element of tuple is a valid argument for indexing rows. If yes, then we assume it is the argument for rows.
            array_test[arg[0],]
            valid_row_indexing = True
        except IndexError:
            valid_row_indexing = False
    return valid_row_indexing

def getNrows(self, useSel=None):
    """ get number of rows of the Db when using or not a selection"""
    if useSel is None:
        useSel = self.useSel
    nrows = self.getSampleNumber(useSel)
    return nrows

def getdbitem(self,arg):
    """
    Extract data from a Db. Use Db[arg]

    Parameters
    ----------
    arg is (rows, columns) or columns
    rows : (optional) int, list of int, or slice. Which rows (samples) to extract.
    columns: str or list of str. Names of the variables to extract.
             int, list of int, or slice. Indices of the variables to extract.
             gstlearn.ELoc. Locator of the variables to extract.
             
    Returns
    -------
    numpy.ndarray
        2D array of shape (nrows, nvars) of the extracted data.
    Examples
    --------
    db["var"] or db[:,"var"] extracts the variable named "var"
    db[5:10,(2,3)] extracts the rows 5 to 9 of the variables of index 2 and 3 (array of shape (5,2))
    db[gl.ELoc.Z] extracts all the variables located with Z.
    """
    nrows = getNrows(self)
        
    selec_rows = has_row_selection(self, arg)
    if selec_rows:
        rows = arg[0]
        columns = arg[1]
    else:
        rows = slice(None,None,None) # extract all rows
        columns = arg
    
    # extract columns
    ColNames = findColumnNames(self, columns)
    nbvar = len(ColNames)
    temp = np.array(self.getColumns(ColNames, self.useSel))
    temp = temp.reshape([nbvar,nrows]).T
            
    # extract rows
    temp = temp[rows,]
    if len(temp.shape) == 2:
      if temp.shape[1] == 1:
        return temp[:,0]
        
    return temp
        
# This function will add a set of vectors (as a numpy array) to a db. 
# If some of the names exist, the corresponding variables will be replaced 
# and not added.

def setdbitem(self,name,tab):
    
    # analyze input arguments
    selec_rows = has_row_selection(self, name)   
    if selec_rows:
        rows = name[0]
        columns = name[1]
    else:
        columns = name
    
    # find existing column names
    arr_columns = np.atleast_1d(columns)
    ColNames = findColumnNames(self, columns) #existing names
    
    # analyze input table
    if isinstance(tab, (float, np.floating, int, np.integer, bool, np.bool_)):
        nrows = getNrows(self)
        nvars = len(ColNames) # this means we will only modify existing columns, not create ones
        tab = np.ones((nrows, nvars))*tab
        if selec_rows:
            tab = np.atleast_2d(tab[rows,:])
    else:
        tab = np.copy(np.float64(tab))
        if len(tab.shape) == 1:
            tab = np.atleast_2d(tab).T
        nrows, nvars = tab.shape
    
    # create list of column names to modify and/or create
    if len(ColNames) == nvars: # modify existing variables only
        new_names = ColNames
     
    elif len(arr_columns) == nvars and is_list_type(arr_columns, (str,np.str_)):
        new_names = arr_columns
        
    elif isinstance(columns, (str,np.str_)) and nvars > 1 and len(ColNames)==0: # create new variables from a unique name
        new_names = gl.generateMultipleNames(columns, nvars)
        
    else:
        raise ValueError("Wrong type or length of input ({0}): the input should correspond"
                         " either to a number of existing variables ({1}) equal to the"
                         " number of columns of the table (nvar={2}), or should be a name or "
                         "list of names of length nvar={2} in order to create new variables.".format(columns, len(ColNames), nvars))
            
    # loop on the column names to modify/create each column
    for i,name in enumerate(new_names):
        # check if existing name
        ExistingNames = findColumnNames(self, name)
        if len(ExistingNames) > 1:
            raise ValueError("There is more than one variable name corresponding to '{}' "
                             "in the Db: {}".format(name, ExistingNames))
            
        if selec_rows:
            useSel = self.useSel
            if len(ExistingNames) == 0: # create new variable
                nrows_tot = getNrows(self, useSel)
                tab_i = np.empty(nrows_tot)
                tab_i.fill(np.nan) # NaNs outside of target rows
            elif len(ExistingNames) == 1: # modify existing variable
                tab_i = self[name]
                
            tab_i = np.squeeze(tab_i)
            tab_i[rows,] = tab[:,i]
            
        else:
            useSel = check_nrows(self,nrows)
            tab_i = np.empty(nrows)
            tab_i[:] = tab[:,i]
        
        tab_i[np.isnan(tab_i)] = np.nan
        VectD = np.double(tab_i)
        self.setColumn(VectD, name, gl.ELoc.UNKNOWN, 0, useSel)
        
    return

setattr(gl.Db,"useSel",False)    
    
setattr(gl.Db,"__getitem__",getdbitem)
setattr(gl.Db,"__setitem__",setdbitem)

def Db_toTL(self, flagLocate=False):
  dat = pd.DataFrame(self.getAllColumns().reshape(-1,self.getSampleNumber()).T, 
    columns = self.getAllNames())
    
#  dat = pd.DataFrame(self.getAllColumns().reshape(-1, self.getSampleNumber()), 
#    columns = self.getAllNames())

  if flagLocate:
    for j,i in enumerate(self.getAllNames()):
      dat[i].locator = self.getLocators()[j] 
  return dat

setattr(gl.Db, "toTL", Db_toTL)

def matrix_toTL(self):
  if self.isSparse():
    Acs = self.getCsToTriplet().toTL()
    return Acs
  else:
    Anp = np.array(self.getValues()).reshape(self.getNRows(),self.getNCols())
    return Anp
  return

setattr(gl.MatrixRectangular, "toTL", matrix_toTL)
setattr(gl.MatrixSquareDiagonal, "toTL", matrix_toTL)
setattr(gl.MatrixSquareDiagonalCst, "toTL", matrix_toTL)
setattr(gl.MatrixSquareGeneral, "toTL", matrix_toTL)
setattr(gl.MatrixSquareSymmetric, "toTL", matrix_toTL)
setattr(gl.MatrixSparse, "toTL", matrix_toTL)

# TODO : Replace Triplet_toTL by MatrixSparse_toTL
def Triplet_toTL(self):
  Acs = sc.csc_matrix((np.array(self.values), 
                      (np.array(self.rows), np.array(self.cols))),
                         shape=(self.nrows, self.ncols))
  return Acs

setattr(gl.Triplet, "toTL", Triplet_toTL)

def table_toTL(self):
  Anp = np.array(self.getValues()).reshape(self.getNRows(),self.getNCols())
  #colnames(mat) <- tab$getColumnNames()
  #rownames(mat) <- tab$getRowNames()
  return Anp

setattr(gl.Table, "toTL", table_toTL)
%}

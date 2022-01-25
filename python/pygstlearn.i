// https://blog.mbedded.ninja/programming/languages/python/python-swig-bindings-from-cplusplus/
%module(directors="1") pygstlearn

// https://stackoverflow.com/a/26035360/3952924

%import "doc/documentation.i"

// gstlearn C++ library must be installed (Keep Order !!!!)
%include ../swig/swig_inc.i
%include ../swig/swig_exp.i

// For suppressing SWIG warning due to -keyword option
#pragma SWIG nowarn=511
#pragma SWIG nowarn=506

/// TODO, include numpy.i ?

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
//  PyObject *item = Py_BuildValue(str);
//  PyObject_Print(item, stdout, 0);
}

void warn_f(const char * str)
{
  py::print(str, "end"_a=" ");  // pybind11 print wrapper with no new line
//  PyObject *item = Py_BuildValue(sstr.str().c_str());
//  PyObject_Print(item, stderr, 0);
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

%extend SpacePoint {
  std::string __repr__() {  return $self->toString(); }
}
%extend std::vector<double> {
  std::string __repr__() {  return ut_vector_string(*$self); }
}
%extend std::vector<int> {
  std::string __repr__() {  return ut_ivector_string(*$self); }
}
%extend Db {
  std::string __repr__() {  return $self->toString(); }
}
%extend Vario {
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
%extend Model {
  std::string __repr__() {  return $self->toString(); }
}
%extend NeighUnique {
  std::string __repr__() {  return $self->toString();  }
}
%extend NeighBench {
  std::string __repr__() {  return $self->toString();  }
}
%extend NeighMoving {
  std::string __repr__() {  return $self->toString();  }
}
%extend NeighImage {
  std::string __repr__() {  return $self->toString();  }
}
%extend AnamHermite {
  std::string __repr__() {   return $self->toString();  }
}
%extend Rule {
  std::string __repr__() {  return $self->toString();  }
}
%extend Interval {
  std::string __repr__() {  return $self->toString();  }
}
%extend AMesh {
  std::string __repr__() {  return $self->toString();  }
}
%extend MeshETurbo {
  std::string __repr__() {  return $self->toString();  }
}
%extend NoStatArray {
  std::string __repr__() {  return $self->toString();  }
}
%extend Rotation {
  std::string __repr__() {  return $self->toString();  }
}
%extend ProjMatrix {
  std::string __repr__() {  return $self->toString();  }
}
%extend Polygons {
  std::string __repr__() {  return $self->toString();  }
};



%pythoncode %{
# Override operator [] for the Db class
# Thanks to Nicolas Desassis:
import gstlearn as gl
import numpy as np


# In[2]:


def getitem(self,name):
    """
    Use Db[name]. Returns the arrays of the fields corresponding to input 'name'.

    Parameters
    ----------
    name : (str or list of str) Name or list of names to be extracted;
           (locator: instance of gstlearn.Eloc) Locator to be extracted;
           (indices) Used as for a numpy array: indices of fields to be extracted,
                     it also works with 2D indices to extract specific samples

    Returns
    -------
    numpy.ndarray
        Array of shape (nsamples, nvar) of the extracted data. Array is squeezed to 1D if nvar=1.

    """
    if self.useSel:
        nrows = self.getActiveSampleNumber()
    else:
        nrows = self.getSampleNumber()
        
    if isinstance(name, gl.ELoc):
        temp = np.array(self.getFieldsByLocator(name, self.useSel))
        nbvar = self.getLocatorNumber(name)
        if nbvar == 0:
            temp = np.array([])

    elif isinstance(np.atleast_1d(name)[0], (str,np.str_)):
        name = np.atleast_1d(name)
        nbvar = len(self.getNames(name))
        temp = np.array(self.getFields(name, self.useSel))
        
    else: # indices of arrays
        array = np.array(self.getFields(useSel=self.useSel))
        nbvar = self.getFieldNumber()
        array = np.reshape(array, (nbvar,nrows))
        temp = array[name]
        return temp.T
        
    temp = temp.reshape([nbvar,nrows]).T
    temp[temp == gl.TEST] = None
    return temp

# This function will add a set of vectors (as a numpy array) to a db. If some of the names exist, the
# corresponding variables will be replaced and not added.

def setitem(self,name,tab):
    
    if len(tab.shape) == 1 :
        temptab = np.atleast_2d(tab).T
    else :
        temptab = tab
    
    nrows = tab.shape[0]
    
    if nrows == self.getActiveSampleNumber() :
        useSel = True
    elif nrows == self.getSampleNumber() or self.getSampleNumber()==0:
        useSel = False
    else :
        if self.getActiveSampleNumber() != self.getSampleNumber():
            raise ValueError("Error of dimension. Your number of lines has to be equal to " +
                str(self.getActiveSampleNumber()) + " or " + str(self.getSampleNumber()))
        else :
            raise ValueError("Error of dimension. Your number of lines has to be equal to " +
                  str(self.getActiveSampleNumber()))
            
    if isinstance(name, (str, np.str_)) :
     	names = self.getNames([name])
    
     	if len(names) == 0 :
         	names = [name]
        
     	if len(names) == 1 and temptab.shape[1] > 1:
         	names = gl.generateMultipleNames(name,temptab.shape[1])
    elif isinstance(name, (list,tuple,np.ndarray)) :
     	names = name
    else :
        raise TypeError("Type of name should be in: 'str', 'numpy.str_', 'list', 'tuple', 'numpy.ndarray'")
    
    vectD = gl.VectorDouble()

    for j in range(temptab.shape[1]):
        vectD.resize(0)
        for i in range(nrows):       
            u = np.double(temptab[i,j])
            if u is None : 
                u = gl.TEST
            vectD.push_back(u)
            
        self.setField(vectD,names[j],useSel)
        
setattr(gl.Db,"useSel",False)    
    
setattr(gl.Db,"__getitem__",getitem)

setattr(gl.Db,"__setitem__",setitem)

		   

		   
%}

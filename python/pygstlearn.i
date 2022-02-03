// https://blog.mbedded.ninja/programming/languages/python/python-swig-bindings-from-cplusplus/
// No need of %module keyword when building using cmake UseSWIG
// TODO: restore directors feature
//%module(directors="1") gstlearn

// https://stackoverflow.com/a/26035360/3952924
%import "doc/documentation.i"

// Include C++ library SWIG interface (Keep Order !!!!)
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
%extend DbGrid {
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
            raise ValueError("Error of dimension. Your number of lines has to be equal to " +
                str(db.getActiveSampleNumber()) + " or " + str(db.getSampleNumber()))
        else :
            raise ValueError("Error of dimension. Your number of lines has to be equal to " +
                  str(db.getActiveSampleNumber()))
    return useSel

def findColumnNames(self, columns):
    """Extract names of columns from Db, given different possible types of arguments"""
    if isinstance(columns, str) or is_list_type(columns, (str, np.str_)): #get variable(s) by name
        names = self.getNames(np.atleast_1d(columns))
    
    elif isinstance(columns, gl.ELoc): #get variable(s) by locator
        names = self.getNamesByLocator(columns)
    
    elif is_list_type(columns, gl.ELoc):
        if not(len(columns)) == 1:
            raise ValueError("The input for columns should not be a list of several Locators")
        names = self.getNamesByLocator(columns[0])
     
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
                print(f"Warning: the index {i} is out of bounds with {Nfields}, this index is ignored")
            else:
                names.append(self.getNameByColIdx(int(i)))
        
    else:
        raise ValueError(f"Argument for columns of wrong type: {type(columns)}")
        
    return np.atleast_1d(names)

def getitem(self,arg):
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
    if self.useSel:
        nrows = self.getActiveSampleNumber()
    else:
        nrows = self.getSampleNumber()
    
    if isinstance(arg, tuple) and isinstance(arg[0], (int,slice)): # 2D (rows, columns)
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
    temp = temp[rows]
        
    temp[temp == gl.TEST] = np.nan
    return temp


# This function will add a set of vectors (as a numpy array) to a db. 
# If some of the names exist, the corresponding variables will be replaced 
# and not added.

def setitem(self,name,tab):
    tab = tab.astype(np.float64)
    if len(tab.shape) == 1 :
        tab = np.atleast_2d(tab).T
    nrows, nvars = tab.shape
    tab[np.isnan(tab)] = gl.TEST
    
    if isinstance(name, tuple) and isinstance(name[0], (int,slice)): # 2D (rows, columns)
        selec_rows = True
        rows = name[0]
        columns = name[1]
    else:
        selec_rows = False
        columns = name
            
    arr_columns = np.atleast_1d(columns)
    ColNames = findColumnNames(self, columns) #existing names
        
    if len(ColNames) == nvars: # modify existing variables only
        new_names = ColNames
     
    elif len(arr_columns) == nvars and is_list_type(arr_columns, (str,np.str_)):
        new_names = arr_columns
        
    elif isinstance(columns, (str,np.str_)) and nvars > 1 and len(ColNames)==0: # create new variables from a unique name
        new_names = gl.generateMultipleNames(columns, nvars)
        
    else:
        raise ValueError(f"Wrong type or length of input ({columns}): the input should correspond"
                         f" either to a number of existing variables ({len(ColNames)}) equal to the"
                         f" number of columns of the table (nvar={nvars}), or should be a name or "
                         f"list of names of length nvar={nvars} in order to create new variables.")
    
    for i,name in enumerate(new_names):
        
        ExistingNames = findColumnNames(self, name)
        if len(ExistingNames) > 1:
            raise ValueError(f"There is more than one variable name corresponding to '{name}' "
                             f"in the Db: {ExistingNames}")
            
        if selec_rows:
            useSel = self.useSel
            if len(ExistingNames) == 0: # create new variable
                if useSel:
                    nrows_tot = self.getActiveSampleNumber()
                else:
                    nrows_tot = self.getSampleNumber()
                tab_i = np.ones(nrows_tot)*gl.TEST # NaNs outside of target rows
            elif len(ExistingNames) == 1: # modify existing variable
                tab_i = self[name]
                
            tab_i = np.squeeze(tab_i)
            tab_i[rows] = tab[:,i]
            
        else:
            tab_i = tab[:,i]
            useSel = check_nrows(self,nrows)
            
        VectD = np.double(tab_i)
        self.setColumn(VectD, name, useSel)
        
    return

setattr(gl.Db,"useSel",False)    
    
setattr(gl.Db,"__getitem__",getitem)

setattr(gl.Db,"__setitem__",setitem)

%}

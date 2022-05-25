// No need of %module keyword for python when building using cmake UseSWIG
// TODO: restore directors feature (inheritance in target language)
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
    temp[temp == gl.TEST] = np.nan
    return temp
        
# This function will add a set of vectors (as a numpy array) to a db. 
# If some of the names exist, the corresponding variables will be replaced 
# and not added.

def setitem(self,name,tab):
    
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
        if len(tab.shape) == 1 :
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
                tab_i = np.ones(nrows_tot)*gl.TEST # NaNs outside of target rows
            elif len(ExistingNames) == 1: # modify existing variable
                tab_i = self[name]
                
            tab_i = np.squeeze(tab_i)
            tab_i[rows,] = tab[:,i]
            
        else:
            useSel = check_nrows(self,nrows)
            tab_i = np.empty(nrows)
            tab_i[:] = tab[:,i]
        
        tab_i[np.isnan(tab_i)] = gl.TEST    
        VectD = np.double(tab_i)
        self.setColumn(VectD, name, useSel)
        
    return

setattr(gl.Db,"useSel",False)    
    
setattr(gl.Db,"__getitem__",getitem)

setattr(gl.Db,"__setitem__",setitem)

# Add plot functions as methods of the class
import gstlearn.plot as gp

setattr(gl.Db,"plot", gp.point)
setattr(gl.Db,"plot_correlation", gp.correlation)
setattr(gl.Db,"plot_hist", gp.hist)

setattr(gl.DbGrid,"plot", gp.grid)
setattr(gl.DbGrid,"plot_grids", gp.grids)
setattr(gl.DbGrid,"plot_point", gp.point)
# plot_correlation and plot_hist are already inherited from the parent class Db

setattr(gl.Vario,"plot", gp.vario)
setattr(gl.Vario,"plot_varioElem", gp.varioElem)
setattr(gl.Vario,"plot_varioDir", gp.varioDir)
setattr(gl.Vario,"plot_varmod", gp.varmod)

setattr(gl.Model,"plot", gp.model)

setattr(gl.Rule,"plot", gp.rule)

setattr(gl.Table,"plot", gp.table)

setattr(gl.Polygons,"plot", gp.polygon)
%}

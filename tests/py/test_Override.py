#!/usr/bin/env python
# coding: utf-8

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
            raise ValueError(f"Error of dimension. Your number of lines ({nrows}) has to be equal to " +
                str(db.getActiveSampleNumber()) + " or " + str(db.getSampleNumber()))
        else :
            raise ValueError(f"Error of dimension. Your number of lines ({nrows}) has to be equal to " +
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
                print(f"Warning: the index {i} is out of bounds with {Nfields}, this index is ignored")
            else:
                names.append(self.getNameByColIdx(int(i)))
        
    else:
        raise ValueError(f"Argument for columns of wrong type: {type(columns)}")
        
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
    if useSel:
        nrows = self.getActiveSampleNumber()
    else:
        nrows = self.getSampleNumber()
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
        raise ValueError(f"Wrong type or length of input ({columns}): the input should correspond"
                         f" either to a number of existing variables ({len(ColNames)}) equal to the"
                         f" number of columns of the table (nvar={nvars}), or should be a name or "
                         f"list of names of length nvar={nvars} in order to create new variables.")
            
    # loop on the column names to modify/create each column
    for i,name in enumerate(new_names):
        # check if existing name
        ExistingNames = findColumnNames(self, name)
        if len(ExistingNames) > 1:
            raise ValueError(f"There is more than one variable name corresponding to '{name}' "
                             f"in the Db: {ExistingNames}")
            
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
    
setattr(gl.Db,"__getitem__",getitem)

setattr(gl.Db,"__setitem__",setitem)

# # Example

# In[3]:

a = gl.DbGrid.create([2,2],[1.,1.])
a

# ## Create  a new variable

# In[4]:


np.random.seed(123)
x = np.random.normal(size=4)
a["var1"] = x
u = np.sum(np.abs(a["var1"][:,0] - np.array([-1.0856306 ,  0.99734545,  0.2829785 , -1.50629471])))
if u>1e-7:
    raise Exception("Problem in simple get or add new variable")

# ## Get a variable

# ## Replace an existing variable

# In[5]:


a["var1"] = 2 * x
u = np.sum(np.abs(a["var1"][:,0] - 2*np.array([-1.0856306 ,  0.99734545,  0.2829785 , -1.50629471])))
if u>1e-7:
    raise Exception("Problem in set existing variable")


# ## Use regexp for get

# In[6]:


a["var2"] = 3*x


# In[7]:


u1 = np.sum(np.abs(a["var*"][:,0] - 2*np.array([-1.0856306 ,  0.99734545,  0.2829785 , -1.50629471])))
u2 = np.sum(np.abs(a["var*"][:,1] - 3*np.array([-1.0856306 ,  0.99734545,  0.2829785 , -1.50629471])))
if u1>1e-7 or u2>1e-7:
    raise Exception("Problem in get with regexp")


# ## Use regexp for replacing several variables

# In[8]:


a["var*"]=a["var*"]>0


# In[9]:


u = np.sum(np.abs(a["var*"]-np.array([[0., 0.],
       [1., 1.],
       [1., 1.],
       [0., 0.]])))
if not u==0 :
    raise Exception("Problem in set existing variable with regexp")


# # Create several variables with one name

# In[10]:


a["newvar"] = np.random.normal(size = (4,3))


# In[11]:


s = 0
for i in a.getAllNames():
    if i =='newvar-1':
        s+=1
    if i =='newvar-2':
        s+=1
    if i == 'newvar-3':
        s+=1
if not s==3:
    raise Exception("Problem to create several new variable with one name (names)")


# In[12]:


u = np.sum(np.abs(a["newvar*"]-np.array([[-0.57860025,  1.65143654, -2.42667924],
       [-0.42891263,  1.26593626, -0.8667404 ],
       [-0.67888615, -0.09470897,  1.49138963],
       [-0.638902  , -0.44398196, -0.43435128]])))
if u > 1e-7:
    raise Exception("Problem to create several new variable with one name (values)")


# In[13]:


v = a["newvar*"]
v[0,0]=None


# In[14]:


a["newvar*"] = v


# In[15]:


if not np.isnan(a["newvar*"][0,0]):
    raise Exception("Problem with nan value conversion")


# ## Add tab to a newly created db (and provide several names)

# In[16]:

a = gl.Db()
a[["var1","var2"]] = np.random.normal(size=(12,2))
a

# In[17]:

u = np.sum(np.abs(a[["var1","var2"]]-np.array([[ 2.20593008,  2.18678609],
       [ 1.0040539 ,  0.3861864 ],
       [ 0.73736858,  1.49073203],
       [-0.93583387,  1.17582904],
       [-1.25388067, -0.6377515 ],
       [ 0.9071052 , -1.4286807 ],
       [-0.14006872, -0.8617549 ],
       [-0.25561937, -2.79858911],
       [-1.7715331 , -0.69987723],
       [ 0.92746243, -0.17363568],
       [ 0.00284592,  0.68822271],
       [-0.87953634,  0.28362732]])))
if u>1e-7:
    raise Exception("Problem with creation of several new variables from several names")

# In[18]:
    
a = gl.Db()
a["var"] = np.random.normal(size=(12,5))
a.setLocators(("var-1","var-2"), gl.ELoc.Z)

u = a[::2,gl.ELoc.Z] - np.array([[-0.80536652, -1.72766949],
       [-1.29408532, -1.03878821],
       [-0.77270871,  0.79486267],
       [ 0.46843912, -0.83115498],
       [ 1.25523737, -0.68886898],
       [ 1.15020554, -1.26735205]])
if np.any(u>1e-7):
    raise Exception("Problem with get from Locator and specific rows")

u = a[2:4] - a[[2,3]]
if np.any(u>1e-7):
    raise Exception("Problem with get indices or slice which give different results")

try:
    a["var-1"] = np.zeros((12,3))
except ValueError:
    None
else:
    raise Exception("This should raise an error (number of existing variables"
                    " does not match the number of columns in the table)")
        
a[::2,"var-1"] = np.ones(12)[::2]
a[::2,"var-6"] = np.ones(12)[::2]
a[1::2,gl.ELoc.Z] = np.zeros((6,2))

try:
    a[gl.ELoc.Z] = np.random.rand(12,3)
except ValueError as ve:
    if ve.__str__()[:29] != "Wrong type or length of input":
        raise Exception("Wrong error is returned")
else:
    raise Exception("This should raise an error (number of existing variables"
                    " does not match the number of columns in the table)")


try:
    a[gl.ELoc.X] = np.random.rand(12)
except ValueError as ve:
    if ve.__str__()[:29] != "Wrong type or length of input":
        raise Exception("Wrong error is returned")
else:
    raise Exception("This should raise an error (setting an array to a non-existing locator)")

a[1] = np.random.rand(12)
a[0:3,(1,2)] = np.random.rand(3,2)

a[(1,2)] # row 1 variable 2
a[[1,2]] # variables 1 and 2, all rows
#a[[1,1]] # que faire ?

a[:10,2:4] = np.ones((10,2))


# In[19]

print("Everything is ok")


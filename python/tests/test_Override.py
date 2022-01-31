#!/usr/bin/env python
# coding: utf-8

# In[1]:


import gstlearn as gl
import numpy as np


# In[2]:
def is_tuple_str(mytuple):
    """Check is a tuple (or iterable) contains only strings"""
    all_str = True
    if not(isinstance(mytuple,(tuple, list, np.ndarray))):
        all_str = False
    i = 0
    while all_str and i<len(mytuple):
        if not(isinstance(mytuple[i], (str, np.str_))):
               all_str = False
        i += 1
    return all_str

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
        Array of shape (nsamples, nvar) of the extracted data.

    """
    
    if self.useSel:
        nrows = self.getActiveSampleNumber()
    else:
        nrows = self.getSampleNumber()
    
    if isinstance(name, tuple) and isinstance(name[0], (int,slice)): # 2D (rows, columns)
        rows = name[0]
        columns = name[1]
    else:
        rows = slice(None,None,None) # extract all rows
        columns = name
    
    # extract columns
    if isinstance(columns, str) or is_tuple_str(columns): #get variable(s) by name
        names = np.atleast_1d(columns)
        nbvar = len(self.getNames(columns))
        temp = np.array(self.getFields(names, self.useSel))
    
    elif isinstance(columns, gl.ELoc): #get variable(s) by locator
        temp = np.array(self.getFieldsByLocator(columns, self.useSel))
        nbvar = self.getLocatorNumber(columns)
        
    else: #indices or slice (column indices)
        array = np.array(self.getAllFields(useSel=self.useSel))
        nbvar_tot = self.getFieldNumber()
        array = np.reshape(array, (nbvar_tot,nrows))
        temp = np.atleast_2d(array[columns])
        nbvar = temp.shape[0]
        
    temp = temp.reshape([nbvar,nrows]).T
            
    # extract rows
    temp = temp[rows]
        
    temp[temp == gl.TEST] = None
    return temp



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

print("Everything is ok")


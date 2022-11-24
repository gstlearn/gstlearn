# # Name Convention for Db in gstlearn

# ## Preamble

# This tutorial gives answers to the frequently asked question regarding the Name Convention used for variables in a Data Base (Db) of gstlearn 

import numpy as np
import gstlearn as gl
import os
import sys

# ## Prepare the Environment
# 
# This paragraph defines the Space Dimension for the whole notebook. It also set the name of the Container (and a Prefix) used if Objects are saved as Neutral Files.

ndim = 2
gl.defineDefaultSpace(gl.ESpaceType.RN,ndim)

gl.ASerializable.setContainerName(True)
gl.ASerializable.setPrefixName("DbTest-");

# The following object will enable having a complete view of the column / attribute manipulation. It will be used later in the notebook.

dbfmt = gl.DbStringFormat()
dbfmt.setFlags(flag_locator=True)

# ## Creating a data file
# 
# A Data Base is created for experimentation. It is constructed as a regular Grid (named **grid**). The variable *nech* will contain the number of samples within *grid*. The number of meshes is voluntarily limited. The mesh is square with dimension 1. The origin (lower left corner) is set to (10,20) in order to be able to distinguish coordinates along first and second axes.

grid = gl.DbGrid.create([5,5], [1,1], [10,20])
nech = grid.getSampleNumber()
print("Number of sample =",nech)
grid

# The data base contains 3 fields, created automatically and respectively called *rank*, *x1* and *x2*. Note that the last two fields are considered as coordinates (locator *x*).

# ## Names

# We now add one new field (named *first*) where values are generated randomly (uniform drawn between 0 and 1). Note that, when adding this new field, a value is returned which corresponds to the number of the newly created *attribute*.
# 
# **Important remark: all numerical variables used to identify a field within a Db are considered as indices, i.e. they are numbered starting from 0**

tab = gl.VH.simulateUniform(nech)
iatt1 = grid.addColumns(tab,"first")
print("Attribute corresponding to 'first' =",iatt1)

# We can double-check the attribute information by visiting the current contents of the *grid* Db. We check that the field *first* is the fourth (i.e. attribute #3).

grid

# Let us add a series (3) of fields created simulateneously. They are filled with a constant value equal to 5. We also define a locator assigned to all the newly created variables: they will be considered as data variable (locator = *z*). 
# Note the returned value: it corresponds to the attribute number assigned to the first new variable.

iatt2 = grid.addColumnsByConstant(3,5.,"second",gl.ELoc.Z)
print("Attribute corresponding to the first variable named 'second-x' =",iatt2)
grid

# Note that the newly created fields are automatically named using the provided string (*second*) as a radix: the variables names are "second-1", "second-2" and "second-3".

# Let us now envisage renaming the variable *second-2* into *first*.

grid.setName("second-2","first")
grid


# As the name *first* already exists, the field has been renamed to *first.1* instead.
# 
# We now wish to rename the field *second-3* into *first*.

grid.setName("second-3","first")
grid

# The automatic renaming procedure has been applied (adding ".1") iteratively until names are all different: the field is now called *first.1.1".

# Now that we have demonstrated the uniqueness of the names, are there are ways to designate a field?
# For the next demonstrations, we first recall the current status of the current Db.
# 
# In order to make the next paragrah more demonstrative, we change the contents of several fields

grid.setColumn(gl.VH.simulateUniform(nech),"second-1")
grid.setColumn(gl.VH.simulateUniform(nech),"first.1")
grid.setColumn(gl.VH.simulateUniform(nech),"first.1.1")

grid

# ### By Name

# As an example, we access to the field named *first.1. For short, only the four first values are systematically printed. 

grid.getColumn("first.1")[0:4]

# ### By Column Index

# We recall that the index numbering starts from 0. Therefore field *first-1*  corresponds to the index 5.

grid.getColumnByColIdx(5)[0:4]

# ### By Attribute Index

grid.getColumnByUID(5)[0:4]

# ### By Locator

# We note that the target variable corresponds to the locator *z2* which is the second one (index 1) or the Z-locator type.

grid.getColumnByLocator(gl.ELoc.Z,1)[0:4]

# ## Difference between Column and Attribute

# We need to recall the *attribute*  value returned when adding the fields:
# - *iatt1* (3) when adding the field named *first*
# - *iatt2* (4) when adding the series of 3 fields (originally named after the radix *second*)
# 
# To better understand, we need to ask for the display of the data base with a specific option which describes the current status of the attributes, either unsorted or through an order driven by the locator

grid.display(dbfmt)

# We can see that the 7 existing fields currently correspond to the 7 first columns of the Data Base *grid*. The second display gives the indices of the locators in use (*x* and *z*) and the indices of the attributes corresponding to the ranks of the items for each locator type.

# Things become more interesting if a field is deleted. To avoid any ambiguity, the field is designated by its name (say *x1*)

grid

grid.deleteColumn("x1")
grid

# The previous printout shows the current contents of the data base where the field *x1* has been suppressed.
# Note an important feature of the *locator* notion. For a given locator type (say *x* for coordinates), the locator type is unique and sorted continuously starting from 1.
# Therefore, when we suppressed the variable *x1* (which corresponded to the locator type *x* and locator rank *1*), the variable *x2* is modified: its name and locator type are not changed but the locator rank is update from *2* to *1*.

# We now look at the attributes internal management

grid.display(dbfmt)

# We can see that the list of attributes has not been reduced: the maximum number of positions is still equal to 7. Instead, the rank of the attribute which corresponded to *x1* is now set to -1, to signify that the column is actually missing. The display sorted by locator does not need any additional explanation.
# 
# Let us now retrieve the information of variable *first.1*  as we did before. We start by addressing the variable by name.

grid.getColumn("first.1")[0:4]


# We can similarly address it by its column index (the column has moved to rank 5)

grid.getColumnByColIdx(4)[0:4]

# The magic of the *attribute* notion is that it can still be used **unchanged**

grid.getColumnByUID(5)[0:4]

# Obviously, trying to read the field which corresponds to the field *x1* (that has just been deleted) returns an empty vector.

grid.getColumnByUID(1)

# ## Remark on Space Dimension

# It might be considered as surprising to see that *grid* is considered as a 2-D Grid while there is only **one** coordinate field (locator *x*). In order to avoid any missunderstanding, let us recall this important fact.
# 
# The data base *grid* is organized as a grid and for that sake, it contains a descrption of the grid organization. This organization is used to elaborate the coordinates (for example when calling *getCoordinate()* method). The coordinate vectors must only be considered as decoration: they will not be used in any internal operation.
# 
# As an example this makes particular sense here as the contents of the variable *x2*, despite its locator rank *1* (i.e. index 0) actually contains the **second** coordinate of the samples, as demonstrated in the next line

grid.getColumnByLocator(gl.ELoc.X, 0)

# Note that at any time, the coordinate vectors can be regenerated. To avoid confusion, the newly generated coordinate fields are named using the radix "X" (uppercase). This feature is obviously only available in the case of a grid

grid.generateCoordinates("X")

grid.getColumnByLocator(gl.ELoc.X, 0)

# Similarly, we can generate a field containing the sample rank (similar as the information contained in the Field #1). Here again, we generate a new field containing this rank information: in order to avoid confusion, the new variable is called *RANK* (uppercase). Note that this field does not have any locator attached.

grid.generateRank("RANK")
grid

# ## Conclusion

# As a conclusion:
# 
# - the variables can be used **safely** when designating them by their **name**
# - the variables can be used easily when addressing them using the locator notion (type and index)
# - the use of (column) index is always valid. This index must be defined precisely when using the variable (it must be updated in case of addition or deletion of other variables)
# - the use of attribute is clever... but it must be used by expert who understands the process. It allows using  fix values, independently of the management of other fields
# 
# We also recall that all numbering refer to indices (0 based numbering). This is the case for *(column) index* as well as *locator index* per locator type.

# # Naming Convention

# Usually, the methods which are designed to add variables in a Db can define some characteristics of the newly created variables using a standardized solution, known as the *NamingConvention* facility. 
# 
# To illustrate this facility, we consider the function *toIndicator* of the class *Limits*. This function stores the indicators of a set of classes (defined in *limits*) in the Db file. We concentrate on the variable *X-1*  which has outcomes varying from 10 to 14 (as demonstrated next).

gl.OptCst.define(gl.ECst.NTCOL,-1)
dbfmt = gl.DbStringFormat()
dbfmt.setFlags(flag_array=True)
grid.display(dbfmt)

# We define a set of adjacent classes with bounds ranging from 10 to 14.

limits = gl.Limits([10, 11, 12, 13, 14])
limits

# We first run the function which uses the defaulted *NamingConvention* where the prefix is set to "Indicator". As a result, 4 new variables are added whose names are composed by concatenating *prefix*, *variable name* and the qualifer *Class* followed by the rank of the class.
# 
# The resulting variables are of the form **Indicator.X-1.Class.1**

limits.toIndicator(grid,"X-1",namconv=gl.NamingConvention("Indicator"))
grid

# In this second trial, we decide to suppress the prefix as well as the variable name

limits.toIndicator(grid,"X-1",namconv=gl.NamingConvention("",False))
grid

import gstlearn as gl
import numpy as np
import pandas as pd

#
# This test is meant to check the manipulation of the new DataBase under Python
#

######################################################
# We intantiate the DbData object and load information

data = gl.DbData()
data.printContents("Checking that the DbData is empty")

gl.mestitle(0, "Adding Columns of various types")

print("- Column VD of type Double with role X, filled with [1.0, 2.0, 3.0]")
data.addColumnD("VD", [1.0, 2.0, 3.0], gl.RoleID(gl.ERole.X))

print("- Column VI of type Int, filled with [5, 6, 7]")
data.addColumnI("VI", [5, 6, 7], gl.RoleID(gl.ERole.Z))

print("- Column VS of type String, filled with ['foo', 'bar', 'baz']")
data.addColumnS("VS", ["foo", "bar", "baz"], gl.RoleID(gl.ERole.Z))

print("- Column VB of type Bool, filled with [True, False, True]")
data.addColumnB("VB", [True, False, True])

print("- Column VDS of type Double, with 5 versions and role F, filled with 3.0")
data.addColumnEmptyD("VDS", 0, 5, gl.RoleID(gl.ERole.F), 3.0)

print("- Column VIS of type Int, with 5 versions and role Z, filled with 1")
data.addColumnEmptyI("VIS", 0, 5, gl.RoleID(gl.ERole.Z), 1)
data.printContents()

############################################
# Special section for Categorical variable #

gl.mestitle(0, "Creating the Dictionary for the categorical variable")
dictionary = gl.Dictionary()
keys = [1, 3, 7]
labels = ["red", "blue", "green"]
for i in range(len(keys)):
    dictionary.addCategory(list(keys)[i], list(labels)[i])
ncategory = dictionary.getNCategories()
dictionary.display()
print("Creating the VectorCategory for the categorical variable")
vc = gl.VectorCategory(3, dictionary)
for i in range(3):
    j = i % ncategory
    vc.setCategory(i, keys[j])

data.addColumnC("VC", vc)
data.printContents()

##############################
# Checking the Neutral Files #

gl.mestitle(1, "Saving and recovering a DbData from a Neutral File")
data.printContents("Before Saving in a Neutral File")
print("VD: ", data.getColumnD("VD"))
print("VI: ", data.getColumnI("VI"))
print("VS: ", data.getColumnS("VS"))
print("VB: ", data.getColumnB("VB"))
print("VDS: ", data.getColumnD("VDS"))
print("VIS: ", data.getColumnI("VIS"))
print("VC: ", data.getColumnC("VC"))
data.dumpToNF("test_DataBase.NF")
vc2 = data.getColumnC("VC")

data2 = gl.DbData.createFromNF("test_DataBase.NF", True)
print("VD: ", data2.getColumnD("VD"))
print("VI: ", data2.getColumnI("VI"))
print("VS: ", data2.getColumnS("VS"))
print("VB: ", data2.getColumnB("VB"))
print("VDS: ", data2.getColumnD("VDS"))
print("VIS: ", data2.getColumnI("VIS"))
print("VC: ", data2.getColumnC("VC"))
data2.dumpToNF("test_DataBase.NF")

########################################################
# Particular test for the case of Categorical variable #
gl.mestitle(1, "Checking the use of a Categorical Variable")
# Reading the category of the first sample and store it in 'cat'
cat = data.getValueC("VC", 0)
# Setting the category 'cat' to the second sample
data.setValueC("VC", 1, cat)
# Cheking the result on the dump of the whole vector
print("VC: ", data.getColumnC("VC"))

####################################################################
# In this part, we check the different ways to construct the ColID #

gl.mestitle(0, "Various ways to specify a column (e.g. for retreiving its name):")
print("- by Name:", data.getName("VD"))
print("- by Name and Version:", data.getName(("VDS", 1)))
print("- by Index:", data.getName(4))
print("- by Index and Version:", data.getName((4, 1)))
print("- by ColID:", data.getName(gl.ColID(1)))
print("- by RoleID:", data.getName(gl.RoleID(gl.ERole.F)))
print("- by RoleID and Version:", data.getName((gl.RoleID(gl.ERole.F), 0)))
print("- by Role:", data.getName(gl.ERole.Z))
print("- by Role and Index:", data.getName((gl.ERole.Z, 0)))
# Next line is the Unique way to define a ColID based on a Role and specifying the Index and the Version
print("- by Role and Index and Version:", data.getName((gl.RoleID(gl.ERole.Z, 0), 2)))

#####################################################################
# In this part, we check the different ways to enquiry the DataBase #

gl.mestitle(0, "Checking the presence of the columns in the DbData:")
print("- Column called 'VD': ", data.hasColumn("VD"))
print("- Column called 'VI': ", data.hasColumn("VI"))
print("- Column called 'VS': ", data.hasColumn("VS"))
print("- Column called 'VB': ", data.hasColumn("VB"))
print("- Column called 'VDS': ", data.hasColumn("VDS"))

gl.mestitle(0, "Checking that we can retrieve the whole contents of any column:")
print(" - Column 0:", data.getColumnD(0))
print(" - Column 1:", data.getColumnI(1))
print(" - Column 2:", data.getColumnS(2))
print(" - Column 3:", data.getColumnB(3))
print(" - Column 4:", data.getColumnD(4))

isample = 2
gl.mestitle(0, "Checking that we can manipulate one Target Sample from any column:")
print("Retrieving a Target Sample (", isample, ") from any column")
print("- From column 0 (specified by its role): ", data.getValueD(gl.ERole.X, isample))
print("- From column 1 (specified by its name): ", data.getValueI("VI", isample))
print("- From column 2 (specified by its index): ", data.getValueS(2, isample))

gl.message("\nModifying the value at Target: Col 0: 4.0, Col 1: 8, Col 2: 'mystring'\n")
data.setValueD(0, isample, 4)
data.setValueI(1, isample, 8)
data.setValueS(2, isample, "mystring")

gl.message("\nChecking the new value of the Target Sample\n")
print("Contents of element #", isample, " for Column 0: ", data.getValueD(0, isample))
print("Contents of element #", isample, " for Column 1: ", data.getValueI(1, isample))
print("Contents of element #", isample, " for Column 2: ", data.getValueS(2, isample))

#########################################################
# In this part, we check the volontary misuse of DbData #

gl.mestitle(0, "Misuses of the DbData")
data.printContents("- Initial situation")

print("Adding Column with existing name (VI) but different type (Double)")
data.addColumnD("VI", [10.0, 11.0, 12.0])

print("Adding Column with existing Role (X) but non consecutive index (10)")
data.addColumnD("VD", [101.0, 102.0, 103.0], gl.RoleID(gl.ERole.X, 10))

print("Adding a Column with an already existing Role (X) and existing Index (0)")
data.addColumnD("VD", [101.0, 102.0, 103.0], gl.RoleID(gl.ERole.X, 0))
data.printContents("- Final situation")

##############################
# Removing Columns of DbData #

gl.mestitle(0, "Deleting a Column (VI.1)")
data.printContents("- Initial situation")
data.deleteColumn("VI.1")
data.printContents("- Final situation")

##########################################################
# Playing with multiple versions in a Column of a DbData #

gl.mestitle(0, "Testing the MultiVersion feature of the DbData")
data.addColumnEmptyD("VDS", nsamples=3, nversion=4, valinit=0.0, forbidNA=True)
data.printContents("After adding the 'VDS' column (with 4 versions)")

##################################
# Additional inquiries on DbData #

gl.mestitle(0, "Various inquiries on the DbData")
print("- Number of columns: ", data.getNColumns())
print("- Number of samples: ", data.getNSamples())
print("- Number of Versions in column 'VDS.1': ", data.getNVersions("VDS.1"))

############################
# Testing errors on DbData #

gl.mestitle(0, "Erroneous operations on the DbData")
data.printContents("- Initial situation")

print("\nTrying to delete a non existing column (VI.22)")
data.deleteColumn("VI.22")

print("\nTrying to add the Column (VD_bad_samples) with a different number of samples")
data.addColumnD("VD_bad_samples", [201.0, 202.0, 203.0, 204.0])

print("\nTrying to set a value in a non existing column (VI.22)")
data.setValueD("VI.22", 0, 52.0)

print("\nTrying to modify the contents of the Column 'VD' (double) by string")
data.setValueS("VD", 0, "Invalid String")

print("\nTrying to use a wrong version (5) in Column 'VDS' (4 versions)")
colid = gl.ColID("VDS", 5)
data.getValueD(("VDS", 5), 0)

print("\nTrying to set a value to NA in Column 'VDS.1' (forbidNA = True)")
data.setValueD("VDS.1", 0, gl.TEST)

data.printContents("- Final situation")

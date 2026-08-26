#
# This test is meant to check the manipulation of the new DataBase under Python
#

# Loading the package

suppressWarnings(suppressMessages(library(gstlearn)))

######################################################
# We intantiate the DbData object and load information

data = DbData()
invisible(data$printContents("Checking that the DbData is empty"))

invisible(mestitle(0, "Adding Columns of various types"))
writeLines("- Column 0 of type Double with role X, filled with [1.0, 2.0, 3.0]")
invisible(data$addColumnD("VD", c(1.0, 2.0, 3.0), RoleID(ERole_X())))

writeLines("- Column 1 of type Int, filled with [5, 6, 7]")
invisible(data$addColumnI("VI", c(5, 6, 7), RoleID(ERole_Z())))

writeLines("- Column 2 of type String, filled with ['foo', 'bar', 'baz']")
invisible(data$addColumnS("VS", c("foo", "bar", "baz"), RoleID(ERole_Z())))

writeLines("- Column 3 of type Bool, filled with [True, False, True]")
invisible(data$addColumnB("VB", c(TRUE, FALSE, TRUE)))

writeLines("- Column 4 of type Double, with 5 versions and role F, filled with 3.0")
invisible(data$addColumnEmptyD("VDS", 0, 5, RoleID(ERole_F()), 3.0))
invisible(data$printContents())

# ####################################################################
# # In this part, we check the different ways to construct the ColID #

invisible(mestitle(0, "Various ways to specify a column (e.g. for retreiving its name):"))
writeLines(paste0("- by Name : ", data$getName("VD")))
writeLines(paste0("- by Name and Version : ", data$getName(list("VDS", 1))))
writeLines(paste0("- by Index :  ", data$getName(4)))
writeLines(paste0("- by Index and Version : ", data$getName(list(4, 1))))

invisible(mestitle(0, "Various ways to specify a column (e.g. for retreiving its name):"))
writeLines(paste0("- by Name : ", data$getName("VD")))
writeLines(paste0("- by Name and Version : ", data$getName(list("VDS", 1))))
writeLines(paste0("- by Index : ", data$getName(4)))
writeLines(paste0("- by Index and Version : ", data$getName(list(4, 1))))

writeLines(paste0("- by RoleID and Version : ", data$getName(list(RoleID(ERole_F()), 0))))
writeLines(paste0("- by Role : ", data$getName(ERole_Z()))) # Pb
writeLines(paste0("- by Role and Index : ", data$getName(list(ERole_Z(), 0)))) # Pb
# Next line is the Unique way to define a ColID based on a Role and specifying the Index and the Version
writeLines(paste0("- by Role and Index and Version : ", data$getName(list(RoleID(ERole_Z(), 0), 2))))

#####################################################################
# In this part, we check the different ways to enquiry the DataBase #

invisible(mestitle(0, "Checking the presence of the columns in the DbData:"))
writeLines(paste0("- Column called 'VD': ", data$hasColumn("VD")))
writeLines(paste0("- Column called 'VI': ", data$hasColumn("VI")))
writeLines(paste0("- Column called 'VS': ", data$hasColumn("VS")))
writeLines(paste0("- Column called 'VB': ", data$hasColumn("VB")))
writeLines(paste0("- Column called 'VDS': ", data$hasColumn("VDS")))

invisible(mestitle(0, "Checking that we can retrieve the whole contents of any column:"))
writeLines(paste0(" - Column 0 : ", paste(data$getColumnD(0), collapse = ", ")))
writeLines(paste0(" - Column 1 : ", paste(data$getColumnI(1), collapse = ", ")))
writeLines(paste0(" - Column 2 : ", paste(data$getColumnS(2), collapse = ", ")))
writeLines(paste0(" - Column 3 : ", paste(data$getColumnB(3), collapse = ", ")))
writeLines(paste0(" - Column 4 : ", paste(data$getColumnD(4), collapse = ", ")))

isample = 2
invisible(mestitle(0, "Checking that we can manipulate one Target Sample from any column:"))
writeLines(paste0("Retrieving a Target Sample (", isample, ") from any column"))
writeLines(paste0("- From column 0 (specified by its role): ", data$getValueD(ERole_X(), isample)))
writeLines(paste0("- From column 1 (specified by its name): ", data$getValueI("VI", isample)))
writeLines(paste0("- From column 2 (specified by its index): ", data$getValueS(2, isample)))

invisible(message("\nModifying the value at Target: Col 0: 4.0, Col 1: 8, Col 2: 'mystring'\n"))
invisible(data$setValueD(0, isample, 4))
invisible(data$setValueI(1, isample, 8))
invisible(data$setValueS(2, isample, "mystring"))

invisible(message("\nChecking the new value of the Target Sample\n"))
writeLines(paste0("Contents of element #", isample, " for Column 0: ", data$getValueD(0, isample)))
writeLines(paste0("Contents of element #", isample, " for Column 1: ", data$getValueI(1, isample)))
writeLines(paste0("Contents of element #", isample, " for Column 2: ", data$getValueS(2, isample)))

#########################################################
# In this part, we check the volontary misuse of DbData #

invisible(mestitle(0, "Misuses of the DbData"))
invisible(data$printContents("- Initial situation"))

writeLines("Adding Column with existing name (VI) but different type (Double)")
invisible(data$addColumnD("VI", c(10.0, 11.0, 12.0)))

writeLines("Adding Column with existing Role (X) but non consecutive index (10)")
invisible(data$addColumnD("VD", c(101.0, 102.0, 103.0), RoleID(ERole_X(), 10)))

writeLines("Adding a Column with an already existing Role (X) and existing Index (0)")
invisible(data$addColumnD("VD", c(101.0, 102.0, 103.0), RoleID(ERole_X(), 0)))
invisible(data$printContents("- Final situation"))

##############################
# Removing Columns of DbData #

invisible(mestitle(0, "Deleting a Column (VI.1)"))
invisible(data$printContents("- Initial situation"))
invisible(data$deleteColumn("VI.1"))
invisible(data$printContents("- Final situation"))

##########################################################
# Playing with multiple versions in a Column of a DbData #

invisible(mestitle(0, "Testing the MultiVersion feature of the DbData"))
invisible(data$addColumnEmptyD("VDS", nsamples = 3, nversion = 4, valinit = 0.0, forbidNA = TRUE))
invisible(data$printContents("After adding the 'VDS' column (with 4 versions)"))

##################################
# Additional inquiries on DbData #

invisible(mestitle(0, "Various inquiries on the DbData"))
writeLines(paste0("- Number of columns: ", data$getNColumns()))
writeLines(paste0("- Number of samples: ", data$getNSamples()))
writeLines(paste0("- Number of Versions in column 'VDS.1': ", data$getNVersions("VDS.1")))

############################
# Testing errors on DbData #

invisible(mestitle(0, "Erroneous operations on the DbData"))
invisible(data$printContents("- Initial situation"))

writeLines("\nTrying to delete a non existing column (VI.22)")
invisible(data$deleteColumn("VI.22"))

writeLines("\nTrying to add the Column (VD_bad_samples) with a different number of samples")
invisible(data$addColumnD("VD_bad_samples", c(201.0, 202.0, 203.0, 204.0)))

writeLines("\nTrying to set a value in a non existing column (VI.22)")
invisible(data$setValueD("VI.22", 0, 52.0))

writeLines("\nTrying to modify the contents of the Column 'VD' (double) by string")
invisible(data$setValueS("VD", 0, "Invalid String"))

writeLines("\nTrying to use a wrong version (5) in Column 'VDS' (4 versions)")
data$getValueD(list("VDS", 5), 0)

writeLines("\nTrying to set a value to NA in Column 'VDS.1' (forbidNA = TRUE)")
invisible(data$setValueD("VDS.1", 0, NA))

invisible(data$printContents("- Final situation"))

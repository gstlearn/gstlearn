# Loading the package

suppressWarnings(suppressMessages(library(gstlearn)))

# Create a sample Db
nfacies = 5
db = Db_createFillRandom(100, ncode = nfacies)

# Transform categorical variable into indicators
limits = Limits(seq(-0.5, 4.5), seq(0.5, 5.5))
limits$display()
# Convert characters string in a NamingConvention
err = limits$toIndicator(db, "code", namconv = "MyVar")

# Display the resulting database
db$display()
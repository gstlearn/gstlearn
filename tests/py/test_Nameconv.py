import gstlearn as gl
import numpy as np

# Create a sample Db
nfacies = 5
db = gl.Db.createFillRandom(100, ncode = nfacies)

# Transform categorical variable into indicators
limits = gl.Limits(np.arange(-0.5,5.5), np.arange(0.5,6.5))
limits.display()
# Convert characters string in a NamingConvention
err = limits.toIndicator(db, "code", namconv = "MyVar")

# Display the resulting database
db.display()
# Benchmark for Db in gstlearn

import gstlearn as gl
import numpy as np
import time

ndim = 2
gl.defineDefaultSpace(gl.ESpaceType.RN,ndim)

gl.ASerializable.setContainerName(True)
gl.ASerializable.setPrefixName("BenchDb-");

# Creating a data file
# 
# Construct a database (point cloud) from random positions

nech = 10000000
idx = [x for x in range(nech)]
db = gl.Db.createFillRandom(nech, 2, 0)
print(f"Number of sample = {nech}")

# Test execution time for Db::getIncrements function 
start_time = time.time()
inc = db.getIncrements(idx, idx)
print(inc.shape)
print("#NO_DIFF# --- %s seconds ---" % round((time.time() - start_time), 4))
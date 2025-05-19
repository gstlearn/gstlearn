import gstlearn as gl


dims = [10, 10, 10]
db = gl.DbGrid.create(nx=dims)
db.dumpToH5("dbgrid.h5")

db2 = gl.DbGrid.createFromH5("dbgrid.h5")

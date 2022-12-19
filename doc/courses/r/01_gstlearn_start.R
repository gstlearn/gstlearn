##
## Fichier test
##
library(gstlearn)
library(ggplot2)

db = DbGrid_create(nx=c(100,100))
model = Model_createFromParam(type = ECov_CUBIC(), range = 30)
err = simtub(NULL, db, model, nbtuba=1000)
plot.grid(db,name="Simu",title="Check is successful!")

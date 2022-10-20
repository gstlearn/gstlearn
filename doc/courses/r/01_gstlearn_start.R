##
## Fichier test
##
library(gstlearn)

db = DbGrid_create(nx=c(100,100))
model = Model_createFromParam(type = ECov_CUBIC(), range = 30)
err = simtub(NULL, db, model, NULL, nbsimu=1, seed=13126, nbtuba = 1000)
# on aurait aime faire ... mais ca ne fonctionne pas!
# err = simtub(, db, model, nbtuba=1000)
plot.grid(db,name="Simu",title="Check is successful!")

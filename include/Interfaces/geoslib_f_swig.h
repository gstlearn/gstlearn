#pragma once

#include "gstlearn_export.hpp"
#include "Interfaces/Database.hpp"
#include "Interfaces/interface_d.hpp"
#include "Neigh/ENeigh.hpp"
#include "Basic/Vector.hpp"
#include "geoslib_enum.h"

class Db;
class Model;

GSTLEARN_EXPORT void migrate_grid_to_point2(const Database &Db_grid,
                                            Database &db_point,
                                            const String &name,
                                            int ldmax = 2,
                                            VectorDouble dmax = VectorDouble());

GSTLEARN_EXPORT void migrate_point_to_grid2(const Database &Db_point,
                                            Database &db_grid,
                                            const String &name,
                                            int ldmax = 2,
                                            VectorDouble dmax = VectorDouble());

GSTLEARN_EXPORT void migrate_grid_to_grid2(const Database &db_grid_in,
                                           Database &db_grid_out,
                                           const String &name,
                                           int ldmax = 2,
                                           VectorDouble dmax = VectorDouble());

//Model* model_auto2(const VarioExp& vario,const VectorInt& lst_model, int uc, int idir, int ivar);
//VectorDouble model_evaluate2(Model* model, const VarioExp& vario, int idir, int ivar);
//
GSTLEARN_EXPORT void kriging2(const Database &dbin,
                              Database &dbout,
                              Model *model,
                              ANeighParam *neighparam);

GSTLEARN_EXPORT void mes_error(ES error);
GSTLEARN_EXPORT VectorDouble affiche(Db *db);
GSTLEARN_EXPORT void my_db_print(Db *db);

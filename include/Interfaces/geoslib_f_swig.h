#ifndef GEOSLIB_F_SWIG_H
#define GEOSLIB_F_SWIG_H

#include "Interfaces/Database.hpp"
#include "Interfaces/interface_d.hpp"
#include "Interfaces/VarioExp.hpp"

void migrate_grid_to_point2(const Database& Db_grid, Database& db_point, const std::string& name, int ldmax = 2, VectorDouble dmax = VectorDouble());

void migrate_point_to_grid2(const Database& Db_point, Database& db_grid, const std::string& name, int ldmax = 2, VectorDouble dmax = VectorDouble());

void migrate_grid_to_grid2(const Database& db_grid_in, Database& db_grid_out, const std::string &name, int ldmax = 2, VectorDouble dmax = VectorDouble());

Model* model_auto2(const VarioExp& vario,const VectorInt& lst_model, int uc, int idir, int ivar);

VectorDouble model_evaluate2(Model* model, const VarioExp& vario, int idir, int ivar);

Neigh* neigh_unique(int ndim);
Neigh* neigh_moving(int ndim, int flag_sector, int flag_rotation, int nmini, int nmaxi, int nsect, int nsmax, int radius, VectorDouble Rotation);

Neigh* my_neigh_init(int ndim, int type, int flag_xvalid, int flag_sector,
                    int flag_aniso, int flag_rotation, int flag_continuous,
                    int nmini, int nmaxi, int nsect, int nsmax, int skip,
                    double width, double radius, double dist_count,
                    VectorDouble nbgh_radius, VectorDouble nbgh_rotmat,
                    VectorDouble nbgh_image);

void kriging2(const Database & dbin, Database &dbout, Model* model, Neigh* neigh);

void mes_error(ES error);
VectorDouble affiche(Db* db);
void my_db_print(Db*  db);
#endif

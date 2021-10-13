#include "Interfaces/geoslib_f_swig.h"

#include <iostream>
#include <numeric>

#include "geoslib_d.h"
#include "geoslib_enum.h"
#include "Interfaces/VariableDouble.hpp"

VectorDouble range(int n)
{
  VectorDouble res;
  int i = 0;
  while(i < n)
  {
    res.push_back(i);
    i++;
  }
  return(res);
}

void migrate_grid_to_point2(const Database& db_grid, Database& db_point, const std::string& name, int ldmax, VectorDouble dmax)
{
  std::vector<double> vec(db_point.getNSamples());
  std::string new_name = name + std::string("_migrate");
  Db* db= db_grid.toGeoslib();

  Db*  dbpointgeos = db_point.toGeoslib(); 
  VectorInt cols(1);
  cols[0] = db_grid.nameIdentify(name);
  migrateByAttribute(db, dbpointgeos, cols, ldmax, dmax, false, false);
  db_point.fromGeoslib(dbpointgeos);
}

void migrate_point_to_grid2(const Database &db_point, Database& db_grid, const std::string& name, int ldmax, VectorDouble dmax)
{
  std::vector<double> vec(db_grid.getGridSize());
  std::string new_name = name + std::string("_migrate");
  Db* dbp = db_point.toGeoslib();

  Db* dbg_geos = db_grid.toGeoslib();
  VectorInt cols(1);
  cols[0] = db_point.nameIdentify(name);
  migrateByAttribute(dbp, dbg_geos, cols, ldmax, dmax, false, false);
  db_grid.fromGeoslib(dbg_geos);
}

void migrate_grid_to_grid2(const Database &dbgin, Database& dbgout, const std::string& name, int ldmax, VectorDouble dmax)
{
  std::vector<double> vec(dbgout.getGridSize());
  std::string new_name = name + std::string("_migrate");
  Db * dbin = dbgin.toGeoslib();

  Db* dbg_geos = dbgout.toGeoslib();
  VectorInt cols(1);
  cols[0] = dbgin.nameIdentify(name);;
  migrateByAttribute(dbin, dbg_geos, cols, ldmax, dmax, false, false);
  dbgout.fromGeoslib(dbg_geos);
}

void mes_error(ES error)
{
  if (error == ES_NAME)
    std::cout << "Error : Duplicate name in Database" << std::endl;
  else if (error == ES_INDEX_OOR)
    std::cout << "Error : Index Out of Range" << std::endl;
  else if (error == ES_SIZE_VAR)
    std::cout << "Error : Var of different size"<<std::endl;
}

/**********************************************************************
 ** Create object Model* According to a Variogram
 ** 
 **********************************************************************/
//:WARNING: Working for one Direction  and one variable
//Model* model_auto2(const VarioExp& vario,const VectorInt& lst_model, int uc, int idir, int ivar)
//{
//
//  Model* res;
//  Vario* v = vario.toGeoslib();
//  VectorDouble vars;
//
//  vars = vario.getVariance();
//
//  //info : 0 flag_linked - field,
//  //      NULL means
//  res = model_init(vario.getNDim(), vario.getNVar(), 0, 0, 0., 0, VectorDouble(), vars);
//
//  if (uc != -1)
//    model_add_drift(res,EDrift::UC,0) ;
//
//  for (int i = 0; i < (int)lst_model.size(); i++)
//  {
//    model_add_cova(res, lst_model[i], 0, 0, TEST, TEST, VectorDouble(), VectorDouble(), vars);
//  }
//
//  Option_AutoFit opt_mauto;
//  Option_VarioFit option;
//  Constraints constraints;
//  model_auto_fit(v, res, false, opt_mauto, constraints, option);
//  return (res);
//}
/*
 **Calculate the value of the model for a set of distances
 */
//VectorDouble model_evaluate2(Model* model,const VarioExp& vario,int idir, int ivar)
//{
//
//  VarioDir vdir = vario.getVarioDir(idir); // work on direction wanted
//  // get information to know how many calculation point will be needed
//  int nlag = vdir.getNLag();
//  int lag = vdir.getLag();
//  double nb_calculation_point = nlag * lag;
//
//  // get direction coefficients
//  VectorDouble codir = vdir.getNormDir().getCoord();
//
//  // create vector of increment , in which model will be evaluate (0:nb_calculation_point)
//  VectorDouble incr(nb_calculation_point);
//  std::iota (std::begin(incr), std::end(incr), 0);
//
//  // create vector in which evaluation will be stored
//  VectorDouble res(nb_calculation_point) ;
//
//  model_evaluate(model, 0, 0, 0, 0, 0, 0, 0, 0, 0, lag * nlag, codir,
//                 incr.data(),res.data());
//  return res;
//}

void kriging2(const Database & dbin, Database &dbout, Model* model, Neigh* neigh)
{
  Db* dbin2 = dbin.toGeoslib();
  Db* dbout2 = dbout.toGeoslib();
  kriging(dbin2,dbout2,model,neigh,EKrigOpt::PONCTUAL,1,1,0);
  dbout.fromGeoslib(dbout2);
  dbout.display();
}

Neigh* neigh_unique(int ndim)
{
  return(neigh_init(ndim, ENeigh::UNIQUE, false, false, false, false, false, 1, 1, 1, 1, 1, 1, 1, 0.5,
                    VectorDouble(), VectorDouble(), VectorInt()));
}

Neigh* neigh_moving(int ndim,int flag_sector, int flag_rotation, int nmini,int nmaxi,int nsect, int nsmax,
                    int radius,VectorDouble Rotation)
{
  return(neigh_init(ndim, ENeigh::MOVING, false, flag_sector, false, flag_rotation,false, nmini, nmaxi, nsect, nsmax,
                    1, 1 , radius, 0.5, VectorDouble(), Rotation, VectorInt()));
}

Neigh* my_neigh_init(int ndim, ENeigh type, int flag_xvalid, int flag_sector,
                     int flag_aniso, int flag_rotation, int flag_continuous,
                     int nmini, int nmaxi, int nsect, int nsmax, int skip,
                     double width, double radius, double dist_count,
                     VectorDouble nbgh_radius, VectorDouble nbgh_rotmat,
                     VectorInt nbgh_image)
{
  return(neigh_init(ndim,type,flag_xvalid,flag_sector,flag_aniso,flag_rotation,
        flag_continuous,nmini,nmaxi,nsect,nsmax,skip,width,
        radius,dist_count, nbgh_radius,nbgh_rotmat,nbgh_image));
}

VectorDouble affiche(Db* db)
{
  return (db->getArrays());
}

void my_db_print(Db* db)
{
  db_print(db,1,1,1,1,1);

}

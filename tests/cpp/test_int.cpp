#include "Interfaces/Database.hpp" 
#include "Interfaces/VariableDouble.hpp" 
#include "Interfaces/VariableBool.hpp" 
#include "Interfaces/VariableInt.hpp" 
#include "Interfaces/VariableString.hpp" 
#include "Interfaces/VariableCategorical.hpp" 
#include "Interfaces/Category.hpp" 
#include "Interfaces/Dictionary.hpp" 
#include "Interfaces/Param.hpp"
#include "Interfaces/ParamCSV.hpp"
#include "Space/SpacePoint.hpp"
#include "Space/SpaceRN.hpp"

#include "geoslib_f.h"
#include "Interfaces/geoslib_f_swig.h"

#include <stdio.h>
#include <iostream>

int* set_ranks(int ncol)
{
  int *ranks;
  int    i = -1;
  
  ranks = (int*)malloc(sizeof(int)*(ncol+ 1));
  
  while(++i < ncol)
  {
    ranks[i] = i;
  }
  return(ranks);
}

std::vector<bool> create_selection(int nech)
{
  std::vector<bool> res;
  int i = -1;
  
  while (++i < nech)
   res.push_back((i+1) % 2);
  return(res);
}

VectorDouble create_z(int n)
{
  VectorDouble res;

  int i = 0;
  while (++i< n + 1 )
    res.push_back(i);
  return(res);
}

void ptg()
{
  ASpaceObject::createGlobalSpace(SPACE_RN,2);

  VectorInt    nx{5,5};
  VectorDouble dx{1,1};
  VectorDouble x0{0,0};
  VectorDouble rotation{0}; 
  
  ParamGrid pgrid(nx, x0, dx, rotation, LOAD_BY_COLUMN);
  Database dbg(pgrid);
  
  Database dbp;

  VariableDouble* v_coord1= new VariableDouble("Coord1");
  VectorDouble v_double{0.5,2.5,4};
  v_coord1->setValues(v_double);
  dbp.addVar(v_coord1);
  
  VariableDouble* v_coord2= new VariableDouble("Coord2");
  VectorDouble v_double2{0.5,2.5,2};
  v_coord2->setValues(v_double2);
  dbp.addVar(v_coord2);
  
  VariableDouble* v_z= new VariableDouble("Z1");
  VectorDouble v_doublez=create_z(3);
  v_z->setValues(v_doublez);
  dbp.addVar(v_z);

  dbp.setRole({"Coord1","Coord2"}, ERoles::COORD);
  dbp.setRole({"Z1"}, ERoles::Z);
  migrate_point_to_grid2(dbp, dbg, "Z1", 0, VectorDouble());
  dbg.display();
}

void gtp()
{
  ASpaceObject::createGlobalSpace(SPACE_RN,2);

  VectorInt nx{5,5};
  VectorDouble dx{1,1};
  VectorDouble x0{0,0};
  VectorDouble rotation{0}; 
  

  ParamGrid pgrid(nx,x0,dx,rotation,LOAD_BY_COLUMN);

  Database dbg(pgrid);
  Database dbp;

  VariableDouble* v_coord1= new VariableDouble("Coord1");
  VectorDouble v_double{0.5,2.5,4};
  v_coord1->setValues(v_double);
  dbp.addVar(v_coord1);
  
  VariableDouble* v_coord2= new VariableDouble("Coord2");
  VectorDouble v_double2{0.5,2.5,2};
  v_coord2->setValues(v_double2);
  dbp.addVar(v_coord2);

  dbp.setRole({"Coord1","Coord2"},ERoles::COORD);
  dbp.display(); 
  VectorDouble z = create_z(25);
  dbg.addVar("Z1",z);

  migrate_grid_to_point2(dbg, dbp, "Z1", 0, VectorDouble());
  
  dbp.display();
}

void gtg()
{
  ASpaceObject::createGlobalSpace(SPACE_RN,2);

  VectorInt nx{5,5};
  VectorDouble dx{1,1};
  VectorDouble x0{0,0};
  VectorDouble rotation{0}; 
  
  ParamGrid pgrid(nx,x0,dx,rotation,LOAD_BY_COLUMN);
  ParamGrid pgrid2({2,2},x0,{2.5,2.5},rotation,LOAD_BY_COLUMN);
  Database dbg(pgrid);
  Database dbg2(pgrid2);

  VectorDouble z = create_z(dbg.getGridSize());
  dbg.addVar("Z1",z);

  migrate_grid_to_grid2(dbg, dbg2, "Z1", 0, VectorDouble());
  
  dbg.display();
  dbg2.display();
}

#ifdef _USE_NETCDF
void  serialize()
{
  ASpaceObject::createGlobalSpace(SPACE_RN,2);

  VectorInt nx{5,5};
  VectorDouble dx{2,1};
  VectorDouble x0{3,0};
  VectorDouble rotation{0}; 
  
  ParamGrid pgrid(nx,x0,dx,rotation,LOAD_BY_COLUMN);
  Database dbg(pgrid);
  
  VectorDouble Z = create_z(25);
  dbg.addVar("Z",Z);
  dbg.display();
  dbg.serialize("dbg.nc");
  dbg.deserialize("dbg.nc");
  dbg.display();

  Database dbp;
  VectorDouble z = create_z(3);
  dbp.addVar("Z1",z);
  
  VariableInt* var_int= new VariableInt("int");
  VectorInt val_int{3,3,4};
  var_int->setValues(val_int);
  dbp.addVar(var_int);
  
  VariableBool* var_bool= new VariableBool("bool");
  VectorBool vec_bool{true, false, false};
  var_bool->setValues(vec_bool);
  dbp.addVar(var_bool);

  VariableString* var_string= new VariableString("string");
  VectorString val_string{"abcde","gfdgr","123"};
  var_string->setValues(val_string);
  dbp.addVar(var_string);
 
  dbp.setRole({"Z1","int"},ERoles::COORD);
 // dbp.setRole({"string","bool"},ERoles::Z);
  dbp.serialize("dbp.nc");
  dbp.deserialize("dbp.nc");
  dbp.display();
}
#endif

void vario(const std::string& file)
{
  ASpaceObject::createGlobalSpace(SPACE_RN, 2);

  //read CSV and create Database
  ParamCSV pcsv(file,",",".",true,0);
  Database database(pcsv);
  
  database.setRole({"x1","x2"},ERoles::COORD);
  database.setRole({"Uniform","Simu"},ERoles::Z);
  
  //database.display();

  ParamVario pVario;
 
  //Give paramater of Vario for one direction
  ParamVarioDir pVarioDir;
  pVarioDir.lag = 5;
  pVarioDir.nlag = 10;
  SpacePoint ndir;
  ndir.setCoordFromAngle({20});
  pVarioDir.normDir=ndir;
  pVario.dirs.push_back(pVarioDir);
  
  ParamVarioDir pVarioDir2;
  pVarioDir2.lag = 4;
  pVarioDir2.nlag = 12;
  SpacePoint ndir2;
  ndir2.setCoordFromAngle({70});
  pVarioDir2.normDir=ndir2;
  pVario.dirs.push_back(pVarioDir2);

//  CalcVarioExp calc;
//  calc.setInputData(database);
//  calc.setParamVario(pVario);
//  calc.run();
//  calc.getVarioExp().display();
}


int main(int argc,char **argv)
{
  setup_license("Demonstration");
  constant_define("NTROW",-1);
  if (argc >= 2)
  {
    if (!strcmp(argv[1],"ptg"))
    {
      ptg();
    }
    else if (!strcmp(argv[1],"gtp"))
    {
      gtp();
    }
    else if (!strcmp(argv[1],"gtg"))
    {
      gtg();
    }
    else if (!strcmp(argv[1],"serialize"))
    {
#ifdef _USE_NETCDF
      serialize();
#else
      std::cout << "Error: 'serialize' test not available (NetCDF is missing)" << std::endl;
#endif
    }
    else if (!strcmp(argv[1],"vario"))
    {
      if (argc == 3)
      {
        vario(argv[2]);
      }
      else
      {
        std::cout << "Error: Missing input file name for 'vario' test (try simunif.csv?)" << std::endl;
      }
    }
    else
    {
      std::cout << argv[1] << " is not a valid argument"<<std::endl;
    }
  }
}


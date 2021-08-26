#include "Interfaces/CalcVarioExp.hpp"
#include "Interfaces/Database.hpp"

#include "Space/ASpace.hpp"

#include "geoslib_f.h"

/*********************************************************************
** Default Constructor
*********************************************************************/
CalcVarioExp::CalcVarioExp(const ASpace* space)
: ASpaceObject(space),
  _param(space),
  _data(NULL),
  _res(space)
{
}

/*********************************************************************
** Destructor
*********************************************************************/
CalcVarioExp::~CalcVarioExp()
{
}

/*********************************************************************
** Override method from ASpaceObject , checking if object is well according
** to dimension
*********************************************************************/
bool CalcVarioExp::isConsistent(const ASpace* space) const
{
  return (_param.isConsistent(space) && _data->isConsistent(space));
}


/*********************************************************************
** ParamVario setter
*********************************************************************/
void CalcVarioExp::setParamVario(const ParamVario &p)
{
  _param = p;
}

/*********************************************************************
** Database setter
*********************************************************************/
//:WARNING: i'm not sure about this
void CalcVarioExp::setInputData(const Database &db)
{
  _data = (Database*)&db;
 // _space = db.getSpace(); 
}

/*********************************************************************
** Database setter
*********************************************************************/
VarioExp CalcVarioExp::getVarioExp() const
{
  return(_res);
}

/*********************************************************************
** Run
*********************************************************************/
void CalcVarioExp::run()
{
  if (check())
  {
    Vario* vario = getVario();
    //:DEBUG:
    //variogram_print(vario,1);
   
    _res.fromGeoslib(vario, _param);
  }
}

/********************************************************************
** Create,calcul and return a Vario*(geoslib struct).
**
** Actually use functions from geoslib( variogram_init, variogram_direction_add
** and variogram_calcul
**
********************************************************************/
Vario* CalcVarioExp::getVario() const
{
  Vario* vario;
 
  //:INFO: 0 : ndate,  0 : scaling factor for transitive variogram, NULL : Mean of variables(only for CALCUL_POISSON), NULL: array of date interval
  vario = variogram_init("vg");
  
  //Add all direction
  for (const auto& dir : _param.dirs)
  {
  //: npatot & size?? (paramOut) les stocker dans paramVarioDir??
    int opc = dir.getFirstOperCond(ROLE_CODE);// correspond to role code in db
    int tol = dir.getFirstTol(ROLE_CODE); // correspond to role code in db
    int bench = dir.getFirstTol(ROLE_COORD); // correspond to ROLE_COORD in db
    //:INFO: 0 : Idate 
    variogram_direction_add(vario,
                            dir.nlag,
                            opc, 0, dir.lag, dir.dlag,
                            dir.pencil.angle, bench,
                            dir.pencil.radius,
                            tol, dir.irregularLags,
                            dir.normDir.getCoord(),
                            dir.gridIncr);
  }
  //:WARNING: 0 (first) : covariance Model (ParamVariofit??) 

  vario->compute(_data->toGeoslib(), "vg", VectorDouble(), VectorDouble(),
                   1, (int)_param.rules);

  return(vario);
}

bool CalcVarioExp::check() const
{
 if (!_data->isConsistent(getSpace()))
 {
    std::cout << "Database is not consistent with Dimension" << std::endl;
    return(false);
 }
 if (!_data->roleExist(ROLE_COORD))
 {
    std::cout << "No Role Coord in Database" << std::endl;
    return(false);
 }
 if (!_data->roleExist(ROLE_Z))
 {
    std::cout << "No Role Z in Database" << std::endl;
    return(false);
 }

 // possible method check in ParamVario
 if (!_param.isConsistent(getSpace()))
 {
    std::cout << "Variogram Parameters are not consistent with Dimension" << std::endl;
    return(false);
 }

 for (const auto& dir : _param.dirs)
 {
   if (dir.getNbCond(ROLE_CODE) > 1)
   {
      std::cout << "Only one cond with Role Code authorize at the moment" << std::endl;
      return(false);
   }

   if (dir.getNbCond(ROLE_COORD) > 1)
   {
      std::cout << "Only one cond with Role Coord authorize at the moment" << std::endl;
      return(false);
   }
   
   for (const auto& cond : dir.paramVarioConds)
   {
      if (cond.role == ROLE_COORD)
      {
        if ((unsigned int) (_data->getNVarRole(ROLE_COORD) - 1) != (unsigned int) cond.iRole)
        {
          std::cout << " Condition with ROLE_COORD can only be on last Variable with ROLE_COORD at the moment" << std::endl;
          return(false);
        }
      }
   }
 }
  return(true); 
}

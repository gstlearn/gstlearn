#ifndef DATABASE_HPP
#define DATABASE_HPP

#include <string>
#include <vector>
#include <string.h>

#include "geoslib_f.h"

#include "Space/ASpaceObject.hpp" // ISA
#include "Interfaces/interface_d.hpp"

#include "Interfaces/ParamGrid.hpp"

class AVariable;
class ASpace;
class ParamCSV;

/*****************************************************************************
 * A Database is a representation of an excel sheet
 * 
 ****************************************************************************/
class Database : public ASpaceObject
{
public:
  Database(const ASpace* space = nullptr);
  Database(const ParamCSV& pcsv,   const ASpace* space = nullptr);
  Database(const ParamGrid& pgrid, const ASpace* space = nullptr);// TODO : Database from ParamGrid only valid in 2D from now
  Database(const Database &db);
  virtual ~Database();
  Database& operator=(const Database& ref);

  virtual bool isConsistent(const ASpace* space) const override;

  virtual String toString(int level = 0) const override;

  unsigned int         getNSamples() const;
  unsigned int         getNVars() const;
  VectorString         getNames() const;
  unsigned int         getNVarRole(Roles role) const;
  Roles                getRole(int ivar) const;
  int                  getIVar(const String& name) const;
  std::pair<Roles,int> getRoleAndIRole(const String& name) const;
  String               getName(Roles role, int i_role) const;
  VectorDouble         getValuesByName(const String& name);
  AVariable*           getVariable(int ivar);

  const ParamGrid& getParamGrid() const { return _pgrid; }

  ES setName(int iatt, const String& name);
  ES setRole(const VectorString& names, Roles role);
  
  ES addVar(const String& name,const VectorDouble& val);
  ES delVar(const String& name);
  ES eraseRole(Roles role);
  ES select(const String& name, const VectorBool& sel);
  
  void printRoles();
  virtual void display_old() const;
  Db* toGeoslib() const;
  void fromGeoslib(Db* db);
  int nameIdentify(const String& name) const;
  int getGridSize() const;
  void reset();
  bool roleExist(Roles role);

#ifdef _USE_NETCDF
  bool serialize(const String& str) const;
  bool deserialize(const String& str);
#endif
  
  ES addVar( AVariable* var); //:WARNING:: can't use this function in python

private:
  ES    addVar(AVariable* var, Roles role);
  int   getNbCoord();
  bool  nameExist(const String& name) const;
  ES    indexOOR(int i) const;

  void  changeKey(std::multimap<Roles,String>::iterator it, Roles new_role);
  std::multimap<Roles,String>::const_iterator getItRole(const String& name) const;

private:
  ParamGrid   _pgrid;
  bool        _isGrid;
  std::vector<AVariable*> _vars;
  std::multimap<Roles,String> _roles; //order have an importance

};

#endif

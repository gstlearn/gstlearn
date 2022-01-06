/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Space/ASpaceObject.hpp"
#include "Interfaces/interface_d.hpp"
#include "Interfaces/ParamGrid.hpp"
#include "geoslib_enum.h"

#include <string>
#include <vector>
#include <string>

class Db;
class AVariable;
class ASpace;
class ParamCSV;

/*****************************************************************************
 * A Database is a representation of an excel sheet
 * 
 ****************************************************************************/
class GSTLEARN_EXPORT Database : public ASpaceObject
{
public:
  Database(const ASpace* space = nullptr);
  Database(const ParamCSV& pcsv,   const ASpace* space = nullptr);
  Database(const ParamGrid& pgrid, const ASpace* space = nullptr);// TODO : Database from ParamGrid only valid in 2D from now
  Database(const Database &db);
  virtual ~Database();
  Database& operator=(const Database& ref);

  virtual bool isConsistent(const ASpace* space) const override;

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  unsigned int         getNSamples() const;
  unsigned int         getNVars() const;
  VectorString         getNames() const;
  unsigned int         getNVarRole(ERoles role) const;
  ERoles                getRole(int ivar) const;
  int                  getIVar(const String& name) const;
  std::pair<ERoles,int> getRoleAndIRole(const String& name) const;
  String               getName(ERoles role, int i_role) const;
  VectorDouble         getValuesByName(const String& name);
  AVariable*           getVariable(int ivar);

  const ParamGrid& getParamGrid() const { return _pgrid; }

  ES setName(int iatt, const String& name);
  ES setRole(const VectorString& names, ERoles role);
  
  ES addVar(const String& name,const VectorDouble& val);
  ES delVar(const String& name);
  ES eraseRole(ERoles role);
  ES select(const String& name, const VectorBool& sel);
  
  void printRoles();
  virtual void display_old() const;
  Db* toGeoslib() const;
  void fromGeoslib(Db* db);
  int nameIdentify(const String& name) const;
  int getGridSize() const;
  void reset();
  bool roleExist(ERoles role);

#ifdef _USE_NETCDF
  bool serialize(const String& str) const;
  bool deserialize(const String& str);
#endif
  
  ES addVar( AVariable* var); //:WARNING:: can't use this function in python

private:
  ES    addVar(AVariable* var, ERoles role);
  int   getNbCoord();
  bool  nameExist(const String& name) const;
  ES    indexOOR(int i) const;

  void  changeKey(std::multimap<ERoles,String>::iterator it, ERoles new_role);
  std::multimap<ERoles,String>::const_iterator getItRole(const String& name) const;

private:
  ParamGrid                    _pgrid;
  bool                         _isGrid;
  std::vector<AVariable*>      _vars;
  std::multimap<ERoles,String> _roles; //order have an importance

};

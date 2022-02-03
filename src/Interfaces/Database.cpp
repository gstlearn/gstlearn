#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "Interfaces/interface_d.hpp"
#include "Interfaces/Database.hpp"
#include "Interfaces/AVariable.hpp"
#include "Interfaces/VariableDouble.hpp"
#include "Interfaces/VariableBool.hpp"
#include "Interfaces/VariableInt.hpp"
#include "Interfaces/VariableString.hpp"
#include "Interfaces/Param.hpp"
#include "Interfaces/ParamCSV.hpp"
#include "Space/ASpace.hpp"
#include "Db/ELoadBy.hpp"
#include "Db/ELoc.hpp"
#include "Db/Db.hpp"

#include <cstddef>
#include <algorithm>

ENUM_DEFINE(ENUM_ROLES)

ENUM_DEFINE(ENUM_CALC_RULES)

Database::Database(const ASpace* space)
    : ASpaceObject(space),
      _pgrid(),
      _isGrid(false),
      _vars(),
      _roles()
{
}

/**
 *  Constructor of Database from ParamCSV.
 *  Create a Db* and transfer data from Db* to _avars.
 *
 *  @param[in] pcsv   Object containing information about how to read csv
 *  @param[in] space  Contextual space pointer
 *
 *  @remark   Each Variable is add as VariableDouble
 *            roles are set to ERoles::NOROLE for each variable
 */

Database::Database(const ParamCSV &pcsv, const ASpace* space)
    : ASpaceObject(space),
      _pgrid(),
      _isGrid(false),
      _vars(),
      _roles()
{
  Db* db;
  int ncol_arg;
  int nrow_arg;
  VectorString names;
  VectorDouble tab;
  int skipnline = pcsv.getSkipNLines();
  String filename = pcsv.getFilePath();
  char sepchar = pcsv.getSeparatorChar();
  char decchar = pcsv.getDecimalChar();

  //read content from csv
  if (csv_table_read(filename, 1, pcsv.getUseHeader(), skipnline,
                     sepchar, decchar, "MISS", -1, -1,
                     &ncol_arg, &nrow_arg, names, tab) != 0)
  {
    // Error, create empty database
    std::cout << "Cannot open csv file!" << std::endl;
    return;
  }

  db = db_create_point(nrow_arg, ncol_arg, ELoadBy::SAMPLE, 0, tab);

  // fill attribute _vars
  int i = 0;
  VectorDouble v_val(db->getSampleNumber());

  while (i < ncol_arg)
  {
    db_vector_get_att(db, i, v_val.data());
    VariableDouble* new_var = new VariableDouble("");
    new_var->setValues(v_val);
    _vars.push_back(new_var);
    _roles.insert(std::pair<ERoles, String>(ERoles::NOROLE, names[i]));
    (*_vars[i]).setName(names[i]);
    i++;
  }

  db_delete(db);
}

/**
 * Constructor of Database from ParamGrid
 *
 *  @param[in] pgrid Grid parameters
 *  @param[in] space Contextual space pointer
 */
Database::Database(const ParamGrid &pgrid, const ASpace* space)
    : ASpaceObject(space),
      _pgrid(pgrid),
      _isGrid(true),
      _vars(),
      _roles()
{
}

/*
 * Copy Constructor.
 * Use Design Pattern Clonage to create Variable with good type
 */
Database::Database(const Database& src)
    : ASpaceObject(src),
      _pgrid(src._pgrid),
      _isGrid(src._isGrid),
      _vars(),
      _roles(src._roles)
{
  for (const auto& v : src._vars)
  {
    _vars.push_back(v->clone());
  }
}

Database::~Database()
{
  unsigned int i = 0;
  while (i < _vars.size())
  {
    delete _vars[i];
    i++;
  }
  _vars.clear();
}

/**
 * Assignment operator
 * Use Design Pattern Clonage to create Variable with good type
 */
Database& Database::operator=(const Database& ref)
{
  if (this != &ref)
  {
    _pgrid = ref._pgrid;
    _isGrid = ref._isGrid;
    _roles = ref._roles;
    for (auto old_v : _vars)
    {
      delete (old_v);
    }
    for (const auto& new_v : ref._vars)
    {
      _vars.push_back(new_v->clone());
      ;
    }
  }
  return (*this);
}

/**
 *  Overload Function From SpaceObject.
 *  Verify if Database is consistent with current space
 */
bool Database::isConsistent(const ASpace* space) const
{
  return (getNVarRole(ERoles::COORD) == space->getNDim());
}

String Database::toString(const AStringFormat* strfmt) const
{
  return toGeoslib()->toString(strfmt);
}

/**
 *  Return the number of rows in the Database
 */
unsigned int Database::getNSamples() const
{
  if (getNVars() > 0)
  {
    return (_vars[0]->getNValues());
  }
  return (0);
}

/**
 *  Return the number of Column in the Database
 */
unsigned int Database::getNVars() const
{
  return (_vars.size());
}

/**
 * Return a Vector containing names of variables
 */
VectorString Database::getNames() const
{
  String name;
  VectorString names;
  int ncol = getNVars();
  int i = 0;

  while (i < ncol)
  {
    name = _vars[i]->getName();
    names.push_back(name);
    i++;
  }
  return (names);
}

/**
 * Return Number of Variable with a given Role
 *
 *  @param[in] role_ref : ENUM_ROLE
 */
unsigned int Database::getNVarRole(ERoles role_ref) const
{
  return (_roles.count(role_ref));
}

/**
 * Return Variable in the i'th position
 *
 *  @param[in] ivar : index of variable
 */
AVariable* Database::getVariable(int ivar)
{
  return (_vars[ivar]);
}

/**
 * Return Role of a variable given index of a variable
 *
 *  @param[in] ivar : indew of variable
 */
ERoles Database::getRole(int ivar) const
{
  for (auto r : _roles)
  {
    if (r.second.compare(_vars[ivar]->getName()) == 0) return (r.first);
  }
  return (ERoles::NOROLE);
}

/**
 * Return index of a variable given his name, if name does not exist in database
 *  -1 is return
 *
 *  @param[in] name: name of variable
 */
int Database::getIVar(const String& name) const
{
  int i = 0;
  for (const auto& v : _vars)
  {
    if (v->getName().compare(name) == 0) return (i);
    i++;
  }
  return (-1);
}

/**
 * Return a pair containing Role and index of the variable inside this role
 * given a variable name
 *
 *  @param[in] name: name of variable
 */
std::pair<ERoles, int> Database::getRoleAndIRole(const String& name) const
{
  std::pair<ERoles, int> res;

  int ivar = getIVar(name);
  ERoles role = getRole(ivar);

  std::pair<std::multimap<ERoles, String>::const_iterator,
      std::multimap<ERoles, String>::const_iterator> ret;
  ret = _roles.equal_range(role);
  int i = 0;
  for (std::multimap<ERoles, String>::const_iterator it = ret.first;
      it != ret.second; ++it)
  {
    if (it->second.compare(name) == 0) res = std::make_pair(role, i);
    i++;
  }
  return (res);
}

/**
 * Return a Variable name given his Role and his rank in the Role
 *
 *  @param[in] role: role of variable
 *  @param[in] i_role: Rank of the role
 */
String Database::getName(ERoles role, int i_role) const
{
  std::multimap<ERoles, String>::const_iterator res;
  res = _roles.find(role);
  int i = 0;
  while (i < i_role)
  {
    i++;
    res++;
  }
  return (res->second);
}

/**
 *  Add a new VariableDouble into the Database
 *
 *  @param[in]   name : name of the Column
 *  @param[in]   values: VectorDouble containing the values according to the type
 */
ES Database::addVar(const String& name, const VectorDouble& values)
{
  if (nameExist(name))
  {
    return (ES_NAME);
  }
  VariableDouble* v;
  v = new VariableDouble(name);
  v->setValues(values);
  return (addVar(v));
}

/**
 *
 *  Add a new Variable into the Database
 *
 *  @param[in]   var : variable to add.
 *
 *  \remark      role of the variable as set to ERoles::NOROLE
 */
/*:WARNING:: WHEN call to This function in python, cause a segfault in Db destructor:
 due to "delete(var)"*/

ES Database::addVar(AVariable* var)
{
  if (!_vars.empty())
  {
    unsigned int size = _vars[0]->getNValues();
    if (var->getNValues() != size)
    {
      return (ES_SIZE_VAR);
    }
  }
  //std::cout << "Adding variable of type:" << typeid(*var).name() << std::endl;
  _vars.push_back(var);
  _roles.insert(std::pair<ERoles, String>(ERoles::NOROLE, var->getName()));
  return (ES_NOERROR);
}

/**
 *  Add a new Variable into the Database with a specific Role
 *  private function , only use in Deserialization

 **/
ES Database::addVar(AVariable* var, ERoles role)
{
  ES es = addVar(var);

  if (es == ES_NOERROR)
    _roles.insert(std::pair<ERoles, String>(role, var->getName()));

  return (es);
}

/**
 *  set the name of a column giving is attribut
 *
 *  @param[in] iatt : indice of the attribute
 *  @param[in] name : name to give to the column
 */
//:TODO:WARNING: We have to change name in multimap too
ES Database::setName(int iatt, const String& name)
{
  //check unicity of name
  if (nameExist(name))
  {
    return (ES_NAME);
  }
  if (indexOOR(iatt))
  {
    return (ES_INDEX_OOR);
  }
  _vars[iatt]->setName(name);
  return (ES_NOERROR);
}

/**
 *  Delete a column from the Database
 *
 *  @param[in]  name : Name of the Column that will be deleted
 *
 */
ES Database::delVar(const String& name)
{
  int iatt = nameIdentify(name);

  if (iatt != -1)
  {
    delete _vars[iatt];
    _vars.erase(_vars.begin() + iatt);
    _roles.erase(getItRole(name));
    return (ES_NOERROR);
  }
  return (ES_NAME);
}

/**
 *  Set a Role for multiple variable
 *  Order is important.
 *
 *  @param[in]  names : names of the variable that will be affected the Role
 *  @param[in]  role  : ENUM ERoles , role given to the columns
 */
ES Database::setRole(const VectorString& names, ERoles role)
{
  int iatt;
  eraseRole(role);
  if (_isGrid && role == ERoles::COORD)
  {
    return (ES_PERMISSION_GRID_COORD);
  }
  // only set Role to Variable found
  for (auto name : names)
  {
    iatt = nameIdentify(name);
    if (iatt != -1)
    {
      std::multimap<ERoles, String>::iterator it = _roles.begin();
      while (it != _roles.end())
      {
        if (name.compare(it->second) == 0)
        {
          changeKey(it, role);
        }
        it++;
      }
    }
    else
    {
      std::cout << "Variable: " << name
                << " doesn't exist in database, it is ignored!" << std::endl;
    }
  }

  return (ES_NOERROR);
}

/**
 * Find all the Variable that have the given role  ans set their role to
 * ERoles::NOROLE
 *
 *  @param[in]  role  : ENUM ERoles , role to be erase
 *
 *  \remark     If Database is a grid and the role is ERoles::COORD, an error
 *              occur: Forbidden to  change or erase the ERoles::COORD for a grid
 */

ES Database::eraseRole(ERoles role)
{

  if (_isGrid && role == ERoles::COORD)
  {
    return (ES_PERMISSION_GRID_COORD);
  }
  std::pair<std::multimap<ERoles, String>::iterator,
      std::multimap<ERoles, String>::iterator> ret;
  ret = _roles.equal_range(role);
  for (std::multimap<ERoles, String>::iterator it = ret.first; it != ret.second;
      ++it)
  {
    String name = it->second;
    _roles.erase(it);
    _roles.insert(std::pair<ERoles, String>(ERoles::NOROLE, name));
  }
  return (ES_NOERROR);
}

/**
 *  Create a new VariableBool , which will have role ERoles::SEL (selection)
 *
 *  @param[in]  name  : name given to the column created
 *  @param sel : Selection
 */

ES Database::select(const String& name, const VectorBool& /*sel*/)
{
  VariableBool* var_bool = new VariableBool(name);
  VectorString names = { name };
  addVar(var_bool);
  setRole(names, ERoles::SEL);
  return (ES_NOERROR);
}

/**
 *  return the pair given a variable Name
 *
 *  @param[in]  name  : name of the Variable
 */

std::multimap<ERoles, String>::const_iterator Database::getItRole(const String& name) const
{
  for (std::multimap<ERoles, String>::const_iterator it = _roles.begin();
      it != _roles.end(); ++it)
  {
    if (it->second.compare(name) == 0) return (it);
  }
  return (_roles.end());
}

/**
 *  Use For Debug, call function db_print from geoslib
 */

void Database::display_old() const
{

  db_print(toGeoslib(), 1, 1, 1, 1, 1);

}

/**
 * Function that clean the object.
 *
 */
void Database::reset()
{
  _pgrid.reset();
  _isGrid = false;
  unsigned int i = 0;
  while (i < _vars.size())
  {
    delete _vars[i];
    i++;
  }
  _vars.clear();
  _roles.clear();
}

/**
 * Check if at least One Variable as the Given Role
 *
 */
bool Database::roleExist(ERoles ref)
{
  if (_roles.find(ref) != _roles.end()) return (true);
  return (false);
}

/*
 *  Check if a name allready exist in the Database
 *
 *  @param[in] name : name to be checked
 */
bool Database::nameExist(const String& test_name) const
{
  VectorString names = getNames();

  for (const auto& name : names)
  {
    if (!test_name.compare(name))
    {
      return (true);
    }
  }
  return (false);
}

/*
 * Tools function to check if the index correspond to an existing Variable
 *
 *  @param[in]  i  index to check
 */
ES Database::indexOOR(int ivar) const
{
  if (ivar < (int) _vars.size() && ivar >= 0)
  {
    return (ES_NOERROR);
  }
  else
  {
    // mes_error(ES_INDEX_OOR);
    return (ES_INDEX_OOR);
  }
}

/**
 * Function returning indice of a Variable given it's name, if there is no
 * Correspondance : -1 will be return.
 *
 *  @param[in] name  name of the variable to search
 */
int Database::nameIdentify(const String &name) const
{
  int i = 0;

  while (i < (int) _vars.size())
  {
    if (!name.compare(_vars[i]->getName()))
    {
      return (i);
    }
    i++;
  }
  return (-1);
}

/**
 * tool function returning number of grid nodes.
 */
int Database::getGridSize() const
{
  if (_isGrid)
  {
    int res = 1;
    for (const auto& n : _pgrid.getNx())
    {
      res = res * n;
    }
    return (res);
  }
  return (-1);
}

/**
 * Create and return a  Db*(geoslib struct) from a Database object
 *
 */
//:Tricky
Db* Database::toGeoslib() const
{
  Db* res;
  std::vector<AVariable*>::const_iterator var;
  VectorDouble vec;
  VectorDouble v;
  int irole, jrole, iatt = 0;
  std::vector<VectorDouble> vec_role(ERoles::getSize());

  // Concatenate all column in one vector
  for (const auto& var : _vars)
  {
    vec = var->getValues();
    v.insert(v.end(), vec.begin(), vec.end());
  }
  //create Db
  if (!_isGrid)
  {
    int nrow = _vars[0]->getValues().size();
    res = db_create_point(nrow, _vars.size(), ELoadBy::COLUMN, 1, v);
  }
  else
  {
    res = db_create_grid(0, 2, 0, ELoadBy::COLUMN, 1, _pgrid.getNx(),
                         _pgrid.getX0(), _pgrid.getDx());
    //:TODO:Add VariableAuto with generator on the fly
  }
  // add name
  int ivar = 0;
  while (ivar < (int) _vars.size())
  {
    if (_isGrid)
    {
      int new_att;
      new_att = res->addColumnsByConstant(1, 0);
      res->setColumnByUID(_vars[ivar]->getValues(),new_att);
    }
    db_name_set(res, iatt, _vars[ivar]->getName().c_str());
    iatt++;
    ivar++;
  }

  /*
   create a vector<double> for each possible Role.
   those vector contains indices of variables which have this role
   remember vec_role is a vector<vector<double>>
   */

  iatt = 0;
  for (auto r : _roles)
  {
    if (r.first != ERoles::NOROLE) vec_role[r.first.getValue()].push_back(getIVar(r.second));
  }

  // assign those role in the create db
  irole = 0;
  while (irole < (int)(ERoles::getSize()))
  {
    if (!vec_role[irole].empty())
    {
      jrole = 0;
      while (jrole < (int) vec_role[irole].size())
      {
        res->setLocatorByUID(vec_role[irole][jrole], ELoc::fromValue(irole), jrole+1);
        jrole++;
      }
    }
    irole++;
  }
  return (res);
}

//:TODO : case grid
void Database::fromGeoslib(DbGrid* db)
{
  int i = 0;

  if (db->isGrid())
  {
    _isGrid = true;
    _pgrid.fromGeoslib(db->getGrid());
  }
  int nvar = db->getColumnNumber();
  while (i < nvar)
  {
    VectorDouble val;
    for (int iech = 0; iech < db->getSampleNumber(); iech++)
      val[iech] = db->getArray(iech,i);
    VariableDouble* new_var = new VariableDouble(db->getNameByColIdx(i));
    new_var->setValues(val);
    addVar(new_var);
    i++;
  }
  VectorString lst_names = getNames();
  int j = 0;
  //fill v_name with name of variable for each possible ERoles.
  auto it = ELoc::getIterator();
  while (it.hasNext())
  {
    if (*it != ELoc::UNKNOWN)
    {
      VectorString names_role;
      int k = 0;
      while (k < db->getFromLocatorNumber(*it))
      {
        String name = lst_names[db->getUIDByLocator(*it,k)];
        names_role.push_back(name);
        k++;
      }
      setRole(names_role, ERoles::fromValue(j));
    }
    it.toNext();
  }
}

/**
 * function to print the content of the multimap
 */
void Database::printRoles()
{
  reset();
  for (const auto& role : _roles)
  {
    std::cout << role.first.getKey() << "    " << role.second << std::endl;
  }
}

/**
 * private
 * Function to change the role of one Variable.(change the key in the multimap)
 */
void Database::changeKey(std::multimap<ERoles, String>::iterator it,
                         ERoles new_role)
{
  std::multimap<ERoles, String>::iterator itCopy = it;
  _roles.erase(itCopy);
  _roles.insert(std::pair<ERoles, String>(new_role, it->second));
}

/**
 * Return a vector containing values for one Variable identify by it's name
 */
VectorDouble Database::getValuesByName(const String& name)
{
  for (const auto& v : _vars)
  {
    if (!v->getName().compare(name)) return (v->getValues());
  }
  std::cout << "Name: " << name << " doesn't exist!" << std::endl;
  return (VectorDouble());
}

#ifdef _USE_NETCDF

#include <netcdf>

using namespace netCDF;
using namespace netCDF::exceptions;

static bool sortbysec(const std::pair<int,int> &a, const std::pair<int,int> &b)
{
  return (a.second < b.second);
}

/**
 * Serialize a Database in a netCDF File
 *
 * @param[in] file_name  name of the created netCDF file
 *
 * \remark actually working For Database and grid (pgrid as argument of DB)
 */
bool Database::serialize(const String& file_name) const
{

  try
  {
    NcFile sfc(file_name, NcFile::replace);

    std::vector<NcDim> dims;
    std::vector<NcVar> variable;
    if (_isGrid)
    {
      // Get size of each dimension
      VectorInt nx = _pgrid.getNx();
      int ndim = nx.size();
      std::vector<int> dim_size(ndim);
      for (int i = 0; i < ndim; i++)
      {
        dim_size[i] = nx[i];
      }

      for (int i = 0; i < ndim; i++)
      {
        String name = "dim_grid" + std::to_string(i);
        NcDim dim = sfc.addDim(name,dim_size[i]);
        //Define the dimensions
        dims.push_back(dim);

        //Define Coordinate netCDF variables
        variable.push_back(sfc.addVar(name, ncDouble, dims[i]));
        /// TODO : Add Serialize / Deserialize to ParamGrid class
        variable[i].putAtt("nx",ncInt,_pgrid.getNx()[i]);
        variable[i].putAtt("dx",ncDouble,_pgrid.getDx()[i]);
        variable[i].putAtt("x0",ncDouble,_pgrid.getX0()[i]);
        variable[i].putAtt("rot",ncDouble,_pgrid.getRotation()[i]);
        //write the coordinate variable Data
        /// TODO : Really needed ?
        variable[i].putVar(_pgrid.getValues(i).data());
      }
    }
    else
    {
      NcDim mydim= sfc.addDim("nSamples",getNSamples());
      dims.push_back(mydim);
    }

    //add each variable to the netcdf file  and create an attribute role with corresponding role
    for (int i = 0; i < (int)getNVars(); i++)
    {
      NcVar var = _vars[i]->serialize(sfc, dims);
      std::pair<ERoles,int> pair = getRoleAndIRole(_vars[i]->getName());
      var.putAtt("role",ncInt ,(int)pair.first.getValue());
      var.putAtt("irole",ncInt, pair.second);
    }
    return(true);
  }
  catch(NcException& e)
  {
    e.what();
    return (false);
  }
  return(true);
}

/**
 * Create a Database from a netCDF File
 *
 * @param[in] filename  name of the netCDF file to be read
 *
 * \remark  actually working For Database and grid (pgrid as argument of DB)
 *
 */
//:TODO: When Deserialize , variable will be sort alphabetically by name in Database
// That's due to the way Netcdf store attribute by attribute
bool Database::deserialize(const String& filename)
{
  reset();

  try
  {
    NcFile datafile(filename, NcFile::read);
    VectorInt vec_nx;
    VectorDouble vec_dx;
    VectorDouble vec_x0;
    VectorDouble vec_rot;
    std::vector<std::vector<std::pair<int,int>>> info_roles(ERoles::getSize()); // used after the treatment of all variable
    // vector of size ROLEMAX. each element concern a role
    // for each role, there is a vector containing the index of the variable and the index of role rank.
    auto vars = datafile.getVars();
    int ivar = 0;
    for(auto pvar : vars)
    {
      String vname = pvar.first;
      NcVar var = pvar.second;

      if (vname.find("dim_grid") != std::string::npos)
      {
        //recup data from netcdf variablesi
        _isGrid = true;
        NcVarAtt att;

        att = var.getAtt("nx");
        int nx;
        att.getValues(&nx);
        vec_nx.push_back(nx);

        att = var.getAtt("dx");
        double dx;
        att.getValues(&dx);
        vec_dx.push_back(dx);

        att = var.getAtt("x0");
        double x0;
        att.getValues(&x0);
        vec_x0.push_back(x0);

        att = var.getAtt("rot");
        double rot;
        att.getValues(&rot);
        vec_rot.push_back(rot);
      }
      else
      {
        //get type attribute
        NcVarAtt att_type;
        String type;
        att_type = var.getAtt("type");
        att_type.getValues(type);
        // get role attribute
        NcVarAtt att_role;
        int role;
        att_role = var.getAtt("role");
        att_role.getValues(&role);
        //get irole attribute and fill info_role
        NcVarAtt att_irole;
        int irole;
        att_irole = var.getAtt("irole");
        att_irole.getValues(&irole);
        std::cout <<role<< std::endl;
        if (role != ERoles::NOROLE.getValue())
        info_roles[role].push_back(std::make_pair(ivar,irole));

        AVariable* v = AVariable::createVariable(type, vname);
        v->deserialize(datafile,var);
        addVar(v, ERoles::fromValue(role));
        ivar++;
      }
    }
    // set Role in good order in Database, process role after role
    int role = 0;
    for (auto info_role : info_roles)
    {
      if(info_role.size() > 1)
      {
        // sort variable with current role, sort is done by second element of the pair: irole.
        std::sort(info_role.begin(),info_role.end(),sortbysec);

        // Create list of name for setRole method
        VectorString list_name;
        for (const auto& p: info_role)
        list_name.push_back(_vars[p.first]->getName());
        //setRole
        setRole(list_name, ERoles::fromValue(role));
      }
    }
    // Create ParamGrid
    if (_isGrid)
    {
      VectorDouble rotation
      { 0};
      ParamGrid pgrid(vec_nx, vec_x0, vec_dx, vec_rot ,ELoadBy::COLUMN);
      _pgrid=pgrid;
    }
  }
  catch(NcException& e)
  {
    e.what();
  }
  return(true);
}

#endif

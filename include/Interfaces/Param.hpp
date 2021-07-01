#ifndef PARAM_HPP
#define PARAM_HPP

#include "Interfaces/AParam.hpp" // ISA
#include "Space/SpacePoint.hpp"  // HASA
#include "Space/SpaceShape.hpp"  // HASA
#include "Interfaces/interface_d.hpp" // USE enum (CalculRules)
#include <vector>          // USE

/**********************************************************************
 ** Class defining the Condition to accept a Pair of Point (refering "code" notion)
 ** This Class is only a storage Class
 **
 ** Contain the following attribute:
 **
 **    bool flagCode   : perform a Selection based on the code
 **    int  opercode   : Code selection Option
 **    double tolCode  : Tolerance on code
 **
 ** Default: Everything is set to 0
 ************************************************************************/
//:WARNING: ParamVarioCond must be attached to a variable
class ParamVarioCond: public AParam
{
public:
  ParamVarioCond()
      : operCond(0),
        tol(0),
        role(ROLE_NOROLE),
        iRole(0)
  {
  }
  virtual ~ParamVarioCond()
  {
  }

  ParamVarioCond(const ParamVarioCond &pvc)
      : operCond(pvc.operCond),
        tol(pvc.tol),
        role(pvc.role),
        iRole(pvc.iRole)
  {
  }

  ParamVarioCond& operator=(const ParamVarioCond& pvc)
  {
    if (this != &pvc)
    {
      operCond = pvc.operCond;
      tol = pvc.tol;
      role = pvc.role;
      iRole = pvc.iRole;
    }
    return (*this);
  }

  // To check consistency of ParamVarioCond, we need a Database.
  virtual bool checkConsistence() const override
  {
    return (true);
  }

  void setOperCode(int op_cond)
  {
    operCond = op_cond;
  }
  void setTol(double to)
  {
    tol = to;
  }
  void setRole(Roles to)
  {
    role = to;
  }
  void setIRole(int to)
  {
    iRole = to;
  }
  /*********************************************************************
   ** Description of the object overriding operator <<
   **********************************************************************/
#ifndef SWIG
  friend std::ostream& operator<<(std::ostream& stream,
                                  const ParamVarioCond pvc)
  {
    stream << " operCond : " << pvc.operCond << std::endl;
    stream << " tolCode : " << pvc.tol << std::endl;
    return (stream);
  }
#endif

  int operCond;
  double tol;
  Roles role;
  int iRole;
};

/**********************************************************************
 ** Class defining the parameter of a variogram for one direction
 ** This Class is only a storage Class
 **
 ** Contain the following attribute:
 **
 **    bool flagGrid       : if calculation are performed on a regular grid
 **    SpacePoint gridIncr : (only if flagGrid) define increment along Grid axes
 **    SpacePoint normDir  : Calculation normalized direction
 **    bool flagRegular    : True if calculation are performed on regular lags
 **    double lag          : value of the regular lag
 **    double dlag         : Tolerance on the lag (proportion of the lag value)
 **    int nlag            : number of lags
 **    VectorDouble irregularLags : Array of values for irregular lags
 **    ParamSpaceShape paramVarioShapes:  Shape criterion for accepting / rejecting pairs
 **
 **:WARNING: For now SpaceShapes is just  a pencil
 **:WARNING: For now just one ParamVarioCond
 **********************************************************************/
class ParamVarioDir: public ASpaceObject
{
public:
  ParamVarioDir(const ASpace* space = nullptr)
      : ASpaceObject(space),
        flagGrid(false),
        gridIncr(space),
        normDir(space),
        flagRegular(true),
        lag(1),
        dlag(0.5),
        nlag(10),
        irregularLags(),
        pencil()
  {
  }

  ParamVarioDir(const ParamVarioDir &pvd)
      : ASpaceObject(pvd),
        flagGrid(pvd.flagGrid),
        gridIncr(pvd.gridIncr),
        normDir(pvd.normDir),
        flagRegular(pvd.flagRegular),
        lag(pvd.lag),
        dlag(pvd.dlag),
        nlag(pvd.nlag),
        irregularLags(pvd.irregularLags),
        pencil(pvd.pencil),
        paramVarioConds(pvd.paramVarioConds)
  {
  }

  ParamVarioDir& operator=(const ParamVarioDir& ref)
  {
    if (this != &ref)
    {
      ASpaceObject::operator=(ref);
      flagGrid = ref.flagGrid;
      gridIncr = ref.gridIncr;
      normDir = ref.normDir;
      flagRegular = ref.flagRegular;
      lag = ref.lag;
      dlag = ref.dlag;
      nlag = ref.nlag;
      irregularLags = ref.irregularLags;
      paramVarioConds = ref.paramVarioConds;
    }
    return (*this);
  }

  virtual ~ParamVarioDir()
  {
  }

  virtual bool isConsistent(const ASpace* space) const override
  {
    return (gridIncr.isConsistent(space) && normDir.isConsistent(space));
  }

  /********************************************************************
   ** Return operCond From first occurence of a Role in ParaVarioConds
   *********************************************************************/
  int getFirstOperCond(Roles role) const
  {
    for (const auto& pvc : paramVarioConds)
    {
      if (pvc.role == role)
      {
        return (pvc.operCond);
      }
    }
    return (0);
  }

  /********************************************************************
   ** Return tol for the first occurence of a Role in ParaVarioConds
   *********************************************************************/
  int getFirstTol(Roles role) const
  {
    for (const auto& pvc : paramVarioConds)
    {
      if (pvc.role == role)
      {
        return (pvc.tol);
      }
    }
    return (0);
  }

  /********************************************************************
   ** Return nb  ParaVarioCond with a specific Role
   *********************************************************************/
  int getNbCond(Roles role) const
  {
    int res = 0;
    for (const auto& pvc : paramVarioConds)
    {
      if (pvc.role == role)
      {
        res++;
      }
    }
    return (res);
  }

  /*********************************************************************
   ** Description of the object overriding operator <<
   **********************************************************************/
#ifndef SWIG
  friend std::ostream& operator<<(std::ostream& stream, const ParamVarioDir pvd)
  {
    stream << "(bool)flagGrid : " << pvd.flagGrid << std::endl;

    stream << "(bool)flagRegular : " << pvd.flagRegular << std::endl;
    stream << " lag : " << pvd.lag << std::endl;
    stream << " dlag: " << pvd.dlag << std::endl;
    stream << " nlag : " << pvd.nlag << std::endl;
    for (const auto& pcond : pvd.paramVarioConds)
    {
      stream << " ParamVarioCond " << std::endl << pcond << std::endl;
    }
    return (stream);
  }
#endif

  bool flagGrid;
  SpacePoint gridIncr;
  SpacePoint normDir;
  bool flagRegular;
  double lag;
  double dlag;
  int nlag;
  VectorDouble irregularLags;
  Pencil pencil;
  //std::vector<ASpaceShape> spaceShapes; //Warning: intersection of Shape not used yet
  std::vector<ParamVarioCond> paramVarioConds;
};

/**********************************************************************
 ** Class defining The Parameter needed for a Variogram Calculation
 ** This Class is only a storage Class
 **
 ** Contain the following attribute:
 **
 **  std::vector<ParamVarioDir>  dirs  : Each Direction calculation Parameter
 **  CalculRules                 rules : Rule of calcul(Lag or Sample)
 **
 ** Default: No direction, User had to create at least one
 **          Variogram , by Lag.
 ************************************************************************/
class ParamVario: public ASpaceObject
{
public:
  ParamVario(const ASpace* space = nullptr)
      : ASpaceObject(space),
        dirs(),
        rules(CALCUL_BY_LAG),
        useWeight(false)
  {
  }
  ParamVario(const ParamVario &p_vario)
      : ASpaceObject(p_vario),
        dirs(p_vario.dirs),
        rules(p_vario.rules),
        useWeight(p_vario.useWeight)
  {
  }
  ParamVario& operator=(const ParamVario& ref)
  {
    if (this != &ref)
    {
      ASpaceObject::operator=(ref);
      dirs = ref.dirs;
      rules = ref.rules;
      useWeight = ref.useWeight;
    }
    return (*this);
  }
  virtual ~ParamVario()
  {
  }

  virtual bool isConsistent(const ASpace* space) const override
  {
    for (const auto& dir : dirs)
    {
      if (!dir.isConsistent(space)) return (false);
    }
    return (true);
  }

  /*********************************************************************
   ** Description of the object overriding operator <<
   **********************************************************************/
#ifndef SWIG
  friend std::ostream& operator<<(std::ostream& stream, const ParamVario pv)
  {
    stream << std::endl << "Number of Dir : " << pv.dirs.size() << std::endl
           << std::endl;
    int idir = 1;
    for (const auto& dir : pv.dirs)
    {
      stream << "Dir n°" << idir << std::endl;
      stream << "------------------------------------------" << std::endl;
      stream << dir;
      stream << "------------------------------------------" << std::endl;
      idir++;
    }
    return (stream);
  }
#endif

  void addDir(const ParamVarioDir& dir)
  {
    dirs.push_back(dir);
  }
  int getNdir()
  {
    return (dirs.size());
  }

  std::vector<ParamVarioDir> dirs;
  CalculRules rules;
  bool useWeight; //:WARNING: unused
};

#endif


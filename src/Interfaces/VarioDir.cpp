#include "Interfaces/VarioDir.hpp"

/**
 *  Default Constructor
 */
VarioDir::VarioDir(Dir* dir, ParamVarioDir pdir, int nvar)
  : _paramVarioDir(pdir)
{
  int j = 0;
  while (j < dir->getLagTotalNumber())
  {
    double sw = dir->getSw(j);
    double hh = dir->getHh(j);
    VectorDouble gg = dir->getGg();
    VarioValue val(sw, hh, gg);
    _lagValues.push_back(val);
    j++;
  }
}

/**
 * Copy constructor
 */
VarioDir::VarioDir(const VarioDir& ref)
  : _lagValues(ref._lagValues),
    _paramVarioDir(ref._paramVarioDir)
{
}

/**
 * Assignment Operator
 */
VarioDir& VarioDir::operator=(const VarioDir& ref)
{
  if (this != &ref)
  {
    _lagValues = ref._lagValues;
    _paramVarioDir = ref._paramVarioDir;
  }
  return (*this);
}

/**
 * Destructor
 */
VarioDir::~VarioDir()
{
}

void VarioDir::display_old() const
{
  std::cout << std::setw(20) << "Sw";
  std::cout << std::setw(20) << "Dist";
  std::cout << std::setw(20) << "Gg";
  std::cout << std::endl;
  for (const auto& val : _lagValues)
  {
    val.display();
    std::cout << std::endl;
  }
}

/**
 * get nlag
 */
int VarioDir::getNLag() const
{
  return (_paramVarioDir.nlag);
}

/**
 * get lag
 */
int VarioDir::getLag() const
{
  return (_paramVarioDir.lag);
}

/**
 * get normdir
 */
SpacePoint VarioDir::getNormDir() const
{
  return (_paramVarioDir.normDir);
}

/**
 * Create Dir* from VarioDir object, fill each variable of Dir* one by one.
 */
Dir* VarioDir::getDir(int nvar, bool flag_asym) const
{
  Dir* dir;
  dir = new(Dir);
  int nb = nvar * (nvar + 1) / 2;

  double bench = 0;
  for (const auto& pvc : _paramVarioDir.paramVarioConds)
  {
    if (pvc.role == ROLE_COORD)
    {
      bench = pvc.tol;
    }
  }

  int ndim = 2;  // TODO this value should come from elsewhere
  dir->init(ndim, _paramVarioDir.nlag, flag_asym,
            _paramVarioDir.getFirstOperCond(ROLE_CODE), 0, _paramVarioDir.lag,
            _paramVarioDir.dlag, _paramVarioDir.pencil.angle, bench,
            _paramVarioDir.pencil.radius, _paramVarioDir.getFirstTol(ROLE_CODE),
            _paramVarioDir.irregularLags, _paramVarioDir.normDir.getCoord(),
            _paramVarioDir.gridIncr.getCoord());

  // Load the variogram contents

  int j = 0;
  while (j < _paramVarioDir.nlag)
  {
    VectorDouble tmp_gg = _lagValues[j].getGg();
    int k = 0;
    while (k < nb)
    {
      int jj = j + k * _paramVarioDir.nlag;
      dir->setSw(jj,_lagValues[j].getSw());
      dir->setHh(jj,_lagValues[j].getDist());
      dir->setGg(jj, tmp_gg[k]);
      dir->setUtilize(jj, 1);
      k++;
    }
    j++;
  }
  return (dir);
}

/*Recup a Vector containing all GG for one direction*/

VectorDouble VarioDir::getGgs(int i)
{
  VectorDouble res;
  for (const auto& val : _lagValues)
  {
    res.push_back(val.getGg()[i]);
  }
  return (res);
}

VectorDouble VarioDir::getDists()
{
  VectorDouble res;
  for (const auto& val : _lagValues)
  {
    res.push_back(val.getDist());
  }
  return (res);
}

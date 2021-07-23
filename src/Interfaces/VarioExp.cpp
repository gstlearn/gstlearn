#include "Interfaces/VarioExp.hpp"
#include "Variogram/Vario.hpp"
#include <iostream>

/*********************************************************************
** Default Constructor
*********************************************************************/
VarioExp::VarioExp(const ASpace* space)
: ASpaceObject(space),
  _ndir(0),
  _nvar(0),
  _variance(),
  _varioDirs(),
  _paramVario(space)
{

}

/*********************************************************************
** Copy Constructor
*********************************************************************/
VarioExp::VarioExp(const VarioExp& ref)
: ASpaceObject(ref),
  _ndir(ref._ndir),
  _nvar(ref._nvar),
  _variance(ref._variance),
  _varioDirs(ref._varioDirs),
  _paramVario(ref._paramVario)
{

}
/*********************************************************************
** Assignment Operator
*********************************************************************/
VarioExp& VarioExp::operator=( const VarioExp& ref)
{
  if(this != &ref)
  {
    ASpaceObject::operator=(ref);
    _ndir = ref._ndir;
    _nvar = ref._nvar;
    _variance = ref._variance;
    _varioDirs= ref._varioDirs;
    _paramVario = ref._paramVario;

  }
  return (*this);
}

/*********************************************************************
**Destructor
*********************************************************************/
VarioExp::~VarioExp()
{

}


/*********************************************************************
**Get a VarioDir given his index od direction
*********************************************************************/
VarioDir VarioExp::getVarioDir(int i) const
{
  return(_varioDirs[i]);
}

/*********************************************************************
**Method overriden from ASpaceObject, check if object is consistent with the dimension
*********************************************************************/
bool VarioExp::isConsistent(const ASpace* space) const
{
  return _paramVario.isConsistent(space);
}


/*********************************************************************
** set VarioExp using a Vario*
*********************************************************************/
void VarioExp::fromGeoslib(Vario* vario, const ParamVario& pvario)
{
 _nvar = vario->getVariableNumber();
 _paramVario = pvario; 
 
 /* recup variance of variables*/
 int nb = (_nvar * (_nvar + 1)) / 2;
 int i = 0;
 while (i < nb)
 {
  _variance.push_back( vario->getVars(i));
  i++;
 }
 /* create A VarioDir for each dir in Dir** geoslib*/
// for (unsigned int i = 0; i < _paramVario.dirs.size(); i++)
//  {
//    Dir dir = vario->getDirs(i);
//    VarioDir new_dir(dir, _paramVario.dirs[i], _nvar);
//    _varioDirs.push_back(new_dir);
//
//    i++;
// }
 _ndir = i;
}

/********************************************************************
** Display
********************************************************************/
void VarioExp::display_old() const
{
  int i = 1;
  for (const auto& dir : _varioDirs)
  {
    std::cout<< "Dir n° " << i << std::endl;
    dir.display();
    i++;
  }
}

/*******************************************************************
**Create Vario* from VarioExp
*******************************************************************/
Vario* VarioExp::toGeoslib() const
{
  Vario* res;
  res = new (Vario);
  res->setCalculName("vg");

  int i = 0;
  while (i < _ndir)
  {
    Dir res_dir = Dir();
    res_dir.init(getNDim(),getNDim(), res->getFlagAsym(), 0, 0, 1., 0., 0., 0.,
                  0., 0., VectorDouble(), VectorDouble(), VectorInt());
    res->addDirs(res_dir);
    // DR: Je ne vois pas ou on donne les informations pour le calcul
    i++;
  }
  return (res);
}

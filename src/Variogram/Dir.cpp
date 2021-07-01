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
#include "Variogram/Dir.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"
#include "geoslib_f.h"

Dir::Dir(int ndim, int npas, double dpas, double toldis, double tolang)
    : _ndim(ndim),
      _nvar(0),
      _flagAsym(0),
      _nPas(npas),
      _optionCode(0),
      _idate(0),
      _dPas(dpas),
      _bench(TEST),
      _cylRad(TEST),
      _tolDist(toldis),
      _tolAngle(90.),
      _tolCode(0.),
      _breaks(),
      _codir(),
      _grincr(),
      _sw(),
      _gg(),
      _hh(),
      _utilize()
{
  _completeDefinition();
}

Dir::Dir(int ndim,
         int nvar,
         int npas,
         int flagAsym,
         double dpas,
         const VectorDouble& gg,
         const VectorDouble& hh,
         const VectorDouble& sw)
    : _ndim(ndim),
      _nvar(0),
      _flagAsym(0),
      _nPas(npas),
      _optionCode(0),
      _idate(0),
      _dPas(dpas),
      _bench(TEST),
      _cylRad(TEST),
      _tolDist(0.5),
      _tolAngle(90.),
      _tolCode(0.),
      _breaks(),
      _codir(),
      _grincr(),
      _sw(),
      _gg(),
      _hh(),
      _utilize()
{
  _completeDefinition();
  resize(nvar, flagAsym);

  if ((int) gg.size() != getSize())
    throw("Wrong dimension for Input 'gg'");
  _gg = gg;
  if (! sw.empty() && (int) sw.size() != getSize())
    throw("Wrong dimension for Input 'sw'");
  _sw = sw;
  if (! hh.empty() && (int) hh.size() != getSize())
    throw("Wrong dimension for Input 'hh'");
  _hh = hh;
}

Dir::Dir(const Dir& r)
    : _ndim(r._ndim),
      _nvar(r._nvar),
      _flagAsym(r._flagAsym),
      _nPas(r._nPas),
      _optionCode(r._optionCode),
      _idate(r._idate),
      _dPas(r._dPas),
      _bench(r._bench),
      _cylRad(r._cylRad),
      _tolDist(r._tolDist),
      _tolAngle(r._tolAngle),
      _tolCode(r._tolCode),
      _breaks(r._breaks),
      _codir(r._codir),
      _grincr(r._grincr),
      _sw(r._sw),
      _gg(r._gg),
      _hh(r._hh),
      _utilize(r._utilize)
{
  _completeDefinition();
}

Dir& Dir::operator=(const Dir& r)
{
  if (this != &r)
  {
    _ndim = r._ndim;
    _nvar = r._nvar;
    _flagAsym = r._flagAsym;
    _nPas = r._nPas;
    _optionCode = r._optionCode;
    _idate = r._idate;
    _dPas = r._dPas;
    _bench = r._bench;
    _cylRad = r._cylRad;
    _tolDist = r._tolDist;
    _tolAngle = r._tolAngle;
    _tolCode = r._tolCode;
    _breaks = r._breaks;
    _codir = r._codir;
    _grincr = r._grincr;
    _sw = r._sw;
    _gg = r._gg;
    _hh = r._hh;
    _utilize = r._utilize;

    _completeDefinition();
  }
  return *this;
}

Dir::~Dir()
{
}

void Dir::init(int ndim,
               int npas,
               double dpas,
               double toldis,
               double tolang,
               int flag_asym,
               int opt_code,
               int idate,
               double bench,
               double cylrad,
               double tolcode,
               VectorDouble breaks,
               VectorDouble codir,
               VectorDouble grincr)
{
  _ndim = ndim;
  _flagAsym = flag_asym;
  _nPas = npas;
  _optionCode = opt_code;
  _idate = idate;
  _dPas = dpas;
  _tolDist = toldis;
  _tolAngle = tolang;
  _bench = bench;
  _cylRad = cylrad;
  _tolCode = tolcode;
  _breaks = breaks;
  _codir  = codir;
  _grincr = grincr;

  _completeDefinition();
}

void Dir::_completeDefinition()
{
  if (! _breaks.empty())
  {
    if (_breaks.size() < 2) _breaks.clear();
  }
  if (_codir.empty())
  {
    _codir.resize(_ndim,0.);
    _codir[0] = 1.;
  }
  if (_grincr.empty())
  {
    _grincr.resize(_ndim,0.);
    _grincr[0] = 1;
  }
}

bool Dir::_isLagValid(int ilag) const
{
  if (ilag < 0 || ilag >= getLagNumber())
  {
    mesArg("Lag Index",ilag,getLagNumber());
    return false;
  }
  return true;
}

bool Dir::_isAddressValid(int iad) const
{
  if (iad < 0 || iad >= getSize())
  {
    mesArg("Variogram Internal Address",iad,getSize());
    return false;
  }
  return true;
}

bool Dir::_isVariableValid(int ivar) const
{
  if (ivar < 0 || ivar >= getVariableNumber())
  {
    mesArg("Variable Index",ivar,getVariableNumber());
    return false;
  }
  return true;
}

bool Dir::_isDimensionValid(int idim) const
{
  if (idim < 0 || idim >= getDimensionNumber())
  {
    mesArg("Space Dimension",idim,getDimensionNumber());
    return false;
  }
  return true;
}

void Dir::resize(int nvar, int flagAsym)
{
  _nvar = nvar;
  _flagAsym = flagAsym;

  int size = getSize();

  _sw.resize(size,0.);
  _gg.resize(size,0.);
  _hh.resize(size,0.);
  _utilize.resize(size,1);
}

double Dir::getHmax(int ivar, int jvar) const
{
  VectorDouble hh = getHh(ivar,jvar);
  return ut_vector_max(hh);
}

int Dir::getLagTotalNumber() const
{
  return ((_flagAsym) ? 2 * _nPas + 1 : _nPas);
}

int Dir::getSize() const
{
  return (getLagTotalNumber() * _nvar * (_nvar + 1) / 2);
}

bool Dir::isCalculated() const
{
  if (_sw.size() <= 0)
  {
    messerr("You must run the Calculations beforehand");
    return false;
  }
  return true;
}

void Dir::copy(const Dir& dir)
{
  if (getLagTotalNumber() != dir.getLagTotalNumber()) return;
  _sw = dir._sw;
  _gg = dir._gg;
  _hh = dir._hh;
}

void Dir::clean(void)
{
  for (int i=0; i<getSize(); i++)
  {
    _sw[i] = 0;
    _hh[i] = 0.;
    _gg[i] = 0.;
    _utilize[i] = 1.;
  }
}

int Dir::getAddress(int ivar, int jvar, int ipas, bool flag_abs, int sens) const
{
   int rank;

   /* Get the order of the variables */

   if (ivar > jvar)
     rank = ivar * (ivar+1)/2 + jvar;
   else
     rank = jvar * (jvar+1)/2 + ivar;

   /* Get the position in the array */

   int iad = 0;
   if (! _flagAsym)
   {
     iad = ipas;
   }
   else
   {
     if (flag_abs)
     {
       iad = ipas;
     }
     else
     {
       switch (sens)
       {
         case 1:
           iad = _nPas + ipas + 1;
           break;

         case -1:
           iad = _nPas - ipas - 1;
           break;

         case 0:
           iad = _nPas;
           break;
       }
     }
   }
   iad += rank * getLagTotalNumber();
   if (iad < 0 || iad > getSize())
     throw("Debordement");
   return(iad);
 }

double Dir::getGg(int iad) const
{
  if (! isCalculated()) return TEST;
  if (! _isAddressValid(iad)) return TEST;
  return _gg[iad];
}

double Dir::getGg(int ivar, int jvar, int ipas) const
{
  if (!  isCalculated()) return TEST;
  if (! _isLagValid(ipas)) return TEST;
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  return _gg[getAddress(ivar,jvar,ipas,true,0)];
}

VectorDouble Dir::getGg(int ivar, int jvar) const
{
  VectorDouble gg;
  if (! _isVariableValid(ivar)) return gg;
  if (! _isVariableValid(jvar)) return gg;
  int npas = getNPas();
  for (int ipas = 0 ; ipas < npas; ipas++)
  {
    int iad = getAddress(ivar,jvar,ipas,true,0);
    if (_sw[iad] > 0.) gg.push_back(_gg[iad]);
  }
  return gg;
}

double Dir::getHh(int iad) const
{
  if (! isCalculated()) return TEST;
  if (! _isAddressValid(iad)) return TEST;
  return _hh[iad];
}

double Dir::getHh(int ivar, int jvar, int ipas) const
{
  if (!  isCalculated()) return TEST;
  if (! _isLagValid(ipas)) return TEST;
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  return _hh[getAddress(ivar,jvar,ipas,true,0)];
}

VectorDouble Dir::getHh(int ivar, int jvar) const
{
  VectorDouble hh;
  if (! _isVariableValid(ivar)) return hh;
  if (! _isVariableValid(jvar)) return hh;
  int npas = getNPas();
  for (int ipas = 0 ; ipas < npas; ipas++)
  {
    int iad = getAddress(ivar,jvar,ipas,true,0);
    if (_sw[iad] > 0.) hh.push_back(_hh[iad]);
  }
  return hh;
}

double Dir::getSw(int iad) const
{
  if (! isCalculated()) return TEST;
  if (! _isAddressValid(iad)) return TEST;
  return _sw[iad];
}

double Dir::getSw(int ivar, int jvar, int ipas) const
{
  if (!  isCalculated()) return TEST;
  if (! _isLagValid(ipas)) return TEST;
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  return _sw[getAddress(ivar,jvar,ipas,true,0)];
}

VectorDouble Dir::getSw(int ivar, int jvar) const
{
  VectorDouble sw;
  if (! _isVariableValid(ivar)) return sw;
  if (! _isVariableValid(jvar)) return sw;
  int npas = getNPas();
  for (int ipas = 0 ; ipas < npas; ipas++)
  {
    int iad = getAddress(ivar,jvar,ipas,true,0);
    if (_sw[iad] > 0.) sw.push_back(_sw[iad]);
  }
  return sw;
}

void Dir::setSw(int iad, double sw)
 {
   if (! isCalculated()) return;
   if (! _isAddressValid(iad)) return;
   _sw[iad] = sw;
 }

 void Dir::setHh(int iad, double hh)
 {
   if (! isCalculated()) return;
   if (! _isAddressValid(iad)) return;
   _hh[iad] = hh;
 }

 void Dir::setGg(int iad, double gg)
 {
   if (! isCalculated()) return;
   if (! _isAddressValid(iad)) return;
   _gg[iad] = gg;
 }

 void Dir::updSw(int iad, double sw)
 {
   if (! isCalculated()) return;
   if (! _isAddressValid(iad)) return;
   _sw[iad] += sw;
 }

 void Dir::updHh(int iad, double hh)
 {
   if (! isCalculated()) return;
   if (! _isAddressValid(iad)) return;
   _hh[iad] += hh;
 }

 void Dir::updGg(int iad, double gg)
 {
   if (! isCalculated()) return;
   if (! _isAddressValid(iad)) return;
   _gg[iad] += gg;
 }

void Dir::setSw(int ivar, int jvar, int ipas, double sw)
{
  if (! isCalculated()) return;
  if (! _isVariableValid(ivar)) return;
  if (! _isVariableValid(jvar)) return;
  if (! _isLagValid(ipas)) return;
  _sw[getAddress(ivar,jvar,ipas,true,0)] = sw;
}

void Dir::setHh(int ivar, int jvar, int ipas, double hh)
{
  if (! isCalculated()) return;
  if (! _isVariableValid(ivar)) return;
  if (! _isVariableValid(jvar)) return;
  if (! _isLagValid(ipas)) return;
  _hh[getAddress(ivar,jvar,ipas,true,0)] = hh;
}

void Dir::setGg(int ivar, int jvar, int ipas, double gg)
{
  if (! isCalculated()) return;
  if (! _isVariableValid(ivar)) return;
  if (! _isVariableValid(jvar)) return;
  if (! _isLagValid(ipas)) return;
  _gg[getAddress(ivar,jvar,ipas,true,0)] = gg;
}

double Dir::getGrincr(int idim) const
{
  if (_grincr.empty()) return 0.;
  if (! _isDimensionValid(idim)) return 0.;
  return _grincr[idim];
}

void Dir::setUtilize(int iad, double val)
{
  if (! isCalculated()) return;
  if (! _isAddressValid(iad)) return;
  _utilize[iad] = val;
}

String Dir::toString(int level) const
{
  std::stringstream sstr;

  int nvar = getVariableNumber();
  sstr << "Number of lags              = " << getNPas() << std::endl;
  int ndim = getDimensionNumber();

  sstr << toVector("Direction coefficients      = ", _codir);
  if (ndim > 1)
  {
    VectorDouble angles(ndim);
    (void) ut_angles_from_codir(ndim,1,_codir,angles);
    sstr << toVector("Direction angles (degrees)  = ", angles);
  }
  // TODO : Display Grid Incr if so

  if (! FFFF(_tolAngle))
    sstr << "Tolerance on direction      = " << toDouble(_tolAngle)
         << " (degrees)" << std::endl;
  if (! FFFF(_bench)  && _bench > 0.)
    sstr << "Slice bench                 = " << toDouble(_bench) << std::endl;
  if (! FFFF(_cylRad) && _cylRad > 0.)
    sstr << "Slice radius                = " << toDouble(_cylRad) << std::endl;

  if (getLagRegular())
  {
    sstr << "Calculation lag             = " << toDouble(getDPas()) << std::endl;
    sstr << "Tolerance on distance       = " << toDouble(100. * getTolDist())
         << " (\% of the lag value)" << std::endl;
  }
  else
  {
    sstr << "Calculation intervals       = " << std::endl;
    for (int i = 0; i < getBreakNumber(); i++)
    {
      sstr << " - Interval " << i + 1 << " = ["
           << toInterval(getBreaks(i), getBreaks(i + 1)) << "]" << std::endl;
    }
  }

  /* Selection on the 'code' */

  if (getOptionCode() == 1)
    sstr << "Selection if Codes are close enough (" << getTolCode() << ")"
         << std::endl;
  if (getOptionCode() == 2)
    sstr << "Selection if Codes are different" << std::endl;

  /* Print the variogram contents */

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++)
    {
      sstr << std::endl;
      if (ivar == jvar)
        sstr << "For variable " << ivar + 1 << std::endl;
      else
        sstr << "For variables " << ivar + 1 << " and " << jvar + 1 << std::endl;
      sstr << toStr("Rank");
      sstr << toStr("Npairs");
      sstr << toStr("Distance");
      sstr << toStr("Value");
      sstr << std::endl;

      for (int i = 0; i < getLagTotalNumber(); i++)
      {
        int j = getAddress(ivar, jvar, i, true, 0);
        if (_sw[j] <= 0) continue;
        int rank = (! getFlagAsym()) ? i : i - getNPas();
        sstr << toInt(rank);
        sstr << toDouble(_sw[j]);
        sstr << toDouble(_hh[j]);
        sstr << toDouble(_gg[j]);
        sstr << std::endl;
      }
    }
  return sstr.str();
}

std::vector<Dir> generateMultipleDirs(int ndim,
                                      int ndir,
                                      int npas,
                                      double dpas,
                                      double toldis)
{
  VectorDouble angles = VectorDouble(1);
  VectorDouble codir  = VectorDouble(ndim);
  std::vector<Dir> dirs;
  for (int idir = 0; idir < ndir; idir++)
  {
    angles[0] = 180. * (double) idir / (double) ndir;
    (void) ut_angles_to_codir(ndim, 1, angles,codir);
    Dir dir = Dir(ndim, npas, dpas, toldis);
    dir.setTolAngle(90. / (double) ndir);
    dir.setCodir(codir);
    dirs.push_back(dir);
  }
  return dirs;
}


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
#include "Variogram/DirC.hpp"
#include "Db/Db.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"
#include "geoslib_f.h"

DirC::DirC(int nvar, const DirParam dirparam)
    : _dirparam(dirparam),
      _nvar(0),
      _sw(),
      _gg(),
      _hh(),
      _utilize()
{
  _internalResize(nvar);
}

DirC::DirC(int nvar,
           const DirParam dirparam,
           const VectorDouble& gg,
           const VectorDouble& hh,
           const VectorDouble& sw)
    : _dirparam(dirparam),
      _nvar(0),
      _sw(),
      _gg(),
      _hh(),
      _utilize()
{
  _internalResize(nvar);

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

DirC::DirC(const DirC& r)
    : _dirparam(r._dirparam),
      _nvar(r._nvar),
      _sw(r._sw),
      _gg(r._gg),
      _hh(r._hh),
      _utilize(r._utilize)
{
}

DirC& DirC::operator=(const DirC& r)
{
  if (this != &r)
  {
    _dirparam = r._dirparam;
    _nvar = r._nvar;
    _sw = r._sw;
    _gg = r._gg;
    _hh = r._hh;
    _utilize = r._utilize;
  }
  return *this;
}

DirC::~DirC()
{
}

bool DirC::_isAddressValid(int iad) const
{
  if (iad < 0 || iad >= getSize())
  {
    mesArg("Variogram Internal Address",iad,getSize());
    return false;
  }
  return true;
}

bool DirC::_isLagValid(int ilag) const
{
  return _dirparam.isLagValid(ilag);
}

bool DirC::_isDimensionValid(int idim) const
{
  return _dirparam.isDimensionValid(idim);
}

bool DirC::_isVariableValid(int ivar) const
{
  if (ivar < 0 || ivar >= getVariableNumber())
  {
    mesArg("Variable Index",ivar,getVariableNumber());
    return false;
  }
  return true;
}

void DirC::_internalResize(int nvar)
{
  _nvar = nvar;
  int size = getSize();

  _sw.resize(size,0.);
  _gg.resize(size,0.);
  _hh.resize(size,0.);
  _utilize.resize(size,1);
}

double DirC::getHmax(int ivar, int jvar) const
{
  VectorDouble hh = getHh(ivar,jvar);
  return ut_vector_max(hh);
}

double DirC::getGmax(int ivar, int jvar, bool flagAbs) const
{
  VectorDouble gg = getGg(ivar,jvar);
  return ut_vector_max(gg, flagAbs);
}

int DirC::getSize() const
{
  return (_dirparam.getLagTotalNumber() * _nvar * (_nvar + 1) / 2);
}

bool DirC::isCalculated() const
{
  if (_sw.size() <= 0)
  {
    messerr("You must run the Calculations beforehand");
    return false;
  }
  return true;
}

void DirC::copy(const DirC& DirC)
{
  if (getSize() != DirC.getSize()) return;
  _sw = DirC._sw;
  _gg = DirC._gg;
  _hh = DirC._hh;
}

void DirC::clean(void)
{
  for (int i=0; i<getSize(); i++)
  {
    _sw[i] = 0;
    _hh[i] = 0.;
    _gg[i] = 0.;
    _utilize[i] = 1.;
  }
}

int DirC::getAddress(int ivar, int jvar, int ipas, bool flag_abs, int sens) const
{
   int rank;

   /* Get the order of the variables */

   if (ivar > jvar)
     rank = ivar * (ivar+1)/2 + jvar;
   else
     rank = jvar * (jvar+1)/2 + ivar;

   /* Get the position in the array */

   int iad = 0;
   if (! _dirparam.getFlagAsym())
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
       int npas = _dirparam.getNPas();
       switch (sens)
       {
         case 1:
           iad = npas + ipas + 1;
           break;

         case -1:
           iad = npas - ipas - 1;
           break;

         case 0:
           iad = npas;
           break;
       }
     }
   }
   iad += rank * _dirparam.getLagTotalNumber();
   if (iad < 0 || iad > getSize())
     throw("Debordement");
   return(iad);
 }

double DirC::getGg(int iad) const
{
  if (! isCalculated()) return TEST;
  if (! _isAddressValid(iad)) return TEST;
  return _gg[iad];
}

double DirC::getGg(int ivar, int jvar, int ipas) const
{
  if (!  isCalculated()) return TEST;
  if (! _isLagValid(ipas)) return TEST;
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  return _gg[getAddress(ivar,jvar,ipas,true,0)];
}

VectorDouble DirC::getGg(int ivar, int jvar) const
{
  VectorDouble gg;
  if (! _isVariableValid(ivar)) return gg;
  if (! _isVariableValid(jvar)) return gg;
  for (int ipas = 0 ; ipas < _dirparam.getNPas(); ipas++)
  {
    int iad = getAddress(ivar,jvar,ipas,true,0);
    if (_sw[iad] > 0.) gg.push_back(_gg[iad]);
  }
  return gg;
}

double DirC::getHh(int iad) const
{
  if (! isCalculated()) return TEST;
  if (! _isAddressValid(iad)) return TEST;
  return _hh[iad];
}

double DirC::getHh(int ivar, int jvar, int ipas) const
{
  if (!  isCalculated()) return TEST;
  if (! _isLagValid(ipas)) return TEST;
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  return _hh[getAddress(ivar,jvar,ipas,true,0)];
}

VectorDouble DirC::getHh(int ivar, int jvar) const
{
  VectorDouble hh;
  if (! _isVariableValid(ivar)) return hh;
  if (! _isVariableValid(jvar)) return hh;
  for (int ipas = 0 ; ipas < _dirparam.getNPas(); ipas++)
  {
    int iad = getAddress(ivar,jvar,ipas,true,0);
    if (_sw[iad] > 0.) hh.push_back(_hh[iad]);
  }
  return hh;
}

double DirC::getSw(int iad) const
{
  if (! isCalculated()) return TEST;
  if (! _isAddressValid(iad)) return TEST;
  return _sw[iad];
}

double DirC::getSw(int ivar, int jvar, int ipas) const
{
  if (!  isCalculated()) return TEST;
  if (! _isLagValid(ipas)) return TEST;
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  return _sw[getAddress(ivar,jvar,ipas,true,0)];
}

VectorDouble DirC::getSw(int ivar, int jvar) const
{
  VectorDouble sw;
  if (! _isVariableValid(ivar)) return sw;
  if (! _isVariableValid(jvar)) return sw;
  for (int ipas = 0 ; ipas < _dirparam.getNPas(); ipas++)
  {
    int iad = getAddress(ivar,jvar,ipas,true,0);
    if (_sw[iad] > 0.) sw.push_back(_sw[iad]);
  }
  return sw;
}

int DirC::getCenter(int ivar, int jvar) const
{
  int center = ITEST;
  if (!_isVariableValid(ivar)) return center;
  if (!_isVariableValid(jvar)) return center;

  if (! _dirparam.getFlagAsym()) return -1;

  int nval = 0;
  for (int i = 0; i < _dirparam.getLagTotalNumber(); i++)
  {
    double ww = getSw(ivar, jvar, i);
    if (i == _dirparam.getNPas()) return (nval + 1);
    if (ww <= 0) continue;
  }
  return center;
}

void DirC::setSw(int iad, double sw)
 {
   if (! isCalculated()) return;
   if (! _isAddressValid(iad)) return;
   _sw[iad] = sw;
 }

 void DirC::setHh(int iad, double hh)
 {
   if (! isCalculated()) return;
   if (! _isAddressValid(iad)) return;
   _hh[iad] = hh;
 }

 void DirC::setGg(int iad, double gg)
 {
   if (! isCalculated()) return;
   if (! _isAddressValid(iad)) return;
   _gg[iad] = gg;
 }

 void DirC::updSw(int iad, double sw)
 {
   if (! isCalculated()) return;
   if (! _isAddressValid(iad)) return;
   _sw[iad] += sw;
 }

 void DirC::updHh(int iad, double hh)
 {
   if (! isCalculated()) return;
   if (! _isAddressValid(iad)) return;
   _hh[iad] += hh;
 }

 void DirC::updGg(int iad, double gg)
 {
   if (! isCalculated()) return;
   if (! _isAddressValid(iad)) return;
   _gg[iad] += gg;
 }

void DirC::setSw(int ivar, int jvar, int ipas, double sw)
{
  if (! isCalculated()) return;
  if (! _isVariableValid(ivar)) return;
  if (! _isVariableValid(jvar)) return;
  if (! _isLagValid(ipas)) return;
  _sw[getAddress(ivar,jvar,ipas,true,0)] = sw;
}

void DirC::setHh(int ivar, int jvar, int ipas, double hh)
{
  if (! isCalculated()) return;
  if (! _isVariableValid(ivar)) return;
  if (! _isVariableValid(jvar)) return;
  if (! _isLagValid(ipas)) return;
  _hh[getAddress(ivar,jvar,ipas,true,0)] = hh;
}

void DirC::setGg(int ivar, int jvar, int ipas, double gg)
{
  if (! isCalculated()) return;
  if (! _isVariableValid(ivar)) return;
  if (! _isVariableValid(jvar)) return;
  if (! _isLagValid(ipas)) return;
  _gg[getAddress(ivar,jvar,ipas,true,0)] = gg;
}

void DirC::setUtilize(int iad, double val)
{
  if (! isCalculated()) return;
  if (! _isAddressValid(iad)) return;
  _utilize[iad] = val;
}

String DirC::toString(int level) const
{
  std::stringstream sstr;

  sstr << _dirparam.toString(level);

  /* Print the variogram contents */

  int nvar = getVariableNumber();
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

      for (int i = 0; i < _dirparam.getLagTotalNumber(); i++)
      {
        int j = getAddress(ivar, jvar, i, true, 0);
        if (_sw[j] <= 0) continue;
        int rank = (! _dirparam.getFlagAsym()) ? i : i - _dirparam.getNPas();
        sstr << toInt(rank);
        sstr << toDouble(_sw[j]);
        sstr << toDouble(_hh[j]);
        sstr << toDouble(_gg[j]);
        sstr << std::endl;
      }
    }
  return sstr.str();
}

/**
 * Set the Central Value (distance zero) for an asymmetric function
 * @param nech : Number of valid samples (will provide an information on the number of pairs)
 * @param rho  : Correlation between pairs of (different) variables
 */
void DirC::patchCenter(int nech, double rho)
{
  if (! _dirparam.getFlagAsym()) return;
  for (int ivar=0; ivar<_nvar; ivar++)
    for (int jvar=0; jvar<=ivar; jvar++)
    {
      // Get the central address
      int iad = getAddress(ivar,jvar,0,false,0);
      setSw(iad,nech);
      setHh(iad,0.);
      if (ivar == jvar)
        setGg(iad,1.);
      else
        setGg(iad,rho);
    }
}

int DirC::fill(int nvar,
              const VectorDouble& sw,
              const VectorDouble& gg,
              const VectorDouble& hh)
{
  _internalResize(nvar);
  int size = getSize();
  if (size != (int) sw.size() ||
      size != (int) hh.size() ||
      size != (int) gg.size())
  {
    messerr("The argument do not have correct dimension");
    return 1;
  }
  for (int i=0; i<size; i++)
  {
    setSw(i, sw[i]);
    setHh(i, hh[i]);
    setGg(i, gg[i]);
  }
  return 0;
}

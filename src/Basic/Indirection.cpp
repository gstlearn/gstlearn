/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Basic/Indirection.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AStringable.hpp"

Indirection::Indirection(int mode)
    : AStringable(),
      _defined(false),
      _mode(mode),
      _nabs(0),
      _nrel(0),
      _vecRToA(),
      _vecAToR(),
      _mapAToR()
{
}

Indirection::Indirection(const Indirection &m)
    : AStringable(m),
      _defined(m._defined),
      _mode(m._mode),
      _nabs(m._nabs),
      _nrel(m._nrel),
      _vecRToA(m._vecRToA),
      _vecAToR(m._vecAToR),
      _mapAToR(m._mapAToR)
{
}

Indirection& Indirection::operator=(const Indirection &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _defined = m._defined;
    _mode = m._mode;
    _nabs = m._nabs;
    _nrel = m._nrel;
    _vecRToA = m._vecRToA;
    _vecAToR = m._vecAToR;
    _mapAToR = m._mapAToR;
  }
  return *this;
}

Indirection::~Indirection()
{
}

String Indirection::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  if (_mode == 0)
    sstr << "- Information (AToR) is stored as Array of Integers" << std::endl;
  else
    sstr << "- Information (AToR) is stored as a Map" << std::endl;
  sstr << "- Number of absolute positions = " <<  _nabs << std::endl;
  sstr << "- Number of active positions   = " <<  _nrel << std::endl;

  return sstr.str();
}

void Indirection::_resetMap()
{
  _mapAToR.clear();
  _vecAToR.clear();
  _vecRToA.clear();
}

void Indirection::setMode(int mode)
{
  if (mode == _mode) return;
  _resetMap();
  _mode = mode;
}

/**
 * Build the needed information from Selection array
 * A sample is active if its 'sel' value is equal to 1
 * @param sel Vector giving the status of all samples (Dimension: absolute)
 */
void Indirection::buildFromSel(const VectorDouble& sel)
{
  _resetMap();
  _nabs = (int) sel.size();
  if (_mode == 0) _vecAToR.resize(_nabs,-1);

  int irel = 0;
  for (int iabs = 0; iabs < _nabs; iabs++)
  {
    if (! sel[iabs]) continue;
    if (_mode == 0)
      _vecAToR[iabs] = irel;
    else
      _mapAToR[iabs] = irel;
    _vecRToA.push_back(iabs);
    irel++;
  }
  _nrel = irel;
  _defined = true;
}

void Indirection::buildFromMap(const std::map<int, int> &map, int nabs)
{
  _resetMap();
  _nabs = nabs;
  _nrel = (int) map.size();

  if (_mode == 0)
    _vecAToR.resize(_nabs, -1);
  else
    _mapAToR = map;
  _vecRToA.resize(_nrel, -1);

  for (auto it = map.begin(); it != map.end(); it++)
  {
    _vecRToA[it->second]  = it->first;
    if (_mode == 0) _vecAToR[it->first] = it->second;
  }
  _defined = true;
}

void Indirection::buildFromRankRInA(const VectorInt& rels, int nabs)
{
  _resetMap();
  _nabs = nabs;
  _nrel = (int) rels.size();

  if (_mode == 0) _vecAToR.resize(_nabs,-1);

  int iabs;
  for (int irel = 0; irel < _nrel; irel++)
  {
    iabs = rels[irel];
    if (_mode == 0)
      _vecAToR[iabs] = irel;
    else
      _mapAToR[iabs] = irel;
  }
  _vecRToA = rels;
  _defined = true;
}

int Indirection::getAToR(int iabs) const
{
  if (_mode == 0)
  {
    return _getArrayAToR(iabs);
  }
  else
  {
    return _getMapAToR(iabs);
  }
}

int Indirection::getRToA(int irel) const
{
  if (_vecRToA.empty()) return irel;
  if (! _isValidRel(irel)) return ITEST;
  return _vecRToA[irel];
}

bool Indirection::_isValidAbs(int iabs) const
{
  if (iabs < 0 || iabs >= getAbsSize())
  {
    mesArg("Absolute Rank", iabs, getAbsSize());
    return false;
  }
  return true;
}

bool Indirection::_isValidRel(int irel) const
{
  if (irel < 0 || irel >= getRelSize())
  {
    mesArg("Relative Rank", irel, getRelSize());
    return false;
  }
  return true;
}

int Indirection::_getArrayAToR(int iabs) const
{
  if (_vecAToR.empty()) return iabs;
  if (! _isValidAbs(iabs)) return ITEST;
  return _vecAToR[iabs];
}

/**
 * Returns the rank of the relative grid node from its absolute index using the Map
 * @param iabs Absolute rank of the grid node
 * @return Rank of the corresponding active (relative) grid node (or -1 is not found)
 */
int Indirection::_getMapAToR(int iabs) const
{
  if (_mapAToR.empty()) return iabs;
  if (_mapAToR.find(iabs) == _mapAToR.end())
    return -1;
  else
  {
    int irel = _mapAToR.find(iabs)->second;
    return irel;
  }
}

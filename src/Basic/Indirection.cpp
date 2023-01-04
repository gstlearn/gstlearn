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
#include "Basic/Indirection.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AStringable.hpp"

Indirection::Indirection(int mode)
    : _mode(mode),
      _nabs(0),
      _nrel(0),
      _AToR(),
      _RToA(),
      _map()
{
}

Indirection::Indirection(const Indirection &m)
    : _mode(m._mode),
      _nabs(m._nabs),
      _nrel(m._nrel),
      _AToR(m._AToR),
      _RToA(m._RToA),
      _map(m._map)
{
}

Indirection& Indirection::operator=(const Indirection &m)
{
  if (this != &m)
  {
    _mode = m._mode;
    _nabs = m._nabs;
    _nrel = m._nrel;
    _AToR = m._AToR;
    _RToA = m._RToA;
    _map  = m._map;
  }
  return *this;
}

Indirection::~Indirection()
{
}

void Indirection::buildFromSel(const VectorDouble& sel, bool verbose)
{
  _map.clear();
  _AToR.clear();
  _RToA.clear();

  if (_mode == 0)
    _buildArrays(sel, verbose);
  else
    _buildMap(sel, verbose);
}

/**
 * Build the list of absolute indices for the only active samples
 * A sample is active if its 'sel' value is equal to 1
 * @param sel Vector giving the status of all samples (Dimension: absolute)
 * @return
 */
void Indirection::_buildArrays(const VectorDouble& sel, bool verbose)
{
  _nabs = (int) sel.size();
  _AToR.resize(_nabs,-1);
  int irel = 0;
  for (int iabs = 0; iabs < _nabs; iabs++)
  {
    if (! sel[iabs]) continue;
    _AToR[iabs] = irel;
    _RToA.push_back(iabs);
    irel++;
  }
  _nrel = irel;

  _mode = 0;
  if (verbose) _printInfo();
}

void Indirection::buildFromMap(const std::map<int, int> &map,
                               int nabs,
                               bool verbose)
{
  _map.clear();
  _AToR.clear();
  _RToA.clear();

  if (_mode == 0)
    _buildArraysFromMap(map, nabs, verbose);
  else
    _setMap(map);
}

void Indirection::_buildArraysFromMap(const std::map<int, int> &map,
                                     int nabs,
                                     bool verbose)
{
  _nabs = nabs;
  _nrel = (int) map.size();
  _AToR.resize(_nabs, -1);
  _RToA.resize(_nrel, -1);

  for (auto it = map.begin(); it != map.end(); it++)
  {
    _AToR[it->first] = it->second;
    _RToA[it->second] = it->first;
  }

  _mode = 0;
  if (verbose) _printInfo();
}

void Indirection::_setMap(const std::map<int, int> &map)
{
  _map = map;
  _mode = 1;
}

void Indirection::_printInfo() const
{
  message("Indexing between Absolute to Relative\n");
  if (_mode == 0)
    message("- Using information stored in Arrays of Integers\n");
  else
    message("- Using information stored in a Map\n");
  message("- Number of absolute positions = %d\n", _nabs);
  message("- Number of active positions   = %d\n", _nrel);
}

/**
 * Creates the map such that MAP[iabs] = iact.
 * A sample is active if its 'sel' value is equal to 1
 * @param sel Vector giving the status of all samples (Dimension: absolute)
 * @param verbose Verbose flag
 */
void Indirection::_buildMap(const VectorDouble& sel, bool verbose)
{
  _map.clear();
  _AToR.clear();
  _RToA.clear();

  _nabs = (int) sel.size();
  int irel   = 0;
  for (int iabs = 0; iabs < _nabs; iabs++)
  {
    if (sel[iabs] == 0) continue;
    _map[iabs] = irel++;
  }

  // Optional control printout
  _mode = 1;
  if (verbose) _printInfo();
}

int Indirection::getAToR(int iabs) const
{
  if (_mode == 0)
  {
    return _getArrayAbsoluteToRelative(iabs);
  }
  else
  {
    return _getMapAbsoluteToRelative(iabs);
  }
}

int Indirection::getRtoA(int irel) const
{
  if (_mode == 0)
  {
    return _getArrayRelativeToAbsolute(irel);
  }
  else
  {
    return _getMapRelativeToAbsolute(irel);
  }
}

bool Indirection::_isValidAbs(int iabs) const
{
  if (iabs < 0 || iabs >= _getAbsSize())
  {
    mesArg("Absolute Rank", iabs, _getAbsSize());
    return false;
  }
  return true;
}

bool Indirection::_isValidRel(int irel) const
{
  if (irel < 0 || irel >= _getRelSize())
  {
    mesArg("Relative Rank", irel, _getRelSize());
    return false;
  }
  return true;
}

int Indirection::_getArrayAbsoluteToRelative(int iabs) const
{
  if (_AToR.empty()) return iabs;
  if (! _isValidAbs(iabs)) return ITEST;
  return _AToR[iabs];
}

int Indirection::_getArrayRelativeToAbsolute(int irel) const
{
  if (_RToA.empty()) return irel;
  if (! _isValidRel(irel)) return ITEST;
  return _RToA[irel];
}

/**
 * Returns the rank of the relative grid node from its absolute index using the Map
 * @param iabs Absolute rank of the grid node
 * @return Rank of the corresponding active (relative) grid node (or -1 is not found)
 */
int Indirection::_getMapAbsoluteToRelative(int iabs) const
{
  if (_map.empty()) return iabs;
  if (_map.find(iabs) == _map.end())
    return -1;
  else
    return _map.find(iabs)->second;
}

int Indirection::_getMapRelativeToAbsolute(int irel) const
{
  if (_map.empty()) return irel;
  auto it = _map.begin();
  std::advance(it, irel);
  return it->first;
}

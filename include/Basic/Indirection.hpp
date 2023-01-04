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
#include "geoslib_define.h"

#include "Basic/VectorNumT.hpp"

#include <map>


class GSTLEARN_EXPORT Indirection {
public:
  Indirection(int mode = 0);
  ~Indirection();
  Indirection(const Indirection &m);
  Indirection& operator=(const Indirection &m);

  void buildFromSel(const VectorDouble& sel, bool verbose = false);
  void buildFromMap(const std::map<int, int> &map,
                    int nabs,
                    bool verbose = false);
  int getAToR(int iabs) const;
  int getRtoA(int irel) const;

private:
  void _buildArrays(const VectorDouble& sel, bool verbose = false);
  void _buildMap(const VectorDouble &sel, bool verbose = false);
  void _buildArraysFromMap(const std::map<int, int> &map,
                           int nabs,
                           bool verbose = false);
  void _setMap(const std::map<int, int> &map);
  int _getMapAbsoluteToRelative(int iabs) const;
  int _getMapRelativeToAbsolute(int irel) const;
  int _getArrayAbsoluteToRelative(int iabs) const;
  int _getArrayRelativeToAbsolute(int irel) const;
  int _getAbsSize() const { return _nabs; }
  int _getRelSize() const { return _nrel; }
  bool _isValidAbs(int iabs) const;
  bool _isValidRel(int irel) const;
  void _printInfo() const;

private:
  int  _mode;     // 0 by array; 1 by MAP
  int  _nabs;     // Number of absolute elements
  int  _nrel;     // Number of relative elements
  VectorInt _AToR;
  VectorInt _RToA;
  std::map<int, int> _map;
};

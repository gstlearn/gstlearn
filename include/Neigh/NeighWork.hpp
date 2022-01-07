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
#include "Neigh/Neigh.hpp"
#include "geoslib_define.h"

class Neigh;
class Db;

class GSTLEARN_EXPORT NeighWork
{
public:
  NeighWork(const Neigh* neigh,
            const Db* dbin,
            bool flag_var_nocheck = false,
            bool flag_simu = false);
  NeighWork(const NeighWork& r);
  NeighWork& operator=(const NeighWork& r);
  virtual ~NeighWork();

  void initialize(const Neigh* neigh,
                  const Db* dbin,
                  bool flag_var_nocheck = false,
                  bool flag_simu = false);
  void clear();
  VectorInt select(Db *dbout, int iech_out);
  bool isInitialized() const { return _flagInitialized; }

private:
  void _unique(Db *dbout, int iech_out, VectorInt& ranks);
  void _bench(Db *dbout, int iech_out, VectorInt& ranks);
  int  _moving(Db *dbout, int iech_out, VectorInt& ranks, double eps = EPSILON9);
  bool _discardUndefined(int iech);
  int  _xvalid(Db *dbout, int iech_in, int iech_out, double eps = EPSILON9);
  int  _movingSectorDefine(double dx, double dy);
  void _movingSectorNsmax(int nsel, VectorInt& ranks);
  void _movingSelect(int nsel, VectorInt& ranks);
  void _display(const VectorInt& ranks);
  double _movingDist(Db *dbout, int iech_in, int iech_out);

private:
  const Neigh* _neigh;
  const Db* _dbin;
  bool _flagInitialized;
  mutable VectorInt _nbghInd;
  mutable VectorInt _nbghIsect;
  mutable VectorInt _nbghNsect;
  mutable VectorDouble _nbghX1;
  mutable VectorDouble _nbghX2;
  mutable VectorDouble _nbghDst;
  bool _flagVariableNoCheck;
  bool _flagSimu;
};

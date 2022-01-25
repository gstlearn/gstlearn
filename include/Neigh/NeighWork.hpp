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
#include "Neigh/ANeighParam.hpp"
#include "geoslib_define.h"

class Neightobedeleted;
class Db;

class GSTLEARN_EXPORT NeighWork
{
public:
  NeighWork(const Db* dbin = nullptr,
            const ANeighParam* neighparam = nullptr,
            bool flag_simu = false);
  NeighWork(const NeighWork& r);
  NeighWork& operator=(const NeighWork& r);
  virtual ~NeighWork();

  void initialize(const Db* dbin,
                  const ANeighParam* neighparam,
                  bool flag_simu = false);
  void clear();
  VectorInt select(Db *dbout,
                   int iech_out,
                   const VectorInt& rankColCok = VectorInt(),
                   bool verbose = false);
  bool isUnchanged() const { return _flagIsUnchanged; }
  VectorDouble summary(Db *dbout,
                       int iech_out,
                       const VectorInt& rankColCok = VectorInt());

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
  void _checkUnchanged(const Db* dbout, int iech_out, const VectorInt& ranks);
  void _clearMemory();
  void _resetFromMemory(bool flagSame, VectorInt& ranks, bool verbose);
  bool _isSameTarget(const Db* dbout,
                     int iech_out,
                     VectorInt& ranks,
                     bool verbose = false);
  bool _isSameTargetBench(const Db* dbout,
                          int iech_out,
                          VectorInt& ranks,
                          bool verbose = false);
  bool _isSameTargetUnique(const Db* dbout,
                           int iech_out,
                           VectorInt& ranks,
                           bool verbose = false);
  void _updateColCok(const VectorInt& rankColCok, VectorInt& ranks);

private:
  const Db* _dbin;
  const ANeighParam* _neighParam;
  bool _flagInitialized;
  bool _flagIsUnchanged;
  mutable VectorInt    _movingInd;
  mutable VectorInt    _movingIsect;
  mutable VectorInt    _movingNsect;
  mutable VectorDouble _movingX1;
  mutable VectorDouble _movingX2;
  mutable VectorDouble _movingDst;
  bool _flagSimu;

  // Following parameters are only kept for optimization
  const Db* _dbout;
  int _iechOut;
  mutable VectorInt _nbghMemo;
};

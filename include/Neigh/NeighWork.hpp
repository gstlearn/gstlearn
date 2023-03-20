/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Neigh/ANeighParam.hpp"
#include "geoslib_define.h"

class Db;

class GSTLEARN_EXPORT NeighWork
{
public:
  NeighWork(const Db* dbin = nullptr,
            const ANeighParam* neighparam = nullptr);
  NeighWork(const NeighWork& r);
  NeighWork& operator=(const NeighWork& r);
  virtual ~NeighWork();

  void initialize(const Db* dbin,
                  const ANeighParam* neighparam);
  void clear();
  VectorInt select(Db *dbout,
                   int iech_out,
                   const VectorInt& rankColCok = VectorInt(),
                   bool verbose = false);
  bool isUnchanged() const { return _flagIsUnchanged; }
  void setIsChanged();
  VectorDouble summary(Db *dbout,
                       int iech_out,
                       const VectorInt& rankColCok = VectorInt());
  void setFlagSimu(bool flagSimu) { _flagSimu = flagSimu; }

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
  bool _belongsToCell(Db* dbout, int iech, int iech_out);
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
  void _updateColCok(const VectorInt& rankColCok, VectorInt& ranks, int iech_out);;
  bool _hiddenByFault(Db* dbout, int iech, int iech_out) const;

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

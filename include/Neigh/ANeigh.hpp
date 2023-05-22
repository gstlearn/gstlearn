/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Enum/ENeigh.hpp"

#include "Basic/ASerializable.hpp"
#include "Space/ASpaceObject.hpp"
#include "geoslib_define.h"

class Db;

class GSTLEARN_EXPORT ANeigh:  public ASpaceObject, public ASerializable
{
public:
  ANeigh(const Db *dbin = nullptr,
         const Db *dbout = nullptr,
         const ASpace* space = nullptr);
  ANeigh(const ANeigh& r);
  ANeigh& operator=(const ANeigh& r);
  virtual ~ANeigh();

  /// ASpaceObject Interface
  virtual bool isConsistent(const ASpace* space) const override { return true; }

  // AStringable Interface overriding
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface for ANeigh
  virtual int attach(const Db *dbin, const Db *dbout);
  virtual VectorInt getNeigh(int iech_out) = 0;
  virtual int getMaxSampleNumber(const Db* db) const = 0;
  virtual bool hasChanged(int iech_out) const { return true; }
  virtual VectorDouble summary(int iech_out) { return VectorDouble(); }
  virtual ENeigh getType() const { return ENeigh::fromKey("UNKNOWN"); }
  virtual bool getFlagContinuous() const { return false; }

  VectorInt select(int iech_out);
  bool isUnchanged() const { return _flagIsUnchanged; }
  void setIsChanged();
  void reset();

  bool getFlagXvalid() const { return _flagXvalid; }
  bool getFlagKFold() const { return _flagKFold; }

  void setFlagXvalid(bool flagXvalid) { _flagXvalid = flagXvalid; }
  void setFlagKFold(bool flagKFold)   { _flagKFold = flagKFold; }
  void setFlagSimu(bool flagSimu)     { _flagSimu = flagSimu; }
  void setRankColCok(const VectorInt &rankColCok) { _rankColCok = rankColCok; }

protected:
  bool _isNbghMemoEmpty() const { return _nbghMemo.empty(); }
  void _neighCompress(VectorInt& ranks);
  void _display(const VectorInt& ranks);
  bool _discardUndefined(int iech);
  int  _xvalid(int iech_in, int iech_out, double eps = EPSILON9);
  bool _isDimensionValid(int idim) const;

  // Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "ANeigh"; }

private:
  bool _isSameTarget(int iech_out);
  void _checkUnchanged(int iech_out, const VectorInt &ranks);
  void _updateColCok(VectorInt& ranks, int iech_out);

protected:
  const Db* _dbin;
  const Db* _dbout;
  const DbGrid* _dbgrid; // Equivalent to dbout, defined only for grid

  VectorInt _rankColCok;
  int  _iechMemo;
  bool _flagSimu;
  bool _flagXvalid;              /* True to suppress the target */
  bool _flagKFold;               /* True to perform a KFold Cross-validation */

private:
  bool _flagIsUnchanged;
  mutable VectorInt  _nbghMemo;
};

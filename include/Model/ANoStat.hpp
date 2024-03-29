/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Model/CovParamId.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/AStringable.hpp"
#include "Db/Db.hpp"
#include "Mesh/AMesh.hpp"

#include <vector>

class Model;
class ACov;

/**
 * \brief
 * NoStat is defined as a series of codes. Each code follows the following syntax with three parts:
 *
 *     G<igrf>M<icov>R<idir>
 *
 * Example:
 *
 *     G1M0R2: Range in the third direction, for the first covariance of the second GRF
 *     M1A: First anisotropy direction angle, for the second covariance of the first GRF
 *
 * Notes:
 * - All numbering is 0-based
 * - First two parts can be omitted: we consider G0 and M0
 * - The keyword for the last part is chosen among:
 *
 *             R : Range
 *             A : Angle
 *             P : Third parameter
 *             V : Sill
 *             S : Scale
 *             T : Tapering range
 *             C : Velocity (advection)
 *             I : Rotation angle for Sphere
 *             H : Anisotropy matrix terms
 */
class GSTLEARN_EXPORT ANoStat : public AStringable, public ICloneable
{
public:
  ANoStat();
  ANoStat(const VectorString& codes);
  ANoStat(const ANoStat &m);
  ANoStat& operator= (const ANoStat &m);
  virtual ~ANoStat();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  bool isNotEmpty() const { return !_items.empty(); }
  bool isDefinedByCov(int icov = -1, int igrf = -1) const;
  bool isDefinedByType(const EConsElem& type, int igrf = -1) const;
  bool isDefinedByCovType(const EConsElem& type, int icov = -1, int igrf = -1) const;
  bool isDefined(const EConsElem &type,
                 int icov = 0,
                 int iv1 = -1,
                 int iv2 = -1,
                 int igrf = -1) const;
  bool isDefinedforAnisotropy(int icov = -1, int igrf = -1) const;
  bool isDefinedforRotation(int icov = -1, int igrf = -1) const;

  virtual double getValue(const EConsElem &type,
                          int icas,
                          int rank,
                          int icov = 0,
                          int iv1 = -1,
                          int iv2 = -1,
                          int igrf = -1) const = 0;
  virtual double getValueByParam(int ipar, int icas, int rank) const = 0;

  virtual int  attachToMesh(const AMesh* mesh, bool verbose = false) const;
  virtual void detachFromMesh() const;
  virtual int  attachToDb(Db* db, int icas, bool verbose = false) const;
  virtual void detachFromDb(Db* db, int icas) const;

  int  manageInfo(int mode, Db *dbin, Db *dbout);

  int  addNoStatElem(int igrf, int icov, const EConsElem& type, int iv1, int iv2);
  int  addNoStatElemByItem(const CovParamId& item);
  int  addNoStatElems(const VectorString& codes);
  void deleteNoStatElem(int ipar);
  void deleteAllNoStatElem();

  int getRank(const EConsElem &type,
              int icov,
              int iv1 = -1,
              int iv2 = -1,
              int igrf = -1) const;
  int getIGrf(int ipar) const { return _items[ipar].getIGrf(); }
  int getICov(int ipar) const { return _items[ipar].getICov(); }
  const EConsElem& getType(int ipar) const { return _items[ipar].getType(); }
  int getIV1 (int ipar) const { return _items[ipar].getIV1(); }
  int getIV2 (int ipar) const { return _items[ipar].getIV2(); }
  int getNoStatElemNumber() const { return static_cast<int>(_items.size()); }
  const std::vector<CovParamId>& getNoStats() const { return _items; }
  const CovParamId getNoStat(int ipar) const { return _items[ipar]; }

  int attachModel(const Model* model);

  bool matchIGrf(int ipar, int igrf0) const { return _items[ipar].matchIGrf(igrf0); }
  bool matchICov(int ipar, int icov0) const { return _items[ipar].matchICov(icov0); }
  bool matchType(int ipar, const EConsElem& type0) const { return _items[ipar].matchType(type0); }
  bool matchIV1(int ipar, int iv10) const { return _items[ipar].matchIV1(iv10); }
  bool matchIV2(int ipar, int iv20) const { return _items[ipar].matchIV2(iv20); }

  const std::vector<CovParamId>& getAllItems() const { return _items; }
  const CovParamId getItems(int ipar) const { return _items[ipar]; }

  bool getInfoFromDb(int ipar,
                     int icas1,
                     int iech1,
                     int icas2,
                     int iech2,
                     double *val1,
                     double *val2) const;

protected:
  void _setAmesh(const AMesh* amesh) const { _amesh = amesh; }
  void _setDbin(const Db* dbin) const { _dbin = dbin; }
  void _setDbout(const Db* dbout) const { _dbout = dbout; }
  bool _isValid(int icas, int rank) const;

private:
  int _understandCode(const String& code,
                      int *igrf,
                      int *icov,
                      EConsElem *type,
                      int *iv1,
                      int *iv2);
  void _updateFromModel(const Model* model);
  bool _checkConsistency() const;

private:
  std::vector<CovParamId> _items;

protected:
  // The following arguments are stored as pointer to ease communication
  // Their list is established as large as possible (even if all of them are not actually used)
  mutable const AMesh* _amesh;
  mutable const Db*    _dbin;
  mutable const Db*    _dbout;
};

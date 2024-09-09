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
#include "Covariances/ParamId.hpp"
#include "Basic/AStringable.hpp"
#include "Db/Db.hpp"
#include "Mesh/AMesh.hpp"

#include <vector>

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
 *
 * - In the Multivariate case, use "V1-2" for the sill of the cross-variogram between variables 1 and 2
 *   Pay attention: "V1" is equivalent to "V1-1" and "V2" is equivalent to "V2-1" (not to "V2-2").
 */
class GSTLEARN_EXPORT ANoStatCov : public AStringable, public ICloneable
{
public:
  ANoStatCov();
  ANoStatCov(const VectorString& codes);
  ANoStatCov(const ANoStatCov &m);
  ANoStatCov& operator= (const ANoStatCov &m);
  virtual ~ANoStatCov();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  bool isNotEmpty() const { return !_items.empty(); }
  bool isDefinedByCov() const;
  bool isDefinedByType(const EConsElem& type) const;
  bool isDefined(const EConsElem &type,
                 int iv1 = -1,
                 int iv2 = -1) const;
  bool isDefinedForVariance() const;
  bool isDefinedforAnisotropy() const;
  bool isDefinedforRotation() const;

  virtual double getValue(const EConsElem &type,
                          int icas,
                          int rank,
                          int iv1 = -1,
                          int iv2 = -1) const = 0;
  virtual double getValueByParam(int ipar, int icas, int rank) const = 0;

  virtual int  attachToMesh(const AMesh* mesh, bool center = true, bool verbose = false) const;
  virtual void detachFromMesh() const;
  virtual int  attachToDb(Db* db, int icas, bool verbose = false) const;
  virtual void detachFromDb(Db* db, int icas) const;

  int  manageInfo(int mode, Db *dbin, Db *dbout) const;

  int  addNoStatElem(const EConsElem& type, int iv1, int iv2);
  int  addNoStatElemByItem(const ParamId& item);
  int  addNoStatElems(const VectorString& codes);
  void deleteNoStatElem(int ipar);
  void deleteAllNoStatElem();

  int getRank(const EConsElem &type,
              int iv1 = -1,
              int iv2 = -1) const;
  const EConsElem& getType(int ipar) const { return _items[ipar].getType(); }
  int getIV1 (int ipar) const { return _items[ipar].getIV1(); }
  int getIV2 (int ipar) const { return _items[ipar].getIV2(); }
  int getNoStatElemNumber() const { return static_cast<int>(_items.size()); }
  const std::vector<ParamId>& getNoStats() const { return _items; }
  ParamId getNoStat(int ipar) const { return _items[ipar]; }

  int attachCova(const ACov* cova);

  bool matchType(int ipar, const EConsElem& type0) const { return _items[ipar].matchType(type0); }
  bool matchIV1(int ipar, int iv10) const { return _items[ipar].matchIV1(iv10); }
  bool matchIV2(int ipar, int iv20) const { return _items[ipar].matchIV2(iv20); }

  const std::vector<ParamId>& getAllItems() const { return _items; }
  ParamId getItems(int ipar) const { return _items[ipar]; }

  bool getInfoFromDb(int ipar,
                     int icas1,
                     int iech1,
                     int icas2,
                     int iech2,
                     double *val1,
                     double *val2) const;

  static void checkCode(const String& code);

protected:
  void _setAmesh(const AMesh* amesh) const { _amesh = amesh; }
  void _setDbin(const Db* dbin) const { _dbin = dbin; }
  void _setDbout(const Db* dbout) const { _dbout = dbout; }
  bool _isValid(int icas, int rank) const;
  const Db* _getDbin() const { return _dbin; }
  const Db* _getDbout() const;
  const AMesh* _getAMesh() const { return _amesh; }

private:
  static int _understandCode(const String& code,
                             EConsElem* type,
                             int* iv1,
                             int* iv2);
  bool _checkConsistency() const;

private:
  std::vector<ParamId> _items;
  mutable const AMesh* _amesh;
  mutable const Db*    _dbin;
  mutable const Db*    _dbout;
};

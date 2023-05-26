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
#include "geoslib_define.h"

#include "Enum/ENeigh.hpp"

#include "Neigh/ANeigh.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/Utilities.hpp"
#include "Geometry/ABiPointCheck.hpp"
#include "Geometry/BiPointCheckDistance.hpp"
#include "Space/SpacePoint.hpp"

class Db;
class Faults;

class GSTLEARN_EXPORT NeighMoving: public ANeigh
{
public:
  NeighMoving(bool flag_xvalid = false,
              int nmaxi = 1000,
              double radius = TEST,
              int nmini = 1,
              int nsect = 1,
              int nsmax = ITEST,
              VectorDouble coeffs = VectorDouble(),
              VectorDouble angles = VectorDouble(),
              bool forceWithinCell = false,
              double distcont = TEST,
              const Faults *faults = nullptr,
              const ASpace* space = nullptr);
  NeighMoving(const NeighMoving& r);
  NeighMoving& operator=(const NeighMoving& r);
  virtual ~NeighMoving();

  /// Interface for ANeigh
  virtual int attach(const Db *dbin, const Db *dbout = nullptr) override;
  virtual VectorInt getNeigh(int iech_out) override;
  virtual bool hasChanged(int iech_out) const override;
  virtual VectorDouble summary(int iech_out) override;
  virtual int getMaxSampleNumber(const Db* db) const override;
  virtual ENeigh getType() const override { return ENeigh::fromKey("MOVING"); }
  virtual bool getFlagContinuous() const override { return (! FFFF(_distCont) && _distCont < 1.); }

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static NeighMoving* create(bool flag_xvalid = false,
                             int nmaxi = 1000,
                             double radius = TEST,
                             int nmini = 1,
                             int nsect = 1,
                             int nsmax = ITEST,
                             VectorDouble coeffs = VectorDouble(),
                             VectorDouble angles = VectorDouble(),
                             bool forceWithinCell = false,
                             double distcont = TEST,
                             const Faults* faults = nullptr,
                             const ASpace* space = nullptr);
  static NeighMoving* createFromNF(const String& neutralFilename, bool verbose = true);

  void addBiPointCheck(const ABiPointCheck* abpc);

  bool getFlagSector() const;
  int getNMaxi() const { return _nMaxi; }
  int getNMini() const { return _nMini; }
  int getNSect() const { return _nSect; }
  int getNSMax() const { return _nSMax; }
  double getDistCont() const { return _distCont; }
  bool   getForceWithinCell() const { return _forceWithinCell; }
  const BiPointCheckDistance* getBiPtDist() const { return _biPtDist; }

  void setNMaxi(int nmaxi) { _nMaxi = nmaxi; }
  void setNMini(int nmini) { _nMini = nmini; }
  void setNSect(int nsect) { _nSect = nsect; }
  void setNSMax(int nsmax) { _nSMax = nsmax; }
  void setDistCont(double distCont) { _distCont = distCont; }
  void setForceWithinCell(bool forceWithinCell) { _forceWithinCell = forceWithinCell; }

  VectorVectorDouble getEllipsoid(const VectorDouble& target, int count = 360) const;
  VectorVectorDouble getSectors(const VectorDouble& target) const;
  VectorVectorDouble getZoomLimits(const VectorDouble& target, double percent=20) const;

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "NeighMoving"; }

private:
  int  _getBiPtsNumber() const { return (int) _bipts.size(); }
  int  _moving(int iech_out, VectorInt& ranks, double eps = EPSILON9);
  bool _belongsToCell(int iech, int iech_out);
  int  _movingSectorDefine(double dx, double dy);
  void _movingSectorNsmax(int nsel, VectorInt& ranks);
  void _movingSelect(int nsel, VectorInt& ranks);
  double _getRadius() const { return _biPtDist->getRadius(); }

private:
  int _nMini;                    /* Minimum number of points in neigh. */
  int _nMaxi;                    /* Maximum number of points in neigh. */
  int _nSect;                    /* Number of 2-D angular sectors */
  int _nSMax;                    /* Maximum number of points per 2-D sector */
  bool _forceWithinCell;         /* Select all samples within a Block */
  double _distCont;              /* Distance for continuous neighborhood */

  BiPointCheckDistance* _biPtDist;
  std::vector<const ABiPointCheck*> _bipts;

  mutable VectorInt    _movingInd;
  mutable VectorInt    _movingIsect;
  mutable VectorInt    _movingNsect;
  mutable VectorDouble _movingDst;

  mutable SpacePoint   _P1;
  mutable SpacePoint   _P2;
};

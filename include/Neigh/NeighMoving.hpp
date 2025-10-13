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

#include "Basic/ICloneable.hpp"
#include "Space/ASpace.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Enum/ENeigh.hpp"

#include "Basic/Utilities.hpp"
#include "Geometry/ABiTargetCheck.hpp"
#include "Geometry/BiTargetCheckDistance.hpp"
#include "Neigh/ANeigh.hpp"
#include "Space/SpaceTarget.hpp"

namespace gstlrn
{
class Db;

/**
 * \brief
 * Moving Neighborhood definition.
 *
 * The Neighborhood is usually meant to select a sub-population from the input Data Base,
 * containing the active samples close to the target.
 *
 * The Moving Neighborhood selects these active samples based on a series of criteria
 * (the corresponding parameters are given between parentheses), such as:
 * - the selected samples should belong to a circle (ellipse) centered on the target sample
 * (circle radius, ellipse orientation and extensions)
 * - the minimum and maximum number of selected samples
 * - the previous circle can be subdivided into angular sectors: the selected samples are
 * taken regularly per sector (maximum number of samples per sector)
 *
 * The neighborhood also offers the possibility to suppress any sample which would be too close to (or coincide with)
 * the target: this is the cross-validation option.
 */
class GSTLEARN_EXPORT NeighMoving: public ANeigh
{
public:
  NeighMoving(bool flag_xvalid             = false,
              Id nmaxi                     = 1000,
              double radius                = TEST,
              Id nmini                     = 1,
              Id nsect                     = 1,
              Id nsmax                     = ITEST,
              const VectorDouble& coeffs   = VectorDouble(),
              const VectorDouble& angles   = VectorDouble(),
              bool useBallTree             = false,
              Id leaf_size                 = 10,
              const ASpaceSharedPtr& space = ASpaceSharedPtr());
  NeighMoving(const NeighMoving& r);
  NeighMoving& operator=(const NeighMoving& r);
  virtual ~NeighMoving();

  IMPLEMENT_CLONING(NeighMoving)

  /// Interface for ANeigh
  Id attach(const Db* dbin, const Db* dbout = nullptr) override;
  void getNeigh(Id iech_out, VectorInt& ranks) override;
  bool hasChanged(Id iech_out) const override;
  VectorDouble summary(Id iech_out) override;
  Id getNSampleMax(const Db* db) const override;
  ENeigh getType() const override { return ENeigh::fromKey("MOVING"); }
  bool getFlagContinuous() const override { return (!FFFF(_distCont) && _distCont < 1.); }

  /// Interface for AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  static NeighMoving* create(bool flag_xvalid             = false,
                             Id nmaxi                     = 1000,
                             double radius                = TEST,
                             Id nmini                     = 1,
                             Id nsect                     = 1,
                             Id nsmax                     = ITEST,
                             const VectorDouble& coeffs   = VectorDouble(),
                             const VectorDouble& angles   = VectorDouble(),
                             bool useBallTree             = false,
                             Id leaf_size                 = 10,
                             const ASpaceSharedPtr& space = ASpaceSharedPtr());
  static NeighMoving* createFromNF(const String& NFFilename, bool verbose = true);

  void addBiTargetCheck(ABiTargetCheck* abpc);

  bool getFlagSector() const;
  Id getNMaxi() const { return _nMaxi; }
  Id getNMini() const { return _nMini; }
  Id getNSect() const { return _nSect; }
  Id getNSMax() const { return _nSMax; }
  double getDistCont() const { return _distCont; }
  const BiTargetCheckDistance* getBiPtDist() const { return _biPtDist; }
  bool getFlagAniso() const { return _biPtDist->getFlagAniso(); }
  bool getFlagRotation() const { return _biPtDist->getFlagRotation(); }
  double getRadius() const { return _biPtDist->getRadius(); }
  double getMaxRadius() const;
  const VectorDouble& getAnisoRotMats() const { return _biPtDist->getAnisoRotMats(); }
  const VectorDouble& getAnisoCoeffs() const { return _biPtDist->getAnisoCoeffs(); }
  double getAnisoCoeff(Id i) const { return _biPtDist->getAnisoCoeff(i); }
  const std::vector<ABiTargetCheck*>& getBipts() const { return _bipts; }
  const ABiTargetCheck* getBipts(Id rank) const { return _bipts[rank]; }

  void setNMaxi(Id nmaxi) { _nMaxi = nmaxi; }
  void setNMini(Id nmini) { _nMini = nmini; }
  void setNSect(Id nsect) { _nSect = nsect; }
  void setNSMax(Id nsmax) { _nSMax = nsmax; }
  void setDistCont(double distCont) { _distCont = distCont; }

  VectorVectorDouble getEllipsoid(const VectorDouble& target, Id count = 360) const;
  VectorVectorDouble getSectors(const VectorDouble& target) const;
  VectorVectorDouble getZoomLimits(const VectorDouble& target, double percent = 20) const;

protected:
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "NeighMoving"; }

private:
  Id _getNBiPts() const { return static_cast<Id>(_bipts.size()); }
  Id _moving(Id iech_out, VectorInt& ranks, double eps = EPSILON9);
  Id _movingSectorDefine(double dx, double dy) const;
  void _movingSectorNsmax(Id nsel, VectorInt& ranks);
  void _movingSelect(Id nsel, VectorInt& ranks);
  double _getRadius() const { return _biPtDist->getRadius(); }
  bool _getAnisotropyElements(double* rx,
                              double* ry,
                              double* theta,
                              double* cosp,
                              double* sinp) const;

private:
  Id _nMini;        /* Minimum number of points in neigh. */
  Id _nMaxi;        /* Maximum number of points in neigh. */
  Id _nSect;        /* Number of 2-D angular sectors */
  Id _nSMax;        /* Maximum number of points per 2-D sector */
  double _distCont; /* Distance for continuous neighborhood */

  BiTargetCheckDistance* _biPtDist;
  std::vector<ABiTargetCheck*> _bipts;

  mutable VectorInt _movingInd;
  mutable VectorInt _movingIsect;
  mutable VectorInt _movingNsect;
  mutable VectorDouble _movingDst;

  mutable const DbGrid* _dbgrid; // Pointer not to be deleted
  mutable SpaceTarget _T1;
  mutable SpaceTarget _T2;
};
} // namespace gstlrn
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

#include "Space/ASpace.hpp"
#include "Tree/Ball.hpp"
#include "gstlearn_export.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Db/Db.hpp"
#include "Faults/Faults.hpp"
#include "Space/SpaceTarget.hpp"
#include "Variogram/DirParam.hpp"

namespace gstlrn
{
  class Model;
  class Db;

#ifndef SWIG
  /**
   * Generic pair iterator over a Db for variogram computations.
   * Supports both spatial BallTree neighborhood queries and standard 1D sorted scans.
   *
   * @param db               Pointer to the Db structure containing samples
   * @param space            Shared pointer to the spatial context
   * @param idir             Rank of the current calculation direction
   * @param dirparam         Directional parameter defining distance/angle limits
   * @param hasDate          Boolean indicating if temporal checks require full pair scanning
   * @param keepPair         Lambda predicate: (idir, T1, T2, &dist) -> bool
   * @param processPair      Lambda action on valid pair: (iech, jech, ilag, dist) -> void
   * @param processOuterEnd  Optional lambda action executed at the end of the outer loop: (iech, weight) -> void
   */
  template<
    typename KeepPairFunc,
    typename ProcessPairFunc,
    typename ProcessOuterEndFunc = std::nullptr_t>
  void loopOnPairs(
    const Db* db,
    const ASpaceSharedPtr& space,
    Id idir,
    const DirParam& dirparam,
    bool hasDate,
    KeepPairFunc&& keepPair,
    ProcessPairFunc&& processPair,
    ProcessOuterEndFunc&& processOuterEnd = nullptr);
#endif // SWIG

  /**
   * \brief
   * Class containing the definition of the criteria for calculating the Spatial (and Temporal) Characteristics
   * from samples contained in a Db.
   *
   * These criteria consist in:
   * - some criteria based on the **dates**: this information will is used for calculating the Temporal Characteristics
   * - a collection of definitions of **Calculation Directions** for Spatial Characteristics.
   * For more information on a Direction definition, please refer to DirParam.hpp
   *
   * Note that this class also stores a pointer to any Faults definition, if to be used during the
   * calculation of the Spatial Characteristics.
   */
  class GSTLEARN_EXPORT VarioParam: public AStringable, public ICloneable
  {
  public:
    VarioParam(
      double scale = 0.,
      const VectorDouble& dates = VectorDouble(),
      const Faults* faults = nullptr);
    VarioParam(
      const VarioParam& VarioParam,
      const VectorInt& dircols,
      const Faults* faults = nullptr);
    VarioParam(const VarioParam& r);
    VarioParam& operator=(const VarioParam& r);
    virtual ~VarioParam();

  public:
    /// ICloneable interface
    IMPLEMENT_CLONING(VarioParam)

    /// AStringable Interface
    String toString(const AStringFormat* strfmt = nullptr) const override;

    /// Shortcuts
    static VarioParam* createOmniDirection(
      Id nlag = 10,
      double dlag = 1.,
      double toldis = 0.5,
      Id opt_code = 0,
      Id idate = 0,
      double bench = TEST,
      double cylrad = TEST,
      double tolcode = 0.,
      const VectorDouble& breaks = VectorDouble(),
      double scale = 0.,
      const VectorDouble& dates = VectorDouble(),
      const VectorDouble& benchdir = VectorDouble(),
      const ASpaceSharedPtr& space = ASpaceSharedPtr());
    static VarioParam* createMultiple(
      Id ndir,
      Id nlag = 10,
      double dlag = 1.,
      double toldis = 0.5,
      double angref = 0.,
      double scale = 0.,
      const VectorDouble& dates = VectorDouble(),
      const ASpaceSharedPtr& space = ASpaceSharedPtr());
    static VarioParam* createMultipleFromGrid(
      const DbGrid* dbgrid,
      Id nlag,
      double scale = 0.,
      const VectorDouble& dates = VectorDouble(),
      const ASpaceSharedPtr& space = ASpaceSharedPtr(),
      Id ndimax = 0);
    static VarioParam* createFromSpaceDimension(
      Id nlag = 10,
      double dlag = 1.,
      double toldis = 0.5,
      double tolang = 45.,
      double scale = 0.,
      const VectorDouble& dates = VectorDouble(),
      const ASpaceSharedPtr& space = ASpaceSharedPtr());
    static VarioParam* createSeveral2D(
      const VectorDouble& angles,
      Id nlag = 10,
      double dlag = 1.,
      double toldis = 0.5,
      double tolang = TEST,
      double scale = 0.,
      const VectorDouble& dates = VectorDouble(),
      const ASpaceSharedPtr& space = ASpaceSharedPtr());
    static VarioParam* createForDb(
      Id ndir = 1,
      Id nlag = 10,
      double dlag = 1.,
      const VectorDouble& angles = VectorDouble(),
      double toldis = 0.5,
      double tolang = TEST,
      const ASpaceSharedPtr& space = ASpaceSharedPtr());
    static VarioParam* createForGrid(
      const DbGrid* dbgrid,
      bool flagAllDirections = true,
      Id nlag = 10,
      const VectorVectorInt& dirincr = VectorVectorInt(),
      const ASpaceSharedPtr& space = ASpaceSharedPtr());

    void addDir(const DirParam& dirparam);
    void addMultiDirs(const std::vector<DirParam>& dirparams);
    void delDir(Id rank);
    void delAllDirs();

    ASpaceSharedPtr getSpace() const { return _dirparams[0].getSpace(); }

    double getScale() const { return _scale; }

    Id getNDate() const { return static_cast<Id>(_dates.size() / 2); }

    Id getNDir() const { return static_cast<Id>(_dirparams.size()); }

    const VectorDouble& getDates() const { return _dates; }

    double getDate(Id idate, Id icas) const;
    Id getNLag(Id idir) const;
    VectorDouble getCodirs(Id idir = 0) const;

    const std::vector<DirParam>& getDirParams() const { return _dirparams; }

    const DirParam& getDirParam(Id idir) const { return _dirparams[idir]; }

    Id getNDim() const;
    bool isDefinedForGrid() const;

    Id hasDate() const
    {
      return (
        getNDate() > 0 && (_dates[0] > MINIMUM_BIG || _dates[1] < MAXIMUM_BIG));
    }

    bool isDateUsed(const Db* db1, const Db* db2 = nullptr) const;

    void setScale(double scale) { _scale = scale; }

    void setDates(const VectorDouble& dates) { _dates = dates; }

    void setDPas(Id idir, const DbGrid* db);
    void setGrincr(Id idir, const VectorInt& grincr);

    String toStringMain(const AStringFormat* strfmt = nullptr) const;

    const Faults* getFaults() const { return _faults; }

    bool hasFaults() const { return _faults != nullptr; }

    void addFaults(const Faults* faults) { _faults = faults; }

  private:
    Id _getAddress(Id ivar, Id jvar) const;
    bool _isVariableValid(Id ivar) const;
    bool _isDirectionValid(Id idir) const;
    bool _isBivariableValid(Id i) const;
    bool _isDateValid(Id idate) const;
    void _initMeans();
    void _initVars();
    VectorDouble _getDirectionInterval(Id idir) const;
    bool _validDefinedFromGrid(const DirParam& dirparam) const;

  private:
    double _scale;
    VectorDouble _dates;
    std::vector<DirParam> _dirparams;
    const Faults* _faults; // Pointer copy (not to be deleted)
  };

  GSTLEARN_EXPORT Db*
    buildDbFromVarioParam(Db* db, const VarioParam& varioparam);

  // -----------------------------------------------------------------------------
  // Template Implementation
  // -----------------------------------------------------------------------------

#ifndef SWIG

  template<
    typename KeepPairFunc,
    typename ProcessPairFunc,
    typename ProcessOuterEndFunc>
  void loopOnPairs(
    const Db* db,
    const ASpaceSharedPtr& space,
    Id idir,
    const DirParam& dirparam,
    bool hasDate,
    KeepPairFunc&& keepPair,
    ProcessPairFunc&& processPair,
    ProcessOuterEndFunc&& processOuterEnd)
  {
    // Local flag the use of Ball Tree sort (to speed up processing)
    bool flagBall = true;

    SpaceTarget T1(space, false);
    SpaceTarget T2(space, false);

    Id nech = db->getNSample();
    double maxdist = dirparam.getMaximumDistance();
    bool hasSel = db->hasLocVariable(ELoc::SEL);
    bool hasWeight = db->hasLocVariable(ELoc::W);

    // Common O(N) pre-filtering: build sorted validRanks vector AND boolean active mask
    VectorInt validRanks;
    validRanks.reserve(nech);
    VectorBool activeSample(nech, false);

    VectorInt rindex = db->getSortArray();

    for (Id iiech = 0; iiech < nech; iiech++)
    {
      Id iech = rindex[iiech];
      if (hasSel && !db->isActive(iech)) continue;
      if (hasWeight && FFFF(db->getWeight(iech))) continue;

      validRanks.push_back(iech);
      activeSample[iech] = true;
    }

    Id nvalid = static_cast<Id>(validRanks.size());
    double dist = 0.;

    // -------------------------------------------------------------------------
    // OPTION 1: Spatial Partitioning via BallTree (ONLY if !hasDate && flagBall)
    // -------------------------------------------------------------------------
    if (flagBall && !hasDate)
    {
      Ball ball(db, nullptr, 10, true, 1, false);
      VectorDouble coords(db->getNDim());
      VectorInt neighbors;

      for (Id i = 0; i < nvalid; i++)
      {
        Id iech = validRanks[i];

        db->getSampleAsSTInPlace(iech, T1);
        db->getCoordinatesInPlace(coords, iech);

        // Query spatial neighbors within maxdist bounding sphere
        ball.queryRadiusInPlace(coords, maxdist, neighbors);

        for (Id jech: neighbors)
        {
          // Enforce unique unordered pairs (j > i)
          if (jech <= iech) continue;

          // Fast O(1) check using the active boolean mask
          if (!activeSample[jech]) continue;

          db->getSampleAsSTInPlace(jech, T2);

          if (!keepPair(idir, T1, T2, &dist)) continue;

          Id ilag = dirparam.getLagRank(dist);
          if (isNA(ilag)) continue;

          processPair(iech, jech, ilag, dist);
        }

        if constexpr (!std::is_same_v<
                        std::decay_t<ProcessOuterEndFunc>, std::nullptr_t>)
        {
          double w1 = hasWeight ? db->getWeight(iech) : 1.0;
          processOuterEnd(iech, w1);
        }
      }
    }
    // -------------------------------------------------------------------------
    // OPTION 2: Standard 1D Sorted Array Scan (Legacy / Date Fallback)
    // -------------------------------------------------------------------------
    else
    {
      for (Id i = 0; i < nvalid; i++)
      {
        Id iech = validRanks[i];
        db->getSampleAsSTInPlace(iech, T1);

        // Triangular scan if no dates, full rectangular scan if dates are active
        Id jstart = hasDate ? 0 : i + 1;
        for (Id j = jstart; j < nvalid; j++)
        {
          Id jech = validRanks[j];

          // Early break along the 1D main projection axis
          if (db->getIncrement1D(jech, iech) > maxdist) break;

          db->getSampleAsSTInPlace(jech, T2);

          if (!keepPair(idir, T1, T2, &dist)) continue;

          Id ilag = dirparam.getLagRank(dist);
          if (isNA(ilag)) continue;

          processPair(iech, jech, ilag, dist);
        }

        if constexpr (!std::is_same_v<
                        std::decay_t<ProcessOuterEndFunc>, std::nullptr_t>)
        {
          double w1 = hasWeight ? db->getWeight(iech) : 1.0;
          processOuterEnd(iech, w1);
        }
      }
    }
  }
#endif // SWIG
} // namespace gstlrn

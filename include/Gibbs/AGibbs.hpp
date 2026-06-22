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

#include "../Matrix/Table.hpp"
#include "Basic/AStringable.hpp"

namespace gstlrn
{
  class Db;

  class GSTLEARN_EXPORT AGibbs: public AStringable
  {
  public:
    AGibbs();
    AGibbs(Db* db);
    AGibbs(
      Db* db,
      Id npgs,
      Id nvar,
      Id nburn,
      Id niter,
      Id seed,
      Id flag_order,
      bool flag_decay);
    AGibbs(const AGibbs& r);
    AGibbs& operator=(const AGibbs& r);
    virtual ~AGibbs();

    /// Interface for AStringable
    String toString(const AStringFormat* strfmt = nullptr) const override;

    /// Interface for AGibbs
    virtual Id calculInitialize(VectorVectorDouble& y, Id isimu, Id ipgs) = 0;
    virtual void update(VectorVectorDouble& y, Id isimu, Id ipgs, Id iter) = 0;
    virtual Id covmatAlloc(bool verbose, bool verboseTimer = false) = 0;
    virtual double getSimulate(
      VectorVectorDouble& y,
      double yk,
      double sk,
      Id icase,
      Id ipgs,
      Id ivar,
      Id iact,
      Id iter) = 0;
    virtual Id checkGibbs(const VectorVectorDouble& y, Id isimu, Id ipgs) = 0;

    virtual void cleanup() {}

    void init(
      Id npgs,
      Id nvar,
      Id nburn,
      Id niter,
      Id seed = 3241,
      Id flag_order = 0,
      bool flag_decay = true);
    Id run(
      VectorVectorDouble& y,
      Id ipgs0 = 0,
      Id isimu0 = 0,
      bool verboseTimer = false,
      bool flagCheck = false);

    Id getNvar() const { return _nvar; }

    void setNvar(Id nvar) { _nvar = nvar; }

    Id getNpgs() const { return _npgs; }

    void setNpgs(Id npgs) { _npgs = npgs; }

    Id getNburn() const { return _nburn; }

    void setNburn(Id nburn) { _nburn = nburn; }

    Id getNiter() const { return _niter; }

    void setNiter(Id niter) { _niter = niter; }

    Id getFlagOrder() const { return _flagOrder; }

    void setFlagOrder(Id flagOrder) { _flagOrder = flagOrder; }

    bool getOptionStats() const { return _optionStats; }

    void setOptionStats(Id option_stats) { _optionStats = option_stats; }

    Db* getDb() const { return _db; }

    VectorVectorDouble allocY() const;
    void storeResult(const VectorVectorDouble& y, Id isimu, Id ipgs);
    Id getNSample() const;
    Id getSampleRank(Id i) const;
    Id getRank(Id ipgs, Id ivar) const;

  protected:
    Id _getDimension() const;
    Id _getSampleRankNumber() const;
    void _statsInit();
    bool _isConstraintTight(Id icase, Id iact, double* value) const;
    void _updateStats(
      const VectorVectorDouble& y,
      Id ipgs,
      Id jter,
      double amort = 0.9);
    void _getBoundsDecay(Id iter, double* vmin, double* vmax) const;
    Id
      _boundsCheck(Id ipgs, Id ivar, Id iact, double* vmin, double* vmax) const;
    void _printInequalities(
      Id iact,
      Id ivar,
      double simval,
      double vmin,
      double vmax) const;
    Id _getNRowStats() const;
    Id _getNColStats() const;
    Id _getColRankStats(Id ipgs, Id ivar, Id mode) const;
    void _displayCurrentVector(
      bool flag_init,
      const VectorVectorDouble& y,
      Id isimu,
      Id ipgs) const;

    const VectorInt& _getRanks() const { return _ranks; }

  private:
    VectorInt _calculateSampleRanks() const;
    Id _getRelativeRank(Id iech);

  private:
    Id _npgs;
    Id _nvar; // or NGRF
    Id _nburn;
    Id _niter;
    Id _flagOrder; // order relationship of the constraints
    //   1 if the ascending order must be honored
    //  -1 if the descending order must be honored
    //   0 if no order relationship must be honored
    bool _flagDecay;
    Id _optionStats; // 0: no storage; 1: printout; 2: Neutral File

    VectorInt _ranks; // Internal array use to store indices of active samples
    // Pointer to the reference Db (only stored for efficiency)
    Db* _db;
    // Optional Table used to store performance statistics (see _optionStats)
    Table _stats;
  };
} // namespace gstlrn

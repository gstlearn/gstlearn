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


class Db;

class GSTLEARN_EXPORT AGibbs : public AStringable
{
public:
  AGibbs();
  AGibbs(Db* db);
  AGibbs(Db* db,
         int npgs,
         int nvar,
         int nburn,
         int niter,
         int seed,
         int flag_order,
         bool flag_decay);
  AGibbs(const AGibbs &r);
  AGibbs& operator=(const AGibbs &r);
  virtual ~AGibbs();

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface for AGibbs
  virtual int calculInitialize(VectorVectorDouble &y, int isimu, int ipgs) = 0;
  virtual void update(VectorVectorDouble& y, int isimu, int ipgs, int iter) = 0;
  virtual int covmatAlloc(bool verbose, bool verboseTimer = false) = 0;
  virtual double getSimulate(VectorVectorDouble& y,
                             double yk,
                             double sk,
                             int icase,
                             int ipgs,
                             int ivar,
                             int iact,
                             int iter) = 0;
  virtual int checkGibbs(const VectorVectorDouble& y, int isimu, int ipgs) = 0;
  virtual void cleanup() { }

  void init(int npgs,
            int nvar,
            int nburn,
            int niter,
            int seed = 3241,
            int flag_order = 0,
            bool flag_decay= true);
  int run(VectorVectorDouble &y,
          int ipgs0 = 0,
          int isimu0 = 0,
          bool verboseTimer = false,
          bool flagCheck = false);

  int getNvar() const { return _nvar; }
  void setNvar(int nvar) { _nvar = nvar; }
  int getNpgs() const { return _npgs; }
  void setNpgs(int npgs) { _npgs = npgs; }
  int getNburn() const { return _nburn; }
  void setNburn(int nburn) { _nburn = nburn; }
  int getNiter() const { return _niter; }
  void setNiter(int niter) { _niter = niter; }
  int getFlagOrder() const { return _flagOrder; }
  void setFlagOrder(int flagOrder) { _flagOrder = flagOrder; }
  bool getOptionStats() const { return _optionStats; }
  void setOptionStats(int option_stats) { _optionStats = option_stats; }
  Db* getDb() const { return _db; }

  VectorVectorDouble allocY() const;
  void storeResult(const VectorVectorDouble& y, int isimu, int ipgs);
  int getSampleNumber() const;
  int getSampleRank(int i) const;
  int getRank(int ipgs, int ivar) const;

protected:
  int  _getDimension() const;
  int  _getSampleRankNumber() const;
  void _statsInit();
  bool _isConstraintTight(int icase, int iact, double* value) const;
  void _updateStats(const VectorVectorDouble &y,
                    int ipgs,
                    int jter,
                    double amort = 0.9);
  void _getBoundsDecay(int iter, double *vmin, double *vmax) const;
  int  _boundsCheck(int ipgs, int ivar, int iact, double *vmin, double *vmax) const;
  void _printInequalities(int iact,
                          int ivar,
                          double simval,
                          double vmin,
                          double vmax) const;
  int _getRowNumberStats() const;
  int _getColNumberStats() const;
  int _getColRankStats(int ipgs, int ivar, int mode) const;
  void _displayCurrentVector(bool flag_init,
                             const VectorVectorDouble& y,
                             int isimu,
                             int ipgs) const;
  const VectorInt& _getRanks() const { return _ranks; }

private:
  VectorInt _calculateSampleRanks() const;
  int  _getRelativeRank(int iech);

private:
  int _npgs;
  int _nvar; // or NGRF
  int _nburn;
  int _niter;
  int _flagOrder; // order relationship of the constraints
  //   1 if the ascending order must be honored
  //  -1 if the descending order must be honored
  //   0 if no order relationship must be honored
  bool _flagDecay;
  int  _optionStats; // 0: no storage; 1: printout; 2: Neutral File

  VectorInt _ranks; // Internal array use to store indices of active samples
  // Pointer to the reference Db (only stored for efficiency)
  Db*   _db;
  // Optional Table used to store performance statistics (see _optionStats)
  Table _stats;
};

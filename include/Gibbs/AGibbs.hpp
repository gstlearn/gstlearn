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

#include "Basic/AStringable.hpp"
#include "Basic/Table.hpp"

class Db;

class AGibbs
{
public:
  AGibbs();
  AGibbs(Db* db);
  AGibbs(Db* db,
         int npgs,
         int nvar,
         int nburn,
         int niter,
         int flag_order,
         bool flag_multi_mono,
         bool flag_decay);
  AGibbs(const AGibbs &r);
  AGibbs& operator=(const AGibbs &r);
  virtual ~AGibbs();

  virtual int calculInitialize(VectorVectorDouble& y,
                               int isimu,
                               int ipgs,
                               bool verbose) = 0;
  virtual void update(VectorVectorDouble& y, int isimu, int ipgs, int iter) = 0;
  virtual int covmatAlloc(bool verbose) = 0;
  virtual double getSimulate(VectorVectorDouble& y,
                             double yk,
                             double sk,
                             int ipgs,
                             int ivar,
                             int iact,
                             int iter) = 0;
  virtual int checkGibbs(const VectorVectorDouble& y, int isimu, int ipgs) = 0;

  virtual int run(VectorVectorDouble& y, int ipgs=0, int isimu=0, bool verbose = false, bool flagCheck = false);

  void init(int npgs,
            int nvar,
            int nburn,
            int niter,
            int flag_order,
            bool flag_multi_mono,
            bool flag_decay);
  void getBoundsDecay(int iter, double *vmin, double *vmax) const;
  void print(bool flag_init,
             const VectorVectorDouble& y,
             int isimu,
             int ipgs) const;

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
  int getDimension() const;
  int getRank(int ipgs, int ivar) const;
  VectorVectorDouble allocY() const;
  void storeResult(const VectorVectorDouble& y, int isimu, int ipgs);


  int getSampleRankNumber() const;
  int getSampleRank(int i) const;
  VectorInt calculateSampleRanks() const;
  void updateStats(const VectorVectorDouble& y,
                   int ipgs,
                   int iter,
                   double amort = 0.9);
  bool isConstraintTight(int ipgs, int ivar, int iact, double* value) const;
  void statsInit();

  bool getFlagDecay() const { return _flagDecay; }
  int  getRelativeRank(int iech);

protected:
  int _boundsCheck(int ipgs, int ivar, int iact, double *vmin, double *vmax);
  void _printInequalities(int iact,
                          int ivar,
                          int nfois,
                          int flag_cv,
                          double simval,
                          double vmin,
                          double vmax) const;
  int _getRowNumberStats() const;
  int _getColNumberStats() const;
  int _getColRankStats(int ipgs, int ivar, int mode) const;


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

  // Pointer to the reference Db (stored for efficiency)
  Db*       _db;

  // Optional Table used to store performance statistics (see _optionStats)
  Table     _stats;
};

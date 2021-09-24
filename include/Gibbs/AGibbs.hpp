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
#include "Stats/StatTable.hpp"
#include "geoslib_enum.h"

class Db;
class Model;

class AGibbs
{
public:
  AGibbs();
  AGibbs(Db* db,Model* model);
  AGibbs(Db* db,
         Model* model,
         int npgs,
         int nvar,
         int nburn,
         int niter,
         int flag_order,
         bool flag_multi_mono,
         bool flag_decay,
         double rho,
         double eps = EPSILON3);
  AGibbs(const AGibbs &r);
  AGibbs& operator=(const AGibbs &r);
  virtual ~AGibbs();

  virtual int calculInitialize(VectorVectorDouble& y,
                               int isimu,
                               int ipgs,
                               bool verbose);
  virtual void update(VectorVectorDouble& y,
                      int isimu,
                      int ipgs,
                      int iter) = 0;
  virtual int covmatAlloc(bool verbose) = 0;

  void init(int npgs, int nvar, int nburn, int niter,
            int flag_order, bool flag_multi_mono, bool flag_decay,
            double rho, double eps);

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
  double getRho() const { return _rho; }
  void setRho(double rho) { _rho = rho; }
  double getSqr() const { return _sqr; }
  void setSqr(double sqr) { _sqr = sqr; }
  double getEps() const { return _eps; }
  void setEps(double eps) { _eps = eps; }
  bool getFlagCategory() const { return _flagCategory; }
  void setFlagCategory(bool flagCategory) { _flagCategory = flagCategory; }
  int getFlagOrder() const { return _flagOrder; }
  void setFlagOrder(int flagOrder) { _flagOrder = flagOrder; }
  bool isFlagStats() const { return _flagStats; }
  void setFlagStats(bool flagStats) { _flagStats = flagStats; }

  int checkGibbs(const VectorVectorDouble& y, int isimu, int ipgs);

  Db* getDb() const { return _db; }
  Model* getModel() const { return _model; }
  int getDimension() const;
  int getRank(int ipgs, int ivar) const;
  VectorVectorDouble allocY() const;
  void storeResult(const VectorVectorDouble& y, int isimu, int ipgs);

  double getSimulate(VectorVectorDouble& y,
                     double yk,
                     double sk,
                     int iact,
                     int ipgs,
                     int ivar,
                     int iter);
  int getSampleRankNumber() const;
  int getSampleRank(int i) const;
  VectorInt calculateSampleRanks() const;
  void updateStats(const VectorVectorDouble& y, int ipgs, int niter);

protected:
  int  _boundsCheck(int iech0, int ipgs, int ivar, double *vmin, double *vmax);
  void _printInequalities(int iact,
                          int ivar,
                          int nfois,
                          int flag_cv,
                          double simval,
                          double vmin,
                          double vmax) const;

private:
  int _npgs;
  int _nvar;
  int _nburn;
  int _niter;
  int _flagOrder; // order relationship of the constraints
  //   1 if the ascending order must be honored
  //  -1 if the descending order must be honored
  //   0 if no order relationship must be honored
  bool _flagCategory; // true for categorical; false for continuous
  bool _flagMultiMono;
  bool _flagDecay;
  bool _flagStats;
  double _rho;
  double _sqr;
  double _eps;
  VectorInt _ranks;
  Db*       _db;
  Model*    _model;
  StatTable _stats;
};

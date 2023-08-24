/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Simulation/ACalcSimulation.hpp"
#include "Simulation/TurningDirection.hpp"
#include "Model/Model.hpp"
#include "Basic/VectorNumT.hpp"

#include "geoslib_define.h"

class Model;
class ANeigh;

class GSTLEARN_EXPORT CalcSimuTurningBands : public ACalcSimulation
{
public:
  CalcSimuTurningBands(int nbsimu = 0,
                       int nbtuba = 0,
                       bool flag_check = false,
                       int seed = 4324324);
  CalcSimuTurningBands(const CalcSimuTurningBands& r) = delete;
  CalcSimuTurningBands& operator=(const CalcSimuTurningBands& r) = delete;
  virtual ~CalcSimuTurningBands();

  int getNBtuba() const { return _nbtuba; }
  void setNBtuba(int nbtuba) { _nbtuba = nbtuba; }
  int getNDirs() const { return (int) _codirs.size(); }

  int simulate(Db *dbin,
               Db *dbout,
               Model* model,
               ANeigh *neigh,
               int icase,
               int flag_bayes = false,
               const VectorDouble& dmean = VectorDouble(),
               const VectorDouble& dcov = VectorDouble(),
               bool flag_pgs = false,
               bool flag_gibbs = false,
               bool flag_dgm = false);
  int simulatePotential(Db *dbiso,
                        Db *dbgrd,
                        Db *dbtgt,
                        Db *dbout,
                        Model* model,
                        double delta);

  static bool isTurningBandsWorkable(const Model *model);

  const VectorDouble& getBayesCov() const { return _bayesCov; }
  void setBayesCov(const VectorDouble &dcov) { _bayesCov = dcov; }
  const VectorDouble& getBayesMean() const { return _bayesMean; }
  void setBayesMean(const VectorDouble &dmean) { _bayesMean = dmean; }
  bool isFlagCheck() const { return _flagCheck; }
  void setFlagCheck(bool flag_check) { _flagCheck = flag_check; }
  bool isFlagBayes() const { return _flagBayes; }
  void setFlagBayes(bool flag_bayes) { _flagBayes = flag_bayes; }
  void setFlagDgm(bool flag_dgm) { _flagDGM = flag_dgm; }
  bool isFlagGibbs() const { return _flagGibbs; }
  void setFlagGibbs(bool flag_gibbs) { _flagGibbs = flag_gibbs; }
  bool isFlagPgs() const { return _flagPGS; }
  void setFlagPgs(bool flag_pgs) { _flagPGS = flag_pgs; }
  int getIcase() const { return _icase; }
  void setIcase(int icase) { _icase = icase; }
  int getNbtuba() const { return _nbtuba; }
  void setNbtuba(int nbtuba) { _nbtuba = nbtuba; }

private:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;

  bool _resize();
  void _simulatePoint(Db *db, const VectorDouble& aic, int icase, int shift);
  void _simulateGrid(DbGrid *db, const VectorDouble& aic, int icase, int shift);;
  void _simulateNugget(Db *db, const VectorDouble& aic, int icase);
  void _simulateGradient(Db *dbgrd, const VectorDouble& aic, double delta);
  void _simulateTangent(Db *dbtgt, const VectorDouble& aic, double delta);
  void _meanCorrect(Db *dbout, int icase);
  void _difference(Db *dbin,
                   Model* model,
                   int icase,
                   bool flag_pgs = false,
                   bool flag_gibbs = false,
                   bool flag_dgm = false);
  void _updateData2ToTarget(Db *dbin,
                            Db *dbout,
                            int icase,
                            bool flag_pgs = false,
                            bool flag_dgm = false);
  void _checkGaussianData2Grid(Db *dbin, Db *dbout, Model *model) const;

  void _setCodirAng(int ibs, int idir, double value) { _codirs[ibs].setAng(idir, value); }
  void _setCodirTmin(int ibs, double value) { _codirs[ibs].setTmin(value); }
  void _setCodirTmax(int ibs, double value) { _codirs[ibs].setTmax(value); }
  void _setCodirScale(int ibs, double value) { _codirs[ibs].setScale(value); }
  void _setCodirT00(int ibs, double value) { _codirs[ibs].setT00(value); }
  void _setCodirDXP(int ibs, double value) { _codirs[ibs].setDXP(value); }
  void _setCodirDYP(int ibs, double value) { _codirs[ibs].setDYP(value); }
  void _setCodirDZP(int ibs, double value) { _codirs[ibs].setDZP(value); }

  VectorDouble _getCodirAng(int ibs) const { return _codirs[ibs].getAng(); }
  double _getCodirAng(int ibs, int idir) const { return _codirs[ibs].getAng(idir); }
  double _getCodirScale(int ibs) { return _codirs[ibs].getScale(); }
  double _getCodirT00(int ibs) const { return _codirs[ibs].getT00(); }
  double _getCodirDXP(int ibs) const { return _codirs[ibs].getDXP(); }
  double _getCodirDYP(int ibs) const { return _codirs[ibs].getDYP(); }
  double _getCodirDZP(int ibs) const { return _codirs[ibs].getDZP(); }
  double _getCodirTmin(int ibs) const { return _codirs[ibs].getTmin(); }
  double _getCodirTmax(int ibs) const { return _codirs[ibs].getTmax(); }

  int  _getAddressBand(int ivar, int is, int ib, int isimu);
  void _setSeedBand(int ivar, int is, int ib, int isimu, int seed);
  int  _getSeedBand(int ivar, int is, int ib, int isimu);

  void _rotateDirections(double a[3], double theta);
  int  _generateDirections(const Db* dbout);
  void _minmax(const Db *db);
  void _setDensity();
  ECov _particularCase(const ECov &type, double param);
  int  _initializeSeedBands();
  VectorDouble _createAIC();
  double _getAIC(const VectorDouble& aic, int icov, int ivar, int jvar);

  double _computeScale(double alpha, double scale);
  VectorDouble _migration(double tmin,
                          double tmax,
                          double scale,
                          double eps = EPSILON5);
  VectorDouble _dilution(double tmin, double tmax, double mesh, double *start);
  void _spectral(const ECov &type,
                 double scale,
                 double param,
                 double *omega,
                 double *phi);
  double _computeScaleKB(double param, double scale);
  void _power1D(int ib,
                double scale,
                double alpha,
                double *omega,
                double *phi,
                double *theta_3,
                double *correc0);
  void _spline1D(int ib,
                 double scale,
                 int k,
                 double *omega,
                 double *phi,
                 double *xi_3,
                 double *correc0);
  void _irfProcess(const ECov &type,
                   const VectorDouble& t,
                   VectorDouble& v0,
                   VectorDouble& v1,
                   VectorDouble& v2);
  int _rankInPoisson(int def_rank, double t0, const VectorDouble& t);
  int _rankRegular(double t0, double tdeb, double scale);
  double _irfCorrec(const ECov &type, double theta1, double scale);
  double _irfProcessSample(const ECov &type,
                           int nt0,
                           double t0,
                           const VectorDouble& t,
                           const VectorDouble& v0,
                           const VectorDouble& v1,
                           const VectorDouble& v2);
  void _getOmegaPhi(int ibs,
                    double omega,
                    double phi,
                    double* cxp,
                    double* sxp,
                    double* cyp,
                    double* syp,
                    double* czp,
                    double* szp,
                    double* c0z,
                    double* s0z);

private:
  int  _nbtuba;
  int  _iattOut;
  int  _icase;
  bool _flagCheck;
  bool _flagBayes;
  bool _flagPGS;
  bool _flagGibbs;
  bool _flagDGM;
  VectorString _nameCoord;
  VectorDouble _bayesMean;
  VectorDouble _bayesCov;
  int _npointSimulated;
  double _field;
  double _theta;
  VectorInt _seedBands;
  std::vector<TurningDirection> _codirs;
};

GSTLEARN_EXPORT int simtub(Db *dbin = nullptr,
                           Db *dbout = nullptr,
                           Model *model = nullptr,
                           ANeigh *neigh = nullptr,
                           int nbsimu = 1,
                           int seed = 43431,
                           int nbtuba = 100,
                           bool flag_dgm = false,
                           bool flag_check = false,
                           const NamingConvention &namconv = NamingConvention("Simu"));
GSTLEARN_EXPORT int simbayes(Db *dbin,
                             Db *dbout,
                             Model *model,
                             ANeigh *neigh,
                             int nbsimu = 1,
                             int seed = 132141,
                             const VectorDouble& dmean = VectorDouble(),
                             const VectorDouble& dcov = VectorDouble(),
                             int nbtuba = 100,
                             bool flag_check = false,
                             const NamingConvention& namconv = NamingConvention("SimBayes"));
GSTLEARN_EXPORT bool isCovValidForTurningBands(const ECov& type);

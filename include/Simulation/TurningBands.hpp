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

#include "gstlearn_export.hpp"
#include "Simulation/ASimulation.hpp"
#include "Simulation/TurningDirection.hpp"
#include "Model/Model.hpp"
#include "Basic/Vector.hpp"

#include "geoslib_define.h"

class Model;

class GSTLEARN_EXPORT TurningBands : public ASimulation {

public:
  TurningBands(int nbsimu = 0,
               int nbtuba = 0,
               const Model* model = nullptr,
               int seed = 4324324);
  TurningBands(const TurningBands& r);
  TurningBands& operator=(const TurningBands& r);
  virtual ~TurningBands();

  int getNBtuba() const { return _nbtuba; }
  void setNBtuba(int nbtuba) { _nbtuba = nbtuba; }

  int simulate(Db *dbin,
               Db *dbout,
               ANeighParam *neighparam,
               int icase,
               int flag_bayes = false,
               const VectorDouble& dmean = VectorDouble(),
               const VectorDouble& dcov = VectorDouble(),
               bool flag_pgs = false,
               bool flag_gibbs = false,
               bool flag_dgm = false,
               double r_coeff = 0.);
  int simulatePotential(Db *dbiso,
                        Db *dbgrd,
                        Db *dbtgt,
                        Db *dbout,
                        double delta);
  void checkGaussianData2Grid(Db *dbin, Db *dbout, Model *model) const;
  static bool isTurningBandsWorkable(const Model *model);

private:
  void _simulatePoint(Db *db, const VectorDouble& aic, int icase, int shift);
  void _simulateGrid(DbGrid *db, const VectorDouble& aic, int icase, int shift);;
  void _simulateNugget(Db *db, const VectorDouble& aic, int icase);
  void _simulateGradient(Db *dbgrd, const VectorDouble& aic, double delta);
  void _simulateTangent(Db *dbtgt, const VectorDouble& aic, double delta);
  void _meanCorrect(Db *dbout, int icase);
  void _difference(Db *dbin,
                   int icase,
                   bool flag_pgs = false,
                   bool flag_gibbs = false,
                   bool flag_dgm = false,
                   double r_coeff = 0.);
  void _updateData2ToTarget(Db *dbin,
                            Db *dbout,
                            int icase,
                            bool flag_pgs = false,
                            bool flag_dgm = false);

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

  int _getNCova() const { return _model->getCovaNumber(); }
  int _getNVar() const { return _model->getVariableNumber(); }
  int _getNBands() const { return (int) _codirs.size(); }
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
  int _nbtuba;
  int _npointSimulated;
  double _field;
  double _theta;
  VectorInt _seedBands;
  std::vector<TurningDirection> _codirs;
  const Model* _model;
};

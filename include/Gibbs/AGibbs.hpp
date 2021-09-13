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
#include "geoslib_enum.h"

class Db;
class Model;

class AGibbs
{
public:
  AGibbs();
  AGibbs(int npgs,
         int ngrf,
         int nbsimu,
         int nburn,
         int niter,
         double rho,
         double eps = EPSILON3);
  AGibbs(const AGibbs &r);
  AGibbs& operator=(const AGibbs &r);
  virtual ~AGibbs();

  void init(int npgs, int ngrf, int nbsimu, int nburn, int niter, double rho, double eps);

  int getNgrf() const { return _ngrf; }
  void setNgrf(int ngrf) { _ngrf = ngrf; }
  int getNpgs() const { return _npgs; }
  void setNpgs(int npgs) { _npgs = npgs; }
  int getNbsimu() const { return _nbsimu; }
  void setNbsimu(int nbsimu) { _nbsimu = nbsimu; }
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

  int getRank(int ipgs, int igrf) const;
  int checkGibbs(Db *db, Model *model, int isimu, int ipgs, int igrf);

protected:
  int _checkMandatoryAttribute(const String& method,
                               Db *db,
                               ENUM_LOCS locatorType);
  int _boundsCheck(Db *db,
                   int iech0,
                   double data,
                   double *vmin,
                   double *vmax,
                   int iech,
                   double value,
                   double vemin,
                   double vemax);
  int _correctBoundsOrder(int flag_category,
                          int flag_order,
                          Db *db,
                          int iech0,
                          int ivar,
                          int icase,
                          int nvar,
                          double *vlmin_arg,
                          double *vlmax_arg);
  void _printInequalities(Db *db,
                          int ifirst,
                          int iech,
                          int ivar,
                          int nfois,
                          int flag_cv,
                          double simval,
                          double vmin,
                          double vmax,
                          double mean,
                          double delta);
  void _gibbsInitPrint(const char *title,
                       Db *dbin,
                       int nvar,
                       int nbsimu,
                       int isimu,
                       int icase);
  void _gibbsIterPrint(const char *title,
                       Db *dbin,
                       int nvar,
                       int isimu,
                       int niter,
                       int icase);

private:
  int _npgs;
  int _ngrf;
  int _nbsimu;
  int _nburn;
  int _niter;
  double _rho;
  double _sqr;
  double _eps;
};

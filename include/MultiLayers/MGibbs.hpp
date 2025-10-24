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

#include "Basic/AStringable.hpp"
#include "Variogram/Vario.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
class M2D_Environ
{
public:
  Id flag_ed;
  Id iatt_fd;
  Id iatt_fg;
  double zmean;
  double zstdv;
  double zeps;
  double zmini;
  double zmaxi;
  double dmini;
  double dmaxi;
  double ystdv;
  VectorDouble dcoef;
};

class Db;
class Model;
class MatrixSquare;
class Vario;
class M2D_Environ;

class GSTLEARN_EXPORT MGibbs: public AStringable
{
public:
  MGibbs(Db* dbin, Db* dbout, Model* model);
  MGibbs(const MGibbs& m);
  MGibbs& operator=(const MGibbs& m);
  virtual ~MGibbs();

  /// Interface to AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  void setOptions(Id nlayer,
                  Id niter,
                  Id seed,
                  Id nbsimu,
                  Id icol_pinch,
                  bool flag_ed,
                  bool flag_drift,
                  bool flag_ce,
                  bool flag_cstd,
                  bool verbose);
  Id run();

private:
  bool _checkValid();
  bool _checkPinchout();
  static void st_m2denv_manage(double ystdv, M2D_Environ& m2denv);
  void st_m2d_stats_init(M2D_Environ& m2denv, double percent = 0.05);
  Id st_m2d_drift_manage(M2D_Environ& m2denv, Id* iatt_f);
  void st_print_db_constraints(const char* title,
                               Db* db,
                               const VectorDouble& ydat) const;
  Db* st_m2d_create_constraints(M2D_Environ& m2denv, Id* number_hard);
  void st_m2d_drift_save(M2D_Environ& m2denv, double* gwork);
  void st_define_locators(M2D_Environ& m2denv, Db* db, Id nvar) const;
  Id st_m2d_initial_elevations(M2D_Environ& m2denv,
                               Db* dbc,
                               VectorDouble& work) const;
  void st_m2d_stats_updt(M2D_Environ& m2denv, Db* dbc) const;
  Id st_m2d_drift_fitting(M2D_Environ& m2denv,
                          Db* dbc,
                          Id number_hard) const;
  Id st_m2d_drift_inc_manage(M2D_Environ& m2denv,
                             Id mode,
                             Db* dbc);
  void st_m2d_set_M(M2D_Environ& m2denv,
                    Id icol_pinch,
                    Db* db,
                    Id iatt) const;
  static double st_m2d_get_drift(M2D_Environ& m2denv,
                                 Db* db,
                                 Id ilayer0,
                                 Id iech0);
  static double st_m2d_external_drift_increment(M2D_Environ& m2denv,
                                                Db* db,
                                                Id ilayer0,
                                                Id iech0);
  static double st_m2d_get_M(M2D_Environ& m2denv,
                             Db* db,
                             Id type,
                             Id ilayer,
                             Id iech);
  static double st_m2d_get_S(M2D_Environ& m2denv,
                             Db* db,
                             Id type,
                             Id ilayer,
                             Id iech);
  static double st_m2d_draw_elevation(M2D_Environ& m2denv,
                                      Id ilayer,
                                      double lower,
                                      double upper);
  Id st_record_sample(M2D_Environ& m2denv,
                      Db* db,
                      Id iech,
                      Id ndim,
                      Id natt,
                      Id bypass,
                      Id* number_arg,
                      double* tab) const;
  static double st_m2d_draw_gaussian(M2D_Environ& m2denv,
                                     Db* dbc,
                                     bool verbose,
                                     Id iter,
                                     Id ilayer,
                                     Id iech,
                                     double Zval,
                                     double Zcum,
                                     double Zmin,
                                     double Zmax,
                                     double Ymean,
                                     double Ysigma);
  void st_convert_Z2Y(M2D_Environ& m2denv,
                      Db* dbc,
                      Id type,
                      Id iech,
                      VectorDouble& tab) const;
  void st_convert_Y2Z(M2D_Environ& m2denv,
                      Db* db,
                      Id type,
                      Id iech,
                      VectorDouble& tab) const;
  Id st_global_gibbs(M2D_Environ& m2denv,
                     Db* dbc,
                     Id iter,
                     double sigma,
                     VectorDouble& ymean,
                     VectorDouble& ydat,
                     VectorDouble& work);
  void st_print_sample(const char* title,
                       M2D_Environ& m2denv,
                       Db* dbc,
                       Id iech,
                       VectorDouble& work) const;
  Id st_check_gibbs_data(const char* title,
                         M2D_Environ& m2denv,
                         Db* dbc,
                         VectorDouble& ydat,
                         VectorDouble& work);
  void st_m2d_vector_extract(M2D_Environ& m2denv,
                             Db* dbc,
                             VectorDouble& ydat,
                             VectorDouble& work);
  static void st_m2d_print_environ(const char* title, M2D_Environ& m2denv);

private:
  Db* _dbin;
  Db* _dbout;
  Model* _model;
  Id _ndim;
  Id _nechin;
  Id _nechout;
  Id _nlayer;
  Id _niter;
  Id _seed;
  Id _nbsimu;
  Id _icolPinch;
  bool _flagED;
  bool _flagDrift;
  bool _flagCE;
  bool _flagCStd;
  bool _verbose;
};

GSTLEARN_EXPORT Id MLayers_spde(Db* dbin,
                                Db* dbout,
                                Model* model,
                                Id nlayer,
                                Id niter,
                                Id seed,
                                Id nbsimu,
                                Id icol_pinch,
                                bool flag_ed    = false,
                                bool flag_drift = false,
                                bool flag_ce    = false,
                                bool flag_cstd  = false,
                                bool verbose    = false);

} // namespace gstlrn

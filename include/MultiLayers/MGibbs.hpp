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
#include "LinearOp/CholeskySparse.hpp"
#include "Variogram/Vario.hpp"
#include "geoslib_d.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
class M2D_Environ
{
public:
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

class QChol
{
public:
  MatrixSparse* Q;
  CholeskySparse* chol;
};

class QSimu
{
public:
  QChol* QCtt;
  QChol* QCtd;
};

class SPDE_Matelem
{
public:
  VectorDouble Lambda;
  MatrixSparse* S;
  MatrixSparse* Aproj;
  QChol* QC;
  std::vector<QChol*> QCov;
  VectorDouble Isill;
  VectorDouble Csill;
  QSimu* qsimu;
  Cheb_Elem* s_cheb;
  AMesh* amesh;
};

class Db;
class Model;
class MatrixSquare;
class MatrixSparse;
class Vario;
class M2D_Environ;
class Cheb_Elem;
class AMesh;
class EConsElem;
class CovAniso;

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
  static void _m2denv_manage(double ystdv, M2D_Environ& m2denv);
  void _stats_init(M2D_Environ& m2denv, double percent = 0.05);
  Id _drift_manage(M2D_Environ& m2denv, Id* iatt_f, double percent = 0.05);
  void _print_db_constraints(const char* title,
                             Db* db,
                             const VectorDouble& ydat) const;
  Db* _create_constraints(Id* number_hard);
  void _drift_save(M2D_Environ& m2denv, double* gwork);
  void _define_locators(Db* db, Id nvar) const;
  Id _initial_elevations(M2D_Environ& m2denv,
                         Db* dbc,
                         VectorDouble& work,
                         Id njter_max = 20) const;
  void _stats_updt(M2D_Environ& m2denv, Db* dbc, double percent = 0.05) const;
  Id _drift_fitting(M2D_Environ& m2denv,
                    Db* dbc,
                    Id number_hard) const;
  Id _drift_inc_manage(M2D_Environ& m2denv,
                       Id mode,
                       Db* dbc);
  void _set_M(M2D_Environ& m2denv,
              Id icol_pinch,
              Db* db,
              Id iatt) const;
  double _get_drift(M2D_Environ& m2denv,
                    Db* db,
                    Id ilayer0,
                    Id iech0) const;
  static double _external_drift_increment(M2D_Environ& m2denv,
                                          Db* db,
                                          Id ilayer0,
                                          Id iech0);
  static double _get_M(M2D_Environ& m2denv,
                       Db* db,
                       Id type,
                       Id ilayer,
                       Id iech);
  static double _get_S(M2D_Environ& m2denv,
                       Db* db,
                       Id type,
                       Id ilayer,
                       Id iech);
  static double _draw_elevation(M2D_Environ& m2denv,
                                Id ilayer,
                                double lower,
                                double upper);
  Id _record_sample(Db* db,
                    Id iech,
                    Id natt,
                    Id bypass,
                    Id* number_arg,
                    double* tab) const;
  static double _draw_gaussian(M2D_Environ& m2denv,
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
  void _convert_Z2Y(M2D_Environ& m2denv,
                    Db* dbc,
                    Id type,
                    Id iech,
                    VectorDouble& tab) const;
  void _convert_Y2Z(M2D_Environ& m2denv,
                    Db* db,
                    Id type,
                    Id iech,
                    VectorDouble& tab) const;
  Id _global_gibbs(M2D_Environ& m2denv,
                   Db* dbc,
                   Id iter,
                   double sigma,
                   VectorDouble& ymean,
                   VectorDouble& ydat,
                   VectorDouble& work);
  void _print_sample(const char* title,
                     M2D_Environ& m2denv,
                     Db* dbc,
                     Id iech,
                     VectorDouble& work) const;
  Id _check_gibbs_data(const char* title,
                       M2D_Environ& m2denv,
                       Db* dbc,
                       VectorDouble& ydat,
                       VectorDouble& work);
  void _vector_extract(M2D_Environ& m2denv,
                       Db* dbc,
                       VectorDouble& ydat,
                       VectorDouble& work);
  void _print_environ(const char* title, M2D_Environ& m2denv) const;
  static Id _check_validity_MS(Db* db,
                               Id ilayer,
                               Id iech,
                               bool flag_positive,
                               bool flag_verbose,
                               double M,
                               double S,
                               double eps = 1.e-3);
  static Id _kriging_cholesky(QChol* QC,
                              VectorDouble& rhs,
                              VectorDouble& work,
                              VectorDouble& z);
  static Id _qchol_cholesky(bool verbose, QChol* QC);
  static Id _simulate_cholesky(QChol* QC, VectorDouble& work, VectorDouble& zsnc);
  static QChol* _derive_Qc(double s2, QChol* Qc, SPDE_Matelem& Matelem);
  static SPDE_Matelem& _get_current_matelem(Id icov);
  static void _matelem_print(Id icov);
  double _get_isill(Id icov, Id ivar, Id jvar) const;
  QSimu* _qsimu_manage(Id mode, QSimu* qsimu);
  static Id _get_nvertex(Id icov);
  Id _fill_Isill(void) const;
  Id _fill_Csill(void) const;
  Id _spde_prepar(Db* dbin,
                  Db* dbout,
                  const VectorDouble& gext,
                  SPDE_Option& s_option);
  Id _fill_Bhetero(Db* dbin, Db* dbout) const;
  MatrixSparse* _extract_Q1_hetero(Id row_var,
                                   Id col_var,
                                   Id row_oper,
                                   Id col_oper,
                                   Id* nrows,
                                   Id* ncols);
  Id _build_QCov(SPDE_Matelem& Matelem);
  Id _spde_build_matrices(bool verbose);
  void _matelem_manage(Id mode);
  Id _spde_check(const Db* dbin,
                 const Db* dbout,
                 Model* model1,
                 bool verbose,
                 const VectorDouble& gext,
                 bool mesh_dbin,
                 bool mesh_dbout,
                 bool flag_advanced,
                 bool flag_est,
                 bool flag_std,
                 bool flag_gibbs,
                 bool flag_modif);
  static SPDE_Option _spde_option_alloc(void);
  static Id _get_rank(Id ivar, Id jvar);
  static MatrixSparse* _extract_Q1_nugget(Id row_var,
                                          Id col_var,
                                          Id* nrows,
                                          Id* ncols);
  static double _get_nugget_sill(Id ivar, Id jvar);
  static void _environ_print(const Db* dbout, const VectorDouble& gext);
  static void _chol_invert(QChol* qctt,
                           VectorDouble& xcr,
                           const VectorDouble& rhs,
                           const VectorDouble& work);
  static void _chol_simulate(QChol* qctt,
                             VectorDouble& simu,
                             const VectorDouble& work);
  Cheb_Elem* _spde_cheb_manage(Id mode,
                               bool verbose,
                               double power,
                               const VectorDouble& blin,
                               MatrixSparse* S,
                               Cheb_Elem* cheb_old);
  static Id _spde_posterior();
  static void _print_concatenate_interval(const char* title,
                                          double lower,
                                          double upper,
                                          Id tail);
  static void _print_constraints_per_point(Id ilayer,
                                           Id iech,
                                           double value,
                                           double drift,
                                           double vgaus,
                                           double lower,
                                           double upper);
  static void _stats_gaus(const char* title,
                          Id nlayer,
                          Id nech,
                          double* ydat);
  Id _active_sample(Db* db, Id nlayer, Id iech, Id bypass) const;
  static void _print_details(Db* dbc, Id nech, Id ilayer);
  Id _migrate_pinch_to_point(Db* dbout, Db* dbc) const;
  static Id _spde_external_copy(SPDE_Matelem& matelem, Id icov0);
  static Id _is_external_AQ_defined(Id icov0);
  static AMesh* _create_meshes(Db* dbin,
                               Db* dbout,
                               const VectorDouble& gext,
                               SPDE_Option& s_option);
  static Id _chebychev_calculate_coeffs(Cheb_Elem* cheb_elem,
                                        bool verbose,
                                        const VectorDouble& blin);
  static double _chebychev_function(double x,
                                    double power,
                                    const VectorDouble& blin);
  static Id _build_Q(SPDE_Matelem& Matelem);
  static MatrixSparse* _spde_build_Q(MatrixSparse* S,
                                     const VectorDouble& Lambda,
                                     Id nblin,
                                     double* blin);
  static VectorDouble _spde_fill_Lambda(AMesh* amesh,
                                        const VectorDouble& TildeC);
  static VectorDouble _spde_fill_TildeC(AMesh* amesh, const double* units);
  MatrixSparse* _spde_fill_S(AMesh* amesh, const double* units);
  static void _tangent_calculate(double center[3],
                                 const double srot[2],
                                 double axes[2][3]);
  static void _project_plane(double center[3],
                             double axes[2][3],
                             double xyz[3],
                             double coeff[2]);
  static void _triangle_center(AMesh* amesh,
                               Id ncorner,
                               Id imesh,
                               double center[3],
                               double xyz[3][3]);
  static VectorInt _get_vertex_ranks(AMesh* amesh, Db* dbin, Db* dbout);
  Id _fill_Bnugget(Db* dbin) const;
  VectorDouble _spde_get_mesh_dimension(AMesh* amesh) const;
  static Id _identify_nostat_param(const EConsElem& type0,
                                   Id icov0 = -1,
                                   Id ivar0 = -1,
                                   Id jvar0 = -1);
  Id _check_model(const Db* dbin, const Db* dbout, Model* model) const;
  static void _convert_exponential2matern(CovAniso* cova);
  void _calcul_update(void);
  void _calcul_init() const;
  void _compute_hh() const;
  void _compute_blin(void) const;
  void _compute_correc(void) const;
  double _spde_compute_correc(double param) const;
  void _print_all(const char* title) const;
  static double _get_sill_total(Id ivar, Id jvar);
  static double _get_cova_range(void);
  static Id _get_ncova(void);
  static double _get_cova_param(void);
  static void _set_title(Id flag_igrf, Id flag_icov, Id rank, const char* title);
  static void _set_filnug(Id flag_filnug);
  static double _get_cova_sill(Id ivar, Id jvar);
  static QChol* _extract_QC_from_Q(const char* title,
                                   QChol* QC_in,
                                   Id row_auth,
                                   Id col_auth);
  static QChol* _qchol_manage(Id mode, QChol* QC);
  static MatrixSparse* _extract_Q_from_Q(MatrixSparse* Q_in, Id row_auth, Id col_auth);
  static void _qchol_filter(const char* title, Id auth);
  static void _print_status(Id auth);
  static Id _get_filnug(void);
  static void _set_model(Model* model);
  static CovAniso* _get_cova(void);
  static CovAniso* _get_nugget(void);
  static bool _is_model_nugget(void);
  static void _environ_init(void);
  static Model* _get_model(void);

private:
  Db* _dbin;
  Db* _dbout;
  Model* _model;
  Id _ndim;
  Id _nvar;
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

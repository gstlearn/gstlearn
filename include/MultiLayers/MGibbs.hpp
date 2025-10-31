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

class SPDE_Matelem
{
public:
  VectorDouble Lambda;
  MatrixSquare Isill;
  VectorDouble Csill;
  MatrixSparse* S;
  MatrixSparse* Aproj;
  QChol* QC;
  QChol* QCtt;
  QChol* QCtd;
  std::vector<QChol*> QCov;
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
  bool isValid() const { return _hasValidEnvironment; }

private:
  bool _checkArguments();
  bool _checkValidAuxiliary() const;
  bool _checkValidOption();
  bool _checkValidPinchout();
  bool _checkSPDE();
  bool _checkGibbsData(const char* title,
                       M2D_Environ& m2denv,
                       VectorDouble& ydat,
                       VectorDouble& work);
  static bool _checkValidityMS(Db* db,
                               Id ilayer,
                               Id iech,
                               bool flag_positive,
                               bool flag_verbose,
                               double M,
                               double S,
                               double eps = 1.e-3);

  Id _getNcovaWithoutNugget(void) const;
  static Id _getNvertex(Id icov);
  static double _getIsill(Id icov, Id ivar, Id jvar);
  double _getNuggetSill(Id ivar, Id jvar);
  double _getSillTotal(Id ivar, Id jvar) const;
  double _getCovaRange(void) const;
  double _getCovaParam(void) const;
  static Id _getRank(Id ivar, Id jvar);
  double _getCovaSill(Id ivar, Id jvar) const;
  static bool _isFilterNugget(void);
  CovAniso* _getCurrentCova(void) const;
  CovAniso* _getNugget(void) const;

  double _getDrift(M2D_Environ& m2denv, Db* db, Id ilayer0, Id iech0) const;
  static double _getM(M2D_Environ& m2denv, Db* db, Id type, Id ilayer, Id iech);
  static double _getS(M2D_Environ& m2denv, Db* db, Id type, Id ilayer, Id iech);
  VectorInt _getVertexRanks(AMesh* amesh) const;
  static SPDE_Matelem& _getCurrentMatelem(Id icov);
  VectorDouble _getMeshDimension(AMesh* amesh) const;

  void _statsInit(M2D_Environ& m2denv, double percent = 0.05);
  void _statsUpdate(M2D_Environ& m2denv, double percent = 0.05) const;

  void _manageMatelem(Id mode);
  static void _manageM2denv(double ystdv, M2D_Environ& m2denv);
  Id _manageDrift(M2D_Environ& m2denv, Id* iatt_f, double percent = 0.05);
  Id _manageDriftIncrement(M2D_Environ& m2denv, Id mode);
  static QChol* _manageQchol(Id mode, QChol* QC);
  Id _manageCheb(SPDE_Matelem& Matelem,
                 double power,
                 const VectorDouble& blin,
                 MatrixSparse* S);

  void _printEnviron(const char* title, M2D_Environ& m2denv) const;
  static void _printMatelem(Id icov);
  void _printDbConstraints(const char* title,
                           Db* db,
                           const VectorDouble& ydat,
                           Id nprint = 10) const;
  void _printSample(const char* title, Id iech, VectorDouble& work) const;
  static void _printConcatenateInterval(const char* title,
                                        double lower,
                                        double upper,
                                        Id tail);
  static void _printConstraintsPerPoint(Id ilayer,
                                        Id iech,
                                        double value,
                                        double drift,
                                        double vgaus,
                                        double lower,
                                        double upper);
  void _printDetails(Id nech, Id ilayer) const;
  void _printAll(const char* title) const;
  static void _printQcholFilter(const char* title, Id auth);
  static void _printStatus(Id auth);

  void _driftSave(M2D_Environ& m2denv, VectorDouble& gwork);
  Id _driftFitting(M2D_Environ& m2denv, Id number_hard) const;
  static double _incrementED(M2D_Environ& m2denv, Db* db, Id ilayer0, Id iech0);

  Id _createConstraints(Id* number_hard);
  void _defineLocators(Db* db, Id nvar) const;
  Id _initializeElevations(M2D_Environ& m2denv,
                           VectorDouble& work,
                           Id njter_max = 20) const;

  static double _drawElevation(M2D_Environ& m2denv, double lower, double upper);

  double _drawGaussian(M2D_Environ& m2denv,
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
  bool _recordSample(Db* db,
                     Id iech,
                     Id natt,
                     Id& number,
                     VectorDouble& tab) const;
  void _converZ2Y(M2D_Environ& m2denv,
                  Id type,
                  Id iech,
                  VectorDouble& tab) const;
  void _convertY2Z(M2D_Environ& m2denv,
                   Db* db,
                   Id type,
                   Id iech,
                   VectorDouble& tab) const;
  void _convertExponentialToMatern(CovAniso* cova) const;

  Id _gibbs(M2D_Environ& m2denv,
            Id iter,
            double sigma,
            VectorDouble& ymean,
            VectorDouble& ydat,
            VectorDouble& work);

  Id _prepar();
  Id _perform(M2D_Environ& m2denv, Id iatt_out, VectorDouble& gwork);

  Id _fillIsill(SPDE_Matelem& Matelem) const;
  Id _fillCsill(SPDE_Matelem& Matelem, Id icov) const;
  Id _fillAproj(SPDE_Matelem& Matelem) const;
  Id _fillBhetero() const;
  Id _fill_Bnugget() const;
  VectorDouble _fillLambda(AMesh* amesh, const VectorDouble& TildeC);
  static VectorDouble _fillTildeC(AMesh* amesh, const VectorDouble& units);
  MatrixSparse* _fillS(AMesh* amesh, const VectorDouble& units);

  void _extractVector(M2D_Environ& m2denv,
                      VectorDouble& ydat,
                      VectorDouble& lwork);
  static MatrixSparse* _extracQ1Hetero(Id row_var,
                                       Id col_var,
                                       Id row_oper,
                                       Id col_oper,
                                       Id* nrows,
                                       Id* ncols);
  static MatrixSparse* _extractQ1Nugget(Id row_var,
                                        Id col_var,
                                        Id* nrows,
                                        Id* ncols);
  QChol* _extractQCfromQ(const char* title,
                         QChol* QC_in,
                         Id row_auth,
                         Id col_auth) const;
  static MatrixSparse* _extractQfromQ(MatrixSparse* Q_in, Id row_auth, Id col_auth);

  Id _buildQ(SPDE_Matelem& Matelem) const;
  Id _buildQCov(SPDE_Matelem& Matelem);
  Id _buildMatrices(SPDE_Matelem& Matelem);

  static void _cholInvert(QChol* qctt,
                          VectorDouble& xcr,
                          const VectorDouble& rhs);
  static void _cholSimulate(QChol* qctt,
                            VectorDouble& simu,
                            const VectorDouble& work);
  static Id _qcholCholesky(bool verbose, QChol* QC);
  static QChol* _deriveQC(double s2, QChol* Qc, SPDE_Matelem& Matelem);

  Id _migratePinchToPoint() const;
  Id _createMeshes(SPDE_Matelem& Matelem);
  Id _chebychevCoefficients(Cheb_Elem* cheb_elem, const VectorDouble& blin) const;
  static double _chebychevFunction(double x,
                                   double power,
                                   const VectorDouble& blin);

  static void _calculateTangent(double center[3],
                                const double srot[2],
                                double axes[2][3]);
  static void _calculateProjectPlane(double center[3],
                                     double axes[2][3],
                                     double xyz[3],
                                     double coeff[2]);
  static void _calculateTriangleCenter(AMesh* amesh,
                                       Id ncorner,
                                       Id imesh,
                                       double center[3],
                                       double xyz[3][3]);

  static Id _identifyNostatParam(const EConsElem& type0,
                                 Id icov0 = -1,
                                 Id ivar0 = -1,
                                 Id jvar0 = -1);

  void _calculInitialize() const;
  void _calculUpdate(void);

  void _computeHH() const;
  void _computeBLin(void) const;
  void _computeCorrec(void) const;

  static void _setTitle(bool flag_icov, Id rank, const char* title);
  static void _setFilterNugget(Id flag_filnug);
  void _setM(M2D_Environ& m2denv,
             Id icol_pinch,
             Db* db,
             Id iatt) const;

private:
  bool _hasValidEnvironment;
  Db* _dbin;
  Db* _dbout;
  Db* _dbc;
  Model* _model;
  Id _ndim;
  Id _nvar;
  Id _nechin;
  Id _nechout;
  Id _nechc;
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
                                Id nlayer       = 1,
                                Id niter        = 10,
                                Id nbsimu       = 3,
                                Id icol_pinch   = -1,
                                bool flag_ed    = false,
                                bool flag_drift = false,
                                bool flag_ce    = false,
                                bool flag_cstd  = false,
                                Id seed         = 13674,
                                bool verbose    = false);

} // namespace gstlrn

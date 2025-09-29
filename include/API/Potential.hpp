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

#include "Basic/VectorNumT.hpp"
#include "Db/Db.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
class Db;
class DbGrid;
class NeighUnique;
class ModelGeneric;

class GSTLEARN_EXPORT Potential
{
public:
  Potential(Db* dbiso,
            Db* dbgrd,
            Db* dbtgt,
            ModelGeneric* model,
            double nugget_grd = 0.,
            double nugget_tgt = 0.);
  Potential(const Potential& r)            = delete;
  Potential& operator=(const Potential& r) = delete;
  virtual ~Potential();

  void setVerbose(bool value) { _verbose = value; }
  Id kriging(DbGrid* dbout,
             bool flag_pot       = true,
             bool flag_grad      = false,
             bool flag_trans     = false,
             bool flag_save_data = false,
             Id option_part      = 0);
  Id simulate(DbGrid* dbout,
              double dist_tempere = TEST,
              bool flag_trans     = false,
              Id seed             = 13443,
              Id nbsimu           = 1,
              Id nbtuba           = 100);
  Id xvalid(bool flag_dist_conv);

private:
  Id GRX(Id i) const { return (_ndim < 1) ? -1 : i + _ngrd * 0; }
  Id GRY(Id i) const { return (_ndim < 2) ? -1 : i + _ngrd * 1; }
  Id GRZ(Id i) const { return (_ndim < 3) ? -1 : i + _ngrd * 2; }
  Id TGT(Id i) const { return _startTgt + i; }
  Id ISC(Id ic, Id i) const { return _startIso + _ptrPerLayer[ic] + (i) - (ic)-1; }
  Id IAD_GRD(Id ig) const { return _rankGrd[ig]; }
  void set_IAD_GRD(Id ig, Id value) { _rankGrd[ig] = value; }
  Id IAD_TGT(Id it) const { return _rankTgt[it]; }
  void set_IAD_TGT(Id it, Id value) { _rankTgt[it] = value; }
  double IAD_ISO(Id ic, Id i) const { return _rankIso[_ptrPerLayer[ic] + (i)]; }
  double TGT_COO(Id it, Id i) const { return _dbtgt->getCoordinate(IAD_TGT(it), i); }
  double TGT_VAL(Id it, Id idim) const { return (idim >= _ndim) ? 0. : _dbtgt->getLocVariable(ELoc::TGTE, IAD_TGT(it), idim); }
  double GRD_COO(Id ig, Id idim) const { return (idim >= _ndim) ? 0. : _dbgrd->getCoordinate(IAD_GRD(ig), idim); }
  double GRD_VAL(Id ig, Id idim) const { return (idim >= _ndim) ? 0. : _dbgrd->getLocVariable(ELoc::G, IAD_GRD(ig), idim); }
  double ISO_COO(Id ic, Id j, Id idim) const { return (idim >= _ndim) ? 0. : _dbiso->getCoordinate(IAD_ISO(ic, j), idim); }
  Id EXT(Id iext) const { return _startExt + (iext); }

  bool _isEnvironmentValid(DbGrid* dbout, Id nring = 1);
  void _environmentManage(bool flag_pot,
                          bool flag_grad,
                          bool flag_trans,
                          Id option_part = 0);
  Id _updateIsopot();
  Id _updateGradient();
  Id _updateTangent();
  Id _updateModel();
  Id _updateFinal();
  void _saveResultData(Db* db,
                       Id nvar,
                       double value,
                       const ELoc& loctype_pot,
                       const ELoc& loctype_grad,
                       VectorInt& uid_pot,
                       VectorInt& uid_grad) const;
  void _calculateCovs(ModelGeneric* model,
                      bool flag_grad,
                      double x1,
                      double y1,
                      double z1,
                      double x2,
                      double y2,
                      double z2,
                      double& covar,
                      VectorDouble& covGp,
                      VectorDouble& covGG) const;
  static void _setLHS(MatrixSymmetric& lhs, Id i, Id j, double value);
  static void _setRHS(MatrixDense& rhs, Id i, Id isol, double value);
  static double _getLHS(MatrixSymmetric& lhs, Id i, Id j);
  double _setMatUV(double ux,
                   double uy,
                   double uz,
                   double vx,
                   double vy,
                   double vz) const;
  double _setMatUAV(const VectorDouble& a,
                    double ux,
                    double uy,
                    double uz,
                    double vx,
                    double vy,
                    double vz) const;
  Id _buildLHS(Db* dbout, MatrixSymmetric& lhs);
  void _buildRHS(bool flag_grad,
                 DbGrid* dbgrid,
                 VectorDouble& coor,
                 MatrixDense& rhs);
  void _blankPartRHS(MatrixDense& rhs) const;
  void _fillDual(VectorDouble& zval);
  void _fillDualSimulation(Id nbsimu, MatrixDense& zvals);
  void _calculatePoint(bool flag_grad,
                       DbGrid* dbgrid,
                       VectorDouble& zdual,
                       MatrixDense& rhs,
                       Db* db_target,
                       Id iech0,
                       VectorDouble& result);
  void _convertPotentialToLayer(Id isimu,
                                const double* potval,
                                VectorDouble& result) const;
  void _estimateResult(bool flag_grad,
                       DbGrid* dbout,
                       double refpot,
                       VectorDouble& zdual,
                       MatrixDense& rhs,
                       double* potval);
  void _estimateData(DbGrid* dbout,
                     double refpot,
                     VectorDouble& zdual,
                     MatrixDense& rhs,
                     Db* db_target,
                     VectorInt& uid_pot,
                     VectorInt& uid_grad);
  void _convertDistance(Id ic0,
                        Id j0,
                        VectorDouble& zval,
                        MatrixSymmetric& lhs_orig_arg,
                        MatrixDense& rhs_arg,
                        double* dist_euc,
                        double* dist_geo);
  void _xvalidCalculate(MatrixSymmetric& lhs,
                        bool flag_dist_conv,
                        VectorDouble& zval,
                        MatrixSymmetric& lhs_orig,
                        MatrixDense& rhs);
  static void _tempere(DbGrid* dbout,
                       Id iech,
                       double dist_tempere,
                       double reskrige,
                       VectorDouble& result);
  void _simcond(double dist_tempere,
                bool flag_trans,
                Id nbsimu,
                DbGrid* dbout,
                double refpot,
                double* potsim,
                VectorDouble& zdual,
                MatrixDense& zduals,
                MatrixDense& rhs);
  void _printResult(Id isimu, double* result, double tgtval) const;
  void _checkData(DbGrid* dbgrid,
                  Id isimu,
                  Id nbsimu,
                  double refpot,
                  VectorDouble& zdual,
                  MatrixDense& rhs);
  static Id _distanceToIsoline(DbGrid* dbout);
  double _evaluateRefPot(DbGrid* dbgrid,
                         VectorDouble& zdual,
                         MatrixDense& rhs);
  void _evaluatePotential(DbGrid* dbgrid,
                          double refpot,
                          Id isimu,
                          Id nbsimu,
                          VectorDouble& zdual,
                          MatrixDense& rhs,
                          double* potval);

  Id _extdriftCreateDb(DbGrid* dbout);
  Id _extdriftCreateModel();
  MatrixDense _extdriftBuildRHS();
  Id _extdriftEstablish(DbGrid* dbout);
  Id _extdriftNeigh(DbGrid* dbgrid);
  Id _extdriftEval(double x0,
                   double y0,
                   double z0,
                   Db* db,
                   double* extval,
                   VectorDouble& extgrd);

private:
  Db* _dbiso;
  Db* _dbgrd;
  Db* _dbtgt;
  DbGrid* _dbExt;
  ModelGeneric* _model;
  ModelGeneric* _modelExt;

  Id _ndim;    /* Space dimension */
  Id _niso;    /* Number of Iso-potential information */
  Id _nlayers; /* Number of Iso-potential values */
  Id _ngrd;    /* Number of gradient information */
  Id _ntgt;    /* Number of tangent information */
  Id _next;    /* Number of external drifts */
  Id _nequa;   /* Number of equations in the System */
  Id _order;   /* Order of the drift */
  Id _nring;
  Id _nfull;
  Id _sizeIso;            /* Matrix size linked to iso-potential */
  Id _sizeGrd;            /* Matrix size linked to gradient */
  Id _sizeTgt;            /* Matrix size linked to tangent */
  Id _sizeDrf;            /* Matrix size linked to Drift functions */
  Id _sizeExt;            /* Matrix size linked to External Drifts */
  Id _startIso;           /* Address of the first iso-potential */
  Id _startGrd;           /* Address of the first gradient */
  Id _startTgt;           /* Address of the first tangent */
  Id _startDrf;           /* Address of the first drift */
  Id _startExt;           /* Address of the first external drift */
  double _nugget_grd;     /* Nugget effect for gradients */
  double _nugget_tgt;     /* Nugget effect for tangents */
  double _rangeExt;       /* Range for external drift */
  VectorInt _nbPerLayer;  /* Array of counts of samples per layer */
  VectorInt _ptrPerLayer; /* Array of ptr per layer */
  VectorInt _rankIso;     /* Array of ranks for Iso-potential */
  VectorInt _rankGrd;     /* Array of ranks for Gradients */
  VectorInt _rankTgt;     /* Array of ranks for Tangents */
  Id _optionPart;
  bool _flagPot;
  bool _flagGrad;
  bool _flagTrans;
  bool _verbose;
  VectorInt _indg;
  VectorInt _indg0;
  VectorDouble _dataExt; // Dimension: nech
  MatrixDense _wgtExt;   // Dimension: nech * 4
  mutable SpacePoint _p1;
  mutable SpacePoint _p2;
};

GSTLEARN_EXPORT Id krigingPotential(Db* dbiso,
                                    Db* dbgrd,
                                    Db* dbtgt,
                                    DbGrid* dbout,
                                    ModelGeneric* model,
                                    double nugget_grd   = 0.,
                                    double nugget_tgt   = 0.,
                                    bool flag_pot       = true,
                                    bool flag_grad      = false,
                                    bool flag_trans     = false,
                                    bool flag_save_data = false,
                                    Id option_part      = 0,
                                    bool verbose        = false);
GSTLEARN_EXPORT Id simulatePotential(Db* dbiso,
                                     Db* dbgrd,
                                     Db* dbtgt,
                                     DbGrid* dbout,
                                     ModelGeneric* model,
                                     double nugget_grd   = 0.,
                                     double nugget_tgt   = 0.,
                                     double dist_tempere = TEST,
                                     bool flag_trans     = false,
                                     Id seed             = 135674,
                                     Id nbsimu           = 1,
                                     Id nbtuba           = 100,
                                     bool verbose        = false);
GSTLEARN_EXPORT Id xvalidPotential(Db* dbiso,
                                   Db* dbgrd,
                                   Db* dbtgt,
                                   ModelGeneric* model,
                                   double nugget_grd   = 0.,
                                   double nugget_tgt   = 0.,
                                   bool flag_dist_conv = false,
                                   bool verbose        = false);
} // namespace gstlrn
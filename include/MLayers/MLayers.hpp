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
#include "Enum/ELoc.hpp"
#include "Variogram/Vario.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
  class Db;
  class Model;
  class MatrixSquare;
  class Vario;

  class GSTLEARN_EXPORT MLayers: public AStringable
  {
  public:
    MLayers();
    MLayers(const MLayers& m);
    MLayers& operator=(const MLayers& m);
    virtual ~MLayers();

    /// Interface to AStringable
    String toString(const AStringFormat* strfmt = nullptr) const override;

    void setGeneral(Db* dbin,
                    DbGrid* dbout,
                    Model* model,
                    bool flag_std,
                    bool flag_same,
                    bool flag_vel,
                    bool flag_cumul,
                    bool flag_ext,
                    bool flag_Z,
                    bool flag_bayes,
                    bool match_time,
                    Id irf_rank,
                    const String& namerefd,
                    const String& namereft,
                    const String& namerefb);
    bool checkValid();
    Id kriging(bool verbose);
    Id vario(Vario* vario, bool verbose);
    Id calculatePrior();

  private:
    Id _getNDrift() const;
    Id _locateSampleInDbout(Id iech) const;
    void _isValidLayer(const char* string, Id ilayer0) const;
    Id _getPropsResult(Id iech, Id ilayer0, VectorDouble& props) const;
    Id _getPropsData(Id iech, Id ilayer0, VectorDouble& props) const;
    double _getDriftResult(Id iech, Id ilayer0);
    double _getDriftData(Id iech, Id ilayer0);
    void _getC00(const VectorDouble& prop1,
                 MatrixSquare& covtab,
                 VectorDouble& c00);
    double _getCIJ(Id ilayer,
                   const VectorDouble& prop1,
                   Id jlayer,
                   const VectorDouble& prop2,
                   const VectorDouble& dd,
                   MatrixSquare& covtab);
    double _getCI0(Id ilayer,
                   const VectorDouble& prop1,
                   Id jlayer,
                   const VectorDouble& dd,
                   MatrixSquare& covtab);
    Id _fillDrift(const VectorDouble& coor,
                  double propval,
                  double drext,
                  Id* ipos_loc,
                  VectorDouble& b) const;
    Id _fillDriftData(VectorInt& seltab,
                      VectorDouble& prop1,
                      VectorDouble& drift,
                      MatrixDense& fftab);
    Id _computeLhsOne(VectorInt& seltab,
                      Id iech0,
                      Id ilayer0,
                      VectorDouble& coor,
                      VectorDouble& prop0,
                      VectorDouble& prop2,
                      MatrixSquare& covtab,
                      VectorDouble& b);
    Id _computeLHS(VectorInt& seltab,
                   VectorDouble& prop1,
                   VectorDouble& prop2,
                   MatrixSquare& covtab,
                   MatrixSquare& a,
                   MatrixSquare& acov);
    Id _computeRHS(VectorDouble& coor,
                   VectorInt& seltab,
                   Id iechout,
                   Id ilayer0,
                   VectorDouble& prop0,
                   VectorDouble& prop2,
                   MatrixSquare& covtab,
                   VectorDouble& b);
    void _getDataVector(VectorInt& seltab, VectorDouble& zval);
    Id _subtractOptimalDrift(bool verbose,
                             VectorInt& seltab,
                             VectorDouble& zval);
    Id _getCloseSample(Id iech0,
                       const VectorDouble& coor,
                       double eps = EPSILON5);
    Id _prepareCollocation(Id iechout,
                           VectorDouble& coor,
                           VectorInt& seltab,
                           MatrixSquare* a,
                           VectorDouble& zval,
                           VectorDouble& prop1,
                           VectorDouble& prop2,
                           MatrixSquare& covtab,
                           VectorDouble& b2,
                           VectorDouble& baux,
                           double* ratio);
    void _estimateRegular(double c00,
                          MatrixSquare* a,
                          VectorDouble& b,
                          VectorDouble& dual,
                          VectorDouble& wgt,
                          double& estim,
                          double& stdev) const;
    void _estimateBayes(double c00,
                        const MatrixSquare* acov,
                        VectorDouble& zval,
                        VectorDouble& b,
                        VectorDouble& wgt,
                        VectorDouble& post_mean,
                        MatrixDense& a0,
                        MatrixSquare& cc,
                        MatrixDense& ss,
                        const MatrixSquare& gs,
                        double& estim,
                        double& stdev) const;
    void _estimate(VectorInt& seltab,
                   MatrixSquare* a,
                   VectorDouble& zval,
                   VectorDouble& dual,
                   VectorDouble& prop1,
                   VectorDouble& prop2,
                   MatrixSquare& covtab,
                   VectorDouble& b,
                   VectorDouble& b2,
                   VectorDouble& baux,
                   VectorDouble& wgt,
                   VectorDouble& c00,
                   MatrixDense& /*fftab*/,
                   MatrixDense& a0,
                   MatrixSquare& cc,
                   MatrixDense& ss,
                   MatrixSquare& gs,
                   VectorDouble& post_mean);
    void _checkAuxiliaryVariables(VectorInt& seltab);
    void _convertResults();
    Id _calculateDriftBayes(bool verbose,
                            MatrixSquare* acov,
                            VectorDouble& zval,
                            MatrixDense& fftab,
                            MatrixDense& a0,
                            MatrixSquare& cc,
                            MatrixDense& ss,
                            MatrixSquare& gs,
                            VectorDouble& post_mean,
                            MatrixSquare& post_vars) const;
    Id _evaluateLag(Vario_Order* vorder,
                    Id ifirst,
                    Id ilast,
                    VectorDouble& zval,
                    Id* nval,
                    double* distsum,
                    VectorInt& stat,
                    VectorDouble& phia,
                    VectorDouble& phib,
                    MatrixSquare& atab,
                    VectorDouble& btab);
    Id _getVarioCHH(Vario_Order* vorder,
                    VectorDouble& zval,
                    Id idir,
                    Vario* vario);
    VectorInt _establishSelection(VectorDouble& props) const;
    static Id _getPrior(Id nech,
                        Id npar,
                        VectorDouble& zval,
                        MatrixDense& fftab,
                        VectorDouble& mean,
                        MatrixSquare& vars);

  private:
    bool _flagSame; /* True if input and output files coincide */
    bool _flagVel; /* True for velocity; False for thickness */
    bool _flagCumul; /* True for cumulating in depth */
    bool _flagExt; /* Use external drift */
    bool _flagZ; /* True if output must be converted into depth */
    bool _flagBayes; /* Use Bayesian drift estimate */
    bool _flagStd; /* Ask for calculation of Std. Error of Estimation */
    Id _colrefD; /* Reference depth map (if >= 0) */
    Id _colrefT; /* Reference time map (if >= 0) */
    Id _colrefB; /* Bottom map (if >= 0) */
    Id _matchTime; /* True if Time provided through External Drift */
    ELoc _ptime; /* Pointer to the Time variables */
    Id _nlayers; /* Number of layers */
    Id _nbfl; /* Number of drift functions */
    Id _nech; /* Number of active samples */
    Id _nechmax; /* Total number of samples */
    Id _npar; /* Nb. of layers * Nb. of drift functions */
    Id _neq; /* Total count of equations */
    Id _irfRank; /* Rank of the IRF */
    Db* _dbin;
    DbGrid* _dbout;
    const Model* _model;
  };

  GSTLEARN_EXPORT Id multilayers_kriging(
    Db* dbin,
    DbGrid* dbout,
    Model* model,
    bool flag_same = false,
    bool flag_Z = true,
    bool flag_vel = false,
    bool flag_cumul = false,
    bool flag_ext = false,
    bool flag_bayes = false,
    bool flag_std = false,
    bool match_time = false,
    Id irf_rank = 0,
    const String& namerefd = String(),
    const String& namereft = String(),
    const String& namerefb = String(),
    bool verbose = false,
    const NamingConvention& namconv = NamingConvention("MLayerKriging"));
  GSTLEARN_EXPORT Id multilayers_vario(Db* dbin,
                                       DbGrid* dbout,
                                       Vario* vario,
                                       bool flag_vel = false,
                                       bool flag_ext = false,
                                       Id irf_rank = 0,
                                       bool match_time = false,
                                       const String& namerefd = String(),
                                       const String& namereft = String(),
                                       bool verbose = false);
  GSTLEARN_EXPORT Id multilayers_getPrior(Db* dbin,
                                          DbGrid* dbout,
                                          Model* model,
                                          bool flag_same = false,
                                          bool flag_vel = false,
                                          bool flag_ext = false,
                                          Id irf_rank = 0,
                                          bool match_time = false,
                                          const String& namerefd = String(),
                                          const String& namereft = String(),
                                          const String& namerefb = String());

} // namespace gstlrn

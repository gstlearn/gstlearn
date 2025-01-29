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

#include "Basic/AFunctional.hpp"
#include "Basic/AStringable.hpp"
#include "Covariances/TabNoStat.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Model/CovInternal.hpp"
#include "gstlearn_export.hpp"
#include "geoslib_define.h"
#include "Basic/NamingConvention.hpp"
#include "Space/ASpaceObject.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovContext.hpp"
#include "Space/SpacePoint.hpp"
#include "Space/ASpace.hpp"

#include <vector>

class Db;
class DbGrid;
class MatrixSquareGeneral;
class MatrixSparse;
class TabNoStat;
class AFunctional;
class CovInternal;

/**
 * \brief
 * Class containing the Covariance part of the Model.
 *
 * It is the uppermost class of the Covariance Tree and is conceived as simple as possible on purpose
 * (in order to let the user defined its own version if necessary): it must simply be able to return its value
 * between two end-point (see eval method).
 *
 * It is mainly implemented in CovAniso.hpp or CovAnisoList.hpp
 */
class GSTLEARN_EXPORT ACov : public ASpaceObject
{
public:
  ACov(const CovContext& ctxt = CovContext());
  ACov(const ACov &r);
  ACov& operator=(const ACov &r);
  virtual ~ACov();

  /// ACov Interface
  virtual int getNVar() const {return _ctxt.getNVar();};
  virtual bool isIndexable() const { return false; }
  virtual bool isNoStat() const { return false; }

  const CovContext& getContext() const { return _ctxt; }
  void setContext(const CovContext& ctxt);
  void updateFromContext() { _updateFromContext(); }
  void copyCovContext(const CovContext& ctxt){ _copyCovContext(ctxt);}
  void initFromContext(){   _initFromContext(); }
  CovContext  getContextCopy() const { return CovContext(_ctxt); }
  /// Calculate the covariance between two variables for 0-distance (stationary case)
  virtual double eval0(int ivar = 0,
                       int jvar = 0,
                       const CovCalcMode* mode = nullptr) const;
  /// Calculate the matrix of covariances for 0-distance (stationary case)
  
  virtual void eval0CovMatBiPointInPlace(MatrixSquareGeneral &mat,
                                 const CovCalcMode *mode = nullptr) const;
  virtual void addEval0CovMatBiPointInPlace(MatrixSquareGeneral &mat,
                                            const CovCalcMode *mode = nullptr) const;
  /// Calculate the covariance between two variables and two points (general case)
  virtual double eval(const SpacePoint& p1,
                      const SpacePoint& p2,
                      int ivar = 0,
                      int jvar = 0,
                      const CovCalcMode* mode = nullptr) const = 0;
  /// Calculate the matrix of covariances between two points (general case)
  virtual void evalCovMatBiPointInPlace(MatrixSquareGeneral &mat,
                                        const SpacePoint &p1,
                                        const SpacePoint &p2,
                                        const CovCalcMode *mode = nullptr) const; 
                                        
  virtual void addEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,
                               const SpacePoint& pwork1, 
                               const SpacePoint& pwork2,
                               const CovCalcMode *mode) const;
                               
  void evalCovKriging(MatrixSquareGeneral &mat,
                      SpacePoint &pwork1,
                      SpacePoint& pout, 
                      const CovCalcMode *mode = nullptr) const;
  virtual double evalCovOnSphere(double alpha,
                                 int degree = 50,
                                 bool flagScaleDistance = false,
                                 const CovCalcMode* mode = nullptr) const
  {
    DECLARE_UNUSED(alpha);
    DECLARE_UNUSED(degree);
    DECLARE_UNUSED(flagScaleDistance);
    DECLARE_UNUSED(mode);
    return TEST;
  }
  virtual VectorDouble evalSpectrumOnSphere(int n,
                                            bool flagNormDistance = false,
                                            bool flagCumul = false) const
  {
    DECLARE_UNUSED(n);
    DECLARE_UNUSED(flagNormDistance);
    DECLARE_UNUSED(flagCumul);
    return VectorDouble();
  }
  virtual double evalSpectrum(const VectorDouble &freq,
                              int ivar,
                              int jvar) const
  {
    DECLARE_UNUSED(freq);
    DECLARE_UNUSED(ivar);
    DECLARE_UNUSED(jvar);
    return TEST;
  }

  virtual void updateCovByPoints(int icas1, int iech1, int icas2, int iech2) const 
  {
    DECLARE_UNUSED(icas1);
    DECLARE_UNUSED(iech1);
    DECLARE_UNUSED(icas2);
    DECLARE_UNUSED(iech2);
  }
  void   attachNoStatDb(const Db* db);
  /////////////////////////////////////////////////////////////////////////////////
  ///

  virtual void optimizationPreProcess(const std::vector<SpacePoint>& p,
                               std::vector<SpacePoint> &p1As) const
  {
    DECLARE_UNUSED(p,p1As)
  }

  void optimizationSetTarget(const SpacePoint &pt) const;
  virtual void optimizationSetTargetByIndex(int iech) const {DECLARE_UNUSED(iech)};
  void optimizationPreProcess(const Db* db) const;
  void optimizationPreProcess(const std::vector<SpacePoint>& p) const;

  void optimizationPostProcess() const;
  virtual bool isOptimEnabled() const {return _isOptimEnabled();}
  virtual MatrixRectangular evalCovMatOptim(const Db* db1,
                                            const Db* db2,
                                            int ivar0               = -1,
                                            int jvar0               = -1,
                                            const VectorInt& nbgh1  = VectorInt(),
                                            const VectorInt& nbgh2  = VectorInt(),
                                            const CovCalcMode* mode = nullptr,
                                            bool cleanOptim         = true) const;
  virtual MatrixRectangular evalCovMatOptimByRanks(const Db* db1,
                                                   const Db* db2,
                                                   const VectorVectorInt& sampleRanks1,
                                                   int ivar0               = -1,
                                                   int jvar0               = -1,
                                                   int iech2               = -1,
                                                   const CovCalcMode* mode = nullptr,
                                                   bool cleanOptim         = true) const;
  virtual MatrixSquareSymmetric evalCovMatSymOptim(const Db* db1,
                                                   const VectorInt& nbgh1  = VectorInt(),
                                                   int ivar0               = -1,
                                                   const CovCalcMode* mode = nullptr,
                                                   bool cleanOptim         = true) const;
  virtual MatrixSquareSymmetric
  evalCovMatSymOptimByRanks(const Db* db1,
                            const VectorVectorInt& sampleRanks1,
                            int ivar0               = -1,
                            const CovCalcMode* mode = nullptr,
                            bool cleanOptim         = true) const;

  VectorDouble eval(const std::vector<SpacePoint>& vec_p1,
                    const std::vector<SpacePoint>& vec_p2,
                    int ivar                = 0,
                    int jvar                = 0,
                    const CovCalcMode* mode = nullptr) const;
  MatrixSquareGeneral eval0Mat(const CovCalcMode* mode = nullptr) const;
  MatrixSquareGeneral evalMat(const SpacePoint& p1,
                              const SpacePoint& p2,
                              const CovCalcMode* mode = nullptr) const;
  MatrixRectangular evalCovMat(const Db* db1,
                               const Db* db2           = nullptr,
                               int ivar0               = -1,
                               int jvar0               = -1,
                               const VectorInt& nbgh1  = VectorInt(),
                               const VectorInt& nbgh2  = VectorInt(),
                               const CovCalcMode* mode = nullptr,
                               bool cleanOptim         = true) const;
  MatrixRectangular evalCovMatByRanks(const Db* db1,
                                      const Db* db2,
                                      const VectorVectorInt& sampleRanks1,
                                      int ivar0               = -1,
                                      int jvar0               = -1,
                                      const int iech2         = 0,
                                      const CovCalcMode* mode = nullptr,
                                      bool cleanOptim         = true) const;
  MatrixSquareSymmetric evalCovMatSym(const Db* db1,
                                      const VectorInt& nbgh1  = VectorInt(),
                                      int ivar0               = -1,
                                      const CovCalcMode* mode = nullptr,
                                      bool cleanOptim         = true) const;
  MatrixSquareSymmetric evalCovMatSymByRanks(const Db* db1,
                                             const VectorVectorInt& sampleRanks1,
                                             int ivar0,
                                             const CovCalcMode* mode,
                                             bool cleanOptim) const;
  double evalIvarIpas(double step,
                      const VectorDouble& dir = VectorDouble(),
                      int ivar                = 0,
                      int jvar                = 0,
                      const CovCalcMode* mode = nullptr) const;
  double evalIvarIpasIncr(const VectorDouble &dincr,
                          int ivar = 0,
                          int jvar = 0,
                          const CovCalcMode* mode = nullptr) const;
  VectorDouble evalIvarNpas(const VectorDouble& vec_step,
                            const VectorDouble& dir = VectorDouble(),
                            int ivar = 0,
                            int jvar = 0,
                            const CovCalcMode* mode = nullptr) const;
  MatrixSquareGeneral evalNvarIpas(double step,
                                   const VectorDouble& dir = VectorDouble(),
                                   const CovCalcMode* mode = nullptr) const;
  MatrixSquareGeneral evalNvarIpasIncr(const VectorDouble &dincr,
                                       const CovCalcMode* mode = nullptr) const;
  double evalIsoIvarIpas(double step,
                         int ivar = 0,
                         int jvar = 0,
                         const CovCalcMode* mode = nullptr) const;
  VectorDouble evalIsoIvarNpas(const VectorDouble& vec_step,
                               int ivar = 0,
                               int jvar = 0,
                               const CovCalcMode* mode = nullptr) const;
  MatrixSquareGeneral evalIsoNvarIpas(double step,
                                      const CovCalcMode* mode = nullptr) const;

  double evalCvv(const VectorDouble& ext,
                 const VectorInt& ndisc,
                 const VectorDouble& angles = VectorDouble(),
                 int ivar = 0,
                 int jvar = 0,
                 const CovCalcMode* mode = nullptr) const;
  double evalCvvShift(const VectorDouble& ext,
                      const VectorInt& ndisc,
                      const VectorDouble& shift,
                      const VectorDouble& angles = VectorDouble(),
                      int ivar = 0,
                      int jvar = 0,
                      const CovCalcMode* mode = nullptr) const;
  MatrixSquareGeneral evalCvvM(const VectorDouble& ext,
                               const VectorInt& ndisc,
                               const VectorDouble& angles = VectorDouble(),
                               const CovCalcMode* mode = nullptr) const;
  double evalCxv(const SpacePoint& p1,
                 const VectorDouble& ext,
                 const VectorInt& ndisc,
                 const VectorDouble& angles = VectorDouble(),
                 const VectorDouble& x0 = VectorDouble(),
                 int ivar = 0,
                 int jvar = 0,
                 const CovCalcMode* mode = nullptr) const;
  double evalCxv(const Db* db,
                 const VectorDouble& ext,
                 const VectorInt& ndisc,
                 const VectorDouble& angles = VectorDouble(),
                 const VectorDouble& x0 = VectorDouble(),
                 int ivar = 0,
                 int jvar = 0,
                 const CovCalcMode* mode = nullptr) const;
  MatrixSquareGeneral evalCxvM(const SpacePoint& p1,
                               const VectorDouble& ext,
                               const VectorInt& ndisc,
                               const VectorDouble& angles = VectorDouble(),
                               const VectorDouble& x0 = VectorDouble(),
                               const CovCalcMode* mode = nullptr) const;
  VectorDouble evalPointToDb(const SpacePoint& p1,
                             const Db* db2,
                             int ivar = 0,
                             int jvar = 0,
                             bool useSel = true,
                             const VectorInt& nbgh2 = VectorInt(),
                             const CovCalcMode* mode = nullptr) const;
  VectorDouble evalPointToDbAsSP(const std::vector<SpacePoint>& p1s,
                                 const SpacePoint& p2,
                                 int ivar = 0,
                                 int jvar = 0,
                                 const CovCalcMode* mode = nullptr) const;
  double evalAverageDbToDb(const Db* db1,
                           const Db* db2,
                           int ivar = 0,
                           int jvar = 0,
                           double eps = 0.,
                           int seed = 434132,
                           const CovCalcMode* mode = nullptr) const;
  double evalAverageIncrToIncr(const VectorVectorDouble& d1,
                               const VectorVectorDouble& d2,
                               int ivar = 0,
                               int jvar = 0,
                               const CovCalcMode* mode = nullptr) const;
  double evalAveragePointToDb(const SpacePoint& p1,
                              const Db* db2,
                              int ivar = 0,
                              int jvar = 0,
                              const CovCalcMode* mode = nullptr) const;
  MatrixSparse* evalCovMatSparse(const Db *db1_arg,
                                    const Db *db2_arg = nullptr,
                                    int ivar0 = -1,
                                    int jvar0 = -1,
                                    const VectorInt &nbgh1 = VectorInt(),
                                    const VectorInt &nbgh2 = VectorInt(),
                                    const CovCalcMode *mode = nullptr,
                                    double eps = EPSILON3) const;
  double extensionVariance(const Db* db,
                           const VectorDouble& ext,
                           const VectorInt& ndisc,
                           const VectorDouble& angles = VectorDouble(),
                           const VectorDouble& x0 = VectorDouble(),
                           int ivar = 0,
                           int jvar = 0) const;
  double samplingDensityVariance(const Db* db,
                                 const VectorDouble& ext,
                                 const VectorInt& ndisc,
                                 const VectorDouble& angles = VectorDouble(),
                                 const VectorDouble& x0 = VectorDouble(),
                                 int ivar = 0,
                                 int jvar = 0) const;
  double specificVolume(const Db *db,
                        double mean,
                        const VectorDouble &ext,
                        const VectorInt &ndisc,
                        const VectorDouble &angles = VectorDouble(),
                        const VectorDouble &x0 = VectorDouble(),
                        int ivar = 0,
                        int jvar = 0) const;
  double coefficientOfVariation(const Db *db,
                                double volume,
                                double mean,
                                const VectorDouble &ext,
                                const VectorInt &ndisc,
                                const VectorDouble &angles = VectorDouble(),
                                const VectorDouble &x0 = VectorDouble(),
                                int ivar = 0,
                                int jvar = 0) const;
  double specificVolumeFromCoV(Db *db,
                               double cov,
                               double mean,
                               const VectorDouble &ext,
                               const VectorInt &ndisc,
                               const VectorDouble &angles = VectorDouble(),
                               const VectorDouble &x0 = VectorDouble(),
                               int ivar = 0,
                               int jvar = 0) const;
  double evaluateOneGeneric(const CovInternal *covint,
                            const VectorDouble &d1 = VectorDouble(),
                            double weight = 1.,
                            const CovCalcMode *mode = nullptr) const;
  double calculateStDev(Db *db1,
                        int iech1,
                        Db *db2,
                        int iech2,
                        bool verbose = false,
                        double factor = 1.,
                        const CovCalcMode *mode = nullptr) const;
  
  VectorDouble evalCovMatV(Db* db1,
                           Db* db2                 = nullptr,
                           int ivar0               = -1,
                           int jvar0               = -1,
                           const VectorInt& nbgh1  = VectorInt(),
                           const VectorInt& nbgh2  = VectorInt(),
                           const CovCalcMode* mode = nullptr) const;
  void evaluateMatInPlace(const CovInternal *covint,
                          const VectorDouble &d1,
                          MatrixSquareGeneral &covtab,
                          bool flag_init = false,
                          double weight = 1.,
                          const CovCalcMode *mode = nullptr) const;
  VectorDouble evaluateFromDb(Db *db,
                              int ivar = 0,
                              int jvar = 0,
                              const CovCalcMode *mode = nullptr) const;
  double evaluateOneIncr(double hh,
                         const VectorDouble &codir = VectorDouble(),
                         int ivar = 0,
                         int jvar = 0,
                         const CovCalcMode *mode = nullptr) const;
  VectorDouble sample(const VectorDouble &h,
                      const VectorDouble &codir = VectorDouble(),
                      int ivar = 0,
                      int jvar = 0,
                      const CovCalcMode* mode = nullptr,
                      const CovInternal* covint = nullptr) const;
  VectorDouble sampleUnitary(const VectorDouble &hh,
                             int ivar = 0,
                             int jvar = 0,
                             VectorDouble codir = VectorDouble(),
                             const CovCalcMode* mode = nullptr);
  VectorDouble envelop(const VectorDouble &hh,
                       int ivar = 0,
                       int jvar = 0,
                       int isign = 1,
                       VectorDouble codir = VectorDouble(),
                       const CovCalcMode* mode = nullptr);
  int buildVmapOnDbGrid(DbGrid *dbgrid, const NamingConvention &namconv = NamingConvention("VMAP")) const;
  double gofToVario(const Vario* vario, bool verbose = true) const;
  static void gofDisplay(double gof,
                         bool byValue                   = true,
                         const VectorDouble& thresholds = {2., 5., 10., 100});
  
  void manage(const Db* db1, const Db* db2) const { _manage(db1, db2); }

  void load(const SpacePoint& p,bool case1) const;

  void loadAndAddEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,const SpacePoint& p1,const SpacePoint&p2,
                                              const CovCalcMode *mode = nullptr) const;

  double loadAndEval(const SpacePoint& p1,
                     const SpacePoint& p2,
                     int ivar,
                     int jvar,
                     const CovCalcMode* mode) const;

  bool checkAndManageNoStatDb(const Db*& db, const String& namecol);

  virtual void updateCovByMesh(int imesh,bool aniso = true) const
  {
    DECLARE_UNUSED(imesh,aniso)
  }
  virtual double getValue(const EConsElem &econs,int iv1,int iv2) const
  {
    DECLARE_UNUSED(econs,iv1,iv2)
    return TEST;
  }
  virtual void makeStationary();
  virtual int makeElemNoStat(const EConsElem &econs, int iv1, int iv2,
                     const AFunctional* func = nullptr, 
                     const Db* db = nullptr,const String& namecol = String());
  void createNoStatTab();
  void informMeshByMesh(const AMesh* amesh) const;
  void informMeshByApex(const AMesh* amesh) const;
  VectorDouble informCoords(const VectorVectorDouble& coords, 
                            const EConsElem& econs,
                            int iv1 = 0, int iv2 = 0) const;
  void informDbIn(const Db* dbin) const;
  void informDbOut(const Db* dbout) const;

  virtual void updateCovByPoints(int icas1, int iech1, int icas2, int iech2)
  {
    DECLARE_UNUSED(icas1);
    DECLARE_UNUSED(iech1);
    DECLARE_UNUSED(icas2);
    DECLARE_UNUSED(iech2);
  }

  int getNDim(int ispace = -1) const { return _ctxt.getNDim(ispace); }

private:
  virtual void _setContext(const CovContext& ctxt) 
  {
    DECLARE_UNUSED(ctxt)
  }
  
  virtual void _manage(const Db* db1,const Db* db2) const 
  {
    DECLARE_UNUSED(db1)
    DECLARE_UNUSED(db2)
  }

  void setNoStatDbIfNecessary(const Db*& db);
private : 
 virtual TabNoStat* _createNoStatTab();

protected:
  void setNVar(int nvar) { _ctxt.setNVar(nvar); }
  virtual void _loadAndAddEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,const SpacePoint& p1,const SpacePoint&p2,
                                              const CovCalcMode *mode = nullptr) const;
  virtual void _optimizationSetTarget(const SpacePoint &pt) const;

  void _setOptimEnabled(bool enabled){ _optimEnabled = enabled;}
  VectorInt _getActiveVariables(int ivar0) const;
  static void _updateCovMatrixSymmetricVerr(const Db* db1,
                                            AMatrix* mat,
                                            const VectorInt& ivars,
                                            const VectorVectorInt& index1);

  virtual void _optimizationPreProcess(const std::vector<SpacePoint>& p) const;
  virtual void _addEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,
                                            const SpacePoint& pwork1, 
                                            const SpacePoint& pwork2,
                                            const CovCalcMode *mode) const;
  double _loadAndEval(const SpacePoint& p1,
                      const SpacePoint& p2,
                      int ivar                = 0,
                      int jvar                = 0,
                      const CovCalcMode* mode = nullptr) const;
  bool _checkDims(int idim, int jdim) const;

  protected:
    virtual void _initFromContext() {};

private:
  virtual void _copyCovContext(const CovContext& ctxt)
  {
    DECLARE_UNUSED(ctxt)
  }

  virtual void _updateFromContext() {};
  virtual void _optimizationPostProcess() const; 
  virtual bool _isOptimEnabled() const {return _optimEnabled;}

  DbGrid* _discretizeBlock(const VectorDouble& ext,
                           const VectorInt& ndisc,
                           const VectorDouble& angles = VectorDouble(),
                           const VectorDouble& x0 = VectorDouble()) const;
  Db* _discretizeBlockRandom(const DbGrid* dbgrid, int seed = 34131) const;
  double _getVolume(const VectorDouble& ext) const;
  

protected:
  CovContext _ctxt;         /* Context */
  bool _optimEnabled;
  mutable bool _isOptimPreProcessed;
  mutable std::vector<SpacePoint> _p1As;
  mutable SpacePoint _p2A;
  const mutable SpacePoint* _pw1;
  const mutable SpacePoint* _pw2;
  
  TabNoStat* _tabNoStat;
};

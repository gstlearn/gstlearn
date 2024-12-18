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
#include "Basic/VectorNumT.hpp"
#include "Covariances/CorAniso.hpp"
#include "Covariances/TabNoStatCovAniso.hpp"
#include "Enum/EConsElem.hpp"
#include "Model/CovInternal.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"
#include "Enum/ECov.hpp"

#include "Basic/ICloneable.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Basic/Tensor.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/ACovFunc.hpp"
#include "Covariances/CovContext.hpp"
#include "Arrays/Array.hpp"
#include "Space/SpacePoint.hpp"
#include <array>
#include <vector>

class Rotation;
class MatrixSquareGeneral;
class MatrixRectangular;
class CovInternal;
/**
 * \brief
 * This class describes an **elementary covariance**.
 *
 * This covariance is described through the following list of parameters:
 * - the covariance **type**: the list of these types is provided in ECov.hpp
 * - the largest set of parameters for any covariance: **range(s)**, **anisotropy angle(s)**, **third parameter**. Some of these parameters
 * do not make sense, depending on the covariance type: e.g. the range for nugget effect, the third parameter for a spherical
 * structure, ...
 * All these parameters are processed and stored as a **tensor** in order to avoid repetitive calculations.
 * - the **sill**. This comes as a square symmetric matrix whose dimension is equal to the number of variables.
 */
class GSTLEARN_EXPORT CovAniso: public ACov, public ICloneable
{
public:
  CovAniso(const ECov& type, const CovContext& ctxt);
  CovAniso(const String& symbol, const CovContext& ctxt);
  CovAniso(const ECov& type,
           double range,
           double param,
           double sill,
           const CovContext& ctxt,
           bool flagRange = true);
  CovAniso(const CovAniso& r);
  CovAniso& operator=(const CovAniso& r);
  virtual ~CovAniso();

  /// ICloneable Interface
  IMPLEMENT_CLONING(CovAniso)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ASpaceObject Interface
  virtual bool isConsistent(const ASpace* space) const override;

  /// ACov Interface
  virtual int getNVariables() const override { return _ctxt.getNVar(); }


  /// ACov Interface
  virtual double eval0(int ivar = 0,
                       int jvar = 0,
                       const CovCalcMode* mode = nullptr) const override;
  virtual double eval(const SpacePoint& p1,
                      const SpacePoint& p2,
                      int ivar = 0,
                      int jvar = 0,
                      const CovCalcMode* mode = nullptr) const override;
 
  virtual void addEval0CovMatBiPointInPlace(MatrixSquareGeneral &mat,
                               const CovCalcMode *mode = nullptr) const override;
 
  virtual double evalCovOnSphere(double alpha,
                                 int degree = 50,
                                 bool flagScaleDistance = true,
                                 const CovCalcMode* mode = nullptr) const override;
  virtual VectorDouble evalSpectrumOnSphere(int n,
                                            bool flagNormDistance = false,
                                            bool flagCumul = false) const override;
  virtual double evalSpectrum(const VectorDouble &freq,
                              int ivar = 0,
                              int jvar = 0) const override;

  virtual double getIntegralRange(int ndisc, double hmax) const;
  virtual String getFormula() const { return _cor.getFormula(); }
  virtual double getBallRadius() const { return TEST; }

  bool isOptimizationInitialized(const Db* db = nullptr) const;
  void _optimizationPreProcess(const std::vector<SpacePoint>& p) const override;
  void optimizationSetTargetByIndex(int iech) const override;

  void evalOptimInPlace(MatrixRectangular& res,
                        const VectorInt& ivars,
                        const VectorVectorInt& index,
                        int ivar2 = 0,
                        int icol = 0,
                        const CovCalcMode *mode = nullptr,
                        bool flagSym = false) const;
  
  bool isValidForTurningBand() const;
  double simulateTurningBand(double t0, TurningBandOperate &operTB) const;
  bool isValidForSpectral() const ;
  MatrixRectangular simulateSpectralOmega(int nb) const;

  static CovAniso* createIsotropic(const CovContext& ctxt,
                                   const ECov& type,
                                   double range,
                                   double sill = 1.,
                                   double param = 1.,
                                   bool flagRange = true);
  static CovAniso* createAnisotropic(const CovContext& ctxt,
                                     const ECov& type,
                                     const VectorDouble& ranges,
                                     double sill = 1.,
                                     double param = 1.,
                                     const VectorDouble& angles = VectorDouble(),
                                     bool flagRange = true);
  static CovAniso* createIsotropicMulti(const CovContext& ctxt,
                                        const ECov& type,
                                        double range,
                                        const MatrixSquareSymmetric& sills,
                                        double param = 1.,
                                        bool flagRange = true);
  static CovAniso* createAnisotropicMulti(const CovContext& ctxt,
                                          const ECov& type,
                                          const VectorDouble& ranges,
                                          const MatrixSquareSymmetric& sills,
                                          double param = 1.,
                                          const VectorDouble& angles = VectorDouble(),
                                          bool flagRange = true);

  void setContext(const CovContext& ctxt);
  void setParam(double param);
  void copyCovContext(const CovContext& ctxt);

  void setSill(double sill); /// Only valid when there is only one variable (in the context)
  void setSill(const MatrixSquareSymmetric& sill);
  void setSill(const VectorDouble& sill);
  void setSill(int ivar, int jvar, double sill);
  void initSill(double value = 0.);

  /// Practical range
  void setRangeIsotropic(double range);
  void setRange(int idim, double range);
  void setRanges(const VectorDouble& ranges);

  void setScale(double scale); /// Make the covariance isotropic
  void setScale(int idim, double scale);
  void setScales(const VectorDouble& scales);

  void setAnisoRotation(const Rotation& rot);
  void setAnisoRotation(const VectorDouble& rot);
  void setAnisoAngles(const VectorDouble& angles);
  void setAnisoAngle(int idim, double angle);

  void setRotationAnglesAndRadius(const VectorDouble& angles = VectorDouble(),
                                  const VectorDouble& ranges = VectorDouble(),
                                  const VectorDouble& scales = VectorDouble());

  const MatrixSquareSymmetric& getSill() const { return _sill; }
  double getSill(int ivar, int jvar) const;
  double getSlope(int ivar, int jvar) const;
  VectorDouble getRanges() const;
  const Rotation& getAnisoRotation() const { return _cor.getAniso().getRotation(); }
  const VectorDouble& getScales() const { return _cor.getAniso().getRadius(); }

  void   setType(const ECov& type);
  double getRange() const;
  double getScale() const;
  bool   getFlagAniso() const { return !isIsotropic(); }
  bool   getFlagRotation() const { return hasRotation(); }
  double getRange(int idim) const { return getRanges()[idim]; }
  double getScale(int idim) const { return getScales()[idim]; }
  VectorDouble getAnisoAngles() const { return _cor.getAniso().getAngles(); }
  const MatrixSquareGeneral& getAnisoRotMat() const { return _cor.getAniso().getMatrixDirect(); }
  const MatrixSquareGeneral& getAnisoInvMat() const { return _cor.getAniso().getMatrixInverse(); }
  VectorDouble getAnisoCoeffs() const;
  double getAnisoAngles(int idim) const { return getAnisoAngles()[idim]; }
  double getAnisoRotMat(int idim, int jdim) const { return _cor.getAniso().getMatrixDirect().getValue(idim,jdim); }
  double getAnisoCoeffs(int idim) const { return getAnisoCoeffs()[idim]; }
  const CovContext& getContext() const { return _ctxt; }
  const ECov& getType() const { return _cor.getType(); }
  double getParam() const;
  double getScadef() const { return _cor.getScadef(); }
  double getParMax() const { return _cor.getParMax(); }
  int    getMaxNDim() const { return _cor.getMaxNDim(); }
  int    getMinOrder() const { return _cor.getMinOrder(); }
  bool   hasInt1D() const { return _cor.hasInt1D(); }
  bool   hasInt2D() const { return _cor.hasInt2D(); }
  int    hasRange() const { return _cor.hasRange(); }
  int    hasParam() const  { return _cor.hasParam(); }
  String getCovName() const { return _cor.getCovName(); }
  bool   isIsotropic() const { return _cor.getAniso().isIsotropic(); }
  bool   isAsymptotic() const { return getScadef() != 1.; }
  bool   hasRotation() const { return _cor.getAniso().hasRotation(); }
  const Tensor& getAniso() const { return _cor.getAniso(); }
  void   setAniso(const Tensor& aniso) { _cor.setAniso(aniso); }
  const ACovFunc* getCova() const { return _cor.getCova(); }
  int    getGradParamNumber() const;
  bool   hasCovDerivative() const { return _cor.hasCovDerivative(); }
  bool   hasCovOnSphere() const { return _cor.hasCovOnSphere(); }
  bool   hasSpectrumOnSphere() const { return _cor.hasSpectrumOnSphere(); }
  bool   hasMarkovCoeffs() const { return _cor.hasMarkovCoeffs(); }
  bool   hasSpectrumOnRn() const { return _cor.hasSpectrumOnRn(); }
  double normalizeOnSphere(int n = 50) const;
  
  //////////////////////// New NoStat methods //////////////////////////
  void   attachNoStatDb(const Db* db);
 
  void   makeRangeNoStatDb( const String &namecol, int idim = 0,              const Db* db = nullptr);
  void   makeScaleNoStatDb( const String &namecol, int idim = 0,              const Db* db = nullptr);
  void   makeAngleNoStatDb( const String &namecol, int idim = 0,              const Db* db = nullptr);
  void   makeSillNoStatDb(  const String &namecol, int ivar = 0, int jvar = 0,const Db* db = nullptr);
  void   makeTensorNoStatDb(const String &namecol, int idim = 0, int jdim = 0,const Db* db = nullptr);
  void   makeParamNoStatDb( const String &namecol,                            const Db* db = nullptr);
  void   makeRangeNoStatFunctional( const AFunctional *func, int idim = 0);
  void   makeScaleNoStatFunctional( const AFunctional *func, int idim = 0);
  void   makeAngleNoStatFunctional( const AFunctional *func, int idim = 0);
  void   makeSillNoStatFunctional(  const AFunctional *func, int ivar = 0, int jvar = 0);
  void   makeTensorNoStatFunctional(const AFunctional *func, int idim = 0, int jdim = 0);
  void   makeParamNoStatFunctional( const AFunctional *func);
  void   makeRangeStationary(int idim = 0);
  void   makeScaleStationary(int idim = 0);
  void   makeAngleStationary(int idim = 0);
  void   makeSillStationary( int ivar = 0, int jvar = 0);
  void   makeTensorStationary(int idim, int jdim);
  void   makeParamStationary();
  void   makeStationary();



  int getNAngles() const {return _cor.getNAngles();}
  int getNRanges() const {return _cor.getNRanges();}
  int getNScales() const {return _cor.getNScales();}
  int getNSills()  const {return _tabNoStat.getNSills();}
  bool isNoStat() const  override
  { 
    return _tabNoStat.isNoStat() || _cor.isNoStat(); 
  };
  bool isNoStatForParam()   const {return _cor.isNoStatForParam();}
  bool isNoStatForTensor()  const {return _cor.isNoStatForTensor();}
  bool isNoStatForAnisotropy() const { return _cor.isNoStatForAnisotropy();}
  bool isNoStatForVariance()   const { return _tabNoStat.isDefinedForVariance();}
  bool isNoStatForRotation()   const { return _cor.isNoStatForRotation();}

  
  VectorDouble evalCovOnSphereVec(const VectorDouble &alpha,
                                  int degree = 50,
                                  bool flagScaleDistance = false,
                                  const CovCalcMode* mode = nullptr) const;
  Array evalCovFFT(const VectorDouble& hmax, int N = 128, int ivar = 0, int jvar = 0) const;
  VectorDouble getMarkovCoeffs() const;
  void setMarkovCoeffs(const VectorDouble& coeffs);
  void setMarkovCoeffsBySquaredPolynomials(VectorDouble coeffs1, VectorDouble coeffs2, double eps = 0);
  void computeMarkovCoeffs();
  double getCorrec() const;
  double getFullCorrec() const;
  int getDimensionNumber() const { return _ctxt.getNDim(); }
  void nostatUpdate(CovInternal *covint);

  CovAniso* createReduce(const VectorInt &validVars) const;

  void informMeshByMesh(const AMesh* amesh) const;
  void informMeshByApex(const AMesh* amesh) const;
  VectorDouble informCoords(const VectorVectorDouble& coords, 
                            const EConsElem& econs,
                            int iv1 = 0, int iv2 = 0) const;
  void informDbIn(const Db* dbin) const;
  void informDbOut(const Db* dbout) const;
  void informMeshByMeshForAnisotropy(const AMesh* amesh) const;
  void informMeshByApexForAnisotropy(const AMesh* amesh) const;
  void informDbInForAnisotropy(const Db* dbin) const;
  void informDbOutForAnisotropy(const Db* dbout) const;
  void informMeshByMeshForSills(const AMesh* amesh) const;
  void informMeshByApexForSills(const AMesh* amesh) const;
  void informDbInForSills(const Db* dbin) const;
  void informDbOutForSills(const Db* dbout) const;

  /// Tell if the use of Optimization is enabled or not

  void updateCovByPoints(int icas1, int iech1, int icas2, int iech2) override;
  void updateCovByMesh(int imesh,bool aniso = true);
  void setOptimEnabled(bool flag) const { _optimEnabled = flag; }
  double getValue(const EConsElem& econs, int iv1, int iv2) const;

protected:
  /// Update internal parameters consistency with the context
  virtual void _addEvalCovMatBiPointInPlace(MatrixSquareGeneral& mat,
                                            const SpacePoint& p1,
                                            const SpacePoint& p2,
                                            const CovCalcMode* mode = nullptr) const override;
  virtual void _updateFromContext();
  virtual void _initFromContext();
  void _optimizationSetTarget(const SpacePoint& pt) const override;

private:
void _optimizationPostProcess() const override; 

bool _isOptimEnabled() const override 
{ 
  return _optimEnabled && !isNoStatForAnisotropy(); 
}
void  _evalOptim(SpacePoint* p1A, SpacePoint* p2A,
                 MatrixSquareGeneral &mat,
                 const CovCalcMode *mode) const;
 void _makeElemNoStat(const EConsElem &econs, int iv1, int iv2,
                      const AFunctional* func = nullptr, 
                      const Db* db = nullptr,const String& namecol = String());

  void _manage(const Db* db1,const Db* db2) const override;

  bool _checkSill(int ivar = 0, int jvar = 0) const;
  bool _checkDims(int idim, int jdim) const;

  void _setNoStatDbIfNecessary(const Db*& db);
  bool _checkAndManageNoStatDb(const Db*& db, const String& namecol);
  bool   _isVariableValid(int ivar) const;
  void   _computeCorrec();
  double _getDetTensor() const;
  void   _optimizationTransformSP(const SpacePoint& ptin, SpacePoint& ptout) const;

private:

  CorAniso _cor;
  CovContext _ctxt;                    /// Context (space, number of variables, ...) // TODO : Really store a copy ?
  mutable MatrixSquareSymmetric _sill; /// Sill matrix (nvar x nvar)
  TabNoStat _tabNoStat;
  mutable bool _optimEnabled;
  // These temporary information is used to speed up processing (optimization functions)
  // They are in a protected section as they may be modified by class hierarchy
};


GSTLEARN_EXPORT double scale2range(const ECov& type, double scale, double param = 1.);
GSTLEARN_EXPORT double range2scale(const ECov& type, double range, double param = 1.);

  
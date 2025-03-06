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
#include "Covariances/CovProportional.hpp"
#include "Model/CovInternal.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"
#include "Enum/ECov.hpp"
#include "Covariances/CorAniso.hpp"
#include "Basic/ICloneable.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Basic/Tensor.hpp"
#include "Covariances/ACovFunc.hpp"
#include "Covariances/CovContext.hpp"
#include "Arrays/Array.hpp"
#include "Space/SpacePoint.hpp"

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
class GSTLEARN_EXPORT CovAniso: public CovProportional
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
  
  /// ACov Interface
  virtual double eval0(int ivar = 0,
                       int jvar = 0,
                       const CovCalcMode* mode = nullptr) const override;
  virtual double eval(const SpacePoint& p1,
                      const SpacePoint& p2,
                      int ivar = 0,
                      int jvar = 0,
                      const CovCalcMode* mode = nullptr) const override;
 
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
  virtual String getFormula() const { return _corAniso->getFormula(); }
  virtual double getBallRadius() const { return TEST; }
 
  bool isValidForTurningBand() const;
  double simulateTurningBand(double t0, TurningBandOperate &operTB) const;
  bool isValidForSpectral() const ;
  MatrixRectangular simulateSpectralOmega(int nb) const;

  static CovAniso* createFromParam(const ECov& type,
                                   double range,
                                   double sill,
                                   double param,
                                   const VectorDouble& ranges,
                                   const MatrixSquareSymmetric& sills,
                                   const VectorDouble& angles,
                                   const ASpaceSharedPtr& space,
                                   bool flagRange);
  static CovAniso* createIsotropic(const CovContext& ctxt,
                                   const ECov& type,
                                   double range,
                                   double sill    = 1.,
                                   double param   = 1.,
                                   bool flagRange = true);
  static CovAniso* createAnisotropic(const CovContext& ctxt,
                                     const ECov& type,
                                     const VectorDouble& ranges,
                                     double sill                = 1.,
                                     double param               = 1.,
                                     const VectorDouble& angles = VectorDouble(),
                                     bool flagRange             = true);
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

  void setParam(double param);

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
  
  const CorAniso* getCorAniso() const { return dynamic_cast<const CorAniso*>(getCor()); } 

  double getSlope(int ivar, int jvar) const;
  VectorDouble getRanges() const;
  const Rotation& getAnisoRotation() const { return _corAniso->getAniso().getRotation(); }
  const VectorDouble& getScales() const { return _corAniso->getAniso().getRadius(); }

  void   setType(const ECov& type);
  double getRange() const;
  double getScale() const;
  bool   getFlagAniso() const { return !isIsotropic(); }
  bool   getFlagRotation() const { return hasRotation(); }
  double getRange(int idim) const { return getRanges()[idim]; }
  double getScale(int idim) const { return getScales()[idim]; }
  VectorDouble getAnisoAngles() const { return _corAniso->getAniso().getAngles(); }
  const MatrixSquareGeneral& getAnisoRotMat() const { return _corAniso->getAniso().getMatrixDirect(); }
  const MatrixSquareGeneral& getAnisoInvMat() const { return _corAniso->getAniso().getMatrixInverse(); }
  VectorDouble getAnisoCoeffs() const;
  double getAnisoAngles(int idim) const { return getAnisoAngles()[idim]; }
  double getAnisoRotMat(int idim, int jdim) const { return _corAniso->getAniso().getMatrixDirect().getValue(idim,jdim); }
  double getAnisoCoeffs(int idim) const { return getAnisoCoeffs()[idim]; }
  const CovContext& getContext() const { return _ctxt; }
  const ECov& getType() const { return _corAniso->getType(); }
  double getParam() const;
  double getScadef() const { return _corAniso->getScadef(); }
  double getParMax() const { return _corAniso->getParMax(); }
  int    getMaxNDim() const { return _corAniso->getMaxNDim(); }
  int    getMinOrder() const { return _corAniso->getMinOrder(); }
  bool   hasInt1D() const { return _corAniso->hasInt1D(); }
  bool   hasInt2D() const { return _corAniso->hasInt2D(); }
  int    hasRange() const { return _corAniso->hasRange(); }
  int    hasParam() const  { return _corAniso->hasParam(); }
  String getCovName() const { return _corAniso->getCovName(); }
  bool   isIsotropic() const { return _corAniso->getAniso().isIsotropic(); }
  bool   isAsymptotic() const { return getScadef() != 1.; }
  bool   hasRotation() const { return _corAniso->getAniso().hasRotation(); }
  const Tensor& getAniso() const { return _corAniso->getAniso(); }
  void   setAniso(const Tensor& aniso) { _corAniso->setAniso(aniso); }
  const ACovFunc* getCorFunc() const { return _corAniso->getCorFunc(); }
  int    getNGradParam() const;
  bool   hasCovDerivative() const { return _corAniso->hasCovDerivative(); }
  bool   hasCovOnSphere() const { return _corAniso->hasCovOnSphere(); }
  bool   hasSpectrumOnSphere() const { return _corAniso->hasSpectrumOnSphere(); }
  bool   hasMarkovCoeffs() const { return _corAniso->hasMarkovCoeffs(); }
  bool   hasSpectrumOnRn() const { return _corAniso->hasSpectrumOnRn(); }
  double normalizeOnSphere(int n = 50) const;
  
  //////////////////////// New NoStat methods //////////////////////////
  


  bool isNoStatForParam()   const {return _corAniso->isNoStatForParam();}
  bool isNoStatForTensor()  const {return _corAniso->isNoStatForTensor();}
  bool isNoStatForAnisotropy() const { return _corAniso->isNoStatForAnisotropy();}
 
  bool isNoStatForRotation()   const { return _corAniso->isNoStatForRotation();}

  void   makeRangeNoStatDb( const String &namecol, int idim = 0,              const Db* db = nullptr);
  void   makeScaleNoStatDb( const String &namecol, int idim = 0,              const Db* db = nullptr);
  void   makeAngleNoStatDb( const String &namecol, int idim = 0,              const Db* db = nullptr);
  
  void   makeTensorNoStatDb(const String &namecol, int idim = 0, int jdim = 0,const Db* db = nullptr);
  void   makeParamNoStatDb( const String &namecol,                            const Db* db = nullptr);
  void   makeRangeNoStatFunctional( const AFunctional *func, int idim = 0);
  void   makeScaleNoStatFunctional( const AFunctional *func, int idim = 0);
  void   makeAngleNoStatFunctional( const AFunctional *func, int idim = 0);
  void   makeTensorNoStatFunctional(const AFunctional *func, int idim = 0, int jdim = 0);
  void   makeParamNoStatFunctional( const AFunctional *func);
  void   makeRangeStationary(int idim = 0);
  void   makeScaleStationary(int idim = 0);
  void   makeAngleStationary(int idim = 0);
  
  void   makeTensorStationary(int idim, int jdim);
  void   makeParamStationary();
  VectorDouble evalCovOnSphereVec(const VectorDouble &alpha,
                                  int degree = 50,
                                  bool flagScaleDistance = false,
                                  const CovCalcMode* mode = nullptr) const;
  Array evalCovFFT(const VectorDouble& hmax, int N = 128, int ivar = 0, int jvar = 0) const;
  VectorDouble getMarkovCoeffs() const;
  void setMarkovCoeffs(const VectorDouble& coeffs);
  void setMarkovCoeffsBySquaredPolynomials(const VectorDouble& coeffs1,
                                           const VectorDouble& coeffs2,
                                           double eps = 0);
  void computeMarkovCoeffs();
  double getCorrec() const;
  double getFullCorrec() const;
  int getNDim() const { return _ctxt.getNDim(); }
  
  CovAniso* createReduce(const VectorInt &validVars) const;

  void informDbInForAnisotropy(const Db* dbin) const;
  void informDbOutForAnisotropy(const Db* dbout) const;

  void setOptimEnabled(bool flag) const
  {
    _optimEnabled = flag;
    _corAniso->setOptimEnabled(flag);
  }

  void informMeshByMeshForAnisotropy(const AMesh* amesh) const;
  void informMeshByApexForAnisotropy(const AMesh* amesh) const;
 
  bool _isOptimEnabled() const override 
  { 
    return _optimEnabled && !isNoStatForAnisotropy(); 
  }
  
  bool isNoStat() const  override
  { 
    return _tabNoStat.isNoStat() || _corAniso->isNoStat(); 
  };
  int getNAngles() const {return _corAniso->getNAngles();}
  int getNRanges() const {return _corAniso->getNRanges();}
  int getNScales() const {return _corAniso->getNScales();}

  void   _computeCorrec();
  double _getDetTensor() const;
  double _getSillValue(int ivar, int jvar, const CovCalcMode* mode) const;

private:
  CorAniso* _corAniso;

  mutable bool _optimEnabled;
  // These temporary information is used to speed up processing (optimization functions)
  // They are in a protected section as they may be modified by class hierarchy
};

GSTLEARN_EXPORT double scale2range(const ECov& type, double scale, double param = 1.);
GSTLEARN_EXPORT double range2scale(const ECov& type, double range, double param = 1.);

  
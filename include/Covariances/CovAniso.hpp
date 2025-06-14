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
#include "Covariances/CovProportional.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"
#include "Enum/ECov.hpp"
#include "Covariances/CorAniso.hpp"
#include "Basic/ICloneable.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Basic/Tensor.hpp"
#include "Covariances/ACovFunc.hpp"
#include "Covariances/CovContext.hpp"
#include "Arrays/Array.hpp"
#include "Space/SpacePoint.hpp"

class Rotation;
class MatrixSquare;
class MatrixDense;
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
  virtual double eval0(int ivar                = 0,
                       int jvar                = 0,
                       const CovCalcMode* mode = nullptr) const override;

  virtual double evalCovOnSphere(double alpha,
                                 int degree              = 50,
                                 bool flagScaleDistance  = true,
                                 const CovCalcMode* mode = nullptr) const override;
  virtual VectorDouble evalSpectrumOnSphere(int n,
                                            bool flagNormDistance = false,
                                            bool flagCumul        = false) const override;
  virtual double evalSpectrum(const VectorDouble& freq,
                              int ivar = 0,
                              int jvar = 0) const override;

  virtual double getIntegralRange(int ndisc, double hmax) const;
  virtual String getFormula() const { return getCorAniso()->getFormula(); }
  virtual double getBallRadius() const { return TEST; }

  static CovAniso* createFromParam(const ECov& type,
                                   double range,
                                   double sill,
                                   double param,
                                   const VectorDouble& ranges,
                                   const MatrixSymmetric& sills,
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
                                        const MatrixSymmetric& sills,
                                        double param   = 1.,
                                        bool flagRange = true);
  static CovAniso* createAnisotropicMulti(const CovContext& ctxt,
                                          const ECov& type,
                                          const VectorDouble& ranges,
                                          const MatrixSymmetric& sills,
                                          double param               = 1.,
                                          const VectorDouble& angles = VectorDouble(),
                                          bool flagRange             = true);

  FORWARD_METHOD_NON_CONST(getCorAniso, setParam)
  FORWARD_METHOD_NON_CONST(getCorAniso, computeMarkovCoeffs)
  FORWARD_METHOD_NON_CONST(getCorAniso, setRangeIsotropic)
  FORWARD_METHOD_NON_CONST(getCorAniso, setRange)
  FORWARD_METHOD_NON_CONST(getCorAniso, setRanges)
  FORWARD_METHOD_NON_CONST(getCorAniso, setScale)
  FORWARD_METHOD_NON_CONST(getCorAniso, setScaleDim)
  FORWARD_METHOD_NON_CONST(getCorAniso, setScales)
  FORWARD_METHOD_NON_CONST(getCorAniso, setAnisoRotationMat)
  FORWARD_METHOD_NON_CONST(getCorAniso, setAnisoRotation)
  FORWARD_METHOD_NON_CONST(getCorAniso, setAnisoAngles)
  FORWARD_METHOD_NON_CONST(getCorAniso, setAnisoAngle)
  FORWARD_METHOD_NON_CONST(getCorAniso, setRotationAnglesAndRadius)
  FORWARD_METHOD_NON_CONST(getCorAniso, setType)

  FORWARD_METHOD_NON_CONST(getCorAniso, makeRangeNoStatDb)
  FORWARD_METHOD_NON_CONST(getCorAniso, makeScaleNoStatDb)
  FORWARD_METHOD_NON_CONST(getCorAniso, makeAngleNoStatDb)

  FORWARD_METHOD_NON_CONST(getCorAniso, makeTensorNoStatDb)
  FORWARD_METHOD_NON_CONST(getCorAniso, makeParamNoStatDb)
  FORWARD_METHOD_NON_CONST(getCorAniso, makeRangeNoStatFunctional)
  FORWARD_METHOD_NON_CONST(getCorAniso, makeScaleNoStatFunctional)
  FORWARD_METHOD_NON_CONST(getCorAniso, makeAngleNoStatFunctional)
  FORWARD_METHOD_NON_CONST(getCorAniso, makeTensorNoStatFunctional)
  FORWARD_METHOD_NON_CONST(getCorAniso, makeParamNoStatFunctional)
  FORWARD_METHOD_NON_CONST(getCorAniso, makeRangeStationary)
  FORWARD_METHOD_NON_CONST(getCorAniso, makeScaleStationary)
  FORWARD_METHOD_NON_CONST(getCorAniso, makeAngleStationary)
  FORWARD_METHOD_NON_CONST(getCorAniso, makeTensorStationary)
  FORWARD_METHOD_NON_CONST(getCorAniso, makeParamStationary)

  FORWARD_METHOD_NON_CONST(getCorAniso, setMarkovCoeffs)
  FORWARD_METHOD_NON_CONST(getCorAniso, setMarkovCoeffsBySquaredPolynomials)

  FORWARD_METHOD_NON_CONST(getCorAniso, informDbInForAnisotropy)
  FORWARD_METHOD_NON_CONST(getCorAniso, informDbOutForAnisotropy)

  FORWARD_METHOD_NON_CONST(getCorAniso, informMeshByMeshForAnisotropy)
  FORWARD_METHOD_NON_CONST(getCorAniso, informMeshByApexForAnisotropy)

  FORWARD_METHOD_NON_CONST(getCorAniso, computeCorrec)

  FORWARD_METHOD(getCorAniso, getScaleIso, TEST)
  FORWARD_METHOD(getCorAniso, getScale, TEST)
  FORWARD_METHOD(getCorAniso, isValidForTurningBand, false)
  FORWARD_METHOD(getCorAniso, simulateTurningBand, TEST)
  FORWARD_METHOD(getCorAniso, isValidForSpectral, false)
  FORWARD_METHOD(getCorAniso, simulateSpectralOmega, MatrixDense())
  FORWARD_METHOD(getCorAniso, getRanges, VectorDouble())
  FORWARD_METHOD(getCorAniso, getScales, VectorDouble())
  FORWARD_METHOD(getCorAniso, getRangeIso, TEST)
  FORWARD_METHOD(getCorAniso, getRange, TEST)
  FORWARD_METHOD(getCorAniso, getAnisoCoeffs, VectorDouble())
  FORWARD_METHOD(getCorAniso, getType, TEST)
  FORWARD_METHOD(getCorAniso, getParam, TEST)
  FORWARD_METHOD(getCorAniso, getScadef, TEST)
  FORWARD_METHOD(getCorAniso, getParMax, TEST)
  FORWARD_METHOD(getCorAniso, getMaxNDim, ITEST)
  FORWARD_METHOD(getCorAniso, getMinOrder, TEST)
  FORWARD_METHOD(getCorAniso, hasInt1D, false)
  FORWARD_METHOD(getCorAniso, hasInt2D, false)
  FORWARD_METHOD(getCorAniso, hasRange, false)
  FORWARD_METHOD(getCorAniso, hasParam, false)
  FORWARD_METHOD(getCorAniso, getNGradParam, ITEST)
  FORWARD_METHOD(getCorAniso, hasCovDerivative, false)
  FORWARD_METHOD(getCorAniso, hasCovOnSphere, false)
  FORWARD_METHOD(getCorAniso, hasSpectrumOnSphere, false)
  FORWARD_METHOD(getCorAniso, hasMarkovCoeffs, false)
  FORWARD_METHOD(getCorAniso, hasSpectrumOnRn, false)
  FORWARD_METHOD(getCorAniso, normalizeOnSphere, false)
  FORWARD_METHOD(getCorAniso, getMarkovCoeffs, VectorDouble())
  FORWARD_METHOD(getCorAniso, getCorrec, false)
  FORWARD_METHOD(getCorAniso, getFullCorrec, false)

  FORWARD_METHOD(getCorAniso, isNoStatForParam, false)
  FORWARD_METHOD(getCorAniso, isNoStatForTensor, false)
  FORWARD_METHOD(getCorAniso, isNoStatForAnisotropy, false)
  FORWARD_METHOD(getCorAniso, isNoStatForRotation, false)

  FORWARD_METHOD(getCorAniso, getNAngles, ITEST)
  FORWARD_METHOD(getCorAniso, getNRanges, ITEST)
  FORWARD_METHOD(getCorAniso, getNScales, ITEST)

  FORWARD_METHOD(getCorAniso, getDetTensor, false)

  double getSlope(int ivar, int jvar) const;
  const Rotation& getAnisoRotation() const { return getCorAniso()->getAniso().getRotation(); }
  const CorAniso* getCorAniso() const { return dynamic_cast<const CorAniso*>(getCor()); }

  bool getFlagAniso() const { return !isIsotropic(); }
  bool getFlagRotation() const { return hasRotation(); }
  const VectorDouble& getAnisoAngles() const { return getCorAniso()->getAniso().getAngles(); }
  const MatrixSquare& getAnisoRotMat() const { return getCorAniso()->getAniso().getMatrixDirect(); }
  const MatrixSquare& getAnisoInvMat() const { return getCorAniso()->getAniso().getMatrixInverse(); }
  double getAnisoAngle(int idim) const { return getAnisoAngles()[idim]; }
  double getAnisoRotMatElement(int idim, int jdim) const { return getCorAniso()->getAniso().getMatrixDirect().getValue(idim, jdim); }
  double getAnisoCoeff(int idim) const { return getAnisoCoeffs()[idim]; }
  const CovContext& getContext() const { return _ctxt; }

  String getCovName() const { return getCorAniso()->getCovName(); }
  bool isIsotropic() const { return getCorAniso()->getAniso().isIsotropic(); }
  bool isAsymptotic() const { return getScadef() != 1.; }
  bool hasRotation() const { return getCorAniso()->getAniso().hasRotation(); }
  const Tensor& getAniso() const { return getCorAniso()->getAniso(); }
  void setAniso(const Tensor& aniso) { getCorAniso()->setAniso(aniso); }
  const ACovFunc* getCorFunc() const { return getCorAniso()->getCorFunc(); }

  VectorDouble evalCovOnSphereVec(const VectorDouble& alpha,
                                  int degree              = 50,
                                  bool flagScaleDistance  = false,
                                  const CovCalcMode* mode = nullptr) const;
  Array evalCovFFT(const VectorDouble& hmax, int N = 128, int ivar = 0, int jvar = 0) const;

  int getNDim() const { return _ctxt.getNDim(); }
  CorAniso* getCorAniso();
  CovAniso* createReduce(const VectorInt& validVars) const;

  bool _isOptimEnabled() const override { return _optimEnabled && !isNoStatForAnisotropy(); }

  std::vector<ParamInfo>& getScalesParam() { return getCorAniso()->getParamInfoScales(); }
  std::vector<ParamInfo>& getAnglesParam() { return getCorAniso()->getParamInfoAngles(); }

  double _getSillValue(int ivar, int jvar, const CovCalcMode* mode) const;

  virtual double _eval(const SpacePoint& p1,
                       const SpacePoint& p2,
                       int ivar                = 0,
                       int jvar                = 0,
                       const CovCalcMode* mode = nullptr) const override;

};

GSTLEARN_EXPORT double scale2range(const ECov& type, double scale, double param = 1.);
GSTLEARN_EXPORT double range2scale(const ECov& type, double range, double param = 1.);

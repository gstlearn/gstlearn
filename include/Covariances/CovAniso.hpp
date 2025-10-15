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

#include "Arrays/Array.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/Tensor.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/ACovFunc.hpp"
#include "Covariances/CorAniso.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/CovProportional.hpp"
#include "Enum/ECov.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
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
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ASpaceObject Interface

  /// ACov Interface
  double eval0(Id ivar                 = 0,
               Id jvar                 = 0,
               const CovCalcMode* mode = nullptr) const override;

  double evalCovOnSphere(double alpha,
                         Id degree               = 50,
                         bool flagScaleDistance  = true,
                         const CovCalcMode* mode = nullptr) const override;
  VectorDouble evalSpectrumOnSphere(Id n,
                                    bool flagNormDistance = false,
                                    bool flagCumul        = false) const override;
  double evalSpectrum(const VectorDouble& freq,
                      Id ivar = 0,
                      Id jvar = 0) const override;

  virtual double getIntegralRange(Id ndisc, double hmax) const;
  virtual String getFormula() const { return getCorAniso()->getFormula(); }
  virtual double getBallRadius() const { return TEST; }

  static CovAniso* createFromParam(const ECov& type,
                                   double range,
                                   double sill                  = 1.,
                                   double param                 = 1,
                                   const VectorDouble& ranges   = VectorDouble(),
                                   const MatrixSymmetric& sills = MatrixSymmetric(),
                                   const VectorDouble& angles   = VectorDouble(),
                                   const ASpaceSharedPtr& space = nullptr,
                                   bool flagRange               = true);
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

  FORWARD_METHOD_NON_CONST(getCorAnisoModify, setParam)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, computeMarkovCoeffs)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, setRangeIsotropic)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, setRange)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, setRanges)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, setScale)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, setScaleDim)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, setScales)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, setAnisoRotationMat)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, setAnisoRotation)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, setAnisoAngles)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, setAnisoAngle)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, setRotationAnglesAndRadius)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, setType)

  FORWARD_METHOD_NON_CONST(getCorAnisoModify, makeRangeNoStatDb)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, makeScaleNoStatDb)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, makeAngleNoStatDb)

  FORWARD_METHOD_NON_CONST(getCorAnisoModify, makeTensorNoStatDb)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, makeParamNoStatDb)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, makeRangeNoStatFunctional)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, makeScaleNoStatFunctional)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, makeAngleNoStatFunctional)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, makeTensorNoStatFunctional)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, makeParamNoStatFunctional)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, makeRangeStationary)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, makeScaleStationary)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, makeAngleStationary)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, makeTensorStationary)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, makeParamStationary)

  FORWARD_METHOD_NON_CONST(getCorAnisoModify, setMarkovCoeffs)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, setMarkovCoeffsBySquaredPolynomials)

  FORWARD_METHOD_NON_CONST(getCorAnisoModify, informDbInForAnisotropy)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, informDbOutForAnisotropy)

  FORWARD_METHOD_NON_CONST(getCorAnisoModify, informMeshByMeshForAnisotropy)
  FORWARD_METHOD_NON_CONST(getCorAnisoModify, informMeshByApexForAnisotropy)

  FORWARD_METHOD_NON_CONST(getCorAnisoModify, computeCorrec)

  FORWARD_METHOD(getCorAniso, getScaleIso, TEST)
  FORWARD_METHOD(getCorAniso, getScale, TEST)
  FORWARD_METHOD(getCorAniso, isValidForTurningBand, false)
  FORWARD_METHOD(getCorAniso, simulateTurningBand, TEST)
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
  FORWARD_METHOD(getCorAniso, getMinOrder, ITEST)
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

  bool isValidForSpectral() const override; //Do not use FORWARD_METHOD because "override" creates compilation issue on Windows
  double getSlope(Id ivar, Id jvar) const;
  const Rotation& getAnisoRotation() const { return getCorAniso()->getAniso().getRotation(); }
  bool getFlagAniso() const { return !isIsotropic(); }
  bool getFlagRotation() const { return hasRotation(); }
  const VectorDouble& getAnisoAngles() const { return getCorAniso()->getAniso().getAngles(); }
  const MatrixSquare& getAnisoRotMat() const { return getCorAniso()->getAniso().getMatrixDirect(); }
  const MatrixSquare& getAnisoInvMat() const { return getCorAniso()->getAniso().getMatrixInverse(); }
  double getAnisoAngle(Id idim) const { return getAnisoAngles()[idim]; }
  double getAnisoRotMatElement(Id idim, Id jdim) const { return getCorAniso()->getAniso().getMatrixDirect().getValue(idim, jdim); }
  double getAnisoCoeff(Id idim) const { return getAnisoCoeffs()[idim]; }
  const CovContext& getContext() const { return _ctxt; }

  String getCovName() const { return getCorAniso()->getCovName(); }
  bool isIsotropic() const { return getCorAniso()->getAniso().isIsotropic(); }
  bool isAsymptotic() const { return getScadef() != 1.; }
  bool hasRotation() const { return getCorAniso()->getAniso().hasRotation(); }
  const Tensor& getAniso() const { return getCorAniso()->getAniso(); }
  void setAniso(const Tensor& aniso) { getCorAnisoModify()->setAniso(aniso); }
  const ACovFunc* getCorFunc() const { return getCorAniso()->getCorFunc(); }

  VectorDouble evalCovOnSphereVec(const VectorDouble& alpha,
                                  Id degree               = 50,
                                  bool flagScaleDistance  = false,
                                  const CovCalcMode* mode = nullptr) const;
  Array evalCovFFT(const VectorDouble& hmax, Id N = 128, Id ivar = 0, Id jvar = 0) const;

  Id getNDim() const { return static_cast<Id>(_ctxt.getNDim()); }
  const CorAniso* getCorAniso() const;
  CorAniso* getCorAnisoModify();
  CovAniso* createReduce(const VectorInt& validVars) const;

  bool _isOptimEnabled() const override { return _optimEnabled && !isNoStatForAnisotropy(); }

  std::vector<ParamInfo>& getScalesParam() { return getCorAnisoModify()->getParamInfoScales(); }
  std::vector<ParamInfo>& getAnglesParam() { return getCorAnisoModify()->getParamInfoAngles(); }

  double _getSillValue(Id ivar, Id jvar, const CovCalcMode* mode) const;

  double _eval(const SpacePoint& p1,
               const SpacePoint& p2,
               Id ivar                 = 0,
               Id jvar                 = 0,
               const CovCalcMode* mode = nullptr) const override;
};

GSTLEARN_EXPORT double scale2range(const ECov& type, double scale, double param = 1.);
GSTLEARN_EXPORT double range2scale(const ECov& type, double range, double param = 1.);
} // namespace gstlrn

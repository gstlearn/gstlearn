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
#include "Basic/ParamInfo.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/TabNoStat.hpp"
#include "Covariances/TabNoStatCovAniso.hpp"
#include "Enum/EConsElem.hpp"
#include "Model/CovInternal.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Enum/ECov.hpp"

#include "Arrays/Array.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/Tensor.hpp"
#include "Covariances/AKernel.hpp"
#include "Covariances/CovContext.hpp"
#include "Space/SpacePoint.hpp"
#include <array>
#include <vector>

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
 */
class GSTLEARN_EXPORT CorAniso: public ACov
{
public:
  CorAniso(const ECov& type, const CovContext& ctxt);
  CorAniso(const String& symbol, const CovContext& ctxt);
  CorAniso(const ECov& type,
           double range,
           double param,
           const CovContext& ctxt,
           bool flagRange = true);
  CorAniso(const CorAniso& r);
  CorAniso& operator=(const CorAniso& r);
  virtual ~CorAniso();

  /// ICloneable Interface
  IMPLEMENT_CLONING(CorAniso)

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ASpaceObject Interface
  bool isConsistent(const ASpace* space) const override;

  /// ACov Interface
  Id getNVar() const override { return 1; }

  /// ACov Interface

  double evalCor(const SpacePoint& p1,
                 const SpacePoint& p2,
                 const CovCalcMode* mode = nullptr,
                 Id ivar                 = 0,
                 Id jvar                 = 0) const; // let ivar and jvar for the future where the
                                     // correlation will be different for multivariate

  double evalCovOnSphere(double alpha,
                         Id degree               = 50,
                         bool flagScaleDistance  = true,
                         const CovCalcMode* mode = nullptr) const override;
  VectorDouble evalSpectrumOnSphere(Id n,
                                    bool flagNormDistance = false,
                                    bool flagCumul        = false) const override;
  double evalSpectrumRatio(const VectorDouble& freq, Id ivar, Id jvar, const ACov* cov0) const override;

    double evalSpectrum(const VectorDouble& freq,
                        Id ivar = 0,
                        Id jvar = 0) const override;

  virtual double getIntegralRange(Id ndisc, double hmax) const;
  virtual String getFormula() const { return _corfunc->getFormula(); }
  virtual double getBallRadius() const { return TEST; }

  static bool isOptimizationInitialized(const std::vector<SpacePoint>& p1As,
                                        const Db* db = nullptr);

  void _optimizationPreProcess(Id mode, const std::vector<SpacePoint>& ps) const override;
  void _optimizationSetTarget(SpacePoint& p) const override;
  void _optimizationPostProcess() const override;

  bool isValidForTurningBand() const;
  double simulateTurningBand(double t0, TurningBandOperate& operTB) const;
  bool isValidForSpectral() const override;
  MatrixDense simulateSpectralOmega(Id nb) const override;

  static CorAniso* createIsotropic(const CovContext& ctxt,
                                   const ECov& type,
                                   double range,
                                   double param   = 1.,
                                   bool flagRange = true);
  static CorAniso* createAnisotropic(const CovContext& ctxt,
                                     const ECov& type,
                                     const VectorDouble& ranges,
                                     double param               = 1.,
                                     const VectorDouble& angles = VectorDouble(),
                                     bool flagRange             = true);
  static CorAniso* createIsotropicMulti(const CovContext& ctxt,
                                        const ECov& type,
                                        double range,
                                        double param   = 1.,
                                        bool flagRange = true);
  static CorAniso* createAnisotropicMulti(const CovContext& ctxt,
                                          const ECov& type,
                                          const VectorDouble& ranges,
                                          double param               = 1.,
                                          const VectorDouble& angles = VectorDouble(),
                                          bool flagRange             = true);

  void setParam(double param, Id ipar = 0);
  void setNoStatFactor(double noStatFactor) { _noStatFactor = noStatFactor; }

  /// Practical range
  void setRangeIsotropic(double range);
  void setRange(Id idim, double range);
  void setRanges(const VectorDouble& ranges);

  void setScale(double scale); /// Make the covariance isotropic
  void setScaleDim(Id idim, double scale);
  void setScales(const VectorDouble& scales);

  void setAnisoRotationMat(const Rotation& rot);
  void setAnisoRotation(const VectorDouble& rot);
  void setAnisoAngles(const VectorDouble& angles);
  void setAnisoAngle(Id idim, double angle);

  void setRotationAnglesAndRadius(const VectorDouble& angles = VectorDouble(),
                                  const VectorDouble& ranges = VectorDouble(),
                                  const VectorDouble& scales = VectorDouble()) const;

  VectorDouble getRanges() const;
  double getRange(Id idim) const { return getRanges()[idim]; }
  double getScale(Id idim) const { return getScales()[idim]; }
  const Rotation& getAnisoRotation() const { return _aniso.getRotation(); }
  const VectorDouble& getScales() const { return _aniso.getRadius(); }

  void setType(const ECov& type);
  double getRangeIso() const;
  double getScaleIso() const;
  bool getFlagAniso() const { return !isIsotropic(); }
  bool getFlagRotation() const { return hasRotation(); }
  VectorDouble getAnisoAngles() const { return _aniso.getAngles(); }
  const MatrixSquare& getAnisoRotMat() const { return _aniso.getMatrixDirect(); }
  const MatrixSquare& getAnisoInvMat() const { return _aniso.getMatrixInverse(); }
  VectorDouble getAnisoCoeffs() const;
  double getAnisoAngle(Id idim) const { return getAnisoAngles()[idim]; }
  double getAnisoRotMat(Id idim, Id jdim) const { return _aniso.getMatrixDirect().getValue(idim, jdim); }
  double getAnisoCoeff(Id idim) const { return getAnisoCoeffs()[idim]; }
  const ECov& getType() const { return _corfunc->getType(); }
  double getParam() const;
  double getScadef() const { return _corfunc->getScadef(); }
  double getParMax() const { return _corfunc->getParMax(); }
  Id getMaxNDim() const { return static_cast<Id>(_corfunc->getMaxNDim()); }
  Id getMinOrder() const { return _corfunc->getMinOrder(); }
  bool hasInt1D() const { return _corfunc->hasInt1D(); }
  bool hasInt2D() const { return _corfunc->hasInt2D(); }
  Id hasRange() const { return _corfunc->hasRange(); }
  Id hasParam() const { return _corfunc->hasParam(); }
  String getCovName() const { return _corfunc->getCovName(); }
  bool isIsotropic() const { return _aniso.isIsotropic(); }
  bool isAsymptotic() const { return getScadef() != 1.; }
  bool hasRotation() const { return _aniso.hasRotation(); }
  const Tensor& getAniso() const { return _aniso; }
  void setAniso(const Tensor& aniso) { _aniso = aniso; }
  const AKernel* getCorFunc() const { return _corfunc; }
  Id getNGradParam() const;
  bool hasCovDerivative() const { return _corfunc->hasCovDerivative(); }
  bool hasCovOnSphere() const { return _corfunc->hasCovOnSphere(); }
  bool hasSpectrumOnSphere() const { return _corfunc->hasSpectrumOnSphere(); }
  bool hasMarkovCoeffs() const { return _corfunc->hasMarkovCoeffs(); }
  bool hasSpectrumOnRn() const { return _corfunc->hasSpectrumOnRn(); }
  double normalizeOnSphere(Id n = 50) const;
  //////////////////////// New NoStat methods //////////////////////////

  void makeRangeNoStatDb(const String& namecol, Id idim = 0, const Db* db = nullptr);
  void makeScaleNoStatDb(const String& namecol, Id idim = 0, const Db* db = nullptr);
  void makeAngleNoStatDb(const String& namecol, Id idim = 0, const Db* db = nullptr);
  void makeTensorNoStatDb(const String& namecol, Id idim = 0, Id jdim = 0, const Db* db = nullptr);
  void makeParamNoStatDb(const String& namecol, const Db* db = nullptr);

  void makeRangeNoStatFunctional(const AFunctional* func, Id idim = 0);
  void makeScaleNoStatFunctional(const AFunctional* func, Id idim = 0);
  void makeAngleNoStatFunctional(const AFunctional* func, Id idim = 0);
  void makeTensorNoStatFunctional(const AFunctional* func, Id idim = 0, Id jdim = 0);
  void makeParamNoStatFunctional(const AFunctional* func);

  void makeRangeStationary(Id idim = 0) const;
  void makeScaleStationary(Id idim = 0) const;
  void makeAngleStationary(Id idim = 0) const;
  void makeTensorStationary(Id idim, Id jdim);
  void makeParamStationary();

  TabNoStatCovAniso* getTabNoStatCovAniso() const { return static_cast<TabNoStatCovAniso*>(_tabNoStat); }

  Id getNAngles() const { return getTabNoStatCovAniso()->getNAngles(); }
  Id getNRanges() const { return getTabNoStatCovAniso()->getNRanges(); }
  Id getNScales() const { return getTabNoStatCovAniso()->getNScales(); }
  bool isNoStatForParam() const { return getTabNoStatCovAniso()->isParam(); }
  bool isNoStatForTensor() const { return getTabNoStatCovAniso()->isDefinedForTensor(); }
  bool isNoStatForAnisotropy() const { return getTabNoStatCovAniso()->isDefinedForAnisotropy(); }
  bool isNoStatForRotation() const { return getTabNoStatCovAniso()->isDefinedForRotation(); }

  VectorDouble evalCovOnSphereVec(const VectorDouble& alpha,
                                  Id degree               = 50,
                                  bool flagScaleDistance  = false,
                                  const CovCalcMode* mode = nullptr) const;
  Array evalCovFFT(const VectorDouble& hmax, Id N = 128, Id ivar = 0, Id jvar = 0) const;
  VectorDouble getMarkovCoeffs() const;
  void setMarkovCoeffs(const VectorDouble& coeffs);

  void setMarkovCoeffsBySquaredPolynomials(VectorDouble coeffs1, VectorDouble coeffs2, double eps = 0);
  void computeMarkovCoeffs();
  double getCorrec() const;
  double getFullCorrec() const;
  void nostatUpdate(CovInternal* covint) const;

  void informMeshByMeshForAnisotropy(const AMesh* amesh) const;
  void informMeshByApexForAnisotropy(const AMesh* amesh) const;
  void informDbInForAnisotropy(const Db* dbin) const;
  void informDbOutForAnisotropy(const Db* dbout) const;

  /// Tell if the use of Optimization is enabled or not

  void updateCovByPoints(Id icas1, Id iech1, Id icas2, Id iech2) const override;
  void updateCovByMesh(Id imesh, bool aniso = true) const override;
  double getValue(const EConsElem& econs, Id iv1, Id iv2) const override;
  void computeCorrec();
  double evalCorFromH(double h, const CovCalcMode* mode) const;
  double getDetTensor() const;

  void optimizationTransformSP(const SpacePoint& ptin, SpacePoint& ptout) const;
  void optimizationTransformSPNew(const SpacePoint& ptin, SpacePoint& ptout) const;
  String toStringParams(const AStringFormat* strfmt = nullptr) const;
  String toStringNoStat(const AStringFormat* strfmt = nullptr, Id i = 0) const;
  void appendParams(ListParams& listParams,
                    std::vector<covmaptype>* gradFuncs = nullptr) override;
  double evalDerivativeBasis(const SpacePoint& p1,
                             const SpacePoint& p2,
                             Id ivar,
                             Id jvar,
                             const CovCalcMode* mode) const;
  void updateCov() override;
  void initParams(const MatrixSymmetric& vars,
                  double href = 1.) override;
  ParamInfo& getParamInfoScale(Id idim) { return _scales[idim]; }
  ParamInfo& getParamInfoAngle(Id idim) { return _angles[idim]; }
  std::vector<ParamInfo>& getParamInfoScales() { return _scales; }
  std::vector<ParamInfo>& getParamInfoAngles() { return _angles; }
  bool getOptimLockIso2d() const { return _optimLockIso2d; }
  void setOptimLockIso2d(bool status) { _optimLockIso2d = status; }
  bool getOptimNoAniso() const { return _optimNoAniso; };
  void setOptimNoAniso(bool status) { _optimNoAniso = status; }
  std::vector<MatrixSquare>& getDRot() const { return _dRot; }

#ifndef SWIG
  Id addEvalCovVecRHSInPlace(vect vect,
                             const VectorInt& index1,
                             Id iech2,
                             const KrigOpt& krigopt,
                             SpacePoint& pin,
                             SpacePoint& pout,
                             VectorDouble& tabwork,
                             double lambda                 = 1.,
                             const ECalcMember& calcMember = ECalcMember::RHS) const override;
#endif

protected:
  /// Update internal parameters consistency with the context
  void _initFromContext() override;

private:
  void _initParamInfo();
  void _handleConstraints();

  bool _isNoStat() const override;
  void _setContext(const CovContext& ctxt) override;
  TabNoStat* _createNoStatTab() override;
  void _copyCovContext(const CovContext& ctxt) override;

  bool _isOptimEnabled() const override
  {
    return _optimEnabled && !isNoStatForAnisotropy();
  }

  void _manage(const Db* db1, const Db* db2) const override;

  bool _checkTensor() const;
  bool _checkRotation() const;
  bool _checkParam() const;

  bool _isVariableValid(Id ivar) const;
  void _updateFromContext() override;

  virtual double _eval(const SpacePoint& p1,
                       const SpacePoint& p2,
                       Id ivar                 = 0,
                       Id jvar                 = 0,
                       const CovCalcMode* mode = nullptr) const override;

private:
  AKernel* _corfunc; /// Basic correlation function
  std::vector<ParamInfo> _scales;
  std::vector<ParamInfo> _angles;
  mutable Tensor _aniso;        /// Anisotropy parameters
  mutable double _noStatFactor; /// Correcting factor for non-stationarity

  bool _optimNoAniso;   // All ranges should be equal
  bool _optimLockIso2d; // Range U and V should be equal
  mutable std::vector<MatrixSquare> _dRot;
  const std::array<EConsElem, 4> _listaniso = {EConsElem::RANGE,
                                               EConsElem::SCALE,
                                               EConsElem::TENSOR,
                                               EConsElem::ANGLE};
  // These temporary information is used to speed up processing (optimization functions)
  // They are in a protected section as they may be modified by class hierarchy
};

} // namespace gstlrn

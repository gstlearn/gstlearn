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

#include <vector>

class Rotation;
class MatrixSquareGeneral;
class MatrixRectangular;

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
  virtual void eval0MatInPlace(MatrixSquareGeneral &mat,
                               const CovCalcMode *mode = nullptr) const override;
  virtual void evalMatInPlace(const SpacePoint &p1,
                              const SpacePoint &p2,
                              MatrixSquareGeneral &mat,
                              const CovCalcMode *mode = nullptr) const override;
  virtual double evalCovOnSphere(double alpha,
                                 int degree = 0,
                                 bool flagScaleDistance = true,
                                 const CovCalcMode* mode = nullptr) const override;
  virtual VectorDouble evalSpectrumOnSphere(int n,
                                            bool flagNormDistance = false,
                                            bool flagCumul = false) const override;
  virtual double evalSpectrum(const VectorDouble &freq,
                              int ivar = 0,
                              int jvar = 0) const override;

  virtual double getIntegralRange(int ndisc, double hmax) const;
  virtual String getFormula() const { return _cova->getFormula(); }
  virtual double getBallRadius() const { return TEST; }

  bool isOptimizationInitialized(const Db* db = nullptr) const;
  void optimizationPreProcess(const Db* db) const;
  void optimizationPostProcess() const;
  void optimizationSetTarget(const SpacePoint& pt) const;
  void optimizationSetTarget(int iech) const;

  void evalOptimInPlace(MatrixRectangular& res,
                        const VectorInt& ivars,
                        const VectorVectorInt& index,
                        int ivar2 = 0,
                        int icol = 0,
                        const CovCalcMode *mode = nullptr,
                        bool flagSym = false) const;
  void evalMatOptimInPlace(int icas1,
                           int iech1,
                           int icas2,
                           int iech2,
                           MatrixSquareGeneral &mat,
                           const CovCalcMode *mode = nullptr) const override;
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
  void setNoStatFactor(double noStatFactor) { _noStatFactor = noStatFactor; }

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
  const Rotation& getAnisoRotation() const { return _aniso.getRotation(); }
  const VectorDouble& getScales() const { return _aniso.getRadius(); }

  void   setType(const ECov& type);
  double getRange() const;
  double getScale() const;
  bool   getFlagAniso() const { return !isIsotropic(); }
  bool   getFlagRotation() const { return hasRotation(); }
  double getRange(int idim) const { return getRanges()[idim]; }
  double getScale(int idim) const { return getScales()[idim]; }
  const VectorDouble getAnisoAngles() const { return _aniso.getAngles(); }
  const MatrixSquareGeneral& getAnisoRotMat() const { return _aniso.getMatrixDirect(); }
  const MatrixSquareGeneral& getAnisoInvMat() const { return _aniso.getMatrixInverse(); }
  const VectorDouble getAnisoCoeffs() const;
  double getAnisoAngles(int idim) const { return getAnisoAngles()[idim]; }
  double getAnisoRotMat(int idim, int jdim) const { return _aniso.getMatrixDirect().getValue(idim,jdim); }
  double getAnisoCoeffs(int idim) const { return getAnisoCoeffs()[idim]; }
  const CovContext& getContext() const { return _ctxt; }
  const ECov& getType() const { return _cova->getType(); }
  double getParam() const;
  double getScadef() const { return _cova->getScadef(); }
  double getParMax() const { return _cova->getParMax(); }
  int    getMaxNDim() const { return _cova->getMaxNDim(); }
  int    getMinOrder() const { return _cova->getMinOrder(); }
  bool   hasInt1D() const { return _cova->hasInt1D(); }
  bool   hasInt2D() const { return _cova->hasInt2D(); }
  int    hasRange() const { return _cova->hasRange(); }
  int    hasParam() const  { return _cova->hasParam(); }
  String getCovName() const { return _cova->getCovName(); }
  bool   isIsotropic() const { return _aniso.isIsotropic(); }
  bool   isAsymptotic() const { return getScadef() != 1.; }
  bool   hasRotation() const { return _aniso.hasRotation(); }
  const Tensor& getAniso() const { return _aniso; }
  void   setAniso(const Tensor& aniso) { _aniso = aniso; }
  const ACovFunc* getCova() const { return _cova; }
  int    getGradParamNumber() const;
  bool   hasCovDerivative() const { return _cova->hasCovDerivative(); }
  bool   hasCovOnSphere() const { return _cova->hasCovOnSphere(); }
  bool   hasSpectrumOnSphere() const { return _cova->hasSpectrumOnSphere(); }
  bool   hasMarkovCoeffs() const { return _cova->hasMarkovCoeffs(); }
  bool   hasSpectrumOnRn() const { return _cova->hasSpectrumOnRn(); }

  VectorDouble evalCovOnSphereVec(const VectorDouble &alpha,
                                  int degree = 0,
                                  bool flagScaleDistance = false,
                                  const CovCalcMode* mode = nullptr) const;
  Array evalCovFFT(const VectorDouble& ext, int N = 128, int ivar = 0, int jvar = 0) const;
  VectorDouble getMarkovCoeffs() const;
  void setMarkovCoeffs(VectorDouble coeffs);

  void setMarkovCoeffsBySquaredPolynomials(VectorDouble coeffs1, VectorDouble coeffs2, double eps = 0);
  void computeMarkovCoeffs();
  double getCorrec() const;
  double getFullCorrec() const;
  int getDimensionNumber() const { return _ctxt.getNDim(); }

  CovAniso* createReduce(const VectorInt &validVars) const;

protected:
  /// Update internal parameters consistency with the context
  virtual void _updateFromContext();
  virtual void _initFromContext();

private:
  bool   _isVariableValid(int ivar) const;
  void   _computeCorrec();
  double _getDetTensor() const;
  void   _optimizationTransformSP(const SpacePoint& ptin, SpacePoint& ptout) const;
  double _evalCovFromH(double h, const CovCalcMode *mode) const;

private:
  CovContext _ctxt;            /// Context (space, number of variables, ...) // TODO : Really store a copy ?
  ACovFunc *_cova;             /// Covariance basic function
  MatrixSquareSymmetric _sill; /// Sill matrix (nvar x nvar)
  Tensor _aniso;               /// Anisotropy parameters
  double _noStatFactor;        /// Correcting factor for non-stationarity
};

GSTLEARN_EXPORT double scale2range(const ECov& type, double scale, double param = 1.);
GSTLEARN_EXPORT double range2scale(const ECov& type, double range, double param = 1.);


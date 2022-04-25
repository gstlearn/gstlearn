/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/Vector.hpp"
#include "Basic/IClonable.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Basic/Tensor.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/ECov.hpp"
#include "Covariances/ACovFunc.hpp"
#include "Covariances/CovContext.hpp"

class Rotation;

class GSTLEARN_EXPORT CovAniso: public ACov, public IClonable
{
public:
  CovAniso(const ECov& type, const CovContext& ctxt);
  CovAniso(const String& symbol, const CovContext& ctxt);
  CovAniso(const ECov& type,
           double range,
           double param,
           double sill,
           const CovContext& ctxt);
  CovAniso(const CovAniso& r);
  CovAniso& operator=(const CovAniso& r);
  virtual ~CovAniso();

  ///////////////////////////////////////////////////
  /// IClonable Interface
  virtual IClonable* clone() const override { return new CovAniso(*this); };
  ///////////////////////////////////////////////////

  ///////////////////////////////////////////////////
  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;
  ///////////////////////////////////////////////////

  ///////////////////////////////////////////////////
  /// ASpaceObject Interface
  virtual bool isConsistent(const ASpace* space) const override;
  ///////////////////////////////////////////////////

  ///////////////////////////////////////////////////
  /// ACov Interface
  virtual int getNVariables() const override { return _ctxt.getNVar(); }
  ///////////////////////////////////////////////////

  ///////////////////////////////////////////////////
  /// ACov Interface
  /**
   * Evaluate the covariance for a pair of variables and a zero distance
   * @param ivar Rank of the first variable
   * @param jvar Rank of the second variable
   * @param mode Reference to the CovCalcMode embedded class
   * @return The covariance value at the origin
   */
  virtual double eval0(int ivar,
                       int jvar,
                       const CovCalcMode& mode = CovCalcMode()) const override;

  /**
   * Evaluate covariance between two points (p1, p2) for two variables (ivar, jvar)
   * @param ivar Rank of the first variable
   * @param jvar Rank of the second variable
   * @param p1   Rank of the first point
   * @param p2   Rank of the second point
   * @param mode Reference to the CovCalcMode embedded class
   * @return The covariance value
   */
  virtual double eval(int ivar,
                      int jvar,
                      const SpacePoint& p1,
                      const SpacePoint& p2,
                      const CovCalcMode& mode = CovCalcMode()) const override;
  ///////////////////////////////////////////////////

  virtual double getIntegralRange(int ndisc, double hmax) const;
  virtual String getFormula() const { return _cova->getFormula(); }
  virtual double getBallRadius() const { return TEST; }

  static CovAniso* createIsotropic(const CovContext& ctxt,
                                   const ECov& type,
                                   double range,
                                   double sill = 1.,
                                   double param = 1.);
  static CovAniso* createAnisotropic(const CovContext& ctxt,
                                     const ECov& type,
                                     const VectorDouble& ranges,
                                     double sill = 1.,
                                     double param = 1.,
                                     const VectorDouble& angles = VectorDouble());
  static CovAniso* createIsotropicMulti(const CovContext& ctxt,
                                        const ECov& type,
                                        double range,
                                        const MatrixSquareGeneral& sills,
                                        double param = 1.);
  static CovAniso* createAnisotropicMulti(const CovContext& ctxt,
                                          const ECov& type,
                                          const VectorDouble& ranges,
                                          const MatrixSquareGeneral& sills,
                                          double param = 1.,
                                          const VectorDouble& angles = VectorDouble());

  void setContext(const CovContext& ctxt);
  void setParam(double param);
  void copyCovContext(const CovContext& ctxt);

  void setSill(double sill); /// Only valid when there is only one variable (in the context)
  void setSill(const MatrixSquareGeneral& sill);
  void setSill(const VectorDouble& sill);
  void setSill(int ivar, int jvar, double sill);
  void initSill(double value = 0.);

  /// Practical range
  void setRange(double range); /// Make the covariance isotropic
  void setRange(int idim, double range);
  void setRanges(const VectorDouble& range);

  void setScale(double scale); /// Make the covariance isotropic
  void setScale(int idim, double scale);
  void setScales(const VectorDouble& scale);

  void setAnisoRotation(const Rotation& rot);
  void setAnisoRotation(const VectorDouble& rot);
  void setAnisoAngles(const VectorDouble& angles);
  void setAnisoAngle(int idim, double angle);

  const MatrixSquareSymmetric& getSill() const { return _sill; }
  double getSill(int ivar, int jvar) const;
  double getSlope(int ivar, int jvar) const;
  VectorDouble getRanges() const;
  const Rotation& getAnisoRotation() const { return _aniso.getRotation(); }
  const VectorDouble& getScales() const { return _aniso.getRadius(); }

  // Look at C2R before touching that
  void   setType(const ECov& type);
  double getRange() const;
  double getTheoretical() const;
  bool   getFlagAniso() const { return !isIsotropic(); }
  bool   getFlagRotation() const { return hasRotation(); }
  double getRange(int idim) const { return getRanges()[idim]; }
  double getScale(int idim) const { return getScales()[idim]; }
  const VectorDouble getAnisoAngles() const { return _aniso.getAngles(); }
  const MatrixSquareGeneral& getAnisoRotMat() const { return _aniso.getMatrixDirect(); }
  const VectorDouble getAnisoRotMatVec() const { return getAnisoRotMat().getValues(); }
  const MatrixSquareGeneral& getAnisoInvMat() const { return _aniso.getMatrixInverse(); }
  const VectorDouble getAnisoInvMatVec() const { return getAnisoInvMat().getValues(); }
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

  static double scale2range(const ECov& type, double scale, double param = 1.);
  static double range2scale(const ECov& type, double range, double param = 1.);

protected:
  /// Update internal parameters consistency with the context
  virtual void _updateFromContext();

private:
  bool   _isVariableValid(int ivar) const;

private:
  CovContext      _ctxt;   /// Context (space, irfDegree, field, ...) // TODO : Really store a copy ?
  ACovFunc*       _cova;   /// Covariance basic function
  MatrixSquareSymmetric    _sill;   /// Sill matrix (nvar x nvar)
  Tensor          _aniso;  /// Anisotropy parameters
};

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

#include "Covariances/CovLMCAnamorphosis.hpp"
#include "Covariances/CovLMCTapering.hpp"
#include "Model/ModelCovList.hpp"
#include "gstlearn_export.hpp"

#include "geoslib_define.h"

#include "Enum/ECalcMember.hpp"
#include "Enum/ECov.hpp"
#include "Enum/EModelProperty.hpp"

#include "Covariances/ACov.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/CovAnisoList.hpp"
#include "Covariances/CovLMGradient.hpp"

#include "Drifts/DriftList.hpp"

#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"
#include "Model/Constraints.hpp"
#include "Covariances/CovAniso.hpp"

#include "Anamorphosis/AAnam.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/ICloneable.hpp"

class Model;
class Db;
class CovLMCTapering;
class CovLMCAnamorphosis;
class CovLMGradient;
class CovInternal;
class MatrixSquareSymmetric;
class CovCalcMode;
class Vario;
class ADrift;
class AnamContinuous;
class AnamHermite;

typedef std::vector<ECov> VectorECov;

/**
 * \brief
 * Class containing the Model Information describing the formal Spatial (or Temporal) Characteristics
 * of the (set of) random variable(s) under study.
 *
 * The Model is essentially a container with two main contents:
 * - the **covariance** part: see ACov.hpp for more information
 * - the **drift** part: see DriftList.hpp for more information
 *
 * The additional member **CovContext** only serves in carrying the following information:
 * - the number of variables: if more than 1, the Model becomes multivariate
 * - the field extension: this information is needed to get a *stationary* version to any covariance
 * - the experimental mean vector and the variance-covariance matrix (used to calibrate the Model)
 */
class GSTLEARN_EXPORT Model : public AStringable, public ASerializable, public ModelCovList
{
public:
  Model(const CovContext& ctxt = CovContext());
  Model(int nvar, int ndim = 2);
  Model(const Model &m);
  Model& operator= (const Model &m);
  virtual ~Model();

public:
  /// ICloneable interface
  IMPLEMENT_CLONING(Model)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

public:
  const CovAnisoList* castInCovAnisoListConst(int icov = -1) const;
  const CovLMCTapering* castInCovLMCTaperingConst() const;
  const CovLMGradient* castInCovLMGradientConst() const;
  const CovLMCAnamorphosis* castInCovLMCAnamorphosisConst() const;

public:
  CovAnisoList* _castInCovAnisoList(int icov = -1);
  CovLMCTapering* _castInCovLMCTapering();
  CovLMGradient* _castInCovLMGradient();
  CovLMCAnamorphosis* _castInCovLMCAnamorphosis();

public:
  int resetFromDb(const Db* db);
  static Model* create(const CovContext& ctxt = CovContext());
  static Model* createFromEnvironment(int nvar, int ndim = 2);
  static Model* createNugget(int nvar, int ndim = 2, double sill = 1.);
  static Model*
  createFromParam(const ECov& type                   = ECov::fromKey("NUGGET"),
                  double range                       = 1.,
                  double sill                        = 1.,
                  double param                       = 1.,
                  const VectorDouble& ranges         = VectorDouble(),
                  const MatrixSquareSymmetric& sills = MatrixSquareSymmetric(),
                  const VectorDouble& angles         = VectorDouble(),
                  const ASpaceSharedPtr& space       = ASpaceSharedPtr(),
                  bool flagRange                     = true);
  static Model* createFromParamOldStyle(const ECov& type             = ECov::fromKey("NUGGET"),
                                        double range                 = 1.,
                                        double sill                  = 1.,
                                        double param                 = 1.,
                                        const VectorDouble& ranges   = VectorDouble(),
                                        const VectorDouble& sills    = VectorDouble(),
                                        const VectorDouble& angles   = VectorDouble(),
                                        const ASpaceSharedPtr& space = ASpaceSharedPtr(),
                                        bool flagRange               = true);
  static Model* createFromDb(const Db* db);
  static Model* createFromNF(const String& neutralFilename,
                             bool verbose = true);
  static Model* createFromVario(Vario* vario,
                                const VectorECov& types        = ECov::fromKeys({"SPHERICAL"}),
                                const Constraints& constraints = Constraints(),
                                const Option_VarioFit& optvar  = Option_VarioFit(),
                                const Option_AutoFit& mauto    = Option_AutoFit(),
                                bool verbose                   = false);
  static Model* createFillRandom(int ndim,
                                 int nvar,
                                 const std::vector<ECov>& types = ECov::fromKeys({"SPHERICAL"}),
                                 double hmax                    = 1,
                                 int order                      = -1,test_gneiting
                                 
                                 int nfex                       = 0,
                                 int seed                       = 13242);
  void setCovAnisoList(const CovAnisoList* covalist);
  void addCov(const CovAniso* cov);
  void addCovFromParam(const ECov& type,
                       double range                       = EPSILON6,
                       double sill                        = 1.,
                       double param                       = 1.,
                       const VectorDouble& ranges         = VectorDouble(),
                       const MatrixSquareSymmetric& sills = MatrixSquareSymmetric(),
                       const VectorDouble& angles         = VectorDouble(),
                       bool flagRange                     = true);
  void addCovFromParamOldStyle(const ECov& type,
                               double range               = EPSILON6,
                               double sill                = 1.,
                               double param               = 1.,
                               const VectorDouble& ranges = VectorDouble(),
                               const VectorDouble& sills  = VectorDouble(),
                               const VectorDouble& angles = VectorDouble(),
                               bool flagRange             = true);

  FORWARD_METHOD(castInCovAnisoListConst, getActiveFactor,ITEST)
  FORWARD_METHOD(castInCovAnisoListConst, getCova)
  FORWARD_METHOD(castInCovAnisoListConst, getNCov,ITEST)
  FORWARD_METHOD(castInCovAnisoListConst, getCovType, ECov::UNKNOWN)
  FORWARD_METHOD(castInCovAnisoListConst, getRange, TEST)
  FORWARD_METHOD(castInCovAnisoListConst, getRanges)
  FORWARD_METHOD(castInCovAnisoListConst, getAngles)
  FORWARD_METHOD(castInCovAnisoListConst, getAnam)
  FORWARD_METHOD(castInCovAnisoListConst, getParam,TEST)
  FORWARD_METHOD(castInCovAnisoListConst, getCovName)
  FORWARD_METHOD(castInCovAnisoListConst, extractCova)
  FORWARD_METHOD(castInCovAnisoListConst, getNGradParam,ITEST)
  FORWARD_METHOD(castInCovAnisoListConst, getMaximumDistance, TEST)
  FORWARD_METHOD(castInCovAnisoListConst, getCovMinIRFOrder, ITEST)
  FORWARD_METHOD(castInCovAnisoListConst, getAnamNClass, ITEST)
  FORWARD_METHOD(castInCovAnisoListConst, hasAnam  , false)
  FORWARD_METHOD(castInCovAnisoListConst, hasNugget, false)
  FORWARD_METHOD(castInCovAnisoListConst, getRankNugget, -1)
  FORWARD_METHOD(castInCovAnisoListConst, getBallRadius, TEST)
  FORWARD_METHOD(castInCovAnisoListConst, hasExternalCov)
  FORWARD_METHOD(castInCovAnisoListConst, isChangeSupportDefined, false)
  FORWARD_METHOD(castInCovAnisoListConst, getAnamHermite) 
  FORWARD_METHOD(castInCovAnisoListConst, evalCovMat)
  FORWARD_METHOD(castInCovAnisoListConst, getCovMode, EModelProperty::NONE)

  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, setActiveFactor)
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, getCova)
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, setSill)
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, setSills)
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, setRangeIsotropic)
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, setMarkovCoeffs)
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, setCovFiltered)
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, setOptimEnabled)
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, normalize)

  FORWARD_METHOD_NON_CONST(_castInCovLMCTapering, setTapeRange)
  FORWARD_METHOD(castInCovLMGradientConst, evalZAndGradients)
  
  int    setAnam(const AAnam* anam, const VectorInt& strcnt = VectorInt());
  int    unsetAnam();
  bool   isFlagGradient() const;
  bool   isFlagGradientNumerical() const;
  bool   isFlagGradientFunctional() const;
  void   switchToGradient();

  ////////////////////////////////////////////////
  /// TODO : to be removed (encapsulation of CovAnisoList)
  const CovAnisoList* getCovAnisoList() const;
  CovAnisoList* getCovAnisoListModify();  
  
  double evalCov(const VectorDouble& incr,
                 int icov = 0,
                 const ECalcMember& member = ECalcMember::fromKey("LHS")) const;

  ////////////////////////////////////////////////
  /// TODO : to be removed (encapsulation of Context)
  void setField(double field);
  /////////////////////////////////////////////////

  Model* duplicate() const;
  Model* createReduce(const VectorInt& validVars) const;

  int getNVar() const;

  int fitFromCovIndices(Vario *vario,
                        const VectorECov &types = ECov::fromKeys({"EXPONENTIAL"}),
                        const Constraints& constraints = Constraints(),
                        const Option_VarioFit& optvar = Option_VarioFit(),
                        const Option_AutoFit& mauto = Option_AutoFit(),
                        bool verbose = false);
  int fit(Vario* vario,
          const VectorECov& types        = ECov::fromKeys({"SPHERICAL"}),
          const Constraints& constraints = Constraints(),
          const Option_VarioFit& optvar = Option_VarioFit(),
          const Option_AutoFit& mauto = Option_AutoFit(),
          bool verbose = false);

  int fitFromVMap(DbGrid *dbmap,
                  const VectorECov &types = ECov::fromKeys({"SPHERICAL"}),
                  const Constraints &constraints = Constraints(),
                  const Option_VarioFit& optvar = Option_VarioFit(),
                  const Option_AutoFit& mauto = Option_AutoFit(),
                  bool verbose = false);
  
  int stabilize(double percent, bool verbose = false);
  int standardize(bool verbose = false);

  static void gofDisplay(double gof,
                         bool byValue                   = true,
                         const VectorDouble& thresholds = {2., 5., 10., 100});
  static VectorECov initCovList(const VectorInt & covranks);

  bool isValid() const;

protected:
  /// Interface to ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "Model"; }

private:
  bool _isValid() const override;
  void _clear();
  void _create();
  void _copyCovContext();
};

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

#include "Anamorphosis/AAnam.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/CovAnisoList.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/CovLMCAnamorphosis.hpp"
#include "Covariances/CovLMCTapering.hpp"
#include "Drifts/DriftList.hpp"
#include "Enum/ECalcMember.hpp"
#include "Enum/ECov.hpp"
#include "Enum/EModelProperty.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Model/Constraints.hpp"
#include "Model/ModelCovList.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{

class Db;
class CovLMCTapering;
class CovLMCAnamorphosis;
class CovInternal;
class CovCalcMode;

class ADrift;

class Model;
class Vario;
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
class GSTLEARN_EXPORT Model: public AStringable, public ASerializable, public ModelCovList
{
public:
  Model(const CovContext& ctxt = CovContext());
  Model(Id nvar, Id ndim = 2);
  Model(const Model& m);
  Model& operator=(const Model& m);
  virtual ~Model();

public:
  /// ICloneable interface
  IMPLEMENT_CLONING(Model)

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

public:
  const CovAnisoList* castInCovAnisoListConst(Id icov = -1) const;
  const CovLMCTapering* castInCovLMCTaperingConst() const;
  const CovLMCAnamorphosis* castInCovLMCAnamorphosisConst() const;

public:
  CovAnisoList* _castInCovAnisoList(Id icov = -1);
  CovLMCTapering* _castInCovLMCTapering();
  CovLMCAnamorphosis* _castInCovLMCAnamorphosis();

public:
  Id resetFromDb(const Db* db);
  static Model* create(const CovContext& ctxt = CovContext());
  static Model* createFromEnvironment(Id nvar, Id ndim = 2);
  static Model* createNugget(Id nvar, Id ndim = 2, double sill = 1.);
  static Model* createFromParam(const ECov& type             = ECov::fromKey("NUGGET"),
                                double range                 = 1.,
                                double sill                  = 1.,
                                double param                 = 1.,
                                const VectorDouble& ranges   = VectorDouble(),
                                const MatrixSymmetric& sills = MatrixSymmetric(),
                                const VectorDouble& angles   = VectorDouble(),
                                const ASpaceSharedPtr& space = ASpaceSharedPtr(),
                                bool flagRange               = true);
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
  static Model* createFromNF(const String& NFFilename, bool verbose = true);
  static Model* createFromVario(Vario* vario,
                                const VectorECov& types        = ECov::fromKeys({"SPHERICAL"}),
                                const Constraints& constraints = Constraints(),
                                const Option_VarioFit& optvar  = Option_VarioFit(),
                                const Option_AutoFit& mauto    = Option_AutoFit(),
                                bool verbose                   = false);
  static Model* createFillRandom(Id ndim,
                                 Id nvar,
                                 const std::vector<ECov>& types = ECov::fromKeys({"SPHERICAL"}),
                                 double hmax                    = 1,
                                 Id order                       = -1,
                                 Id nfex                        = 0,
                                 Id seed                        = 13242);
  void setCovAnisoList(const CovAnisoList* covalist);
  void addCovAniso(const CovAniso& cov);
  void addCovFromParam(const ECov& type,
                       double range                 = EPSILON6,
                       double sill                  = 1.,
                       double param                 = 1.,
                       const VectorDouble& ranges   = VectorDouble(),
                       const MatrixSymmetric& sills = MatrixSymmetric(),
                       const VectorDouble& angles   = VectorDouble(),
                       bool flagRange               = true);
  void addCovFromParamOldStyle(const ECov& type,
                               double range               = EPSILON6,
                               double sill                = 1.,
                               double param               = 1.,
                               const VectorDouble& ranges = VectorDouble(),
                               const VectorDouble& sills  = VectorDouble(),
                               const VectorDouble& angles = VectorDouble(),
                               bool flagRange             = true);

  FORWARD_METHOD(castInCovAnisoListConst, getActiveFactor, ITEST)
  FORWARD_METHOD(castInCovAnisoListConst, getCovAniso)
  FORWARD_METHOD(castInCovAnisoListConst, getNCov, ITEST)
  FORWARD_METHOD(castInCovAnisoListConst, getCovType, ECov::UNKNOWN)
  FORWARD_METHOD(castInCovAnisoListConst, getRange, TEST)
  FORWARD_METHOD(castInCovAnisoListConst, getRanges)
  FORWARD_METHOD(castInCovAnisoListConst, getAngles)
  FORWARD_METHOD(castInCovAnisoListConst, getAnam)
  FORWARD_METHOD(castInCovAnisoListConst, getParam, TEST)
  FORWARD_METHOD(castInCovAnisoListConst, getCovName)
  FORWARD_METHOD(castInCovAnisoListConst, extractCova)
  FORWARD_METHOD(castInCovAnisoListConst, getNGradParam, ITEST)
  FORWARD_METHOD(castInCovAnisoListConst, getMaximumDistance, TEST)
  FORWARD_METHOD(castInCovAnisoListConst, getCovMinIRFOrder, ITEST)
  FORWARD_METHOD(castInCovAnisoListConst, getAnamNClass, ITEST)
  FORWARD_METHOD(castInCovAnisoListConst, hasAnam, false)
  FORWARD_METHOD(castInCovAnisoListConst, hasNugget, false)
  FORWARD_METHOD(castInCovAnisoListConst, getRankNugget, -1)
  FORWARD_METHOD(castInCovAnisoListConst, getBallRadius, TEST)
  FORWARD_METHOD(castInCovAnisoListConst, hasExternalCov)
  FORWARD_METHOD(castInCovAnisoListConst, isChangeSupportDefined, false)
  FORWARD_METHOD(castInCovAnisoListConst, getAnamHermite)
  FORWARD_METHOD(castInCovAnisoListConst, getCovMode, EModelProperty::NONE)

  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, setActiveFactor)
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, getCovAniso)
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, setRangeIsotropic)
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, setMarkovCoeffs)
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, makeRangeNoStatDb);
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, makeScaleNoStatDb);
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, makeAngleNoStatDb);
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, makeTensorNoStatDb);
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, makeParamNoStatDb);
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, makeRangeNoStatFunctional);
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, makeScaleNoStatFunctional);
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, makeAngleNoStatFunctional);
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, makeTensorNoStatFunctional);
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, makeParamNoStatFunctional);
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, makeRangeStationary);
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, makeScaleStationary);
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, makeAngleStationary);
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, makeTensorStationary);
  FORWARD_METHOD_NON_CONST(_castInCovAnisoList, makeParamStationary);

  FORWARD_METHOD_NON_CONST(_castInCovLMCTapering, setTapeRange)

  Id setAnam(const AAnam* anam, const VectorInt& strcnt = VectorInt());
  Id unsetAnam();

  ////////////////////////////////////////////////
  /// TODO : to be removed (encapsulation of CovAnisoList)
  const CovAnisoList* getCovAnisoList() const;
  CovAnisoList* getCovAnisoListModify();

  double evalCovFromIncr(const VectorDouble& incr,
                         Id icov                   = 0,
                         const ECalcMember& member = ECalcMember::fromKey("LHS")) const;

  Model* duplicate() const;
  Model* createReduce(const VectorInt& validVars) const;

  Id getNVar() const;

  Id fitFromCovIndices(Vario* vario,
                       const VectorECov& types        = ECov::fromKeys({"EXPONENTIAL"}),
                       const Constraints& constraints = Constraints(),
                       const Option_VarioFit& optvar  = Option_VarioFit(),
                       const Option_AutoFit& mauto    = Option_AutoFit(),
                       bool verbose                   = false);
  Id fit(Vario* vario,
         const VectorECov& types        = ECov::fromKeys({"SPHERICAL"}),
         const Constraints& constraints = Constraints(),
         const Option_VarioFit& optvar  = Option_VarioFit(),
         const Option_AutoFit& mauto    = Option_AutoFit(),
         bool verbose                   = false);

  Id fitFromVMap(DbGrid* dbmap,
                 const VectorECov& types        = ECov::fromKeys({"SPHERICAL"}),
                 const Constraints& constraints = Constraints(),
                 const Option_VarioFit& optvar  = Option_VarioFit(),
                 const Option_AutoFit& mauto    = Option_AutoFit(),
                 bool verbose                   = false);

  Id stabilize(double percent, bool verbose = false);
  Id standardize(bool verbose = false);

  static void gofDisplay(double gof,
                         bool byValue                   = true,
                         const VectorDouble& thresholds = {2., 5., 10., 100});
  static VectorECov initCovList(const VectorInt& covranks);

  bool isValid() const;

protected:
  /// Interface to ASerializable
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "Model"; }

private:
  bool _isValid() const override;
  void _clear();
  void _create();
  void _copyCovContext();
};
} // namespace gstlrn
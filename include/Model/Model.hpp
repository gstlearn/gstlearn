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
#include "geoslib_define.h"

#include "Enum/ECalcMember.hpp"
#include "Enum/ECov.hpp"
#include "Enum/EConsElem.hpp"
#include "Enum/EModelProperty.hpp"

#include "Covariances/ACov.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/ACovAnisoList.hpp"
#include "Covariances/CovLMGradient.hpp"

#include "Drifts/DriftList.hpp"

#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"
#include "Model/Constraints.hpp"
#include "Model/CovParamId.hpp"
#include "Covariances/CovAniso.hpp"

#include "Anamorphosis/AAnam.hpp"

#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/ICloneable.hpp"

class Model;
class Db;
class CovInternal;
class MatrixSquareSymmetric;
class CovCalcMode;
class Vario;
class ANoStat;
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
class GSTLEARN_EXPORT Model : public AStringable, public ASerializable, public ICloneable
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

  int resetFromDb(const Db* db);

  static Model* create(const CovContext& ctxt = CovContext());
  static Model* createFromEnvironment(int nvar, int ndim = 2);
  static Model* createFromParam(const ECov& type = ECov::fromKey("NUGGET"),
                                double range = 1.,
                                double sill = 1.,
                                double param = 1.,
                                const VectorDouble& ranges = VectorDouble(),
                                const VectorDouble& sills = VectorDouble(),
                                const VectorDouble& angles = VectorDouble(),
                                const ASpace* space = nullptr,
                                bool flagRange = true);
  static Model* createFromDb(const Db* db);
  static Model* createFromNF(const String& neutralFilename, bool verbose = true);

  void   setCovList(const ACovAnisoList* covalist);
  void   addCov(const CovAniso* cov);
  void   addCovFromParam(const ECov& type,
                         double range = 0.,
                         double sill = 1.,
                         double param = 1.,
                         const VectorDouble& ranges = VectorDouble(),
                         const VectorDouble& sills  = VectorDouble(),
                         const VectorDouble& angles = VectorDouble(),
                         bool flagRange = true);
  void   delCova(int icov);
  void   delAllCovas();
  void   setDriftList(const DriftList* driftlist);
  void   setDriftIRF(int order = 0, int nfex = 0);
  void   addDrift(const ADrift* drift);  // TODO: check that the same driftM has not been already defined
  void   setDrifts(const VectorString& driftSymbols);
  void   delDrift(int rank);
  void   delAllDrifts();
  int    setAnam(const AAnam* anam, const VectorInt& strcnt = VectorInt());
  int    unsetAnam();
  bool   isFlagGradient() const;
  bool   isFlagGradientNumerical() const;
  bool   isFlagGradientFunctional() const;
  bool   isFlagLinked() const;
  CovAniso extractCova(int icov) const;
  void   switchToGradient();

  ////////////////////////////////////////////////
  /// TODO : to be removed (encapsulation of ACovAnisoList)
  const ACovAnisoList* getCovAnisoList() const;
  const CovAniso* getCova(unsigned int icov) const;
  CovAniso* getCova(unsigned int icov);
  int getCovaNumber() const;
  const ECov& getCovaType(int icov) const;
  const MatrixSquareSymmetric getSillValues(int icov) const;
  double getSill(int icov, int ivar, int jvar) const;
  double getParam(int icov) const;
  bool isCovaFiltered(int icov) const;
  bool isStationary() const;
  String getCovName(int icov) const;
  int getGradParamNumber(int icov) const;
  double getTotalSill(int ivar, int jvar) const;
  double getBallRadius() const;
  const AnamHermite* getAnamHermite() const;

  double getMaximumDistance() const;
  int getCovaMinIRFOrder() const;
  bool hasAnam() const;
  const AAnam* getAnam() const;
  bool isChangeSupportDefined() const;
  void normalize(double sill);
  bool hasNugget() const;
  VectorInt getActiveCovList() const;
  VectorInt getAllActiveCovList() const;
  bool isAllActiveCovList() const;
  void setTapeRange(double range);

  void setOptimEnabled(bool flagOptim) { _cova->setOptimEnabled(flagOptim); }
  bool isOptimEnabled() const { return _cova->isOptimEnabled(); }

  double eval0(int ivar = 0,
               int jvar = 0,
               const CovCalcMode* mode = nullptr) const
  {
    return _cova->eval0(ivar, jvar, mode);
  }
  MatrixSquareGeneral eval0Nvar(const CovCalcMode* mode = nullptr) const
  {
    return _cova->eval0Mat(mode);
  }
  int isNoStat() const
  {
    return _cova->isNoStat();
  }
  const ANoStat* getNoStat() const
  {
    return _cova->getNoStat();
  }
  ANoStat* getNoStatModify() const
  {
    return _cova->getNoStatModify();
  }
  void eval0MatInPlace(MatrixSquareGeneral &mat,
                       const CovCalcMode *mode = nullptr) const
  {
    _cova->eval0MatInPlace(mat, mode);
  }
  double eval(const SpacePoint& p1,
              const SpacePoint& p2,
              int ivar = 0,
              int jvar = 0,
              const CovCalcMode* mode = nullptr) const
  {
    return _cova->eval(p1, p2, ivar, jvar, mode);
  }
  MatrixSquareGeneral evalNvarIpas(double step,
                                   const VectorDouble& dir = VectorDouble(),
                                   const VectorDouble& center = VectorDouble(),
                                   const CovCalcMode* mode = nullptr) const
  {
    return _cova->evalNvarIpas(step, dir, center, mode);
  }
  MatrixSquareGeneral evalMat(const SpacePoint& p1,
                              const SpacePoint& p2,
                              const CovCalcMode* mode = nullptr) const
  {
    return _cova->evalMat(p1, p2, mode);
  }

  void evalMatInPlace(const SpacePoint &p1,
                      const SpacePoint &p2,
                      MatrixSquareGeneral &mat,
                      const CovCalcMode* mode = nullptr) const
  {
    _cova->evalMatInPlace(p1, p2, mat, mode);
  }

  MatrixSquareGeneral evalNvarIpasIncr(const VectorDouble& dincr,
                                       const CovCalcMode* mode = nullptr) const
  {
    return _cova->evalNvarIpasIncr(dincr, mode);
  }
  VectorDouble evalIvarNpas(const VectorDouble& vec_step,
                            const VectorDouble& dir = VectorDouble(),
                            int ivar = 0,
                            int jvar = 0,
                            const VectorDouble& center = VectorDouble(),
                            const CovCalcMode* mode = nullptr) const
  {
    return _cova->evalIvarNpas(vec_step, dir, ivar, jvar, center, mode);
  }
  double evalIvarIpas(double step,
                      const VectorDouble& dir = VectorDouble(),
                      int ivar = 0,
                      int jvar = 0,
                      const VectorDouble& center = VectorDouble(),
                      const CovCalcMode* mode = nullptr) const
  {
    return _cova->evalIvarIpas(step, dir, ivar, jvar, center, mode);
  }
  double evalCvv(const VectorDouble& ext,
                 const VectorInt& ndisc,
                 const VectorDouble& angles = VectorDouble(),
                 int ivar = 0,
                 int jvar = 0,
                 const CovCalcMode* mode = nullptr) const
  {
    return _cova->evalCvv(ext, ndisc, angles, ivar, jvar, mode);
  }
  double evalCvvShift(const VectorDouble& ext,
                      const VectorInt& ndisc,
                      const VectorDouble& shift,
                      const VectorDouble& angles = VectorDouble(),
                      int ivar = 0,
                      int jvar = 0,
                      const CovCalcMode* mode = nullptr) const
  {
    return _cova->evalCvvShift(ext, ndisc, shift, angles, ivar, jvar, mode);
  }
  MatrixSquareGeneral evalCvvM(const VectorDouble& ext,
                               const VectorInt& ndisc,
                               const VectorDouble& angles = VectorDouble(),
                               const CovCalcMode* mode = nullptr)
  {
    return _cova->evalCvvM(ext, ndisc, angles, mode);
  }
  double evalCxv(const SpacePoint& p1,
                 const VectorDouble& ext,
                 const VectorInt& ndisc,
                 const VectorDouble& angles = VectorDouble(),
                 const VectorDouble& x0 = VectorDouble(),
                 int ivar = 0,
                 int jvar = 0,
                 const CovCalcMode* mode = nullptr)
  {
    return _cova->evalCxv(p1, ext, ndisc, angles, x0, ivar, jvar, mode);
  }
  MatrixSquareGeneral evalCxvM(const SpacePoint& p1,
                               const VectorDouble& ext,
                               const VectorInt& ndisc,
                               const VectorDouble& angles = VectorDouble(),
                               const VectorDouble& x0 = VectorDouble(),
                               const CovCalcMode* mode = nullptr)
  {
    return _cova->evalCxvM(p1, ext, ndisc, angles, x0, mode);
  }
  VectorDouble evalPointToDb(const SpacePoint& p1,
                             const Db* db2,
                             int ivar = 0,
                             int jvar = 0,
                             bool useSel = true,
                             const VectorInt& nbgh2 = VectorInt(),
                             const CovCalcMode* mode = nullptr)
  {
    return _cova->evalPointToDb(p1, db2, ivar, jvar, useSel, nbgh2, mode);
  }
  VectorDouble evalPointToDbAsSP(const std::vector<SpacePoint>& p1s,
                                 const SpacePoint& p2,
                                 int ivar = 0,
                                 int jvar = 0,
                                 const CovCalcMode* mode = nullptr) const
  {
    return _cova->evalPointToDbAsSP(p1s, p2, ivar, jvar, mode);
  }
  double evalAverageDbToDb(const Db* db1,
                           const Db* db2,
                           int ivar = 0,
                           int jvar = 0,
                           const CovCalcMode* mode = nullptr) const
  {
    return _cova->evalAverageDbToDb(db1, db2, ivar, jvar, mode);
  }
  double evalAverageIncrToIncr(const VectorVectorDouble& d1,
                               const VectorVectorDouble& d2,
                               int ivar = 0,
                               int jvar = 0,
                               const CovCalcMode* mode = nullptr) const
  {
    return _cova->evalAverageIncrToIncr(d1, d2, ivar, jvar, mode);
  }

  double evalAveragePointToDb(const SpacePoint& p1,
                              const Db* db2,
                              int ivar = 0,
                              int jvar = 0,
                              const CovCalcMode* mode = nullptr)
  {
    return _cova->evalAveragePointToDb(p1, db2, ivar, jvar, mode);
  }
  MatrixRectangular evalCovMatrix(Db* db1,
                                  Db* db2 = nullptr,
                                  int ivar = 0,
                                  int jvar = 0,
                                  const VectorInt& nbgh1 = VectorInt(),
                                  const VectorInt& nbgh2 = VectorInt(),
                                  const CovCalcMode* mode = nullptr)
  {
    return _cova->evalCovMatrix(db1, db2, ivar, jvar, nbgh1, nbgh2, mode);
  }

  void evalMatOptimInPlace(int icas1,
                           int iech1,
                           int icas2,
                           int iech2,
                           MatrixSquareGeneral &mat,
                           const CovCalcMode *mode = nullptr) const
  {
    _cova->evalMatOptimInPlace(icas1, iech1, icas2, iech2, mat, mode);
  }

  VectorVectorDouble evalCovMatrixOptim(const Db *db1,
                                        const Db *db2 = nullptr,
                                        int ivar = 0,
                                        int jvar = 0,
                                        const CovCalcMode *mode = nullptr);

  double extensionVariance(const Db* db,
                           const VectorDouble& ext,
                           const VectorInt& ndisc,
                           const VectorDouble& angles = VectorDouble(),
                           const VectorDouble& x0 = VectorDouble(),
                           int ivar = 0,
                           int jvar = 0)
  {
    return _cova->extensionVariance(db, ext, ndisc, angles, x0, ivar, jvar);
  }
  double samplingDensityVariance(const Db* db,
                                 const VectorDouble& ext,
                                 const VectorInt& ndisc,
                                 const VectorDouble& angles = VectorDouble(),
                                 const VectorDouble& x0 = VectorDouble(),
                                 int ivar = 0,
                                 int jvar = 0) const
  {
    return _cova->samplingDensityVariance(db, ext, ndisc, angles, x0, ivar, jvar);
  }
  double specificVolume(const Db *db,
                        double mean,
                        const VectorDouble &ext,
                        const VectorInt &ndisc,
                        const VectorDouble &angles = VectorDouble(),
                        const VectorDouble &x0 = VectorDouble(),
                        int ivar = 0,
                        int jvar = 0) const
  {
    return _cova->specificVolume(db, mean, ext, ndisc, angles, x0, ivar, jvar);
  }
  double coefficientOfVariation(const Db *db,
                                double volume,
                                double mean,
                                const VectorDouble &ext,
                                const VectorInt &ndisc,
                                const VectorDouble &angles = VectorDouble(),
                                const VectorDouble &x0 = VectorDouble(),
                                int ivar = 0,
                                int jvar = 0) const
  {
    return _cova->coefficientOfVariation(db, volume, mean, ext, ndisc, angles, x0, ivar, jvar);
  }
  double specificVolumeFromCoV(Db *db,
                               double cov,
                               double mean,
                               const VectorDouble &ext,
                               const VectorInt &ndisc,
                               const VectorDouble &angles = VectorDouble(),
                               const VectorDouble &x0 = VectorDouble(),
                               int ivar = 0,
                               int jvar = 0) const
  {
    return _cova->specificVolumeFromCoV(db, cov, mean, ext, ndisc, angles, x0, ivar, jvar);
  }
  void evalZAndGradients(const SpacePoint& p1,
                         const SpacePoint& p2,
                         double& covVal,
                         VectorDouble& covGp,
                         VectorDouble& covGG,
                         const CovCalcMode* mode = nullptr,
                         bool flagGrad = false) const;
  void evalZAndGradients(const VectorDouble& vec,
                         double& covVal,
                         VectorDouble& covGp,
                         VectorDouble& covGG,
                         const CovCalcMode* mode = nullptr,
                         bool flagGrad = false) const;

  void setSill(int icov, int ivar, int jvar, double value);
  void setCovaFiltered(int icov, bool filtered);
  void updateCovByPoints(int icas1, int iech1, int icas2, int iech2);
  void updateCovByMesh(int imesh);
  void setActiveFactor(int iclass);
  int  getActiveFactor() const;
  int  getAnamNClass() const;
  /////////////////////////////////////////////////

  ////////////////////////////////////////////////
  /// TODO : to be removed (encapsulation of DriftList)
  const DriftList* getDriftList()                  const;
  const ADrift* getDrift(int il)                   const;
  ADrift* getDrift(int il)                              ;
  int getDriftNumber()                             const;
  int getExternalDriftNumber()                     const;
  int getRankFext(int il)                          const;
  const VectorDouble& getDriftCoefs()              const;
  double getDriftCoef(int ivar, int il, int ib)    const;
  int getDriftEquationNumber()                     const;
  bool isDriftFiltered(unsigned int il)            const;
  bool isDriftDefined(const VectorInt &powers, int rank_fex = 0) const;
  bool isDriftDifferentDefined(const VectorInt &powers, int rank_fex = -1) const;
  int getDriftMaxIRFOrder(void) const;

  void resetDriftCoef() { _driftList->resetDriftCoeff(); }
  void setDriftCoef(int ivar, int il, int ib, double coeff)    ;
  void setDriftFiltered(int il, bool filtered)                 ;
  VectorDouble getDriftByColumn(const Db* db, int ib, bool useSel=true);
  VectorVectorDouble getDrifts(const Db* db, bool useSel=true) ;

  double evalDrift(const Db* db,
                   int iech,
                   int il,
                   const ECalcMember& member = ECalcMember::fromKey("LHS")) const;
  VectorDouble evalDriftVec(const Db* db,
                            int iech,
                            const ECalcMember& member = ECalcMember::fromKey("LHS")) const;
  VectorDouble evalDrifts(const Db* db,
                          const VectorDouble& coeffs,
                          int ivar = 0,
                          bool useSel = false) const;
  void evalDriftVecInPlace(const Db* db,
                           int iech,
                           const ECalcMember& member,
                           VectorDouble& drftab) const;
  double evalDriftCoef(const Db* db,
                        int iech,
                        int ivar,
                        const double* coef) const;
  /////////////////////////////////////////////////

  ////////////////////////////////////////////////
  /// TODO : to be removed (encapsulation of Context)
  const CovContext& getContext() const { return _ctxt; }
  const VectorDouble& getMeans() const { return _ctxt.getMean(); }
  double getMean(int ivar) const { return _ctxt.getMean(ivar); }
  const VectorDouble& getCovar0s() const { return _ctxt.getCovar0(); }
  double getCovar0(int ivar, int jvar) const { return _ctxt.getCovar0(ivar,jvar); }
  double getField() const               { return _ctxt.getField(); }
  int getDimensionNumber() const        { return _ctxt.getNDim(); }

  void setMeans(const VectorDouble& mean);
  void setMean(double mean, int ivar=0);
  void setCovar0s(const VectorDouble& covar0);
  void setCovar0(int ivar, int jvar, double covar0);
  void setField(double field);
  /////////////////////////////////////////////////

  /////////////////////////////////////////////////
  /// Shortcut for Non-stationary
  int  addNoStat(const ANoStat* anostat);
  int  getNoStatElemNumber() const;
  int  addNoStatElem(int igrf, int icov, const EConsElem& type, int iv1, int iv2);
  int  addNoStatElems(const VectorString& codes);
  CovParamId getCovParamId(int ipar) const;
  ////////////////////////////////////////////////

  const EModelProperty& getCovMode() const;
  Model* duplicate() const;
  Model* reduce(const VectorInt& validVars) const;

  int getVariableNumber() const
  {
    // TODO/ the strange next line have been commented out.
    // There should be either validated or suppressed
    //    if (isFlagGradient())
    //      return 3; // This strange number of variables is linked to the Gradient calculation
    //    else
    // However, note used for Gradient (Functional type) in Potential
      return _ctxt.getNVar();
  }

  int hasExternalCov() const;

  MatrixSquareSymmetric covMatrixM(Db *db1,
                                   Db *db2 = nullptr,
                                   int ivar = -1,
                                   int jvar = -1,
                                   const CovCalcMode* mode = nullptr);
  VectorDouble covMatrixV(Db *db1,
                          Db *db2 = nullptr,
                          int ivar = 0,
                          int jvar = 0,
                          const CovCalcMode* mode = nullptr);
  void covMatrix(VectorDouble& covmat,
                 Db *db1,
                 Db *db2 = nullptr,
                 int ivar = 0,
                 int jvar = 0,
                 const CovCalcMode* mode = nullptr);
  VectorDouble sample(const VectorDouble& hh,
                      int ivar = 0,
                      int jvar = 0,
                      VectorDouble codir = VectorDouble(),
                      const CovCalcMode* mode = nullptr);
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
  int fitFromCovIndices(Vario *vario,
                        const VectorECov &types = ECov::fromKeys({"EXPONENTIAL"}),
                        const Constraints& constraints = Constraints(),
                        Option_VarioFit optvar = Option_VarioFit(),
                        Option_AutoFit mauto = Option_AutoFit(),
                        bool verbose = false);
  int fit(Vario *vario,
          const VectorECov& types = ECov::fromKeys({"SPHERICAL"}),
          const Constraints& constraints = Constraints(),
          Option_VarioFit optvar = Option_VarioFit(),
          Option_AutoFit mauto = Option_AutoFit(),
          bool verbose = false);

  int fitFromVMap(DbGrid *dbmap,
                  const VectorECov &types = ECov::fromKeys({"SPHERICAL"}),
                  const Constraints &constraints = Constraints(),
                  Option_VarioFit optvar = Option_VarioFit(),
                  Option_AutoFit mauto = Option_AutoFit(),
                  bool verbose = false);

  double gofToVario(const Vario* vario, bool verbose = true);
  void gofDisplay(double gof, bool byValue = true,
                  const VectorDouble& thresholds = {2., 5., 10., 100});
  VectorECov initCovList(const VectorInt & covranks);

  bool isValid() const;

protected:
  /// Interface to ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "Model"; }

  const ACovAnisoList* _castInCovAnisoListConst(int icov = -1) const;
  ACovAnisoList*       _castInCovAnisoList(int icov = -1);

private:
  void _clear();
  void _create();
  void _copyCovContext();

private:
  ACov*          _cova;         /* Generic Covariance structure */
  DriftList*     _driftList;    /* Series of Drift functions */
  CovContext     _ctxt;         /* Context */
};

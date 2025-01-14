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

  int resetFromDb(const Db* db);

  static Model* create(const CovContext& ctxt = CovContext());
  static Model* createFromEnvironment(int nvar, int ndim = 2);
  static Model* createNugget(int nvar, int ndim = 2, double sill = 1.);
  static Model* createFromParam(const ECov& type = ECov::fromKey("NUGGET"),
                                double range     = 1.,
                                double sill      = 1.,
                                double param     = 1.,
                                const VectorDouble& ranges = VectorDouble(),
                                const MatrixSquareSymmetric& sills  = MatrixSquareSymmetric(),
                                const VectorDouble& angles = VectorDouble(),
                                const ASpace* space        = nullptr,
                                bool flagRange             = true);
  static Model* createFromParamOldStyle(const ECov& type = ECov::fromKey("NUGGET"),
                          double range               = 1.,
                          double sill                = 1.,
                          double param               = 1.,
                          const VectorDouble& ranges = VectorDouble(),
                          const VectorDouble& sills  = VectorDouble(),
                          const VectorDouble& angles = VectorDouble(),
                          const ASpace* space        = nullptr,
                          bool flagRange             = true);
  static Model* createFromDb(const Db* db);
  static Model* createFromNF(const String& neutralFilename,
                             bool verbose = true);
  static Model* createFromVario(Vario* vario,
                  const VectorECov& types = ECov::fromKeys({"SPHERICAL"}),
                  const Constraints& constraints = Constraints(),
                  const Option_VarioFit& optvar  = Option_VarioFit(),
                  const Option_AutoFit& mauto    = Option_AutoFit(),
                  bool verbose                   = false);

  void   setCovList(const CovAnisoList* covalist);
  void addCov(const CovAniso* cov);
  void
  addCovFromParam(const ECov& type,
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
 
  void   setDriftList(const DriftList* driftlist);
  void   setDriftIRF(int order = 0, int nfex = 0);
  void   setFlagLinked(bool flagLinked);
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
  void switchToGradient();
  bool   hasDrift() const;

  ////////////////////////////////////////////////
  /// TODO : to be removed (encapsulation of CovAnisoList)
  const CovAnisoList* getCovAnisoList() const;
  CovAnisoList* getCovAnisoListModify();

  const CovAniso* getCova(int icov) const;
  CovAniso* getCova(int icov);
  int getCovaNumber(bool skipNugget = false) const;
  const ECov& getCovaType(int icov) const;
  
  double getRange(int icov) const;
  VectorDouble getRanges(int icov) const;
  double getParam(int icov) const;
  String getCovName(int icov) const;
  int getGradParamNumber(int icov) const;
  

  double getBallRadius() const;
  const AnamHermite* getAnamHermite() const;

  double getMaximumDistance() const;
  int getCovaMinIRFOrder() const;
  bool hasAnam() const;
  const AAnam* getAnam() const;
  bool isChangeSupportDefined() const;
  void normalize(double sill);
  bool hasNugget() const;
  int  getRankNugget() const;

  void setTapeRange(double range);

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

  /**
   * Calculate the Matrix of covariance for zero distance
   * @param mat   Covariance matrix (Dimension: nvar * nvar)
   * @param mode  Calculation Options
   *
   * @remarks: Matrix 'mat' should be dimensioned and initialized beforehand
   */
  void eval0MatInPlace(MatrixSquareGeneral &mat,
                       const CovCalcMode *mode = nullptr) const
  {
    _cova->eval0CovMatBiPointInPlace(mat, mode);
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
                                   const VectorDouble& dir,
                                   const CovCalcMode* mode = nullptr) const
  {
    return _cova->evalNvarIpas(step, dir, mode);
  }
  MatrixSquareGeneral evalMat(const SpacePoint& p1,
                              const SpacePoint& p2,
                              const CovCalcMode* mode = nullptr) const
  {
    return _cova->evalMat(p1, p2, mode);
  }

  /**
   * Calculate the Matrix of covariance between two space points
   * @param p1 Reference of the first space point
   * @param p2 Reference of the second space point
   * @param mat   Covariance matrix (Dimension: nvar * nvar)
   * @param mode  Calculation Options
   *
   * @remarks: Matrix 'mat' should be dimensioned and initialized beforehand
   */
  void evalMatInPlace(const SpacePoint &p1,
                      const SpacePoint &p2,
                      MatrixSquareGeneral &mat,
                      const CovCalcMode* mode = nullptr) const
  {
    _cova->evalCovMatBiPointInPlace(mat, p1, p2, mode);
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
                            const CovCalcMode* mode = nullptr) const
  {
    return _cova->evalIvarNpas(vec_step, dir, ivar, jvar, mode);
  }
  double evalIvarIpas(double step,
                      const VectorDouble& dir = VectorDouble(),
                      int ivar = 0,
                      int jvar = 0,
                      const CovCalcMode* mode = nullptr) const
  {
    return _cova->evalIvarIpas(step, dir, ivar, jvar, mode);
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
                           double eps = 0.,
                           int seed = 434132,
                           const CovCalcMode* mode = nullptr) const
  {
    return _cova->evalAverageDbToDb(db1, db2, ivar, jvar, eps, seed, mode);
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
  /**
   * \defgroup Model Model: Set of classes for processing Model contents
   *
   **/

  /** @addtogroup Model_0 Calculating Covariance Matrix
   * \ingroup Model
   *
   * These functions are meant to calculate the covariance Matrix between two Dbs
   * or between a Db and itself.
   * They take into account the presence of a possible selection
   * They also account for heterotopy (if Z-variables are defined in the Db(s)
   *
   * @param  db1   First Db
   * @param  db2   (Optional second Db)
   * @param  ivar0 Rank of the selected variable in the first Db (-1 for all variables)
   * @param  jvar0 Rank of the selected variable in the second Db (-1 for all variables)
   * @param  nbgh1 Vector of indices of active samples in first Db (optional)
   * @param  nbgh2 Vector of indices of active samples in second Db (optional)
   * @param  mode  CovCalcMode structure
   *
   * @remarks The returned matrix if dimension to nrows * ncols where
   * @remarks each term is the product of the number of active samples
   * @remarks by the number of samples where the variable is defined
   *
   * @note 'dbin' and 'dbout' cannot be made 'const' as they can be updated
   * @note due to the presence of 'nostat'
   *
   * @return A Matrix either in Dense or Sparse format
   *
   *  @{
   */
  VectorDouble evalCovMatrixV(Db *db1,
                              Db *db2 = nullptr,
                              int ivar0 = -1,
                              int jvar0 = -1,
                              const VectorInt &nbgh1 = VectorInt(),
                              const VectorInt &nbgh2 = VectorInt(),
                              const CovCalcMode *mode = nullptr)
  {
    if (_cova == nullptr) return VectorDouble();
    return _cova->evalCovMatrix(db1, db2, ivar0, jvar0, nbgh1, nbgh2, mode).getValues();
  }
  MatrixRectangular evalCovMatrixOptim(const Db *db1,
                                       const Db *db2 = nullptr,
                                       int ivar0 = -1,
                                       int jvar0 = -1,
                                       const VectorInt &nbgh1 = VectorInt(),
                                       const VectorInt &nbgh2 = VectorInt(),
                                       const CovCalcMode *mode = nullptr)
  {
    const CovAnisoList *covalist = _castInCovAnisoListConst();
    if (covalist == nullptr) return MatrixRectangular();
    return covalist->evalCovMatrixOptim(db1, db2, ivar0, jvar0, nbgh1, nbgh2, mode);
  }

  MatrixSquareSymmetric evalCovMatrixSymmetricOptim(const Db *db1,
                                                    int ivar0 = -1,
                                                    const VectorInt &nbgh1 = VectorInt(),
                                                    const CovCalcMode *mode = nullptr)
  {
    const CovAnisoList *covalist = _castInCovAnisoListConst();
    if (covalist == nullptr) return MatrixRectangular();
    return covalist->evalCovMatrixSymmetricOptim(db1, ivar0, nbgh1, mode);
  }
  
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

  double evalCov(const VectorDouble& incr,
                 int icov = 0,
                 const ECalcMember& member = ECalcMember::fromKey("LHS")) const;

  void setSill(int icov, int ivar, int jvar, double value);
  void setRangeIsotropic(int icov, double range);
  void setMarkovCoeffs(int icov, const VectorDouble& coeffs);
  void setCovaFiltered(int icov, bool filtered);
  void setActiveFactor(int iclass);
  int  getActiveFactor() const;
  int  getAnamNClass() const;
  /////////////////////////////////////////////////

  ////////////////////////////////////////////////
  /// TODO : to be removed (encapsulation of DriftList)
  const DriftList* getDriftList()                  const;
  const ADrift* getDrift(int il)                   const;
  int  getDriftNumber()                            const;
  int  getExternalDriftNumber()                    const;
  int  getRankFext(int il)                         const;
  int  getDriftEquationNumber()                    const;
  bool isDriftFiltered(unsigned int il)            const;
  int  getDriftMaxIRFOrder(void)                   const;
  bool isDriftDefined(const VectorInt &powers, int rank_fex = 0) const;
  bool isDriftDifferentDefined(const VectorInt &powers, int rank_fex = -1) const;
  bool isDriftSampleDefined(const Db *db,
                            int ib,
                            int nech,
                            const VectorInt &nbgh,
                            const ELoc &loctype) const;

  void setDriftFiltered(int il, bool filtered)                 ;
  VectorVectorDouble getDrifts(const Db* db, bool useSel=true) ;
  void setBetaHat(const VectorDouble &betaHat);

  double evalDrift(const Db* db,
                   int iech,
                   int il,
                   const ECalcMember& member = ECalcMember::fromKey("LHS")) const;
  double evalDriftValue(const Db *db,
                        int iech,
                        int ivar,
                        int ib,
                        const ECalcMember &member = ECalcMember::fromKey("LHS")) const;
  VectorDouble evalDriftBySample(const Db* db,
                                 int iech,
                                 const ECalcMember& member = ECalcMember::fromKey("LHS")) const;
  void evalDriftBySampleInPlace(const Db *db,
                                int iech,
                                const ECalcMember &member,
                                VectorDouble &drftab) const;


  double evalDriftVarCoef(const Db *db,
                          int iech,
                          int ivar,
                          const VectorDouble &coeffs) const;
  VectorDouble evalDriftVarCoefs(const Db *db,
                                 const VectorDouble &coeffs,
                                 int ivar = 0,
                                 bool useSel = false) const;
  /////////////////////////////////////////////////

  ////////////////////////////////////////////////
  /// TODO : to be removed (encapsulation of Context)
  const CovContext& getContext() const { return _ctxt; }
  const ASpace*     getASpace() const { return _ctxt.getASpace(); }
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

  const EModelProperty& getCovMode() const;
  Model* duplicate() const;
  Model* createReduce(const VectorInt& validVars) const;

  int getVariableNumber() const
  {
    // TODO/ the strange next line have been commented out.
    // There should be either validated or suppressed
    //if (isFlagGradient())
    //      return 3; // This strange number of variables is linked to the Gradient calculation
    //    else
    // However, note used for Gradient (Functional type) in Potential
    int nvar = _cova->getNVariables();
    if (nvar <= 0)
      nvar = _ctxt.getNVar();
    return nvar;
  }

  int hasExternalCov() const;

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
  int buildVmapOnDbGrid(DbGrid *dbgrid, const NamingConvention &namconv = NamingConvention("VMAP")) const;
  int stabilize(double percent, bool verbose = false);
  int standardize(bool verbose = false);

  double gofToVario(const Vario* vario, bool verbose = true);
  static void gofDisplay(double gof,
                         bool byValue                   = true,
                         const VectorDouble& thresholds = {2., 5., 10., 100});
  static VectorECov initCovList(const VectorInt & covranks);

  bool isValid() const;

  VectorDouble sample(const VectorDouble &h,
                      const VectorDouble &codir = VectorDouble(),
                      int ivar = 0,
                      int jvar = 0,
                      const CovCalcMode* mode = nullptr,
                      const CovInternal* covint = nullptr);
  double evaluateOneIncr(double hh,
                         const VectorDouble &codir = VectorDouble(),
                         int ivar = 0,
                         int jvar = 0,
                         const CovCalcMode *mode = nullptr);
  void evaluateMatInPlace(const CovInternal *covint,
                          const VectorDouble &d1,
                          MatrixSquareGeneral &covtab,
                          bool flag_init = false,
                          double weight = 1.,
                          const CovCalcMode *mode = nullptr);
  double evaluateOneGeneric(const CovInternal *covint,
                            const VectorDouble &d1 = VectorDouble(),
                            double weight = 1.,
                            const CovCalcMode *mode = nullptr);
  VectorDouble evaluateFromDb(Db *db,
                              int ivar = 0,
                              int jvar = 0,
                              const CovCalcMode *mode = nullptr);
  double calculateStdev(Db *db1,
                        int iech1,
                        Db *db2,
                        int iech2,
                        bool verbose = false,
                        double factor = 1.,
                        const CovCalcMode *mode = nullptr);

  double computeLogLikelihood(const Db* db, bool verbose = false);

protected:
  /// Interface to ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "Model"; }

  const CovAnisoList* _castInCovAnisoListConst(int icov = -1) const;
  CovAnisoList*       _castInCovAnisoList(int icov = -1);

private:
  void _clear();
  void _create();
  void _copyCovContext();

  MatrixSquareSymmetric _dummy;
};

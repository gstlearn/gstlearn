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
#include "geoslib_define.h"

#include "Enum/ECalcMember.hpp"
#include "Enum/ECov.hpp"
#include "Enum/EDrift.hpp"
#include "Enum/EConsElem.hpp"
#include "Enum/EModelProperty.hpp"

#include "Covariances/CovContext.hpp"
#include "Covariances/ACovAnisoList.hpp"
#include "Covariances/CovLMGradient.hpp"

#include "Drifts/DriftList.hpp"

#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"
#include "Model/Constraints.hpp"
#include "Model/CovParamId.hpp"
#include "Covariances/CovAniso.hpp"

#include "Matrix/MatrixRectangular.hpp"

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
class ADriftElem;

typedef std::vector<ECov> VectorECov;

/// TODO : Create AModel which inherits from ACov ?
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
  void   delCova(int rank);
  void   delAllCovas();
  void   setDriftList(const DriftList* driftlist);
  void   setDriftIRF(int order = 0, int nfex = 0);
  void   addDrift(const ADriftElem* drift);
  void   setDrifts(const VectorString& driftSymbols);
  void   delDrift(int rank);
  void   delAllDrifts();
  int    addNoStat(const ANoStat* anostat);
  int    setAnam(const AAnam* anam, const VectorInt& strcnt = VectorInt());
  bool   isFlagGradient() const;
  bool   isFlagGradientNumerical() const;
  bool   isFlagGradientFunctional() const;
  bool   isFlagLinked() const;
  CovAniso extractCova(int icov) const { return _covaList->extractCova(icov); }
  void   switchToGradient();

  ////////////////////////////////////////////////
  /// TODO : to be removed (encapsulation of ACovAnisoList)
  const ACovAnisoList* getCovAnisoList() const { return _covaList; }
  ACovAnisoList* getCovAnisoList() { return _covaList; } // Needed for dynamic cast
  const CovAniso* getCova(unsigned int icov) const;
  CovAniso* getCova(unsigned int icov);
  int getCovaNumber() const;
  const ECov& getCovaType(int icov) const;
  const MatrixSquareSymmetric& getSill(int icov) const;
  double getSill(int icov, int ivar, int jvar) const;
  double getParam(int icov) const;
  bool isCovaFiltered(int icov) const;
  String getCovName(int icov) const;
  int getGradParamNumber(int icov) const;
  double getTotalSill(int ivar, int jvar) const;
  double getBallRadius() const;
  double getMaximumDistance() const { return _covaList->getMaximumDistance(); }
  int    getMinOrder() const { return _covaList->getMinOrder(); }
  bool   hasAnam() const { return _covaList->hasAnam(); }
  const AAnam* getAnam() { return _covaList->getAnam(); }
  void normalize(double sill) { _covaList->normalize(sill); }
  bool hasNugget() const { return _covaList->hasNugget(); }

  double eval0(int ivar,
               int jvar,
               const CovCalcMode& mode = CovCalcMode()) const
  {
    return _covaList->eval0(ivar, jvar, mode);
  }
  MatrixSquareGeneral eval0Nvar(const CovCalcMode& mode = CovCalcMode()) const
  {
    return _covaList->eval0Nvar(mode);
  }
  double eval(int ivar,
              int jvar,
              const SpacePoint& p1,
              const SpacePoint& p2,
              const CovCalcMode& mode = CovCalcMode()) const
  {
    return _covaList->eval(ivar, jvar, p1, p2, mode);
  }
  MatrixSquareGeneral evalNvarIpas(double step,
                                   const VectorDouble& dir = VectorDouble(),
                                   const VectorDouble& center = VectorDouble(),
                                   const CovCalcMode& mode = CovCalcMode()) const
  {
    return _covaList->evalNvarIpas(step, dir, center, mode);
  }
  MatrixSquareGeneral evalNvarIpas(const VectorDouble& dincr,
                                   const CovCalcMode& mode = CovCalcMode()) const
  {
    return _covaList->evalNvarIpas(dincr, mode);
  }
  VectorDouble evalIvarNpas(int ivar,
                            int jvar,
                            const VectorDouble& vec_step,
                            const VectorDouble& dir = VectorDouble(),
                            const VectorDouble& center = VectorDouble(),
                            const CovCalcMode& mode = CovCalcMode()) const
  {
    return _covaList->evalIvarNpas(ivar, jvar, vec_step, dir, center, mode);
  }
  double evalIvarIpas(int ivar,
                      int jvar,
                      double step,
                      const VectorDouble& dir = VectorDouble(),
                      const VectorDouble& center = VectorDouble(),
                      const CovCalcMode& mode = CovCalcMode()) const
  {
    return _covaList->evalIvarIpas(ivar, jvar, step, dir, center, mode);
  }
  double evalCvv(const VectorDouble& ext,
                 const VectorInt& ndisc,
                 const VectorDouble& angles = VectorDouble(),
                 int ivar = 0,
                 int jvar = 0,
                 const CovCalcMode& mode = CovCalcMode()) const
  {
    return _covaList->evalCvv(ext, ndisc, angles, ivar, jvar, mode);
  }
  double evalCvvShift(const VectorDouble& ext,
                      const VectorInt& ndisc,
                      const VectorDouble& shift,
                      const VectorDouble& angles = VectorDouble(),
                      int ivar = 0,
                      int jvar = 0,
                      const CovCalcMode& mode = CovCalcMode()) const
  {
    return _covaList->evalCvvShift(ext, ndisc, shift, angles, ivar, jvar, mode);
  }
  MatrixSquareGeneral evalCvvM(const VectorDouble& ext,
                               const VectorInt& ndisc,
                               const VectorDouble& angles = VectorDouble(),
                               const CovCalcMode& mode = CovCalcMode())
  {
    return _covaList->evalCvvM(ext, ndisc, angles, mode);
  }
  double evalCxv(const SpacePoint& p1,
                 const VectorDouble& ext,
                 const VectorInt& ndisc,
                 const VectorDouble& angles = VectorDouble(),
                 const VectorDouble& x0 = VectorDouble(),
                 int ivar = 0,
                 int jvar = 0,
                 const CovCalcMode& mode = CovCalcMode())
  {
    return _covaList->evalCxv(p1, ext, ndisc, angles, x0, ivar, jvar, mode);
  }
  MatrixSquareGeneral evalCxvM(const SpacePoint& p1,
                               const VectorDouble& ext,
                               const VectorInt& ndisc,
                               const VectorDouble& angles = VectorDouble(),
                               const VectorDouble& x0 = VectorDouble(),
                               const CovCalcMode& mode = CovCalcMode())
  {
    return _covaList->evalCxvM(p1, ext, ndisc, angles, x0, mode);
  }
  VectorDouble evalPointToDb(const SpacePoint& p1,
                             const Db* db2,
                             int ivar = 0,
                             int jvar = 0,
                             bool useSel = true,
                             const CovCalcMode& mode = CovCalcMode())
  {
    return _covaList->evalPointToDb(p1, db2, ivar, jvar, useSel, mode);
  }
  double evalAverageDbToDb(const Db* db1,
                           const Db* db2,
                           int ivar = 0,
                           int jvar = 0,
                           const CovCalcMode& mode = CovCalcMode())
  {
    return _covaList->evalAverageDbToDb(db1, db2, ivar, jvar, mode);
  }
  double evalAveragePointToDb(const SpacePoint& p1,
                              const Db* db2,
                              int ivar = 0,
                              int jvar = 0,
                              const CovCalcMode& mode = CovCalcMode())
  {
    return _covaList->evalAveragePointToDb(p1, db2, ivar, jvar, mode);
  }
  MatrixRectangular evalCovMatrix(const Db* db1,
                                  const Db* db2 = nullptr,
                                  int ivar = 0,
                                  int jvar = 0,
                                  const CovCalcMode& mode = CovCalcMode())
  {
    return _covaList->evalCovMatrix(db1, db2, ivar, jvar, mode);
  }
  double extensionVariance(const Db* db,
                           const VectorDouble& ext,
                           const VectorInt& ndisc,
                           const VectorDouble& angles = VectorDouble(),
                           const VectorDouble& x0 = VectorDouble(),
                           int ivar = 0,
                           int jvar = 0)
  {
    return _covaList->extensionVariance(db, ext, ndisc, angles, x0, ivar, jvar);
  }
  double samplingDensityVariance(const Db* db,
                                 const VectorDouble& ext,
                                 const VectorInt& ndisc,
                                 const VectorDouble& angles = VectorDouble(),
                                 const VectorDouble& x0 = VectorDouble(),
                                 int ivar = 0,
                                 int jvar = 0) const
  {
    return _covaList->samplingDensityVariance(db, ext, ndisc, angles, x0, ivar, jvar);
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
    return _covaList->specificVolume(db, mean, ext, ndisc, angles, x0, ivar, jvar);
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
    return _covaList->coefficientOfVariation(db, volume, mean, ext, ndisc, angles, x0, ivar, jvar);
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
    return _covaList->specificVolumeFromCoV(db, cov, mean, ext, ndisc, angles, x0, ivar, jvar);
  }
  void evalZAndGradients(const SpacePoint& p1,
                         const SpacePoint& p2,
                         double& covVal,
                         VectorDouble& covGp,
                         VectorDouble& covGG,
                         const CovCalcMode& mode = CovCalcMode(),
                         bool flagGrad = false) const
  {
    CovLMGradient* covgrad = dynamic_cast<CovLMGradient *>(_covaList);
    if (covgrad != nullptr)
      covgrad->evalZAndGradients(p1, p2, covVal, covGp, covGG, mode, flagGrad);
  }
  void evalZAndGradients(const VectorDouble& vec,
                         double& covVal,
                         VectorDouble& covGp,
                         VectorDouble& covGG,
                         const CovCalcMode& mode = CovCalcMode(),
                         bool flagGrad = false) const
  {
    CovLMGradient* covgrad = dynamic_cast<CovLMGradient *>(_covaList);
    if (covgrad != nullptr)
      covgrad->evalZAndGradients(vec, covVal, covGp, covGG, mode, flagGrad);
  }

  void setSill(int icov, int ivar, int jvar, double value);
  void setCovaFiltered(int icov, bool filtered);
  void setActiveFactor(int iclass) { _covaList->setActiveFactor(iclass); }
  int  getActiveFactor() const { return _covaList->getActiveFactor(); }
  int  getAnamNClass() const { return _covaList->getAnamNClass(); }
  /////////////////////////////////////////////////

  ////////////////////////////////////////////////
  /// TODO : to be removed (encapsulation of DriftList)
  const DriftList* getDriftList()                  const;
  const ADriftElem* getDrift(int il)               const;
  ADriftElem* getDrift(int il)                          ;
  int getDriftNumber()                             const;
  int getExternalDriftNumber()                     const;
  const EDrift& getDriftType(int il)               const;
  int getRankFext(int il)                          const;
  const VectorDouble& getCoefDrifts()              const;
  double getCoefDrift(int ivar, int il, int ib)    const;
  int getDriftEquationNumber()                     const;
  bool isDriftFiltered(unsigned int il)            const;
  bool isDriftDefined(const EDrift& type0)         const;
  bool isDriftDifferentDefined(const EDrift& type0) const;
  int getMaximumOrder(void) const { return _driftList->getMaximumOrder(); }

  void setCoefDrift(int ivar, int il, int ib, double coeff)    ;
  void setCoefDriftByRank(int rank, double coeff)              ;
  void setDriftFiltered(int il, bool filtered)                 ;
  VectorDouble getDrift(const Db* db, int ib, bool useSel=true);
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
  double _evalDriftCoef(const Db* db,
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
  void setMean(int ivar, double mean);
  void setCovar0s(const VectorDouble& covar0);
  void setCovar0(int ivar, int jvar, double covar0);
  void setField(double field);
  /////////////////////////////////////////////////

  /////////////////////////////////////////////////
  /// Shortcut for Non-stationary
  int  isNoStat() const;
  const ANoStat* getNoStat() const { return _noStat; }
  int  getNoStatElemNumber() const;
  int  addNoStatElem(int igrf, int icov, const EConsElem& type, int iv1, int iv2);
  int  addNoStatElems(const VectorString& codes);
  int  getNoStatElemIcov(int ipar);
  const EConsElem& getNoStatElemType(int ipar);
  CovParamId getCovParamId(int ipar) const;
  ////////////////////////////////////////////////

  const EModelProperty& getCovMode() const;
  Model* duplicate() const;

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

  void covMatrix(VectorDouble& covmat,
                 Db *db1,
                 Db *db2 = nullptr,
                 int ivar0 = 0,
                 int jvar0 = 0,
                 int flag_norm = 0,
                 int flag_cov = 1);
  VectorDouble sample(double hmax,
                      int nh = 100,
                      int ivar = 0,
                      int jvar = 0,
                      VectorDouble codir = VectorDouble(),
                      int norder = 0,
                      bool asCov = false,
                      bool addZero = false);

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

private:
  void _clear();
  void _create();
  void _copyCovContext();

private:
  /// TODO : Transform to ACov in place of ACovAnisoList (to be put in AModel)
  ACovAnisoList* _covaList;     /* Series of Covariance structures */
  DriftList*     _driftList;    /* Series of Drift functions */
  ANoStat*       _noStat;       /* Description of Non-stationary Model */
  CovContext     _ctxt;         /* Context */
};

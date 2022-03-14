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

#include "Covariances/ECov.hpp"
#include "Covariances/ECalcMember.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/ACovAnisoList.hpp"

#include "Drifts/EDrift.hpp"
#include "Drifts/DriftList.hpp"

#include "Model/EModelProperty.hpp"
#include "Model/EConsElem.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"
#include "Model/Constraints.hpp"
#include "Model/CovParamId.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/IClonable.hpp"

class Model;
class Db;
class Drift;
class CovInternal;
class MatrixSquareSymmetric;
class CovCalcMode;
class Vario;
class CovAniso;
class ANoStat;
class ADriftElem;

class GSTLEARN_EXPORT Model : public AStringable, public ASerializable, public IClonable
{
public:
  Model(const CovContext& ctxt = CovContext());
  Model(const Model &m);
  Model& operator= (const Model &m);
  virtual ~Model();

public:
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;
  virtual IClonable* clone() const override { return new Model(*this); }

  int resetFromDb(const Db* db);

  int dumpToNF(const String& neutralFilename, bool verbose = false) const;
  static Model* create(const CovContext& ctxt = CovContext());
  static Model* createFromDb(const Db* db);
  static Model* createFromNF(const String& neutralFilename, bool verbose = false);

  void   setCovList(const ACovAnisoList* covalist);
  void   addCov(const CovAniso* cov);
  void   delCova(int rank);
  void   delAllCovas();
  void   setDriftList(const DriftList* driftlist);
  void   addDrift(const ADriftElem* drift);
  void   addDrift(const VectorString& driftSymbols);
  void   delDrift(int rank);
  void   delAllDrifts();
  int    addNoStat(const ANoStat* anostat);
  bool   isFlagGradient() const;
  bool   isFlagGradientNumerical() const;
  bool   isFlagGradientFunctional() const;
  bool   isFlagLinked() const;

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
  void setSill(int icov, int ivar, int jvar, double value);
  void setCovaFiltered(int icov, bool filtered);
  double getMaximumDistance() const { return _covaList->getMaximumDistance(); }
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
                   const ECalcMember& member = ECalcMember::LHS) const;
  VectorDouble evalDriftVec(const Db* db,
                            int iech,
                            const ECalcMember& member = ECalcMember::LHS) const;
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
    if (isFlagGradient())
      return 3; // This strange number of variables is linked to the Gradient calculation
    else
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
                      bool addZero = false);

  // TODO : Remove Model::fit duplicate declaration
  int fitFromCovIndices(Vario *vario,
                        const VectorInt& types,
                        bool verbose = false,
                        Option_AutoFit mauto = Option_AutoFit(),
                        const Constraints& constraints = Constraints(),
                        Option_VarioFit optvar = Option_VarioFit());
  int fit(Vario *vario,
          const std::vector<ECov>& types,
          bool verbose = false,
          Option_AutoFit mauto = Option_AutoFit(),
          const Constraints& constraints = Constraints(),
          Option_VarioFit optvar = Option_VarioFit());

  double gofToVario(const Vario* vario);

protected:
  virtual int _deserialize(std::istream& is, bool verbose = false) override;
  virtual int _serialize(std::ostream& os, bool verbose = false) const override;

private:
  void _clear();
  void _create();
  void _copyCovContext();

private:
  ACovAnisoList* _covaList;     /* Series of Covariance structures */
  DriftList*     _driftList;    /* Series of Drift functions */
  ANoStat*       _noStat;       /* Description of Non-stationary Model */
  CovContext     _ctxt;         /* Context */
};

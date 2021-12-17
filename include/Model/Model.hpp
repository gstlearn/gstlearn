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

//Enums
#include "Covariances/ECov.hpp"
#include "Model/EModelProperty.hpp"
#include "Covariances/ECalcMember.hpp"
#include "Model/EConsElem.hpp"
#include "Drifts/EDrift.hpp"

#include "Model/ModTrans.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"
#include "Model/Constraints.hpp"
#include "Model/CovParamId.hpp"
#include "Covariances/CovContext.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/IClonable.hpp"

class Model;
class Db;
class Drift;
class CovInternal;
class MatrixSquareSymmetric;
class CovCalcMode;
class ACovAnisoList;
class Vario;
class CovAniso;
class ANoStat;
class DriftList;
class ADriftElem;

class GSTLEARN_EXPORT Model : public AStringable, public ASerializable, public IClonable
{
public:
  Model(const CovContext& ctxt, bool flagGradient = false, bool flagLinked = false);
  Model(const Db *db, bool flagGradient = false, bool flagLinked = false);
  Model(const String& neutralFileName, bool verbose = false);
  Model(const Model &m);
  Model& operator= (const Model &m);
  virtual ~Model();

public:
  virtual String toString(int level = 0) const override;
  int deSerialize(const String& filename, bool verbose = false) override;
  int serialize(const String& filename, bool verbose = false) const override;
  virtual IClonable* clone() const override { return new Model(*this); }

  /// TODO : to be converted as internal member
  const CovContext& getContext() const { return _ctxt; }

  void   addCova(const CovAniso* cov);
  void   delCova(int rank);
  void   delAllCovas();
  void   addDrift(const ADriftElem* drift);
  void   addDrift(const VectorString& driftSymbols);
  void   delDrift(int rank);
  void   delAllDrifts();
  int    addNoStat(const ANoStat* anostat);
  bool   isFlagGradient() const { return _flagGradient; }
  bool   isFlagLinked() const { return _flagLinked; }

  ////////////////////////////////////////////////
  /// TODO : to be removed (encapsulation of ACovAnisoList)
  const ACovAnisoList* getCovAnisoList() const;
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

  void setSill(int icov, int ivar, int jvar, double value);
  void setCovaFiltered(int icov, bool filtered);
  /////////////////////////////////////////////////

  ////////////////////////////////////////////////
  /// TODO : to be removed (encapsulation of DriftList)
  const DriftList* getDriftList()                  const;
  const ADriftElem* getDrift(int il)               const;
  ADriftElem* getDrift(int il)                          ;
  int getDriftNumber()                             const;
  const EDrift& getDriftType(int il)               const;
  int getRankFext(int il)                          const;
  const VectorDouble& getCoefDrifts()              const;
  double getCoefDrift(int ivar, int il, int ib)    const;
  int getDriftEquationNumber()                     const;
  bool isDriftFiltered(unsigned int il)            const;

  void setCoefDrift(int ivar, int il, int ib, double coeff)    ;
  void setCoefDriftByRank(int rank, double coeff)              ;
  void setDriftFiltered(int il, bool filtered)                 ;
  VectorDouble getDrift(const Db* db, int ib, bool useSel=true);
  VectorVectorDouble getDrifts(const Db* db, bool useSel=true) ;

  double evaluateDrift(const Db* db, int iech, int il,
                       const ECalcMember& member = ECalcMember::LHS) const;
  /////////////////////////////////////////////////

  ////////////////////////////////////////////////
  /// TODO : to be removed (encapsulation of Context)
  const VectorDouble& getMeans() const { return _ctxt.getMean(); }
  double getMean(int ivar) const { return _ctxt.getMean(ivar); }
  const VectorDouble& getCovar0s() const { return _ctxt.getCovar0(); }
  double getCovar0(int ivar, int jvar) const { return _ctxt.getCovar0(ivar,jvar); }
  void setMeans(const VectorDouble& mean) { _ctxt.setMean(mean); }
  void setMean(int ivar, double mean) { _ctxt.setMean(ivar, mean); }
  void setCovar0s(const VectorDouble& covar0) { _ctxt.setCovar0(covar0); }
  void setCovar0(int ivar, int jvar, double covar0) { _ctxt.setCovar0(ivar,jvar,covar0); }
  /////////////////////////////////////////////////

  /////////////////////////////////////////////////
  /// Shortcut for Non-stationary
  int isNoStat() const;
  const ANoStat* getNoStat() const { return _noStat; }
  int  getNoStatElemNumber() const;
  int  addNoStatElem(int igrf, int icov, const EConsElem& type, int iv1, int iv2);
  int  addNoStatElems(const VectorString& codes);
  int  getNoStatElemIcov(int ipar);
  const EConsElem& getNoStatElemType(int ipar);
  CovParamId getCovParamId(int ipar) const;
  ////////////////////////////////////////////////

  double getField() const                       { return _ctxt.getField(); }
  ModTrans& getModTrans()                       { return _modTrans; }
  int getDimensionNumber() const                { return _ctxt.getNDim(); }
  void setField(double field)                   { _ctxt.setField(field); }
  const EModelProperty& getModTransMode() const { return _modTrans.getModTransMode(); }
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
                      int norder = 0);

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

private:
  void _create(bool flagGradient, bool flagLinked);
  void _destroy();

private:
  bool           _flagGradient;
  bool           _flagLinked;
  ACovAnisoList* _covaList;     /* Series of Covariance structures */
  DriftList*     _driftList;    /* Series of Drift functions */
  ModTrans       _modTrans;     /* Covariance Transformation */
  ANoStat*       _noStat;       /* Description of Non-stationary Model */
  CovContext     _ctxt;         /* Context */

public:
  void (*generic_cov_function)(CovInternal *cov_nostat,
                               Model *model,
                               const CovCalcMode& mode,
                               int flag_init,
                               double weight,
                               VectorDouble d1,
                               double *covtab);
};

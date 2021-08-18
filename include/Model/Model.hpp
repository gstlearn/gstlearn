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

#include "Db/Db.hpp"
#include "Model/ModTrans.hpp"
#include "Drifts/ADriftList.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"
#include "Model/Constraints.hpp"
#include "Model/ANoStat.hpp"

#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Covariances/ACovAnisoList.hpp"
#include "Covariances/CovContext.hpp"
#include "Space/SpaceRN.hpp"

class Model;
class Cova;
class Drift;
class ModTrans;
class CovInternal;
class MatrixCSSym;
class CovCalcMode;
class Vario;
class CovAniso;

class Model : public AStringable, public ASerializable
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

  /// TODO : to be converted as internal member
  const CovContext& getContext() const { return _ctxt; }

  void   addCova(const CovAniso* cov);
  void   delCova(int rank);
  void   delAllCovas();
  void   addDrift(const ADriftElem* drift);
  void   addDrift(const VectorString& driftSymbols);
  void   delDrift(int rank);
  void   delAllDrifts();
  int    addNoStat(ANoStat* anostat);
  bool   isFlagGradient() const { return _flagGradient; }
  bool   isFlagLinked() const { return _flagLinked; }

  ////////////////////////////////////////////////
  /// TODO : to be removed (encapsulation of ACovAnisoList)
  const ACovAnisoList* getCovAnisoList()           const { return _covaList; }
  const CovAniso* getCova(unsigned int icov)       const { return _covaList->getCova(icov); }
  CovAniso* getCova(unsigned int icov)                   { return _covaList->getCova(icov); }
  int getCovaNumber()                              const { return _covaList->getCovNumber(); }
  ENUM_COVS getCovaType(int icov)                  const { return _covaList->getType(icov); }
  const MatrixCSGeneral& getSill(int icov)         const { return _covaList->getSill(icov); }
  double getSill(int icov, int ivar, int jvar)     const { return _covaList->getSill(icov, ivar, jvar); }
  double getParam(int icov)                        const { return _covaList->getParam(icov); }
  bool isCovaFiltered(int icov)                    const { return _covaList->isFiltered(icov); }
  String getCovName(int icov)                      const { return _covaList->getCovName(icov); }
  int getGradParamNumber(int icov)                 const { return _covaList->getGradParamNumber(icov); }

  void setSill(int icov, int ivar, int jvar, double value) { _covaList->setSill(icov, ivar, jvar, value); }
  void setCovaFiltered(int icov, bool filtered)            { _covaList->setFiltered(icov, filtered); }
  /////////////////////////////////////////////////

  ////////////////////////////////////////////////
  /// TODO : to be removed (encapsulation of ADriftList)
  /////////////////////////////////////////////////
  const ADriftList* getDriftList()                 const { return _driftList; }
  const ADriftElem* getDrift(int il)               const { return _driftList->getDrift(il); }
  ADriftElem* getDrift(int il)                           { return _driftList->getDrift(il); }
  int getDriftNumber()                             const { return _driftList->getDriftNumber(); }
  ENUM_DRIFTS getDriftType(int il)                 const { return _driftList->getType(il); }
  int getRankFext(int il)                          const { return _driftList->getRankFex(il); }
  const VectorDouble& getCoefDrift()               const { return _driftList->getCoefDrift(); }
  double getCoefDrift(int ivar, int il, int ib)    const { return _driftList->getCoefDrift(ivar, il, ib); }
  int getDriftEquationNumber()                     const { return _driftList->getDriftEquationNumber(); }
  bool isDriftFiltered(unsigned int il)            const { return _driftList->isFiltered(il); }

  void setCoefDrift(int ivar, int il, int ib, double coeff) { _driftList->setCoefDrift(ivar, il, ib, coeff); }
  void setCoefDrift(int rank, double coeff)                 { _driftList->setCoefDrift(rank, coeff); }
  void setDriftFiltered(int il, bool filtered)              { _driftList->setFiltered(il, filtered); }

  double evaluateDrift(const Db* db, int iech, int il, int member = MEMBER_LHS) const;
  /////////////////////////////////////////////////

  ////////////////////////////////////////////////
  /// TODO : to be removed (encapsulation of Context)
  const VectorDouble& getMean() const { return _ctxt.getMean(); }
  double getMean(int ivar) const { return _ctxt.getMean(ivar); }
  const VectorDouble& getCovar0() const { return _ctxt.getCovar0(); }
  double getCovar0(int ivar, int jvar) const { return _ctxt.getCovar0(ivar,jvar); }
  void setMean(const VectorDouble& mean) { _ctxt.setMean(mean); }
  void setMean(int ivar, double mean) { _ctxt.setMean(ivar, mean); }
  void setCovar0(const VectorDouble& covar0) { _ctxt.setCovar0(covar0); }
  void setCovar0(int ivar, int jvar, double covar0) { _ctxt.setCovar0(ivar,jvar,covar0); }
  /////////////////////////////////////////////////

  /////////////////////////////////////////////////
  /// Shortcut for Non-stationary
  int isNoStat() const;
  const ANoStat* getNoStat() { return _noStat; }
  int  getNoStatElemNumber() const;
  void addNoStatElem(int igrf, int icov, ENUM_CONS type, int iv1, int iv2);
  void addNoStatElems(const VectorString& codes);
  int  getNoStatElemIcov(int ipar);
  ENUM_CONS getNoStatElemType(int ipar);
  ConsItem getConsItem(int ipar) const;
  ////////////////////////////////////////////////

  double getField() const            { return _ctxt.getField(); }
  ModTrans& getModTrans()            { return _modTrans; }
  int getDimensionNumber() const     { return _ctxt.getNDim(); }
  void setField(double field)        { _ctxt.setField(field); }
  int getModTransMode() const        { return _modTrans.getModTransMode(); }

  int getVariableNumber() const
  {
    if (isFlagGradient())
      return 3; // This strange number of variables is linked to the Gradient calculation
    else
      return _ctxt.getNVar();
  }

  int hasExternalCov() const;

  VectorDouble sample(double hmax,
                    int nh = 100,
                    int ivar = 0,
                    int jvar = 0,
                    VectorDouble codir = VectorDouble(),
                    int norder = 0);

  int fit(Vario *vario,
          const std::vector<int>& types,
          bool verbose = false,
          Option_AutoFit mauto = Option_AutoFit(),
          const Constraints& constraints = Constraints(),
          Option_VarioFit optvar = Option_VarioFit());
  int fit(Vario *vario,
          const std::vector<ENUM_COVS>& types,
          bool verbose = false,
          Option_AutoFit mauto = Option_AutoFit(),
          const Constraints& constraints = Constraints(),
          Option_VarioFit optvar = Option_VarioFit());

private:
  void _create(bool flagGradient, bool flagLinked);
  void _destroy();

private:
  bool           _flagGradient;
  bool           _flagLinked;
  ACovAnisoList* _covaList;     /* Series of Covariance structures */
  ADriftList*    _driftList;    /* Series of Drift functions */
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

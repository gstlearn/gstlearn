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
#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Db/Db.hpp"
#include "Anamorphosis/AAnam.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Anamorphosis/AnamDiscreteDD.hpp"
#include "Anamorphosis/AnamDiscreteIR.hpp"
#include "Polynomials/Hermite.hpp"
#include "Stats/Selectivity.hpp"
#include "Stats/ESelectivity.hpp"

#include <math.h>

#define NCOLS 7
#define Z_CUT 0
#define T_EST 1
#define Q_EST 2
#define B_EST 3
#define M_EST 4
#define T_STD 5
#define Q_STD 6

Selectivity::Selectivity(int ncut)
    : AStringable(),
      _Zcut(ncut),
      _stats(ncut, NCOLS),
      _zmax(TEST),
      _proba(TEST),
      _flagTonnageCorrect(false),
      _numberQT(),
      _rankQT()
{
  ut_vector_fill(_Zcut, TEST);
  _stats.setColumnNames(_getAllNames());
  _stats.fill(TEST);
}

Selectivity::Selectivity(const VectorDouble &zcuts,
                         double zmax,
                         double proba,
                         bool flag_tonnage_correct)
    : AStringable(),
      _Zcut(zcuts),
      _stats(),
      _zmax(zmax),
      _proba(proba),
      _flagTonnageCorrect(flag_tonnage_correct),
      _numberQT(),
      _rankQT()
{
  _stats.init(getNCuts(), NCOLS);
  _stats.setColumnNames(_getAllNames());
  _stats.fill(TEST);
}

Selectivity::Selectivity(const Selectivity &m)
    : AStringable(m),
      _Zcut(m._Zcut),
      _stats(m._stats),
      _zmax(m._zmax),
      _proba(m._proba),
      _flagTonnageCorrect(m._flagTonnageCorrect),
      _numberQT(m._numberQT),
      _rankQT(m._rankQT)
{
}

Selectivity& Selectivity::operator=(const Selectivity &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _Zcut  = m._Zcut;
    _stats = m._stats;
    _zmax  = m._zmax;
    _proba = m._proba;
    _flagTonnageCorrect = m._flagTonnageCorrect;
    _numberQT = m._numberQT;
    _rankQT = m._rankQT;
  }
  return *this;
}

Selectivity::~Selectivity()
{
}

Selectivity* Selectivity::create(int ncut)
{
  Selectivity* selectivity = new Selectivity(ncut);
  return selectivity;
}

Selectivity* Selectivity::create(const VectorDouble& zcut)
{
  Selectivity* selectivity = new Selectivity(zcut);
  return selectivity;
}

Selectivity* Selectivity::createByCodes(const std::vector<ESelectivity>& codes,
                                        const VectorDouble& zcuts,
                                        bool flag_est,
                                        bool flag_std,
                                        double proba,
                                        bool verbose)
{
  Selectivity* selectivity = new Selectivity(zcuts);
  selectivity->defineRecoveries(codes, flag_est, flag_std, proba, verbose);
  return selectivity;
}

Selectivity* Selectivity::createByKeys(const VectorString& keys,
                                       const VectorDouble& zcuts,
                                       bool flag_est,
                                       bool flag_std,
                                       double zmax,
                                       bool flag_tonnage_correct,
                                       double proba,
                                       bool verbose)
{
  std::vector<ESelectivity> codes;

  for (int i = 0; i < (int) keys.size(); i++)
  {
    ESelectivity code = ESelectivity::fromKey(keys[i]);
    if (code == ESelectivity::UNKNOWN) continue;
    codes.push_back(code);
  }

  Selectivity* selectivity = new Selectivity(zcuts, zmax, flag_tonnage_correct);
  selectivity->defineRecoveries(codes, flag_est, flag_std, proba, verbose);
  return selectivity;
}

void Selectivity::resetCuts(const VectorDouble& zcuts)
{
  _Zcut = zcuts;
  _stats.init(getNCuts(), NCOLS);
  _stats.setColumnNames(_getAllNames());
  _stats.fill(TEST);
}

int Selectivity::calculateFromDb(const Db* db, bool autoCuts)
{
  if (getNCuts() <= 0)
  {
    messerr("You must define 'zcuts'");
    return 1;
  }
  if (db == nullptr)
  {
    messerr("You must provide a valid 'Db'");
    return 1;
  }
  if (db->getVariableNumber() != 1)
  {
    messerr("The 'Db' must contain a SINGLE variable");
    return 1;
  }

  // Extract the array of data and weights
  VectorDouble tab = db->getColumnByLocator(ELoc::Z, 0, true);
  VectorDouble wtab;
  if (db->hasWeight())
    wtab = db->getColumnByLocator(ELoc::W, 0, true);

  return calculateFromArray(tab, wtab, autoCuts);
}

int Selectivity::calculateFromArray(const VectorDouble &tab,
                                    const VectorDouble &weights,
                                    bool autoCuts)
{
  if (getNCuts() <= 0)
  {
    messerr("You must define 'cuts'");
    return 1;
  }
  if (tab.empty())
  {
    messerr("You must provide a valid 'tab'");
    return 1;
  }

  int nech = (int) tab.size();
  VectorDouble wtab = weights;
  if (wtab.empty())
  {
    wtab.resize(nech,1);
  }
  else
  {
    if (nech != (int) wtab.size())
    {
      messerr("Arguments 'tab' and 'weights' should have same dimension");
      return 1;
    }
  }

  // Allocate the returned structure

  int ncut = getNCuts();

  // Define the cutoffs (optional)

  if (autoCuts)
    _defineAutomaticCutoffs(tab);

  // Perform calculations

  double tonref = ut_vector_cumul(wtab);
  for (int icut = 0; icut < ncut; icut++)
  {
    double coupure = getZcut(icut);

    double tonnage = 0.;
    double metal = 0.;
    for (int iech = 0; iech < nech; iech++)
    {
      double x = tab[iech];
      double w = wtab[iech];
      double indic = (x >= coupure);

      tonnage += w * indic;
      metal   += x * w * indic;
    }

    tonnage /= tonref;
    metal   /= tonref;
    double benefit = metal - tonnage * coupure;
    double grade;
    if (tonnage <= 0.)
    {
      grade = TEST;
    }
    else
    {
      grade = metal / tonnage;
    }
    setZcut(icut, coupure);
    setTest(icut, tonnage);
    setQest(icut, metal);
    setBest(icut, benefit);
    setMest(icut, grade);
  }
  return 0;
}

int Selectivity::calculateFromAnam(AAnam* anam)
{
  AnamHermite* anamH = dynamic_cast<AnamHermite*>(anam);
  if (anamH != nullptr)
  {
    anamH->_globalSelectivity(this);
    return 0;
  }

  AnamDiscreteDD* anamDD = dynamic_cast<AnamDiscreteDD*>(anam);
  if (anamDD != nullptr)
  {
    anamDD->_globalSelectivity(this);
    return 0;
  }

  AnamDiscreteIR* anamIR = dynamic_cast<AnamDiscreteIR*>(anam);
  if (anamIR != nullptr)
  {
    anamIR->_globalSelectivity(this);
    return 0;
  }

  messerr("Code not yet implemented for current anamorphosis");
  return 1;
}

const Table Selectivity::eval(const Db *db, bool autoCuts)
{
  (void) calculateFromDb(db, autoCuts);
  return getStats();
}
const Table Selectivity::eval(const VectorDouble &tab,
                              const VectorDouble &weights,
                              bool autoCuts)
{
  (void) calculateFromArray(tab, weights, autoCuts);
  return getStats();
}
const Table Selectivity::eval(AAnam *anam)
{
  (void) calculateFromAnam(anam);
  return getStats();
}

/****************************************************************************/
/*!
 **  Interpolate the Grade-Tonnage curves
 **
 ** \param[in] zcuts    Array of cutoffs
 ** \param[in] selecin  Input Selectivity
 ** \param[in] verbose  Verbose flag
 **
 *****************************************************************************/
Selectivity* Selectivity::createInterpolation(const VectorDouble& zcuts,
                                              const Selectivity& selecin,
                                              bool verbose)
{
  double tval, qval;

  int nclass = selecin.getNCuts();
  int ncuts = (int) zcuts.size();

  Selectivity* selectivity = new Selectivity(ncuts);
  for (int icut = 0; icut < ncuts; icut++)
  {
    double zval = zcuts[icut];
    selectivity->setZcut(icut, zval);

    /* Find interval to which cutoffs belongs */

    int iclass = -1;
    for (int jclass = 0; jclass < nclass - 1 && iclass < 0; jclass++)
    {
      double valmin = MIN(selecin.getZcut(jclass), selecin.getZcut(jclass + 1));
      double valmax = MAX(selecin.getZcut(jclass), selecin.getZcut(jclass + 1));
      if (zval >= valmin && zval <= valmax) iclass = jclass;
    }

    if (iclass >= 0 && iclass < nclass)
    {

      /* Assuming that cutoffs belongs to the interval the class 'iclass' */

      double zi0 = selecin.getZcut(iclass);
      double zi1 = (iclass + 1 > nclass - 1) ? 0. : selecin.getZcut(iclass + 1);
      double ti0 = selecin.getTest(iclass);
      double ti1 = (iclass + 1 > nclass - 1) ? 0. : selecin.getTest(iclass + 1);
      double qi0 = selecin.getQest(iclass);
      double qi1 = selecin.getQest(iclass + 1);
      selectivity->_interpolateInterval(zval, zi0, zi1, ti0, ti1, qi0, qi1, &tval, &qval);
      selectivity->setTest(icut, tval);
      selectivity->setQest(icut, qval);
    }
    else
    {
      selectivity->setTest(icut, 0.);
      selectivity->setQest(icut, 0.);
    }
  }

  selectivity->calculateBenefitAndGrade();
  if (verbose) selectivity->dumpGini();

  return selectivity;
}

void Selectivity::setZcut(int iclass, double zcut)
{
  if (! _isValidCut(iclass)) return;
  _stats.setValue(iclass, Z_CUT, zcut);
  _Zcut[iclass] = zcut;
}
double Selectivity::getZcut(int iclass) const
{
  if (! _isValidCut(iclass)) return(TEST);
  return _Zcut[iclass];
}
void Selectivity::setBest(int iclass, double best)
{
  if (! _isValidCut(iclass)) return;
  _stats.setValue(iclass, B_EST, best);
}
void Selectivity::setMest(int iclass, double mest)
{
  if (! _isValidCut(iclass)) return;
  _stats.setValue(iclass, M_EST, mest);
}

void Selectivity::setQest(int iclass, double qest)
{
  if (! _isValidCut(iclass)) return;
  _stats.setValue(iclass, Q_EST, qest);
}
void Selectivity::setQstd(int iclass, double qstd)
{
  if (! _isValidCut(iclass)) return;
  _stats.setValue(iclass, Q_STD, qstd);
}
void Selectivity::setTest(int iclass, double test)
{
  if (! _isValidCut(iclass)) return;
  _stats.setValue(iclass, T_EST, test);
}
void Selectivity::setTstd(int iclass, double tstd)
{
  if (! _isValidCut(iclass)) return;
  _stats.setValue(iclass, T_STD, tstd);
}
double Selectivity::getBest(int iclass) const
{
  if (! _isValidCut(iclass)) return(TEST);
  return _stats.getValue(iclass, B_EST);
}
double Selectivity::getMest(int iclass) const
{
  if (! _isValidCut(iclass)) return(TEST);
  return _stats.getValue(iclass, M_EST);
}
double Selectivity::getQest(int iclass) const
{
  if (! _isValidCut(iclass)) return(TEST);
  return _stats.getValue(iclass, Q_EST);
}
double Selectivity::getQstd(int iclass) const
{
  if (! _isValidCut(iclass)) return(TEST);
  return _stats.getValue(iclass, Q_STD);
}
double Selectivity::getTest(int iclass) const
{
  if (! _isValidCut(iclass)) return(TEST);
  return _stats.getValue(iclass, T_EST);
}
double Selectivity::getTstd(int iclass) const
{
  if (! _isValidCut(iclass)) return(TEST);
  return _stats.getValue(iclass, T_STD);
}

bool Selectivity::_isValidCut(int iclass) const
{
  if (iclass < 0 || iclass >= getNCuts())
  {
    mesArg("Selectivity Class", iclass, getNCuts());
    return false;
  }
  return true;
}

/*****************************************************************************/
/*!
 **  Analyze the contents of the codes
 **
 ** \param[in]  codes        Array of selectivity codes
 ** \param[in]  flag_est     True for estimation
 ** \param[in]  flag_std     True for st. dev.
 ** \param[in]  proba        Probability value (or TEST)
 ** \param[in]  verbose      Verbose flag
 **
 *****************************************************************************/
void Selectivity::defineRecoveries(const std::vector<ESelectivity>& codes,
                                   bool flag_est,
                                   bool flag_std,
                                   double proba,
                                   bool verbose)
{
  int ncode = (int) codes.size();
  _proba = proba;
  _numberQT.reset(getNQT(), 2);
  _numberQT.fill(0);

  // Optional printout (title)

  if (verbose) mestitle(1, "List of options");

  /* Check for the presence of other codes */

  for (int icode = 0; icode < ncode; icode++)
  {
    const ESelectivity& code = codes[icode];

    int key = code.getValue();
    switch (code.toEnum())
    {
      case ESelectivity::E_UNKNOWN:
        break;

      case ESelectivity::E_Z:
        if (flag_est)
        {
          _numberQT.setValue(key, 0, 1);
          if (verbose) _printQTvars("Average", 1, 1);
        }
        if (flag_std)
        {
          _numberQT.setValue(key, 1, 1);
          if (verbose) _printQTvars("Average", 2, 1);
        }
        break;

      case ESelectivity::E_T:
        if (getNCuts() <= 0) break;
        if (flag_est)
        {
          _numberQT.setValue(key, 0, getNCuts());
          if (verbose) _printQTvars("Tonnage", 1, getNCuts());
        }
        if (flag_std)
        {
          _numberQT.setValue(key, 1, getNCuts());
          if (verbose) _printQTvars("Tonnage", 2, getNCuts());
        }
        break;

      case ESelectivity::E_Q:
        if (getNCuts() <= 0) break;
        if (flag_est)
        {
          _numberQT.setValue(key, 0, getNCuts());
          if (verbose) _printQTvars("Metal Quantity", 1, getNCuts());
        }
        if (flag_std)
        {
          _numberQT.setValue(key, 1, getNCuts());
          if (verbose) _printQTvars("Metal Quantity", 2, getNCuts());
        }
        break;

      case ESelectivity::E_B:
        if (getNCuts() <= 0) break;
        if (flag_est)
        {
          _numberQT.setValue(key, 0, getNCuts());
          if (verbose) _printQTvars("Conventional Benefit", 1, getNCuts());
        }
        break;

      case ESelectivity::E_M:
        if (getNCuts() <= 0) break;
        if (flag_est)
        {
          _numberQT.setValue(key, 1, getNCuts());
          if (verbose) _printQTvars("Average Metal", 1, getNCuts());
        }
        break;

      case ESelectivity::E_PROBA:
        if (getNCuts() <= 0) break;
        if (flag_est)
        {
          _numberQT.setValue(key, 0, getNCuts());
          if (verbose) _printQTvars("Probability", 1, getNCuts());
        }
        break;

      case ESelectivity::E_QUANT:
        if (FFFF(_proba)) break;
        if (flag_est)
        {
          _numberQT.setValue(key, 0, 1);
          if (verbose) _printQTvars("Quantile", 1, 1);
        }
        break;
    }
  }

  /* Count the total number of variables */

  int ntotal = getVariableNumber();
  if (ntotal <= 0)
  {
    messerr("The number of variables calculated is zero");
    messerr("Check the recovery function (the number of cutoffs is %d)", getNCuts());
  }
  else
    _defineVariableRanks();
}

void Selectivity::_defineVariableRanks()
{
  _rankQT.reset(getNQT(), 2);
  _rankQT.fill(-1);

  int rank = 0;
  for (int i = 0; i < getNQT(); i++)
  {
    if (_numberQT.getValue(i, 0) > 0)
    {
      _rankQT.setValue(i, 0, rank);
      rank += _numberQT.getValue(i, 0);
    }
    if (_numberQT.getValue(i, 1) > 0)
    {
      _rankQT.setValue(i, 1, rank);
      rank += _numberQT.getValue(i, 1);
    }
  }
}

bool Selectivity::_isMultiplied(const ESelectivity& code) const
{
  if (code == ESelectivity::UNKNOWN ||
      code == ESelectivity::Z       ||
      code == ESelectivity::PROBA   ||
      code == ESelectivity::QUANT)
    return false;
  else
    return true;
}

VectorString Selectivity::getVariableNames() const
{
  VectorString names;

  for (int i = 0; i < getNQT(); i++)
  {
    if (_numberQT.getValue(i, 0) > 0)
    {
      if (_isMultiplied(i))
        for (int icut = 0; icut < getNCuts(); icut++)
          names.push_back(getVariableName(i, icut, 0));
      else
        names.push_back(getVariableName(i, 0, 0));
    }
    if (_numberQT.getValue(i, 1) > 0)
    {
      if (_isMultiplied(i))
        for (int icut = 0; icut < getNCuts(); icut++)
          names.push_back(getVariableName(i, icut, 1));
      else
        names.push_back(getVariableName(i, 0, 1));
    }
  }
  return names;
}

String Selectivity::getVariableName(const ESelectivity& code, int icut, int mode) const
{
  switch (code.toEnum())
  {
    case ESelectivity::E_UNKNOWN:
      break;

    case ESelectivity::E_Z:
      if (mode == 0)
        return("Z-Est");
      else
        return("Z-Std");
      break;

    case ESelectivity::E_T:
      if (mode == 0)
        return(encodeString("T-Est", _Zcut[icut]));
      else
        return(encodeString("T-Std", _Zcut[icut]));
      break;

    case ESelectivity::E_Q:
      if (mode == 0)
        return(encodeString("Q-Est", _Zcut[icut]));
      else
        return(encodeString("Q-Std", _Zcut[icut]));
      break;

    case ESelectivity::E_B:
      if (mode == 0)
        return(encodeString("B-Est", _Zcut[icut]));
      else
        return(encodeString("B-Std", _Zcut[icut]));
      break;

    case ESelectivity::E_M:
      if (mode == 0)
        return(encodeString("M-Est", _Zcut[icut]));
      else
        return(encodeString("M-Std", _Zcut[icut]));
      break;

    case ESelectivity::E_PROBA:
      if (mode == 0)
        return("Proba-Est");
      else
        return("Proba-Std");
      break;

    case ESelectivity::E_QUANT:
      if (mode == 0)
        return("Quant-Est");
      else
        return("Quant-Std");
      break;
  }

  return String();
}

String Selectivity::getVariableName(int rank0) const
{
  int rank = -1;
  for (int code = 0; code < getNQT(); code++)
  {
    if (_numberQT.getValue(code, 0) > 0)
    {
      if (_isMultiplied(code))
        for (int icut = 0; icut < getNCuts(); icut++)
        {
          rank++;
          if (rank == rank0) return getVariableName(code, icut, 0);
        }
      else
      {
        rank++;
        if (rank == rank0) return getVariableName(code, 0, 0);
      }
    }
    if (_numberQT.getValue(code, 1) > 0)
    {
      if (_isMultiplied(code))
        for (int icut = 0; icut < getNCuts(); icut++)
        {
          rank++;
          if (rank == rank0) return getVariableName(code, icut, 1);
        }
      else
      {
        rank++;
        if (rank == rank0) return getVariableName(code, 0, 1);
      }
    }
  }
  return String();
}

/*****************************************************************************/
/*!
 **  Print the contents of the qtvars structure
 **
 ** \param[in]  title        Title
 ** \param[in]  type         1 for estimation; 2 for stdev
 ** \param[in]  number       Number of cutoffs
 **
 *****************************************************************************/
void Selectivity::_printQTvars(const char *title, int type, int number) const
{
  message("- %s", title);
  if (type == 1)
    message(" (Estimation)");
  else
    message(" (St. Deviation)");
  message(": %d\n", number);
}

bool Selectivity::_isRecoveryDefined() const
{
  if (_numberQT.empty())
  {
    messerr("No recovery function has been defined yet");
    return false;
  }
  return true;
}
int Selectivity::getVariableNumber() const
{
  int ntotal = 0;
  if (! _isRecoveryDefined()) return ntotal;
  for (int i = 0; i < getNQT(); i++)
  {
    ntotal += _numberQT.getValue(i, 0);
    ntotal += _numberQT.getValue(i, 1);
  }
  return ntotal;
}

bool Selectivity::isUsed(const ESelectivity& code) const
{
  if (code == ESelectivity::UNKNOWN) return false;
  if (! _isRecoveryDefined()) return false;
  int key = code.getValue();
  if (_numberQT.getValue(key, 0) > 0) return true;
  if (_numberQT.getValue(key, 1) > 0) return true;
  return false;
}

bool Selectivity::isUsedEst(const ESelectivity& code) const
{
  if (code == ESelectivity::UNKNOWN) return false;
  if (! _isRecoveryDefined()) return false;
  int key = code.getValue();
  if (_numberQT.getValue(key, 0) > 0) return true;
  return false;
}

bool Selectivity::isUsedStD(const ESelectivity& code) const
{
  if (code == ESelectivity::UNKNOWN) return false;
  int key = code.getValue();
  if (_numberQT.getValue(key, 1) > 0) return true;
  return false;
}

bool Selectivity::isNeededT() const
{
  if (isUsed(ESelectivity::T)) return true;
  if (isUsed(ESelectivity::B)) return true;
  if (isUsed(ESelectivity::M)) return true;
  if (isUsed(ESelectivity::PROBA)) return true;
  return false;
}

bool Selectivity::isNeededQ() const
{
  if (isUsed(ESelectivity::Q)) return true;
  if (isUsed(ESelectivity::B)) return true;
  if (isUsed(ESelectivity::M)) return true;
  return false;
}

int Selectivity::getAddressQTEst(const ESelectivity& code, int iptr0, int rank) const
{
  if (code == ESelectivity::UNKNOWN) return -1;
  int key = code.getValue();
  if (rank < 0 || rank >= _numberQT.getValue(key, 0)) return -1;
  return (iptr0 + _rankQT.getValue(key, 0) + rank);
}

int Selectivity::getAddressQTStd(const ESelectivity& code, int iptr0, int rank) const
{
  if (code == ESelectivity::UNKNOWN) return -1;
  if (! _isRecoveryDefined()) return -1;
  int key = code.getValue();
  if (rank < 0 || rank >= _numberQT.getValue(key, 1)) return -1;
  return (iptr0 + _rankQT.getValue(key, 1) + rank);
}

/****************************************************************************/
/*!
 **  From the cutoff, tonnage and metal quantity, derive
 **  the conventional benefit and the average recovered grade
 **
 *****************************************************************************/
void Selectivity::calculateBenefitAndGrade()
{
  int iclass;
  double zval, tval, qval;
  int nclass = getNCuts();

  for (iclass = 0; iclass < nclass; iclass++)
  {
    zval = getZcut(iclass);
    tval = getTest(iclass);
    qval = getQest(iclass);
    setBest(iclass, qval - zval * tval);
    setMest(iclass ,(ABS(tval) < EPSILON6) ? TEST : qval / tval);
  }
}

/****************************************************************************/
/*!
 **  Calculate and print the Gini index
 **
 *****************************************************************************/
void Selectivity::dumpGini()
{
  int nclass = getNCuts();

  double gini = 1.;
  for (int iclass = 0; iclass < nclass - 1; iclass++)
    gini -= ((getTest(iclass)
        - getTest(iclass + 1))
             * (getQest(iclass) + getQest(iclass + 1)));

  message("Gini calculated on %d classes\n", nclass);
  message("Value of the Gini index = %lf\n", gini);
}

/*****************************************************************************/
/*!
 **  Correct the order relationship for Tonnage
 **
 *****************************************************************************/
void Selectivity::correctTonnageOrder()
{
  if (! isFlagTonnageCorrect()) return;
  int nclass = getNCuts();
  VectorDouble ta(nclass);
  VectorDouble tb(nclass);

  for (int iclass = nclass - 1; iclass >= 0; iclass--)
  {
    double auxval = getTest(iclass);
    if (iclass < nclass - 1) auxval = MAX(ta[iclass + 1], auxval);
    ta[iclass] = MIN(1., MAX(0., auxval));
  }

  for (int iclass = 0; iclass < nclass; iclass++)
  {
    double auxval = getTest(iclass);
    if (iclass > 0) auxval = MIN(tb[iclass - 1], auxval);
    tb[iclass] = MAX(0., MIN(1., auxval));
  }

  for (int iclass = 0; iclass < nclass; iclass++)
    setTest(iclass, 0.5 * (ta[iclass] + tb[iclass]));
}

/*****************************************************************************/
/*!
 **  Interpolate the QT within an interval
 **
 ** \param[in]  zval     Cutoff value
 ** \param[in]  zi0      Lower cutoff of the interval
 ** \param[in]  zi1      Upper cutoff of the interval
 ** \param[in]  ti0      Lower tonnage of the interval
 ** \param[in]  ti1      Upper tonnage of the interval
 ** \param[in]  qi0      Lower metal quantity of the interval
 ** \param[in]  qi1      Upper metal quantity of the interval
 ** \param[in]  tol      Tolerance
 **
 ** \param[out] tval     Tonnage for the current cutoff
 ** \param[out] qval     Metal quantity for the current cutoff
 **
 *****************************************************************************/
void Selectivity::_interpolateInterval(double zval,
                                       double zi0,
                                       double zi1,
                                       double ti0,
                                       double ti1,
                                       double qi0,
                                       double qi1,
                                       double *tval,
                                       double *qval,
                                       double tol)
{
  double dzi = zi1 - zi0;
  double dti = ti1 - ti0;
  double zmoy = (qi1 - qi0) / (ti1 - ti0);
  double aa0 = (zi1 - zmoy) / (zmoy - zi0);

  if (ABS(zval - zi0) < tol)
  {
    (*tval) = ti0;
    (*qval) = qi0;
    return;
  }

  if (ABS(zval - zi1) < tol)
  {
    (*tval) = ti1;
    (*qval) = qi1;
    return;
  }

  double u = (zval - zi0) / dzi;
  (*tval) = (u <= 0.) ? ti0 : ti0 + dti * pow(u, 1. / aa0);
  (*qval) = (u <= 0.) ? qi0 :
      qi0 + zi0 * ((*tval) - ti0)
      + dzi * dti * pow(u, 1. + 1. / aa0) / (1. + aa0);
}

/****************************************************************************/
/*!
 **  Store the local results of the recovery
 **
 ** \param[in]  db          Db structure containing the factors (Z-locators)
 ** \param[in]  iech0       Rank of the target sample
 ** \param[in]  iptr        Rank for storing the results
 ** \param[in]  zestim      Estimated grade
 ** \param[in]  zstdev      St. dev.
 **
 *****************************************************************************/
void Selectivity::storeInDb(Db *db,
                            int iech0,
                            int iptr,
                            double zestim,
                            double zstdev)
{
  int ncut = getNCuts();

  /* Store the recovered grade */

  if (isUsedEst(ESelectivity::Z))
    db->setArray(iech0, getAddressQTEst(ESelectivity::Z, iptr), zestim);
  if (isUsedStD(ESelectivity::Z))
    db->setArray(iech0, getAddressQTStd(ESelectivity::Z, iptr), zstdev);

  /* Loop on the cutoff classes */

  for (int icut = 0; icut < ncut; icut++)
  {
    double tval = getTest(icut);
    double qval = getQest(icut);
    double bval = getBest(icut);
    double mval = getMest(icut);
    double tstd = getTstd(icut);
    double qstd = getQstd(icut);

    // Tonnage

    if (isUsedEst(ESelectivity::T))
      db->setArray(iech0,
                   getAddressQTEst(ESelectivity::T, iptr, icut), tval);
    if (isUsedStD(ESelectivity::T))
      db->setArray(iech0,
                   getAddressQTStd(ESelectivity::T, iptr, icut), tstd);

    // Metal Quantity

    if (isUsedEst(ESelectivity::Q))
      db->setArray(iech0,
                   getAddressQTEst(ESelectivity::Q, iptr, icut), qval);
    if (isUsedStD(ESelectivity::Q))
      db->setArray(iech0,
                   getAddressQTStd(ESelectivity::Q, iptr, icut), qstd);

    // Conventional Benefit
    if (isUsedEst(ESelectivity::B))
      db->setArray(iech0,
                   getAddressQTEst(ESelectivity::B, iptr, icut), bval);

    // Average recovered Grade
    if (isUsedEst(ESelectivity::M))
      db->setArray(iech0,
                   getAddressQTEst(ESelectivity::M, iptr, icut), mval);
  }
}

/*****************************************************************************/
/*!
 **  Interpolate the QT curves (Local estimation)
 **
 ** \param[in]  selecin  Selectivity
 **
 *****************************************************************************/
void Selectivity::interpolateSelectivity(const Selectivity* selecin)
{
  double tval, qval;

  double z_max = getZmax();
  int nclass = selecin->getNCuts();
  int ncutmine = getNCuts();
  VectorDouble zz(nclass + 2);
  VectorDouble TT(nclass + 2);
  VectorDouble QQ(nclass + 2);
  VectorDouble zcuts = getZcut();

  /* Load arrays */

  int ncleff = 1;
  TT[0] = QQ[0] = 0.;
  for (int iclass = 0; iclass < nclass; iclass++)
  {
    int jclass = nclass - iclass - 1;
    if (selecin->getTest(jclass) <= TT[ncleff - 1]) continue;
    TT[ncleff] = selecin->getTest(jclass);
    QQ[ncleff] = selecin->getQest(jclass);
    ncleff++;
  }
  zz[0] = z_max;
  for (int iclass = 0; iclass < ncleff - 1; iclass++)
    zz[iclass + 1] = (QQ[iclass + 2] - QQ[iclass]) / (TT[iclass + 2] - TT[iclass]);
  zz[ncleff - 1] = 0.;
  if (FFFF(z_max)) zz[0] = 2 * zz[1];

  for (int icut = 0; icut < ncutmine; icut++)
  {
    double zval = zcuts[icut];
    setZcut(icut, zval);

    /* Find interval to which cutoffs belongs */

    int iclass = -1;
    for (int jclass = 0; jclass < ncleff && iclass < 0; jclass++)
      if ((zval - zz[jclass]) * (zval - zz[jclass + 1]) <= 0) iclass = jclass;

    /* Assuming that cutoffs belongs to the interval the class 'iclass' */

    double zi0 = zz[iclass];
    double zi1 = (iclass + 1 > ncleff - 1) ? 0. : zz[iclass + 1];
    double ti0 = TT[iclass];
    double ti1 = (iclass + 1 > ncleff - 1) ? 0. : TT[iclass + 1];
    double qi0 = QQ[iclass];
    double qi1 = QQ[iclass + 1];
    _interpolateInterval(zval, zi0, zi1, ti0, ti1, qi0, qi1, &tval, &qval);
    setTest(icut, tval);
    setQest(icut, qval);
  }
  return;
}

VectorString Selectivity::_getAllNames() const
{
  VectorString names;
  names.push_back("Z-Cut");
  names.push_back("T-Est");
  names.push_back("Q-Est");
  names.push_back("B-Est");
  names.push_back("M-Est");
  names.push_back("T-Std");
  names.push_back("Q-Std");
  return names;
}

String Selectivity::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  int ncut = getNCuts();
  if (ncut <= 0) return sstr.str();
  sstr << toTitle(0, "Selectivity Curves");
  sstr << _stats.toString();
  return sstr.str();
}

int Selectivity::getNumberQTEst(const ESelectivity& code) const
{
  return _numberQT.getValue(code.getValue(), 0);
}

int Selectivity::getNumberQTStd(const ESelectivity& code) const
{
  return _numberQT.getValue(code.getValue(), 1);
}

const VectorInt Selectivity::getNumberQTEst() const
{
  return _numberQT.getValuesPerColumn(0);
}

const VectorInt Selectivity::getNumberQTStd() const
{
  return _numberQT.getValuesPerColumn(1);
}

void Selectivity::_defineAutomaticCutoffs(const VectorDouble& tab, double eps)
{
  double zmin = ut_vector_min(tab);
  double zmax = ut_vector_max(tab) + eps;
  int ncuts = getNCuts();

  if (ncuts <= 1)
  {
    messerr("Automatic Cutoffs can only be calculated for more than 1 cutoff");
    return;
  }
  for (int icut = 0; icut < ncuts; icut++)
    _Zcut[icut] = zmin + (zmax - zmin) * (double) icut / ((double) (ncuts - 1));
}

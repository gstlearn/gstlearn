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

#include "Stats/Selectivity.hpp"
#include "Stats/ESelectivity.hpp"

Selectivity::Selectivity(int ncut)
    : _nCut(ncut),
      _Zcut(ncut)
{
  ut_vector_fill(_Zcut, TEST);
}

Selectivity::Selectivity(const VectorDouble& zcuts)
    : _nCut(0),
      _Zcut(zcuts)
{
  _nCut = (int) _Zcut.size();
}

Selectivity::Selectivity(const Selectivity &m)
    : _nCut(m._nCut),
      _Zcut(m._Zcut)
{

}

Selectivity& Selectivity::operator=(const Selectivity &m)
{
  if (this != &m)
  {
    _nCut = m._nCut;
    _Zcut = m._Zcut;
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

Selectivity* Selectivity::createByCodes(const std::vector<ESelectivity>& codes,
                                        bool flag_est,
                                        bool flag_std,
                                        bool flag_inter,
                                        int ncut,
                                        double proba,
                                        bool verbose)
{
  Selectivity* selectivity = new Selectivity(ncut);
  selectivity->defineRecoveries(codes, flag_est, flag_std, flag_inter, proba, verbose);
  return selectivity;
}

Selectivity* Selectivity::createByKeys(const VectorString& scodes,
                                       bool flag_est,
                                       bool flag_std,
                                       bool flag_inter,
                                       int ncut,
                                       double proba,
                                       bool verbose)
{
  std::vector<ESelectivity> codes;

  for (int i = 0; i < (int) scodes.size(); i++)
  {
    ESelectivity code = ESelectivity::fromKey(scodes[i]);
    if (code == ESelectivity::UNKNOWN) continue;
    codes.push_back(code);
  }

  Selectivity* selectivity = new Selectivity(ncut);
  selectivity->defineRecoveries(codes, flag_est, flag_std, flag_inter, proba, verbose);
  return selectivity;
}

void Selectivity::setZcut(int icut, double zcut)
{
  if (! _isValid(icut)) return;
  _Zcut[icut] = zcut;
}

double Selectivity::getZcut(int icut) const
{
  if (! _isValid(icut)) return(TEST);
  return _Zcut[icut];
}

bool Selectivity::_isValid(int icut) const
{
  if (icut < 0 || icut >= getNCuts())
  {
    mesArg("Selectivity Class", icut, getNCuts());
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
 ** \param[in]  flag_inter   QT must be interpolated
 ** \param[in]  proba        Probability value (or TEST)
 ** \param[in]  verbose      Verbose flag
 **
 *****************************************************************************/
void Selectivity::defineRecoveries(const std::vector<ESelectivity>& codes,
                                   bool flag_est,
                                   bool flag_std,
                                   bool flag_inter,
                                   double proba,
                                   bool verbose)
{
  int ncode = (int) codes.size();
  if (flag_inter) flag_std = false;
  _numberQTEst.resize(getNQT(), 0);
  _numberQTStd.resize(getNQT(), 0);

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
          _numberQTEst[key] = 1;
          if (verbose) _printQTvars("Average", 1, 1);
        }
        if (flag_std)
        {
          _numberQTStd[key] = 1;
          if (verbose) _printQTvars("Average", 2, 1);
        }
        break;

      case ESelectivity::E_T:
        if (_nCut <= 0) break;
        if (flag_est)
        {
          _numberQTEst[key] =  _nCut;
          if (verbose) _printQTvars("Tonnage", 1, _nCut);
        }
        if (flag_std)
        {
          _numberQTStd[key] = _nCut;
          if (verbose) _printQTvars("Tonnage", 2, _nCut);
        }
        break;

      case ESelectivity::E_Q:
        if (_nCut <= 0) break;
        if (flag_est)
        {
          _numberQTEst[key] = _nCut;
          if (verbose) _printQTvars("Metal Quantity", 1, _nCut);
        }
        if (flag_std)
        {
          _numberQTStd[key] = _nCut;
          if (verbose) _printQTvars("Metal Quantity", 2, _nCut);
        }
        break;

      case ESelectivity::E_B:
        if (_nCut <= 0) break;
        if (flag_est)
        {
          _numberQTEst[key] = _nCut;
          if (verbose) _printQTvars("Conventional Benefit", 1, _nCut);
        }
        break;

      case ESelectivity::E_M:
        if (_nCut <= 0) break;
        if (flag_est)
        {
          _numberQTEst[key] = _nCut;
          if (verbose) _printQTvars("Average Metal", 1, _nCut);
        }
        break;

      case ESelectivity::E_PROBA:
        if (_nCut <= 0) break;
        if (flag_est)
        {
          _numberQTEst[key] = _nCut;
          if (verbose) _printQTvars("Probability", 1, _nCut);
        }
        break;

      case ESelectivity::E_QUANT:
        if (FFFF(proba)) break;
        if (flag_est)
        {
          _numberQTEst[key] = 1;
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
    messerr("Check the recovery function (the number of cutoffs is %d)", _nCut);
  }
  else
    _defineVariableRanks();
}

void Selectivity::_defineVariableRanks()
{
  _rankQTEst.resize(getNQT(),-1);
  _rankQTStd.resize(getNQT(),-1);

  int rank = 0;
  for (int i = 0; i < getNQT(); i++)
  {
    if (_numberQTEst[i] > 0)
    {
      _rankQTEst[i] = rank;
      rank += _numberQTEst[i];
    }
    if (_numberQTStd[i] > 0)
    {
      _rankQTStd[i] = rank;
      rank += _numberQTStd[i];
    }
  }
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

int Selectivity::getVariableNumber() const
{
  int ntotal = 0;
  for (int i = 0; i < getNQT(); i++)
  {
    if (_numberQTEst[i]) ntotal += _numberQTEst[i];
    if (_numberQTStd[i]) ntotal += _numberQTStd[i];
  }
  return ntotal;
}


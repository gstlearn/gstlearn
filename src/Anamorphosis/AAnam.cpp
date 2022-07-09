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
#include "geoslib_enum.h"
#include "geoslib_old_f.h"

#include "Anamorphosis/AAnam.hpp"
#include "Db/Db.hpp"
#include "Basic/AException.hpp"
#include "Stats/Selectivity.hpp"

#include "math.h"

#define QT_EST    0
#define QT_STD    1
#define QT_VARS(i,j)              (qt_vars[(i) + 2 * (j)])
#define QT_FLAG(j)                (QT_VARS(QT_EST,j) > 0 || \
                                   QT_VARS(QT_STD,j) > 0)

AAnam::AAnam()
    : AStringable(),
      ASerializable()
{
}

AAnam::AAnam(const AAnam &m)
    : AStringable(m),
      ASerializable(m)
{
}

AAnam& AAnam::operator=(const AAnam &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    ASerializable::operator=(m);
  }
  return *this;
}

AAnam::~AAnam()
{
}

VectorDouble AAnam::z2factor(double /*z*/, const VectorInt& /*nfact*/) const
{
  messerr("This function is not programmed yet");
  return VectorDouble();
}

/*****************************************************************************/
/*!
 **  Find the coefficient of change of support
 **
 ** \return  Value for the change of support coefficient
 **
 ** \param[in]  cvv      Mean covariance value over a block
 **
 *****************************************************************************/
double AAnam::invertVariance(double cvv)
{
  if (! allowChangeSupport()) return TEST;
  double s0, s1, s2, var0, var1;
  int converge, niter;
  static int niter_max = 1000;

  s1 = 0.;
  var1 = computeVariance(s1);

  /* Dichotomy */

  s2 = 1.;
  converge = niter = 0;
  while (!converge)
  {
    niter++;
    s0 = (s1 + s2) / 2.;
    var0 = computeVariance(s0);
    converge = (ABS(var0 - cvv) < EPSILON8 || niter > niter_max);
    if ((var1 - cvv) * (var0 - cvv) < 0.)
    {
      s2 = s0;
    }
    else
    {
      s1 = s0;
      var1 = var0;
    }
  }
  return (s0);
}

/*****************************************************************************/
/*!
 **  Calculates the block variance
 **
 ** \return Value of the block variance (as a function of support coefficient)
 **
 *****************************************************************************/
double AAnam::computeVariance(double /*sval*/) const
{
  messerr("This function is not programmed yet");
  return TEST;
}

int AAnam::updatePointToBlock(double /*r_coef*/)
{
  messerr("This function is not programmed yet");
  return 1;
}

double AAnam::RawToTransformValue(double /*z*/) const
{
  if (! hasGaussian())
    messerr("This function is not possible");
  else
    messerr("This function is not programmed yet");
  return TEST;
}

double AAnam::TransformToRawValue(double /*y*/) const
{
  if (! hasGaussian())
    messerr("This function is not available");
  else
    messerr("This function is not programmed yet");
  return TEST;
}

/*****************************************************************************/
/*!
 **  Check if a sample must be considered or not
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  iech         Rank of the target sample
 ** \param[in]  cols_est     Array of columns for factor estimation
 ** \param[in]  cols_std     Array of columns for factor st. dev.
 **
 *****************************************************************************/
bool AAnam::_isSampleSkipped(Db *db,
                             int iech,
                             const VectorInt& cols_est,
                             const VectorInt& cols_std)
{
  double value;

  if (!db->isActive(iech)) return true;

  int nb_est = (int) cols_est.size();
  if (!cols_est.empty())
  {
    for (int ivar = 0; ivar < nb_est; ivar++)
    {
      value = db->getArray(iech, cols_est[ivar]);
      if (FFFF(value)) return true;
    }
  }

  int nb_std = (int) cols_std.size();
  if (!cols_std.empty())
  {
    for (int ivar = 0; ivar < nb_std; ivar++)
    {
      value = db->getArray(iech, cols_std[ivar]);
      if (FFFF(value)) return true;
    }
  }
  return false;
}

/****************************************************************************/
/*!
 **  Store the local results of the recovery
 **
 ** \param[in]  db          Db structure containing the factors (Z-locators)
 ** \param[in]  iech0       Rank of the target sample
 ** \param[in]  iptr        Rank for storing the results
 ** \param[in]  codes       Array of codes for stored results
 ** \param[in]  qt_vars     Array with the number of output variables
 ** \param[in]  zestim      Estimated grade
 ** \param[in]  zstdev      St. dev.
 ** \param[in]  calest      Selectivity
 **
 *****************************************************************************/
void AAnam::recoveryLocal(Db *db,
                          int iech0,
                          int iptr,
                          const VectorInt& codes,
                          const VectorInt& qt_vars,
                          double zestim,
                          double zstdev,
                          const Selectivity& calest)
{
  int jptr = iptr;
  int nclass = calest.getNCuts();
  int ncode = (int) codes.size();

  /* Store the recovered grade */

  if (codes[0] == 0)
  {
    if (QT_VARS(QT_EST,ANAM_QT_Z) > 0) db->setArray(iech0, jptr++, zestim);
    if (QT_VARS(QT_STD,ANAM_QT_Z) > 0) db->setArray(iech0, jptr++, zstdev);
  }

  /* Loop on the recovery functions */

  for (int icode = 0; icode < ncode; icode++)
  {

    /* Loop on the cutoff classes */

    for (int iclass = 0; iclass < nclass; iclass++)
    {
      double tval = calest.getTest(iclass);
      double qval = calest.getQest(iclass);
      double bval = calest.getBest(iclass);
      double mval = calest.getMest(iclass);
      double tstd = calest.getTstd(iclass);
      double qstd = calest.getQstd(iclass);

      switch (codes[icode])
      {
        case 1: /* Tonnage */
          if (QT_VARS(QT_EST,ANAM_QT_T) > 0) db->setArray(iech0, jptr++, tval);
          if (QT_VARS(QT_STD,ANAM_QT_T) > 0) db->setArray(iech0, jptr++, tstd);
          break;

        case 2: /* Metal Quantity */
          if (QT_VARS(QT_EST,ANAM_QT_Q) > 0) db->setArray(iech0, jptr++, qval);
          if (QT_VARS(QT_STD,ANAM_QT_Q) > 0) db->setArray(iech0, jptr++, qstd);
          break;

        case 3: /* Conventional Benefit */
          if (QT_VARS(QT_EST,ANAM_QT_B) > 0) db->setArray(iech0, jptr++, bval);
          break;

        case 4: /* Average recovered grade */
          if (QT_VARS(QT_EST,ANAM_QT_M) > 0) db->setArray(iech0, jptr++, mval);
          break;
      }
    }
  }
}

/*****************************************************************************/
/*!
 **  Analyze the contents of the codes
 **
 ** \returns The number of different variables to be calculated
 **
 ** \param[in]  verbose      Verbose flag
 ** \param[in]  codes        Array of codes for stored results
 ** \param[in]  nb_est       Number of columns for factor estimation
 ** \param[in]  nb_std       Number of columns for factor st. dev.
 ** \param[in]  ncut         Number of cutoffs
 ** \param[in]  proba        Probability value
 ** \param[in]  flag_inter   QT must be interpolated
 **
 ** \param[out] qt_vars      Array with the number of output variables
 **
 ** \remark When the number of cutoff is zero, the flag of T, Q, B and M
 ** \remark are set to 0
 ** \remark When QT are interpolated, no variance can be calculated
 **
 ** \remark When the resulting number of variables is zero, an error
 ** \remark message is issued
 **
 *****************************************************************************/
int AAnam::codeAnalyze(bool verbose,
                       const VectorInt& codes,
                       int nb_est,
                       int nb_std,
                       int ncut,
                       double proba,
                       int flag_inter,
                       VectorInt& qt_vars) const
{
  int ncode = (int) codes.size();
  bool flag_est = nb_est > 0;
  bool flag_std = nb_std > 0 && !flag_inter;
  for (int i = 0; i < 2 * ANAM_N_QT; i++) qt_vars[i] = 0;

  // Optional printout (title)

  if (verbose) mestitle(1, "List of options");

  /* Check for the presence of other codes */

  for (int icode = 0; icode < ncode; icode++)
  {
    switch (codes[icode])
    {
      case ANAM_QT_Z:
        if (flag_est)
        {
          QT_VARS(QT_EST,ANAM_QT_Z) = 1;
          if (verbose) _printQTvars("Average", 1, 1);
        }
        if (flag_std)
        {
          QT_VARS(QT_STD,ANAM_QT_Z) = 1;
          if (verbose) _printQTvars("Average", 2, 1);
        }
        break;

      case ANAM_QT_T:
        if (! _isNcutValid(ncut)) return (0);
        if (flag_est)
        {
          QT_VARS(QT_EST,ANAM_QT_T) = ncut;
          if (verbose) _printQTvars("Tonnage", 1, ncut);
        }
        if (flag_std)
        {
          QT_VARS(QT_STD,ANAM_QT_T) = ncut;
          if (verbose) _printQTvars("Tonnage", 2, ncut);
        }
        break;

      case ANAM_QT_Q:
        if (! _isNcutValid(ncut)) return (0);
        if (flag_est)
        {
          QT_VARS(QT_EST,ANAM_QT_Q) = ncut;
          if (verbose) _printQTvars("Metal Quantity", 1, ncut);
        }
        if (flag_std)
        {
          QT_VARS(QT_STD,ANAM_QT_Q) = ncut;
          if (verbose) _printQTvars("Metal Quantity", 2, ncut);
        }
        break;

      case ANAM_QT_B:
        if (! _isNcutValid(ncut)) return (0);
        if (flag_est)
        {
          QT_VARS(QT_EST,ANAM_QT_B) = ncut;
          if (verbose) _printQTvars("Conventional Benefit", 1, ncut);
        }
        break;

      case ANAM_QT_M:
        if (! _isNcutValid(ncut)) return (0);
        if (flag_est)
        {
          QT_VARS(QT_EST,ANAM_QT_M) = ncut;
          if (verbose) _printQTvars("Average Metal", 1, ncut);
        }
        break;

      case ANAM_QT_PROBA:
        if (! _isNcutValid(ncut)) return (0);
        if (flag_est)
        {
          QT_VARS(QT_EST,ANAM_QT_PROBA) = ncut;
          if (verbose) _printQTvars("Probability", 1, ncut);
        }
        break;

      case ANAM_QT_QUANT:
        if (! _isProbaValid(proba)) return (0);
        if (flag_est)
        {
          QT_VARS(QT_EST,ANAM_QT_QUANT) = 1;
          if (verbose) _printQTvars("Quantile", 1, 1);
        }
        break;
    }
  }

  /* Count the total number of variables */

  int ntotal = 0;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < ANAM_N_QT; j++)
      ntotal += QT_VARS(i, j);

  if (ntotal <= 0)
  {
    messerr("The number of variables calculated is zero");
    messerr("Check the recovery function (the number of cutoffs is %d)", ncut);
  }

  return (ntotal);
}

bool AAnam::_isNcutValid(int ncut) const
{
  if (ncut <= 0)
  {
    messerr("The computing option requires Cutoffs to be defiend");
    return false;
  }
  return true;
}

bool AAnam::_isProbaValid(double proba) const
{
  if (FFFF(proba))
  {
    messerr("The computing option requires Proba to be defined");
    return false;
  }
  if (proba < 0 || proba > 1)
  {
    messerr("The computing option requires Proba to lie in [0,1]");
    return false;
  }
  return true;
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
void AAnam::_printQTvars(const char *title, int type, int number) const
{
  message("- %s", title);
  if (type == 1)
    message(" (Estimation)");
  else
    message(" (St. Deviation)");
  message(": %d\n", number);
}

/*****************************************************************************/
/*!
 **  Calculate the factors corresponding to an input data vector
 **
 ** \return  Error return code
 **
 ** \param[in]  db          Db structure
 ** \param[in]  ifacs       Array of factor ranks
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int AAnam::DbZToFactor(Db *db,
                       const VectorInt& ifacs,
                       const NamingConvention& namconv)
{
  if (db == nullptr)
  {
    messerr("You must define the 'db' argument");
    return 1;
  }
  int nvar = db->getVariableNumber();
  if (nvar != 1)
  {
    messerr("This function is only coded for the monovariate Db");
    return 1;
  }
  int nfact = (int) ifacs.size();
  if (nfact <= 0)
  {
    messerr("You must define the list of factors");
    return 1;
  }
  int nmax = getNFactor();
  for (int ifac = 0; ifac < nfact; ifac++)
    if (ifacs[ifac] < 1 || ifacs[ifac] > nmax)
    {
      messerr("Error in the rank of the factor(%d): it should lie in [1,%d]",
              ifacs[ifac], nmax);
      return 1;
    }

  /* Create the factors */

  int iptr = db->addColumnsByConstant(nfact, TEST);
  if (iptr <= 0) return 1;

  // Loop on the samples

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    double zval = db->getVariable(iech, 0);
    if (FFFF(zval)) continue;
    VectorDouble factors = z2factor(zval, ifacs);
    if (factors.empty()) continue;
    for (int ifac = 0; ifac < nfact; ifac++)
      db->setArray(iech, iptr + ifac, factors[ifac]);
  }

  /* Set the error return code */

  namconv.setNamesAndLocators(db, ELoc::Z, nfact, db, iptr, "Factor");

  return 0;
}

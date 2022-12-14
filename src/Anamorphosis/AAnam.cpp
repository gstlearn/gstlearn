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

#include <math.h>

#define QT_EST    0
#define QT_STD    1
#define QT_VARS(i,j)              (qt_vars[(i) + 2 * (j)])
#define QT_FLAG(j)                (QT_VARS(QT_EST,j) > 0 || \
                                   QT_VARS(QT_STD,j) > 0)

AAnam::AAnam()
    : AStringable(),
      ASerializable(),
      _flagFitted(false)
{
}

AAnam::AAnam(const AAnam &m)
    : AStringable(m),
      ASerializable(m),
      _flagFitted(m._flagFitted)
{
}

AAnam& AAnam::operator=(const AAnam &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    ASerializable::operator=(m);
    _flagFitted = m._flagFitted;
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
double AAnam::invertVariance(double cvv) const
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

bool AAnam::_isNcutValid(int ncut) const
{
  if (ncut <= 0)
  {
    messerr("The computing option requires Cutoffs to be defined");
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

VectorDouble AAnam::RawToTransformVec(const VectorDouble& z) const
{
  VectorDouble y = VectorDouble(z.size(), TEST);
  for (int i = 0; i < (int) z.size(); i++)
    y[i] = RawToTransformValue(z[i]);
  return y;
}

VectorDouble AAnam::TransformToRawVec(const VectorDouble& y) const
{
  VectorDouble z = VectorDouble(y.size(), TEST);
  for (int i = 0; i < (int) z.size(); i++)
    z[i] = TransformToRawValue(y[i]);
  return z;
}

int AAnam::fitFromLocator(Db *db, const ELoc& locatorType)
{
  int number = db->getLocatorNumber(locatorType);
  if (number != 1)
  {
    messerr("The number of items for locator(%d) is %d. It should be 1",
            locatorType.getValue(),number);
    return 1;
  }
  VectorDouble tab = db->getColumnByLocator(locatorType,0,true);
  VectorDouble wt;
  if (db->hasWeight())
    wt = db->getColumnByLocator(ELoc::W,0,true);

  if (fitFromArray(tab, wt)) return 1;
  _flagFitted = true;
  return 0;
}

int AAnam::fit(Db *db, const String& name)
{
  VectorDouble tab = db->getColumn(name,true);
  VectorDouble wt;
  if (db->hasWeight())
    wt = db->getColumnByLocator(ELoc::W,0,true);

  if (fitFromArray(tab, wt)) return 1;
  _flagFitted = true;
  return 0;
}

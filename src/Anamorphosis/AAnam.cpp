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

#include "Anamorphosis/AAnam.hpp"
#include "Db/Db.hpp"
#include "Basic/AException.hpp"
#include "Stats/Selectivity.hpp"

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
 ** \param[in]  power    Power of the change of support coefficient
 **
 *****************************************************************************/
double AAnam::calculateR(double cvv, double power)
{
  if (! hasChangeSupport()) return TEST;
  double s0, s1, s2, var0, var1;
  int converge, niter;
  static int niter_max = 1000;

  s1 = 0.;
  var1 = getBlockVariance(s1, power);

  /* Dichotomy */

  s2 = 1.;
  converge = niter = 0;
  while (!converge)
  {
    niter++;
    s0 = (s1 + s2) / 2.;
    var0 = getBlockVariance(s0, power);
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
double AAnam::getBlockVariance(double /*sval*/, double /*power*/) const
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
  int nclass = calest.getNClass();
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


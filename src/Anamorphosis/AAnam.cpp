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
#include "Anamorphosis/AAnam.hpp"
#include "Basic/AException.hpp"

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
 ** \param[in]  sval     Tentative Coefficient of change of support
 ** \param[in]  power    Power of the change of support coefficient
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

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
#include "Covariances/CovGradientFunctional.hpp"

#include "geoslib_f.h"
#include "geoslib_f_private.h"

#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include <math.h>

#define TR(i,j) (Tr[ndim * (i) + (j)])

CovGradientFunctional::CovGradientFunctional(const ENUM_COVS& type,
                                           const CovContext& ctxt)
    : ACovGradient(type, ctxt)
{
}

CovGradientFunctional::CovGradientFunctional(const CovGradientFunctional &r)
    : ACovGradient(r)
{
}

CovGradientFunctional& CovGradientFunctional::operator=(const CovGradientFunctional &r)
{
  if (this != &r)
  {
    ACovGradient::operator =(r);
  }
  return *this;
}

CovGradientFunctional::~CovGradientFunctional()
{
}

IClonable* CovGradientFunctional::clone() const
{
  return new CovGradientFunctional(*this);
}

/**
 * Calculate the square of the transformation matrix which transforms
 * a vector into its isotropic equivalent
 * @param d     Vector giving the distance in initial space
 * @param u     Vector of Gradient
 * @param trttr Matrix (Tr)^t %*% (Tr)
 * @remarks Limited to Space Dimension <= 3
 */
void CovGradientFunctional::_calculateTrTtr(const VectorDouble& d,
                                            VectorDouble& u,
                                            VectorDouble& trttr) const
{
  int ndim = getContext().getNDim();

  VectorDouble Tr(9, 0.);
  VectorDouble h(3,0.);

  ut_vector_fill(u, 0.);
  ut_vector_fill(trttr, 0.);

  // Matrix Tr = diag(coeffs) . R

  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
    {
      if (i >= ndim || j >= ndim) continue;
      TR(i,j) = getAniso().getRotation().getMatrixDirect().getValue(i,j) / getScale(i);
    }

  // Calculate the t(Tr) %*% Tr matrix

  int ecr = 0;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
    {
      double prod = 0.;
      for (int k=0; k<3; k++) prod += TR(k,i) * TR(k,j);
      trttr[ecr++] = prod;
    }

  // Rotate the Vector

  h[0] = Tr[0] * d[0] + Tr[1] * d[1] + Tr[2] * d[2];
  h[1] = Tr[3] * d[0] + Tr[4] * d[1] + Tr[5] * d[2];
  h[2] = Tr[6] * d[0] + Tr[7] * d[1] + Tr[8] * d[2];

  u[0] = Tr[0] * h[0] + Tr[3] * h[1] + Tr[6] * h[2];
  u[1] = Tr[1] * h[0] + Tr[4] * h[1] + Tr[7] * h[2];
  u[2] = Tr[2] * h[0] + Tr[5] * h[1] + Tr[8] * h[2];

  return;
 }

/**
 * Evaluates the covariance and gradient components
 * This function is restricted to the monovariate case
 * This function is limited to the only functions for which the
 * covariance(Point-Gradient) and covariance(Gradient-Gradient)
 * has been coded: i.e. Cubic or Gaussian,
 * in addition to the nugget effect.
 *
 * If 'flag_grad' == 0, then the output array 'covGG' is not filled.
 *
 * @param p1     First point of the Increment
 * @param p2     Second point of the increment
 * @param covVal Covariance value
 * @param covGp  Covariance <G[i](x0+x,y0+y,z0+z), P(x0,y0,z0)> (dim=3)
 * @param covGg  Covariance <G[i](x0+x,y0+y,z0+z), G[j](x0,y0,z0)> (dim=3)
 * @param mode   CovCalcMode structure
 * @param flagGrad true if the Gradient must be calculated
 *
 * @remarks: The returned arguments covVal, covGp and covGG are incremented here.
 * @remarks: They must have been initialized beforehand
 */

void CovGradientFunctional::evalZAndGradients(const SpacePoint& p1,
                                              const SpacePoint& p2,
                                              double& covVal,
                                              VectorDouble& covGp,
                                              VectorDouble& covGg,
                                              const CovCalcMode& mode,
                                              bool flagGrad) const
{
  VectorDouble trttr(9),u(3);

  // Calculate the isotropic distance

  double h = getSpace()->getDistance(p1, p2, getAniso());
  VectorDouble d1 = ut_vector_subtract(p1.getCoord(), p2.getCoord());

  //  Calculate the covariance

  double covar = getSill(0,0) * getCova()->evalCov(h);
  covVal += covar;
  if (getCova()->getType() == COV_NUGGET) return;

  _calculateTrTtr(d1, u, trttr);
  double dcovsr = getSill(0,0) * getCova()->evalCovDerivative(1,h);

  //  Case where distance is null

  if (h < EPSGRAD)
  {
    if (flagGrad)
    {
      for (int i = 0; i < 9; i++)
        covGg[i] -= dcovsr * trttr[i];
    }
  }
  else
  {

    //  Calculate covariance between point and gradient

    for (int i = 0; i < 3; i++)
    {
      covGp[i] += u[i] * dcovsr;
    }

    //  Calculate the covariance between gradient and gradient

    if (flagGrad)
    {
      double d2cov = getSill(0,0) * getCova()->evalCovDerivative(2,h);
      double a = (dcovsr - d2cov) / (h * h);
      if (getAniso().isIsotrop())
      {

        //  Isotropic case

        double b = dcovsr * trttr[0];
        int ecr = 0;
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
          {
            covGg[ecr] += a * u[i] * u[j];
            if (i == j) covGg[ecr] -= b;
            ecr++;
          }
      }
      else
      {

        //  Anisotropic case

        int ecr = 0;
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
          {
            covGg[ecr] += a * u[i] * u[j] - dcovsr * trttr[ecr];
            ecr++;
          }
      }
    }
  }
}

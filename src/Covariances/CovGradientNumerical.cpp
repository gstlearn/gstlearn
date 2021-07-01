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
#include "Covariances/CovGradientNumerical.hpp"

#include "geoslib_f.h"
#include "geoslib_f_private.h"

#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include <math.h>

#include "Covariances/CovFactory.hpp"

#define TR(i,j)                (Tr[(i) * 3 + (j)])

CovGradientNumerical::CovGradientNumerical(const ENUM_COVS& type,
                                           const CovContext& ctxt)
    : ACovGradient(type, ctxt)
{
}

CovGradientNumerical::CovGradientNumerical(const CovGradientNumerical &r)
    : ACovGradient(r)
{
}

CovGradientNumerical& CovGradientNumerical::operator=(const CovGradientNumerical &r)
{
  if (this != &r)
  {
    ACovGradient::operator =(r);
  }
  return *this;
}

CovGradientNumerical::~CovGradientNumerical()
{
}

IClonable* CovGradientNumerical::clone() const
{
  return new CovGradientNumerical(*this);
}

double CovGradientNumerical::_evalZZ(int ivar,
                                     int jvar,
                                     const SpacePoint& p1,
                                     const SpacePoint& p2,
                                     const CovCalcMode& mode) const
{
  return ACovGradient::eval(ivar,jvar,p1,p2,mode);
}

double CovGradientNumerical::_evalZGrad(int ivar,
                                        int jvar,
                                        int idim,
                                        const SpacePoint& p1,
                                        const SpacePoint& p2,
                                        const CovCalcMode& mode) const
{
  SpacePoint paux;
  double radius = getContext().getBallRadius();
  int ndim = getContext().getNDim();
  VectorDouble vec(ndim, 0);

  vec[idim] = radius / 2.;
  paux = p2;
  paux.move(vec);
  double covp0 = eval(ivar, jvar, p1, paux, mode);
  vec[idim] = -radius / 2.;
  paux = p2;
  paux.move(vec);
  double covm0 = eval(ivar, jvar, p1, paux, mode);

  double cov = (covm0 - covp0) / radius;
  return (cov);
}

double CovGradientNumerical::_evalGradGrad(int ivar,
                                           int jvar,
                                           int idim,
                                           int jdim,
                                           const SpacePoint& p1,
                                           const SpacePoint& p2,
                                           const CovCalcMode& mode) const
{
  SpacePoint paux;
  double radius = getContext().getBallRadius();
  int ndim = getContext().getNDim();
  VectorDouble vec(ndim, 0);

  vec[idim] = -radius / 2.;
  vec[jdim] =  radius / 2.;
  paux = p2;
  paux.move(vec);
  double covmp = eval(ivar, jvar, p1, paux, mode);
  vec[idim] = -radius / 2.;
  vec[jdim] = -radius / 2.;
  paux = p2;
  paux.move(vec);
  double covmm = eval(ivar, jvar, p1, paux, mode);
  vec[idim] =  radius / 2.;
  vec[jdim] = -radius / 2.;
  paux = p2;
  paux.move(vec);
  double covpm = eval(ivar, jvar, p1, paux, mode);
  vec[idim] = radius / 2.;
  vec[jdim] = radius / 2.;
  paux = p2;
  paux.move(vec);
  double covpp = eval(ivar, jvar, p1, paux, mode);

  double cov = (covmm + covpp - covmp - covpm) / (radius * radius);
  return(cov);
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
void CovGradientNumerical::evalZAndGradients(const SpacePoint& p1,
                                             const SpacePoint& p2,
                                             double& covVal,
                                             VectorDouble& covGp,
                                             VectorDouble& covGg,
                                             const CovCalcMode& mode,
                                             bool flagGrad) const
{
  //  Calculate the covariance

  double covar = _evalZZ(0,0,p1,p2,mode);
  covVal += covar;
  if (getCova()->getType() == COV_NUGGET) return;

  //  Calculate covariance between point and gradient

  for (int i = 0; i < 3; i++)
  {
    covGp[i] += _evalZGrad(0,0,i,p1,p2,mode);
  }

  //  Calculate the covariance between gradient and gradient

  if (flagGrad)
  {
    int ecr = 0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
      {
        covGg[ecr++] += _evalGradGrad(0, 0, i, j, p1, p2, mode);
      }
  }
}

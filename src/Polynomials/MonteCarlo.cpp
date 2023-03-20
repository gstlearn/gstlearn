/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Polynomials/Hermite.hpp"
#include "Polynomials/MonteCarlo.hpp"
#include "Basic/Law.hpp"

#include <math.h>

/**
 * Calculate: int phi(r*y + u * sqrt(1-r^2)) g(u) du
 *
 * @param yc Cutoff value
 * @param r  Change of support coefficient
 * @param psi Vector of Hermite coefficients
 * @return Vector of returned values for all Hermite coefficients
 */
double integralGaussHermite(double yc, double r, const VectorDouble &psi)
{
  int nbpoly = static_cast<int>(psi.size()) - 1;
  VectorDouble vect = hermitePolynomials(yc, 1., nbpoly);
  double value = hermiteSeries(vect, psi);
  return r * r * value;
}
void normalizeResults(int nbsimu, double &valest)
{
  valest /= (double) nbsimu;
}

void normalizeResults(int nbsimu, double &valest, double &valstd)
{
  valest /= (double) nbsimu;
  valstd = valstd / nbsimu - valest * valest;
  valstd = (valstd > 0.) ? sqrt(valstd) :
                           0.;
}

VectorDouble MCCondExp(VectorDouble krigest,
                       VectorDouble krigstd,
                       const VectorDouble &psi,
                       int nbsimu)
{
  VectorDouble condexp;

  int nech = static_cast<int>(krigest.size());
  condexp.resize(nech, 0.);

  for (int iech = 0; iech < nech; iech++)
  {
    double valest = 0.;
    for (int isimu = 0; isimu < nbsimu; isimu++)
    {
      double y = krigest[iech] + krigstd[iech] * law_gaussian();
      double z = hermiteCondExpElement(y, 0., psi);
      valest += z;
    }
    normalizeResults(nbsimu, valest);
    condexp[iech] = valest;
  }
  return condexp;
}

double MCCondExpElement(double krigest,
                        double krigstd,
                        const VectorDouble &psi,
                        int nbsimu)
{
  double valest = 0.;
  for (int isimu = 0; isimu < nbsimu; isimu++)
  {
    double y = krigest + krigstd * law_gaussian();
    double z = hermiteCondExpElement(y, 0., psi);
    valest += z;
  }
  normalizeResults(nbsimu, valest);
  return valest;
}

VectorDouble MCCondStd(VectorDouble krigest,
                       VectorDouble krigstd,
                       const VectorDouble &psi,
                       int nbsimu)
{
  VectorDouble condstd;

  int nech = static_cast<int>(krigest.size());
  condstd.resize(nech);

  for (int iech = 0; iech < nech; iech++)
  {
    double valest = 0.;
    double valstd = 0.;
    for (int isimu = 0; isimu < nbsimu; isimu++)
    {
      double y = krigest[iech] + krigstd[iech] * law_gaussian();
      double z = hermiteCondExpElement(y, 0., psi);
      valest += z;
      valstd += z * z;
    }
    normalizeResults(nbsimu, valest, valstd);
    condstd[iech] = valstd;
  }
  return condstd;
}

double MCCondStdElement(double krigest,
                        double krigstd,
                        const VectorDouble &psi,
                        int nbsimu)
{
  double valest = 0.;
  double valstd = 0.;
  for (int isimu = 0; isimu < nbsimu; isimu++)
  {
    double y = krigest + krigstd * law_gaussian();
    double z = hermiteCondExpElement(y, 0., psi);
    valest += z;
    valstd += z * z;
  }
  normalizeResults(nbsimu, valest, valstd);
  return valstd;
}

VectorDouble MCIndicator(double yc,
                         VectorDouble krigest,
                         VectorDouble krigstd,
                         int nbsimu)
{
  VectorDouble proba;

  int nech = static_cast<int>(krigest.size());
  proba.resize(nech, 0.);

  for (int iech = 0; iech < nech; iech++)
  {
    double valest = 0.;
    for (int isimu = 0; isimu < nbsimu; isimu++)
    {
      double y = krigest[iech] + krigstd[iech] * law_gaussian();
      if (y > yc)
      {
        valest += 1.;
      }
    }
    normalizeResults(nbsimu, valest);
    proba[iech] = valest;
  }
  return proba;
}

double MCIndicatorElement(double yc, double krigest, double krigstd, int nbsimu)
{
  double proba = 0.;
  for (int isimu = 0; isimu < nbsimu; isimu++)
  {
    double y = krigest + krigstd * law_gaussian();
    if (y > yc)
    {
      proba += 1.;
    }
  }
  normalizeResults(nbsimu, proba);
  return proba;
}

VectorDouble MCIndicatorStd(double yc,
                            VectorDouble krigest,
                            VectorDouble krigstd,
                            int nbsimu)
{
  VectorDouble probstd = MCIndicator(yc, krigest, krigstd, nbsimu);

  int nech = static_cast<int>(krigest.size());

  for (int iech = 0; iech < nech; iech++)
  {
    double proba = probstd[iech];
    probstd[iech] = sqrt(proba * (1. - proba));
  }
  return probstd;
}

double MCIndicatorStdElement(double yc,
                             double krigest,
                             double krigstd,
                             int nbsimu)
{
  double proba = MCIndicatorElement(yc, krigest, krigstd, nbsimu);
  double probstd = sqrt(proba * (1. - proba));
  return probstd;
}

VectorDouble MCMetal(double yc,
                     VectorDouble krigest,
                     VectorDouble krigstd,
                     const VectorDouble &psi,
                     int nbsimu)
{
  VectorDouble metal;

  int nech = static_cast<int>(krigest.size());
  metal.resize(nech, 0.);

  for (int iech = 0; iech < nech; iech++)
  {
    double valest = 0.;
    for (int isimu = 0; isimu < nbsimu; isimu++)
    {
      double y = krigest[iech] + krigstd[iech] * law_gaussian();
      if (y > yc)
      {
        valest += hermiteCondExpElement(y, 0., psi);
      }
    }
    normalizeResults(nbsimu, valest);
    metal[iech] = valest;
  }
  return metal;
}

double MCMetalElement(double yc,
                      double krigest,
                      double krigstd,
                      const VectorDouble &psi,
                      int nbsimu)
{
  double metal = 0.;
  for (int isimu = 0; isimu < nbsimu; isimu++)
  {
    double y = krigest + krigstd * law_gaussian();
    if (y > yc)
    {
      metal += hermiteCondExpElement(y, 0., psi);
    }
  }
  normalizeResults(nbsimu, metal);
  return metal;
}

VectorDouble MCMetalStd(double yc,
                        VectorDouble krigest,
                        VectorDouble krigstd,
                        const VectorDouble &psi,
                        int nbsimu)
{
  VectorDouble metstd;

  int nech = static_cast<int>(krigest.size());
  metstd.resize(nech, 0.);

  for (int iech = 0; iech < nech; iech++)
  {
    double valest = 0.;
    double valstd = 0.;
    for (int isimu = 0; isimu < nbsimu; isimu++)
    {
      double y = krigest[iech] + krigstd[iech] * law_gaussian();
      if (y > yc)
      {
        double z = hermiteCondExpElement(y, 0., psi);
        valest += z;
        valstd += z * z;
      }
    }
    normalizeResults(nbsimu, valest, valstd);
    metstd[iech] = valstd;
  }
  return metstd;
}

double MCMetalStdElement(double yc,
                         double krigest,
                         double krigstd,
                         const VectorDouble &psi,
                         int nbsimu)
{
  double metstd = 0.;
  double metest = 0.;
  for (int isimu = 0; isimu < nbsimu; isimu++)
  {
    double y = krigest + krigstd * law_gaussian();
    if (y > yc)
    {
      double z = hermiteCondExpElement(y, 0., psi);
      metest += z;
      metstd += z * z;
    }
  }
  normalizeResults(nbsimu, metest, metstd);
  return metstd;
}


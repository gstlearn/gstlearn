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
#include "geoslib_old_f.h"

#include "Basic/Vector.hpp"
#include "Db/Db.hpp"
#include "Stats/Classical.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"

#include <math.h>

/****************************************************************************/
/*!
 **  Calculate the quantile which corresponds to a given probability
 **
 ** \return  Quantile value
 **
 ** \param[in]  tab   Array of outcomes per sample
 ** \param[in]  ntab  Number of active values
 ** \param[in]  proba Probability value (between 0 and 1)
 **
 *****************************************************************************/
double _getQuantile(VectorDouble &tab, int ntab, double proba)
{
  int rank;
  double p1, p2, v1, v2, value;

  if (FFFF(proba)) return (TEST);

  ut_sort_double(0, ntab, NULL, tab.data());
  rank = (int) (proba * (double) ntab);
  v1 = tab[rank];

  if (rank < ntab - 1)
  {
    v2 = tab[rank + 1];
    p1 = (double) rank / (double) ntab;
    p2 = (double) (1 + rank) / (double) ntab;
    value = v1 + (proba - p1) * (v2 - v1) / (p2 - p1);
  }
  else
  {
    value = v1;
  }
  return (value);
}

VectorString statsNames(const std::vector<EStatOption>& opers)
{
  VectorString names;
  for (int i = 0; i < (int) opers.size(); i++)
  {
    EStatOption oper = opers[i];
    names.push_back(oper.getKey());
  }
  return names;
}

void dbStatisticsVariables(Db *db,
                           const VectorInt &iatts,
                           const std::vector<EStatOption>& opers,
                           int iptr0,
                           double vmin,
                           double vmax,
                           double proba)
{
  int noper = static_cast<int>(opers.size());
  if (noper <= 0) return;
  int natt = static_cast<int>(iatts.size());
  if (natt <= 0) return;

  /* Loop on the samples */

  VectorDouble local(natt);
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;

    /* Loop on the variables */

    int neff = 0;
    int nperc = 0;
    double mean = 0.;
    double var = 0.;
    double stdv = 0.;
    double sum = 0.;
    double metal = 0.;
    double mini = 1.e30;
    double maxi = -1.e30;
    for (int iatt = 0; iatt < natt; iatt++)
    {
      int jatt = iatts[iatt];
      double value = db->getArray(iech, jatt);
      if (FFFF(value)) continue;

      local[neff] = value;
      neff++;
      mean += value;
      sum += value;
      var += value * value;
      if (value < mini) mini = value;
      if (value > maxi) maxi = value;
      if (!FFFF(vmin) && value < vmin) continue;
      if (!FFFF(vmax) && value > vmax) continue;
      metal += value;
      nperc++;
    }

    // Normalization

    if (neff > 0)
    {
      mean /= neff;
      var = var / neff - mean * mean;
      stdv = (var >= 0) ? sqrt(var) : 0.;
    }

    /* Set the return array */

    for (int i = 0; i < noper; i++)
    {
      double tab = TEST;
      if (neff > 0)
      {
        if (opers[i] == EStatOption::NUM)
          tab = (double) neff;
        else if (opers[i] == EStatOption::MEAN)
          tab = mean;
        else if (opers[i] == EStatOption::VAR)
          tab = var;
        else if (opers[i] == EStatOption::STDV)
          tab = stdv;
        else if (opers[i] == EStatOption::MINI)
          tab = mini;
        else if (opers[i] == EStatOption::MAXI)
          tab = maxi;
        else if (opers[i] == EStatOption::SUM)
          tab = sum;
        else if (opers[i] == EStatOption::PROP)
          tab = (double) nperc / (double) neff;
        else if (opers[i] == EStatOption::QUANT)
          tab = _getQuantile(local, neff, proba);
        else if (opers[i] == EStatOption::T)
          tab = (double) nperc / (double) neff;
        else if (opers[i] == EStatOption::Q)
          tab = metal / (double) neff;
        else if (opers[i] == EStatOption::M)
          tab = (nperc > 0) ? metal / (double) nperc : TEST;
        else if (opers[i] == EStatOption::B)
          tab = (!FFFF(vmin)) ? (metal - vmin) / (double) neff : TEST;
        else
          return;
      }
      db->setArray(iech, iptr0 + i, tab);
    }
  }
}

/****************************************************************************/
/*!
 **  A Matrix containing the monovariate statistics for the set of variables
 **
 ** \return  The vector of statistics organized by variable
 **
 ** \param[in]  db         Db structure
 ** \param[in]  iatts      Vector of attribute ranks
 ** \param[in]  opers      List of the operator ranks
 ** \param[in]  flagIso    Restrain statistics to isotopic samples
 ** \param[in]  proba      Probability value (between 0 and 1)
 ** \param[in]  vmin       Minimum threshold
 ** \param[in]  vmax       Maximum threshold
 **
 *****************************************************************************/
VectorDouble dbStatisticsMono(Db *db,
                              const VectorInt &iatts,
                              const std::vector<EStatOption>& opers,
                              bool flagIso,
                              double proba,
                              double vmin,
                              double vmax)
{
  int noper = static_cast<int>(opers.size());
  int natt = static_cast<int>(iatts.size());
  int nech = db->getSampleNumber();

  // Find the Isotopic samples (optional)

  VectorDouble tab;
  VectorDouble local(nech, 0.);
  VectorInt accept(nech, 0);
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    accept[iech] = 0;
    if (! db->isActive(iech)) continue;
    accept[iech] = 1;
    int nundef = 0;
    for (int iatt = 0; iatt < natt; iatt++)
    {
      double value = db->getArray(iech, iatts[iatt]);
      if (FFFF(value)) nundef++;
    }
    if (flagIso && nundef > 0) accept[iech] = 0;
  }

  /* Loop on the attributes */

  for (int iatt = 0; iatt < natt; iatt++)
  {
    int neff = 0;
    int nperc = 0;
    double mean = 0.;
    double var = 0.;
    double stdv = 0.;
    double sum = 0.;
    double metal = 0.;
    double mini = 1.e30;
    double maxi = -1.e30;

    /* Loop on the samples */

    for (int iech = 0; iech < db->getSampleNumber(); iech++)
    {
      if (accept[iech] == 0) continue;
      double value = db->getArray(iech, iatts[iatt]);

      local[neff] = value;
      neff++;
      mean += value;
      sum += value;
      var += value * value;
      if (value < mini) mini = value;
      if (value > maxi) maxi = value;
      if (!FFFF(vmin) && value < vmin) continue;
      if (!FFFF(vmax) && value > vmax) continue;
      metal += value;
      nperc++;
    }

    // Normalization

    if (neff > 0)
    {
      mean /= neff;
      var = var / neff - mean * mean;
      stdv = (var >= 0) ? sqrt(var) : 0.;
    }

    // Constitute the array to be printed

    for (int i = 0; i < noper; i++)
    {
      if (opers[i] == EStatOption::NUM) tab.push_back((double) neff);
      if (neff > 0)
      {
        if (opers[i] == EStatOption::MEAN) tab.push_back(mean);
        if (opers[i] == EStatOption::VAR)  tab.push_back(var);
        if (opers[i] == EStatOption::STDV) tab.push_back(stdv);
        if (opers[i] == EStatOption::MINI) tab.push_back(mini);
        if (opers[i] == EStatOption::MAXI) tab.push_back(maxi);
        if (opers[i] == EStatOption::SUM)  tab.push_back(sum);
        if (opers[i] == EStatOption::PROP)
          tab.push_back((double) nperc / (double) neff);
        if (opers[i] == EStatOption::QUANT)
          tab.push_back(_getQuantile(local, neff, proba));
        if (opers[i] == EStatOption::T) tab.push_back((double) nperc / (double) neff);
        if (opers[i] == EStatOption::Q) tab.push_back(metal / (double) neff);
        if (opers[i] == EStatOption::M)
          tab.push_back((nperc > 0) ? metal / (double) nperc : TEST);
        if (opers[i] == EStatOption::B)
          tab.push_back((!FFFF(vmin)) ? (metal - vmin) / (double) neff : TEST);
      }
      else
      {
        if (opers[i] == EStatOption::MEAN) tab.push_back(TEST);
        if (opers[i] == EStatOption::VAR) tab.push_back(TEST);
        if (opers[i] == EStatOption::STDV) tab.push_back(TEST);
        if (opers[i] == EStatOption::MINI) tab.push_back(TEST);
        if (opers[i] == EStatOption::MAXI) tab.push_back(TEST);
        if (opers[i] == EStatOption::SUM) tab.push_back(TEST);
        if (opers[i] == EStatOption::PROP) tab.push_back(TEST);
        if (opers[i] == EStatOption::QUANT) tab.push_back(TEST);
        if (opers[i] == EStatOption::T) tab.push_back(TEST);
        if (opers[i] == EStatOption::Q) tab.push_back(TEST);
        if (opers[i] == EStatOption::M) tab.push_back(TEST);
        if (opers[i] == EStatOption::B) tab.push_back(TEST);
      }
    }
  }
  return tab;
}

/****************************************************************************/
/*!
 **  Considering that the Unique variable is a Facies (positive integer)
 **  returns the vector of proportions

 ** \return  The vector of proportions per Facies
 **
 ** \param[in]  db         Db structure
 **
 *****************************************************************************/
VectorDouble dbStatisticsFacies(Db *db)
{
  VectorDouble props;

  if (db->getLocatorNumber(ELoc::Z) != 1)
  {
    messerr(
        "This function requires the number of variables (%d) to be equal to 1",
        db->getLocatorNumber(ELoc::Z));
    return props;
  }
  int nech = db->getSampleNumber();
  int nfac = db->getFaciesNumber();

  // Calculate the proportions

  props.resize(nfac, 0.);
  int neff = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    if (!db->isActiveAndDefined(iech, 0)) continue;
    int ifac = (int) db->getVariable(iech, 0);
    if (ifac <= 0) continue;
    props[ifac - 1] += 1.;
    neff++;
  }

  // Normalization

  if (neff > 0)
  {
    for (int ifac = 0; ifac < nfac; ifac++)
      props[ifac] /= (double) neff;
  }
  return props;
}

/****************************************************************************/
/*!
 **  Considering that the Unique variable is an Indicator (0 or 1)
 **  returns the proportion of 1

 ** \return  The vector of proportions per Facies
 **
 ** \param[in]  db         Db structure
 **
 *****************************************************************************/
double dbStatisticsIndicator(Db *db)
{
  if (db->getLocatorNumber(ELoc::Z) != 1)
  {
    messerr(
        "This function requires the number of variables (%d) to be equal to 1",
        db->getLocatorNumber(ELoc::Z));
    return TEST;
  }

  // Calculate the proportions

  double prop = 0.;
  int neff = 0;
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActiveAndDefined(iech, 0)) continue;
    int ifac = (int) db->getVariable(iech, 0);
    if (ifac == 1) prop += 1.;
    neff++;
  }

  // Normalization

  if (neff > 0) prop = prop / neff;
  return prop;
}

/****************************************************************************/
/*!
 ** Calculates the correlation matrix for a set of variables of a Db
 **
 ** \return  Error Return code
 **
 ** \param[in]  db          Db structure
 ** \param[in]  iatts       Vector of attribute ranks
 ** \param[in]  flagIso    Restrain statistics to isotopic samples
 **
 *****************************************************************************/
VectorDouble dbStatisticsMulti(Db *db, const VectorInt &iatts, bool flagIso)
{
  int natt = static_cast<int>(iatts.size());

  /* Preliminary checks */

  if (natt <= 1)
    messerr("Correlation matrix will not be printed for a single variable");

  /* Core allocation */

  VectorDouble data(natt, 0.);
  VectorDouble mean(natt, 0.);
  VectorDouble var(natt, 0.);
  VectorDouble num(natt, 0.);
  VectorDouble cov(natt * natt, 0.);

  /* Loop on the samples */

  int numiso = 0;
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;

    /* Look for isotopic sample */

    int nundef = 0;
    for (int iatt = 0; iatt < natt; iatt++)
    {
      data[iatt] = db->getArray(iech, iatts[iatt]);
      if (FFFF(data[iatt])) nundef++;
    }
    if (flagIso && nundef > 0) continue;

    /* Calculate the 1-point statistics */

    for (int iatt = 0; iatt < natt; iatt++)
    {
      if (FFFF(data[iatt])) continue;
      num[iatt] += 1.;
      mean[iatt] += data[iatt];
      var[iatt] += data[iatt] * data[iatt];
    }

    if (nundef > 0) continue;
    numiso++;
    int ijatt = 0;
    for (int iatt = 0; iatt < natt; iatt++)
      for (int jatt = 0; jatt < natt; jatt++)
      {
        cov[ijatt] += data[iatt] * data[jatt];
        ijatt++;
      }
  }

  /* Normalization */

  for (int iatt = 0; iatt < natt; iatt++)
  {
    if (num[iatt] > 0)
    {
      mean[iatt] /= num[iatt];
      var[iatt] /= num[iatt];
      var[iatt] -= mean[iatt] * mean[iatt];
      if (var[iatt] <= 0) var[iatt] = 0.;
    }
  }
  if (numiso > 0)
  {
    int ijatt = 0;
    for (int iatt = 0; iatt < natt; iatt++)
      for (int jatt = 0; jatt < natt; jatt++)
      {
        cov[ijatt] /= numiso;
        cov[ijatt] -= mean[iatt] * mean[jatt];
        cov[ijatt] /= sqrt(var[iatt] * var[jatt]);
        ijatt++;
      }
  }
  return cov;
}

/****************************************************************************/
/*!
 **  Print the monovariate statistics between different variables of a Db
 **
 ** \return  Error Return code
 **
 ** \param[in]  stats       Array of statistics (organized by variable)
 ** \param[in]  opers       List of the operator ranks
 ** \param[in]  varnames    List of variables
 ** \param[in]  title       Title for the printout (optional)
 **
 *****************************************************************************/
String statisticsMonoPrint(const VectorDouble &stats,
                           const std::vector<EStatOption>& opers,
                           const VectorString &varnames,
                           const String &title)
{
  int noper = static_cast<int>(opers.size());
  int natt = static_cast<int>(varnames.size());
  std::stringstream sstr;

  // Constitute the vector of row and column names
  VectorString colnames = statsNames(opers);

  // Printout the matrix
  sstr << toMatrix(title, colnames, varnames, false, noper, natt, stats, true);

  return sstr.str();
}

/****************************************************************************/
/*!
 **  Print the multivariate statistics between different variables of a Db
 **
 ** \return  Error Return code
 **
 ** \param[in]  stats       Matrix of variance-covariance
 ** \param[in]  varnames    Vector of variable names
 ** \param[in]  title       Title for the printout (optional)
 **
 *****************************************************************************/
String statisticsMultiPrint(const VectorDouble &stats,
                                            const VectorString &varnames,
                                            const String &title)
{
  int natt = static_cast<int>(varnames.size());
  std::stringstream sstr;

  sstr
      << toMatrix(title, VectorString(), VectorString(), true, natt, natt,
                  stats, true);

  return sstr.str();
}

bool regressionLoad(Db *db1,
                    Db *db2,
                    int iech,
                    int icol0,
                    const VectorInt &icols,
                    int mode,
                    int flagCste,
                    double *value,
                    VectorDouble &x)
{
  int ecr = 0;
  switch (mode)
  {
    case 0:
      *value = db1->getArray(iech, icol0);
      if (flagCste) x[ecr++] = 1.;
      for (int icol = 0; icol < (int) icols.size(); icol++)
        x[ecr++] = db2->getArray(iech, icols[icol]);
      break;

    case 1:
      int nfex = db2->getExternalDriftNumber();
      *value = db1->getVariable(iech, 0);
      if (flagCste) x[ecr++] = 1.;
      for (int i = 0; i < nfex; i++)
        x[ecr++] = db2->getExternalDrift(iech, i);
      break;
  }

  bool flagTest = false;
  for (int i = 0; i < (int) x.size() && !flagTest; i++)
    flagTest = FFFF(x[i]);
  return (FFFF(*value) || flagTest);
}

bool regressionCheck(Db *db1,
                     Db *db2,
                     int mode,
                     int icol0,
                     const VectorInt &icols)
{
  int ncol = (int) icols.size();
  int nfex = db2->getExternalDriftNumber();

  if (db1->getVariableNumber() != 1)
  {
    messerr("This method is restricted to the Monovariate case");
    return false;
  }

  switch (mode)
  {
    case 0:
      if (icol0 < 0 || icol0 >= db1->getColumnNumber())
      {
        messerr("The regression requires a valid target variable");
        return false;
      }
      for (int icol = 0; icol < ncol; icol++)
      {
        if (icols[icol] < 0 || icols[icol] >= db2->getColumnNumber())
        {
          messerr("The regression requires a valid auxiliary variable (#%d)",
                  icol + 1);
          return false;
        }
      }
      break;

    case 1:
      if (nfex <= 0)
      {
        messerr("The multivariate regression is designated");
        messerr("as a function of several drift variables");
        messerr("The Db contains %d drift variables", nfex);
        return false;
      }
      break;
  }
  return true;
}

/****************************************************************************/
/*!
 **  Evaluate the regression
 **
 ** \return  Error return code
 **
 ** \param[in,out]  db1        Db descriptor (for target variable)
 ** \param[in]  db2            Db descriptor (for auxiliary variables)
 ** \param[in]  mode           Type of calculation
 ** \li                        0 : standard multivariate case
 ** \li                        1 : using external drifts
 ** \param[in]  icol0          Rank of the target variable
 ** \param[in]  icols          Vector of ranks of the explanatory variables
 ** \param[in]  flagCste       The constant is added as explanatory variable
 ** \param[in]  verbose        Verbose option
 **
 ** \remark  The flag_mode indicates the type of regression calculation:
 ** \remark  0 : V[icol] as a function of V[icols[i]]
 ** \remark  1 : Z1 as a function of the different Fi's
 **
 *****************************************************************************/
ResRegr regression(Db *db1,
                   Db *db2,
                   int mode,
                   int icol0,
                   const VectorInt& icols,
                   bool flagCste,
                   bool verbose)
{
  ResRegr regr;

  int nfex = db2->getExternalDriftNumber();
  int nech = db1->getSampleNumber();
  int ncol = (int) icols.size();
  int size = 0;
  switch (mode)
  {
    case 0:
      size = ncol;
      if (flagCste) size++;
      break;
    case 1:
      size = nfex;
      if (flagCste) size++;
      break;
  }

  /* Preliminary checks */

  if (! regressionCheck(db1, db2, mode, icol0, icols)) return regr;

  /* Core allocation */

  VectorDouble x(size,0.);
  VectorDouble b(size,0.);
  MatrixSquareSymmetric a(size);

  /* Loop on the samples */

  int number = 0;
  double prod = 0.;
  double mean = 0.;
  double value = 0.;

  for (int iech=0; iech < nech; iech++)
  {
    if (! db1->isActive(iech)) continue;

    /* Get the information for the current sample */

    if (regressionLoad(db1, db2, iech, icol0, icols, mode, flagCste, &value, x))
      continue;

    prod += value * value;
    mean += value;
    number++;

    /* Update the matrices */

    for (int i = 0; i < size; i++)
    {
      b[i] += value * x[i];
      for (int j=0; j<=i; j++)
        a.setValue(i, j, a.getValue(i,j) + x[i] * x[j]);
    }
  }

  if (number <= 0)
  {
    messerr("No sample found where variables are defined");
    return regr;
  }

  /* Solve the regression system */

  int pivot = a.solve(b, x);
  if (pivot > 0)
  {
    messerr("Error during regression calculation: pivot %d is null", pivot);
    return regr;
  }

  // Normalization
  mean /= number;
  regr.count = number;
  regr.nvar = ncol;
  regr.flagCste = flagCste;
  regr.coeffs = x;
  regr.variance = prod / number - mean * mean;

  /* Calculate the residuals */

  for (int i = 0; i < size; i++)
  {
    prod -= 2. * x[i] * b[i];
    for (int j = 0; j < size; j++)
      prod += x[i] * x[j] * a.getValue(i,j);
  }
  regr.varres = prod / number;

  /************/
  /* Printout */
  /************/

  if (verbose)
  {
    mestitle(1, "Linear regression:");
    for (int i = 0; i < size; i++)
      message("Explanatory variable Aux.#%d - Coefficient = %lf\n", i + 1, x[i]);
    message("Variance of Residuals = %lf\n", regr.varres);
  }
  return regr;
}

void regrprint(const ResRegr& regr)
{
  mestitle(1, "Linear Regression");
  message("- Calculated on %d active values\n",regr.count);

  int ecr = 0;
  if (regr.flagCste)
    message("- Constant term           = %lf\n",regr.coeffs[ecr++]);
  for (int ivar = 0; ivar < regr.nvar; ivar++)
    message("- Explanatory Variable #%d = %lf\n", ivar+1, regr.coeffs[ecr++]);

  message("- Initial variance        = %lf\n",regr.variance);
  message("- Variance of residuals   = %lf\n",regr.varres);
}

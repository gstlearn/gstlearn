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
#include "Db/DbGrid.hpp"
#include "Stats/Classical.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "Basic/AStringable.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Model/Model.hpp"

#include <math.h>
#include <string.h>

/****************************************************************************/
/*!
 **  Update the proportions
 **
 ** \param[in]  dbin       Db for the input grid
 ** \param[in]  indg       Array of grid indices
 ** \param[in]  nfacies    Number of facies
 **
 ** \param[out] prop       Array of proportions
 **
 *****************************************************************************/
void _updateProportions(DbGrid *dbin,
                        VectorInt &indg,
                        int nfacies,
                        VectorDouble &prop)
{
  int rank = dbin->getGrid().indiceToRank(indg);
  int ifac = (int) dbin->getVariable(rank, 0);
  if (ifac < 1 || ifac > nfacies) return;
  prop[ifac - 1] += 1.;
}

/****************************************************************************/
/*!
 **  Update the transitions
 **
 ** \param[in]  dbin       DbGrid for the input grid
 ** \param[in]  pos        Rank of the montee axis
 ** \param[in]  indg       Array of grid indices
 ** \param[in]  nfacies    Number of facies
 ** \param[in]  orient     Orientation
 **
 ** \param[out] trans      Array of transitions
 **
 *****************************************************************************/
void _updateTransition(DbGrid *dbin,
                       int pos,
                       VectorInt& indg,
                       int nfacies,
                       int orient,
                       VectorDouble& trans)
{
  int jpos = indg[pos] + orient;
  if (jpos <= 0 || jpos >= dbin->getNX(pos)) return;
  int ifac1 = (int) dbin->getVariable(dbin->getGrid().indiceToRank(indg), 0);
  indg[pos] += orient;
  int ifac2 = (int) dbin->getVariable(dbin->getGrid().indiceToRank(indg), 0);
  indg[pos] -= orient;

  if (ifac1 < 1 || ifac1 > nfacies || ifac2 < 1 || ifac2 > nfacies) return;
  trans[(ifac1 - 1) * nfacies + (ifac2 - 1)] += 1.;
}

/****************************************************************************/
/*!
 **  Scale the proportions and store the proportions
 **
 ** \param[in]  dbout      Db for the output grid
 ** \param[in]  iptr       Writing pointer
 ** \param[in]  iech       Rank of the target sample
 ** \param[in]  nitem      Number of items
 ** \param[in]  tab        Array of cumulative statistics
 **
 *****************************************************************************/
void _scaleAndAffect(Db *dbout,
                     int iptr,
                     int iech,
                     int nitem,
                     VectorDouble &tab)
{
  double value;

  double total = 0.;
  for (int ifac = 0; ifac < nitem; ifac++)
    total += tab[ifac];

  for (int ifac = 0; ifac < nitem; ifac++)
  {
    if (total <= 0.)
      value = TEST;
    else
      value = tab[ifac] / total;
    dbout->setArray(iech, iptr + ifac, value);
  }
}

/****************************************************************************/
/*!
 **  Check the operator name is mentioned within a list
 **
 ** \return  1 if the operator is mentioned; 0 otherwise
 **
 ** \param[in]  opers Array of operators
 ** \param[in]  refe  Reference operator
 **
 ** \remarks If the array 'opers' if empty, any name is considered as valid
 **
 *****************************************************************************/
bool _operExists(const std::vector<EStatOption>& opers,
                 const EStatOption& refe)
{
  int noper = static_cast<int>(opers.size());
  if (noper == 0) return (1);
  for (int i = 0; i < noper; i++)
  {
    if (opers[i] == refe) return (1);
  }
  return (0);
}

bool _regressionCheck(Db *db1,
                      int icol0,
                      const VectorInt &icols,
                      int mode,
                      Db *db2,
                      const Model *model)
{
  int ncol = (int) icols.size();
  int nfex = db2->getExternalDriftNumber();

  switch (mode)
  {
    case 0:
      if (icol0 < 0 || icol0 >= db1->getUIDMaxNumber())
      {
        messerr("The regression requires a valid target variable");
        return false;
      }
      for (int icol = 0; icol < ncol; icol++)
      {
        if (icols[icol] < 0 || icols[icol] >= db2->getUIDMaxNumber())
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

    case 2:
      if (model == nullptr)
      {
        messerr("Model should be defined");
        return false;
      }
      if (model->getDriftNumber() <= 0)
      {
        messerr("The number of Drift equations in the Model should be positive");
        return false;
      }
  }
  return true;
}

bool _regressionLoad(Db *db1,
                     Db *db2,
                     int iech,
                     int icol0,
                     const VectorInt &icols,
                     int mode,
                     int flagCste,
                     const Model *model,
                     double *value,
                     VectorDouble &x)
{
  int nfex = 0;
  int nbfl = 0;

  int ecr  = 0;
  switch (mode)
  {
    case 0:
      *value = db1->getArray(iech, icol0);
      if (flagCste) x[ecr++] = 1.;
      for (int icol = 0; icol < (int) icols.size(); icol++)
        x[ecr++] = db2->getArray(iech, icols[icol]);
      break;

    case 1:
      nfex = db2->getExternalDriftNumber();
      *value = db1->getVariable(iech, 0);
      if (flagCste) x[ecr++] = 1.;
      for (int i = 0; i < nfex; i++)
        x[ecr++] = db2->getExternalDrift(iech, i);
      break;

    case 2:
      nbfl = model->getDriftNumber();
      *value = db1->getVariable(iech, 0);
      for (int i = 0; i < nbfl; i++)
         x[ecr++] = model->evalDrift(db2, iech, i, ECalcMember::LHS);
      break;
  }

  bool flagTest = false;
  for (int i = 0; i < (int) x.size() && !flagTest; i++)
    flagTest = FFFF(x[i]);
  return (FFFF(*value) || flagTest);
}

void _regrprint(const ResRegr& regr)
{
  mestitle(1, "Linear Regression");
  message("- Calculated on %d active values\n",regr.count);

  int ecr = 0;
  int nvar = regr.nvar;
  if (regr.flagCste) nvar--;

  if (regr.flagCste)
    message("- Constant term           = %lf\n",regr.coeffs[ecr++]);
  for (int ivar = 0; ivar < nvar; ivar++)
    message("- Explanatory Variable #%d = %lf\n", ivar+1, regr.coeffs[ecr++]);

  message("- Initial variance        = %lf\n",regr.variance);
  message("- Variance of residuals   = %lf\n",regr.varres);
}

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
                           double proba,
                           double vmin,
                           double vmax)
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

VectorDouble dbStatisticsMono(Db *db,
                              const VectorString& names,
                              const std::vector<EStatOption> &opers,
                              bool flagIso,
                              double proba,
                              double vmin,
                              double vmax)
{
  VectorInt iatts = db->getUIDs(names);
  return dbStatisticsMono(db, iatts, opers, flagIso, proba, vmin, vmax);
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
        if (opers[i] == EStatOption::MEAN)  tab.push_back(TEST);
        if (opers[i] == EStatOption::VAR)   tab.push_back(TEST);
        if (opers[i] == EStatOption::STDV)  tab.push_back(TEST);
        if (opers[i] == EStatOption::MINI)  tab.push_back(TEST);
        if (opers[i] == EStatOption::MAXI)  tab.push_back(TEST);
        if (opers[i] == EStatOption::SUM)   tab.push_back(TEST);
        if (opers[i] == EStatOption::PROP)  tab.push_back(TEST);
        if (opers[i] == EStatOption::QUANT) tab.push_back(TEST);
        if (opers[i] == EStatOption::T)     tab.push_back(TEST);
        if (opers[i] == EStatOption::Q)     tab.push_back(TEST);
        if (opers[i] == EStatOption::M)     tab.push_back(TEST);
        if (opers[i] == EStatOption::B)     tab.push_back(TEST);
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

VectorDouble dbStatisticsMulti(Db *db, const VectorString &names, bool flagIso)
{
  VectorInt iatts = db->getUIDs(names);
  return dbStatisticsMulti(db, iatts, flagIso);
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
 ** \param[in]  names       List of variables
 ** \param[in]  title       Title for the printout (optional)
 **
 *****************************************************************************/
String statisticsMonoPrint(const VectorDouble &stats,
                           const std::vector<EStatOption>& opers,
                           const VectorString &names,
                           const String &title)
{
  int noper = static_cast<int>(opers.size());
  int natt = static_cast<int>(names.size());
  std::stringstream sstr;

  // Constitute the vector of row and column names
  VectorString colnames = statsNames(opers);

  // Printout the matrix
  sstr << toMatrix(title, colnames, names, false, noper, natt, stats, true);

  return sstr.str();
}

/****************************************************************************/
/*!
 **  Print the multivariate statistics between different variables of a Db
 **
 ** \return  Error Return code
 **
 ** \param[in]  stats       Matrix of variance-covariance
 ** \param[in]  names       Vector of variable names
 ** \param[in]  title       Title for the printout (optional)
 **
 *****************************************************************************/
String statisticsMultiPrint(const VectorDouble &stats,
                            const VectorString &names,
                            const String &title)
{
  int natt = static_cast<int>(names.size());
  std::stringstream sstr;

  sstr << toMatrix(title, VectorString(), VectorString(), true, natt, natt,
                   stats, true);

  return sstr.str();
}

ResRegr regression(Db *db1,
                   const String &name0,
                   const VectorString &names,
                   int mode,
                   bool flagCste,
                   Db *db2,
                   const Model *model,
                   bool verbose)
{
  if (db1 == nullptr) return ResRegr();
  if (db2 == nullptr) db2 = db1;

  int icol0 = db1->getUID(name0);
  VectorInt icols;
  if (! names.empty()) icols = db2->getUIDs(names);
  return regressionByUID(db1, icol0, icols, mode, flagCste, db2, model, verbose);
}

/****************************************************************************/
/*!
 **  Evaluate the regression
 **
 ** \return  Error return code
 **
 ** \param[in,out]  db1        Db descriptor (for target variable)

 ** \param[in]  icol0          Rank of the target variable
 ** \param[in]  icols          Vector of ranks of the explanatory variables
 ** \param[in]  mode           Type of calculation
 ** \li                        0 : standard multivariate case
 ** \li                        1 : using external drifts
 ** \li                        2 : using standard drift functions (in 'model')
 ** \param[in]  flagCste       The constant is added as explanatory variable
 ** \param[in]  db2            Db descriptor (for auxiliary variables)
 ** \param[in]  model          Model (only used for Drift functions if mode==2)
 ** \param[in]  verbose        Verbose option
 **
 ** \remark  The flag_mode indicates the type of regression calculation:
 ** \remark  0 : V[icol] as a function of V[icols[i]]
 ** \remark  1 : Z1 as a function of the different Fi's
 **
 *****************************************************************************/
ResRegr regressionByUID(Db *db1,
                        int icol0,
                        const VectorInt &icols,
                        int mode,
                        bool flagCste,
                        Db *db2,
                        const Model *model,
                        bool verbose)
{
  ResRegr regr;

  if (db1 == nullptr) return regr;
  if (db2 == nullptr) db2 = db1;

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
    case 2:
      size = model->getDriftNumber();
      break;
  }

  /* Preliminary checks */

  if (! _regressionCheck(db1, icol0, icols, mode, db2, model)) return regr;

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

    if (_regressionLoad(db1, db2, iech, icol0, icols, mode, flagCste, model,
                        &value, x)) continue;

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
  regr.nvar = size;
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
    _regrprint(regr);
  }
  return regr;
}

/****************************************************************************/
/*!
 **  Evaluate the regression
 **
 ** \return  Error return code
 **
 ** \param[in,out]  db1        Db descriptor (for target variable)
 ** \param[in]  iptr0          Storing address
 ** \param[in]  name0          Name of the target variable
 ** \param[in]  names          Vector of names of the explanatory variables
 ** \param[in]  mode           Type of calculation
 ** \li                        0 : standard multivariate case
 ** \li                        1 : using external drifts
 ** \li                        2 : using standard drift functions (mode==2)
 ** \param[in]  flagCste       The constant is added as explanatory variable]
 ** \param[in]  db2            Db descriptor (for auxiliary variables)
 ** \param[in]  model          Model structure (used for mode==2)
 **
 ** \remark  The flag_mode indicates the type of regression calculation:
 ** \remark  0 : V[icol] as a function of V[icols[i]]
 ** \remark  1 : Z1 as a function of the different Fi's
 **
 ** \remark  The Db1 structure is modified: the column (iptr0) of the Db1
 ** \remark  is added by this function; it contains the value
 ** \remark  of the residuals at each datum (or TEST if the residual has not
 ** \remark  been calculated).
 **
 *****************************************************************************/
int regressionApply(Db *db1,
                    int iptr0,
                    const String& name0,
                    const VectorString& names,
                    int mode,
                    bool flagCste,
                    Db *db2,
                    const Model* model)
{
  ResRegr regr;
  if (db2 == nullptr) db2 = db1;
  int icol0 = db1->getUID(name0);
  VectorInt icols;
  if (! names.empty()) icols = db2->getUIDs(names);

  regr = regressionByUID(db1, icol0, icols, mode, flagCste, db2, model);

  /* Preliminary checks */

  if (! _regressionCheck(db1, icol0, icols, mode, db2, model)) return 1;

  /* Store the regression error at sample points */

  int size = (int) regr.coeffs.size();
  double value = 0;
  VectorDouble x(size);

  for (int iech = 0; iech < db1->getSampleNumber(); iech++)
  {
    if (db1->isActive(iech))
    {
      /* Get the information for the current sample */

      if (_regressionLoad(db1, db2, iech, icol0, icols, mode, flagCste, model,
                          &value, x))
      {
        value = TEST;
      }
      else
      {
        for (int i = 0; i < size; i++)
          value -= x[i] * regr.coeffs[i];
      }
    }
    db1->setArray(iech, iptr0, value);
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Calculates the "montee" from a grid into a 1-D grid
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       Db for the input grid
 ** \param[in]  dbout      Db for the output grid
 ** \param[in]  pos        Rank of the montee axis (starting from 0)
 ** \param[in]  nfacies    Number of facies
 ** \param[in]  radius     Radius of the neighborhood
 **
 *****************************************************************************/
int statisticsProportion(DbGrid *dbin,
                         DbGrid *dbout,
                         int pos,
                         int nfacies,
                         int radius)
{
  int ndim = dbin->getNDim();
  if (ndim != 2 && ndim != 3)
  {
    messerr("This function is limited to 2-D or 3-D input grids");
    return 1;
  }
  if (pos < 0 || pos >= ndim)
  {
    messerr("The rank of the 'montee' axis should lie between 1 and %d", ndim);
    return 1;
  }
  if (dbin->getNX(pos) != dbout->getNX(0) ||
      dbin->getX0(pos) != dbout->getX0(0) ||
      dbin->getDX(pos) != dbout->getDX(0))
  {
    messerr("The 1-D output grid does not match input grid");
    return 1;
  }
  if (!dbin->isVariableNumberComparedTo(1)) return 1;

  /* Core allocation */

  int ngrid = dbin->getNX(pos);
  VectorDouble prop(nfacies, 0.);
  VectorInt indg(ndim);

  /* Create the new variables in the output file */

  int iptr = dbout->addColumnsByConstant(nfacies, TEST);
  if (iptr < 0) return 1;

  /* Loop on the elements of the output grid */

  for (int iech = 0; iech < ngrid; iech++)
  {
    for (int ifac = 0; ifac < nfacies; ifac++) prop[ifac] = 0.;

    if (ndim == 2)
    {
      int aux1 = (pos + 1) % ndim;
      for (int ishift = -radius; ishift <= radius; ishift++)
      {
        indg[pos] = iech + ishift;
        if (indg[pos] < 0 || indg[pos] >= ngrid) continue;
        for (int i1 = 0; i1 < dbin->getNX(aux1); i1++)
        {
          indg[aux1] = i1;
          _updateProportions(dbin, indg, nfacies, prop);
        }
        _scaleAndAffect(dbout, iptr, iech, nfacies, prop);
      }
    }
    else
    {
      int bux1 = (pos + 1) % ndim;
      int bux2 = (pos + 2) % ndim;
      int aux1 = MIN(bux1, bux2);
      int aux2 = MAX(bux1, bux2);
      for (int ishift = -radius; ishift <= radius; ishift++)
      {
        indg[pos] = iech + ishift;
        if (indg[pos] < 0 || indg[pos] >= ngrid) continue;
        for (int i1 = 0; i1 < dbin->getNX(aux1); i1++)
          for (int i2 = 0; i2 < dbin->getNX(aux2); i2++)
          {
            indg[aux1] = i1;
            indg[aux2] = i2;
            _updateProportions(dbin, indg, nfacies, prop);
          }
        _scaleAndAffect(dbout, iptr, iech, nfacies, prop);
      }
    }
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Calculates the transition from a grid into a 1-D grid
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       Db for the input grid
 ** \param[in]  dbout      Db for the output grid
 ** \param[in]  pos        Rank of the montee axis (starting from 0)
 ** \param[in]  nfacies    Number of facies
 ** \param[in]  radius     Radius of the neighborhood
 ** \param[in]  orient     Orientation (+1 or -1)
 **
 *****************************************************************************/
int statisticsTransition(DbGrid *dbin,
                         DbGrid *dbout,
                         int pos,
                         int nfacies,
                         int radius,
                         int orient)
{
  int ndim = dbin->getNDim();
  if (ndim != 2 && ndim != 3)
  {
    messerr("This function is limited to 2-D or 3-D input grids");
    return 1;
  }
  if (pos < 0 || pos >= ndim)
  {
    messerr("The rank of the 'montee' axis should lie between 1 and %d", ndim);
    return 1;
  }
  if (dbin->getNX(pos) != dbout->getNX(0) ||
      dbin->getX0(pos) != dbout->getX0(0) ||
      dbin->getDX(pos) != dbout->getDX(0))
  {
    messerr("The 1-D output grid does not match input grid");
    return 1;
  }
  if (!dbin->isVariableNumberComparedTo(1)) return 1;

  /* Core allocation */

  int ngrid = dbin->getNX(pos);
  int nitem = nfacies * nfacies;
  VectorDouble trans(nitem, 0.);
  VectorInt indg(ndim);

  /* Create the new variables in the output file */

  int iptr = dbout->addColumnsByConstant(nfacies * nfacies, TEST);
  if (iptr < 0) return 1;

  /* Loop on the elements of the output grid */

  for (int iech = 0; iech < ngrid; iech++)
  {
    indg[pos] = iech;
    for (int item = 0; item < nitem; item++) trans[item] = 0.;

    if (ndim == 2)
    {
      int aux1 = (pos + 1) % ndim;
      for (int ishift = -radius; ishift <= radius; ishift++)
      {
        indg[pos] = iech + ishift;
        if (indg[pos] < 0 || indg[pos] >= ngrid) continue;
        for (int i1 = 0; i1 < dbin->getNX(aux1); i1++)
        {
          indg[aux1] = i1;
          _updateTransition(dbin, pos, indg, nfacies, orient, trans);
        }
        _scaleAndAffect(dbout, iptr, iech, nitem, trans);
      }
    }
    else
    {
      int bux1 = (pos + 1) % ndim;
      int bux2 = (pos + 2) % ndim;
      int aux1 = MIN(bux1, bux2);
      int aux2 = MAX(bux1, bux2);
      for (int ishift = -radius; ishift <= radius; ishift++)
      {
        indg[pos] = iech + ishift;
        if (indg[pos] < 0 || indg[pos] >= ngrid) continue;
        for (int i1 = 0; i1 < dbin->getNX(aux1); i1++)
          for (int i2 = 0; i2 < dbin->getNX(aux2); i2++)
          {
            indg[aux1] = i1;
            indg[aux2] = i2;
            _updateTransition(dbin, pos, indg, nfacies, orient, trans);
          }
        _scaleAndAffect(dbout, iptr, iech, nitem, trans);
      }
    }
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Constitute the name of the row
 **
 ** \param[in]  radix       Radix for the different variables (optional)
 ** \param[in]  ncol        Number of variables
 ** \param[in]  icol        Rank of the variable
 ** \param[in]  name        Variables name
 ** \param[in]  string      String array
 **
 *****************************************************************************/
void _getRowname(const String &radix,
                 int ncol,
                 int icol,
                 const String &name,
                 char *string)
{
  if (!radix.empty())
    (void) gslSPrintf(string, "%s-%d", radix.c_str(), icol + 1);
  else if (!name.empty())
    (void) gslSPrintf(string, "%s", name.c_str());
  else if (ncol > 1)
    (void) gslSPrintf(string, "Variable-%d", icol + 1);
  else
    (void) gslSPrintf(string, "Variable");
}

void dbStatisticsPrint(const Db *db,
                    const VectorString &names,
                    const std::vector<EStatOption>& opers,
                    bool flagIso,
                    bool flagCorrel,
                    const String &title,
                    const String &radix)
{
  VectorInt iatts = db->getUIDs(names);
  if (iatts.size() <= 0) return;
  dbStatisticsPrint(db, iatts, opers, flagIso, flagCorrel, title, radix);
}

/****************************************************************************/
/*!
 **  Print the multivariate statistics between different variables of a Db
 **
 ** \param[in]  db          Db structure
 ** \param[in]  iatts_arg   Ranks of the attributes (empty = all)
 ** \param[in]  opers       Array of operators
 ** \param[in]  flagIso     Restrain statistics to isotopic samples
 ** \param[in]  flagCorrel  True if the correlations must be calculated
 ** \param[in]  title       Title for the printout (optional)
 ** \param[in]  radix       Radix for the different variables (optional)

 **
 *****************************************************************************/
void dbStatisticsPrint(const Db *db,
                       const VectorInt &iatts_arg,
                       const std::vector<EStatOption>& opers,
                       bool flagIso,
                       bool flagCorrel,
                       const String &title,
                       const String &radix)
{
  char string[50];

  /* Initializations */

  VectorInt iatts = iatts_arg;
  if (iatts.empty()) iatts = db->getAllUIDs();
  int ncol = static_cast<int>(iatts.size());

  /* Preliminary checks */

  if (flagCorrel && ncol <= 1)
    messerr("Correlation matrix will not be printed for a single variable");

  /* Core allocation */

  VectorDouble data(ncol, 0.);
  VectorDouble mean(ncol, 0.);
  VectorDouble mini(ncol, 0.);
  VectorDouble maxi(ncol, 0.);
  VectorDouble var(ncol,  0.);
  VectorDouble num(ncol,  0.);
  VectorDouble cov;
  if (flagCorrel) cov.resize(ncol * ncol, 0.);

  /* Initializations */

  int numiso = 0;
  int ijcol = 0;
  for (int icol = 0; icol < ncol; icol++)
  {
    mean[icol] = var[icol] = num[icol] = 0.;
    mini[icol] = 1.e30;
    maxi[icol] = -1.e30;
    if (flagCorrel)
      for (int jcol = 0; jcol < ncol; jcol++, ijcol++)
        cov[ijcol] = 0.;
  }

  /* Loop on the samples */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;

    /* Look for isotopic sample */

    int nundef = 0;
    for (int icol = 0; icol < ncol; icol++)
    {
      data[icol] = db->getArray(iech, iatts[icol]);
      if (FFFF(data[icol])) nundef++;
    }
    if (flagIso && nundef > 0) continue;

    /* Calculate the 1-point statistics */

    for (int icol = 0; icol < ncol; icol++)
    {
      if (FFFF(data[icol])) continue;
      num[icol] += 1.;
      mean[icol] += data[icol];
      var[icol] += data[icol] * data[icol];
      if (data[icol] < mini[icol]) mini[icol] = data[icol];
      if (data[icol] > maxi[icol]) maxi[icol] = data[icol];
    }

    if (nundef > 0) continue;
    numiso++;
    if (flagCorrel)
      for (int icol = ijcol = 0; icol < ncol; icol++)
        for (int jcol = 0; jcol < ncol; jcol++, ijcol++)
          cov[ijcol] += data[icol] * data[jcol];
  }

  /* Normalization */

  for (int icol = 0; icol < ncol; icol++)
  {
    if (num[icol] > 0)
    {
      mean[icol] /= num[icol];
      var[icol] /= num[icol];
      var[icol] -= mean[icol] * mean[icol];
      if (var[icol] <= 0) var[icol] = 0.;
    }
  }
  if (numiso > 0 && flagCorrel)
  {
    int ijcol = 0;
    for (int icol = 0; icol < ncol; icol++)
      for (int jcol = 0; jcol < ncol; jcol++, ijcol++)
      {
        cov[ijcol] /= numiso;
        cov[ijcol] -= mean[icol] * mean[jcol];
        cov[ijcol] /= sqrt(var[icol] * var[jcol]);
      }
  }

  /************/
  /* Printout */
  /************/

  if (! title.empty()) mestitle(1, title.c_str());

  /* Calculate the maximum size of the variable */

  int taille = 0;
  for (int icol = 0; icol < ncol; icol++)
  {
    _getRowname(radix, ncol, icol, db_name_get_by_att(db, iatts[icol]), string);
    taille = MAX(taille, (int ) strlen(string));
  }

  /* Print the header of the monovariate statistics */

  tab_print_rowname(" ", taille);
  if (_operExists(opers, EStatOption::NUM))
    tab_prints(NULL, "Number");
  if (_operExists(opers, EStatOption::MINI))
    tab_prints(NULL, "Minimum");
  if (_operExists(opers, EStatOption::MAXI))
    tab_prints(NULL, "Maximum");
  if (_operExists(opers, EStatOption::MEAN))
    tab_prints(NULL, "Mean");
  if (_operExists(opers, EStatOption::STDV))
    tab_prints(NULL, "St. Dev.");
  if (_operExists(opers, EStatOption::VAR))
    tab_prints(NULL, "Variance");
  message("\n");

  /* Print the monovariate statistics */

  for (int icol = 0; icol < ncol; icol++)
  {
    _getRowname(radix, ncol, icol, db_name_get_by_att(db, iatts[icol]), string);
    tab_print_rowname(string, taille);

    if (_operExists(opers, EStatOption::NUM))
      tab_printi(NULL, (int) num[icol]);
    if (num[icol] > 0)
    {
      if (_operExists(opers, EStatOption::MINI))
        tab_printg(NULL, mini[icol]);
      if (_operExists(opers, EStatOption::MAXI))
        tab_printg(NULL, maxi[icol]);
      if (_operExists(opers, EStatOption::MEAN))
        tab_printg(NULL, mean[icol]);
      if (_operExists(opers, EStatOption::STDV))
        tab_printg(NULL, sqrt(var[icol]));
      if (_operExists(opers, EStatOption::VAR))
        tab_printg(NULL, var[icol]);
    }
    else
    {
      if (_operExists(opers, EStatOption::MINI))
        tab_prints(NULL, STRING_NA);
      if (_operExists(opers, EStatOption::MAXI))
        tab_prints(NULL, STRING_NA);
      if (_operExists(opers, EStatOption::MEAN))
        tab_prints(NULL, STRING_NA);
      if (_operExists(opers, EStatOption::STDV))
        tab_prints(NULL, STRING_NA);
      if (_operExists(opers, EStatOption::VAR))
        tab_prints(NULL, STRING_NA);
    }
    message("\n");
  }
  message("\n");

  /* Print the correlation matrix  and count of isotopic samples */

  if (ncol > 1 && numiso > 0 && flagCorrel)
  {
    message("Number of isotopic active samples = %d\n", numiso);
    print_matrix("Correlation matrix", 0, 1, ncol, ncol, NULL, cov.data());
    message("\n");
  }

  return;
}

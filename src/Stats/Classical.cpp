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
#include "Db/Db.hpp"
#include "Stats/Classical.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "geoslib_f.h"

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
double _getQuantile(VectorDouble& tab, int ntab, double proba)
{
  int    rank;
  double p1,p2,v1,v2,value;

  if (FFFF(proba)) return(TEST);

  ut_sort_double(0,ntab,NULL,tab.data());
  rank = (int) (proba * (double) ntab);
  v1   = tab[rank];

  if (rank < ntab - 1)
  {
    v2 = tab[rank+1];
    p1 = (double) rank    / (double) ntab;
    p2 = (double) (1+rank)/ (double) ntab;
    value = v1 + (proba - p1) * (v2 - v1) / (p2 - p1);
  }
  else
  {
    value = v1;
  }
  return(value);
}

void _statList(void)
{
  messerr("List of operators available:");
  messerr("num    : Number of defined values");
  messerr("mean   : Mean over the defined values");
  messerr("var    : Variance over the defined values");
  messerr("stdv   : Standard Deviation over the defined values");
  messerr("mini   : Minimum over the defined values");
  messerr("maxi   : Maximum over the defined values");
  messerr("sum    : Sum over the variables");
  messerr("prop   : Proportion of values within [vmin;vmax]");
  messerr("quant  : Quantile corresponding to given probability");
  messerr("T      : Tonnage within [vmin;vmax]");
  messerr("Q      : Metal quantity within [vmin;vmax]");
  messerr("M      : Recovered mean within [vmin;vmax]");
  messerr("B      : Conventional Benefit within [vmin;vmax]");
}

ENUM_STATS _statIdentify(const String& oper)
{
  if (matchKeyword(oper,"num",false))   return STAT_NUM;
  if (matchKeyword(oper,"mean",false))  return STAT_MEAN;
  if (matchKeyword(oper,"var",false))   return STAT_VAR;
  if (matchKeyword(oper,"stdv",false))  return STAT_STDV;
  if (matchKeyword(oper,"mini",false))  return STAT_MINI;
  if (matchKeyword(oper,"maxi",false))  return STAT_MAXI;
  if (matchKeyword(oper,"sum",false))   return STAT_SUM;
  if (matchKeyword(oper,"prop",false))  return STAT_PROP;
  if (matchKeyword(oper,"quant",false)) return STAT_QUANT;
  if (matchKeyword(oper,"t",false))     return STAT_T;
  if (matchKeyword(oper,"q",false))     return STAT_Q;
  if (matchKeyword(oper,"m",false))     return STAT_M;
  if (matchKeyword(oper,"b",false))     return STAT_B;
  messerr("Invalid operator name (%s)",oper.c_str());
  _statList();
  return STAT_UNKNOWN;
}

String statsName(int ioper)
{
  if (ioper == STAT_NUM)   return "Number";
  if (ioper == STAT_MEAN)  return "Mean";
  if (ioper == STAT_VAR)   return "Variance";
  if (ioper == STAT_STDV)  return "St. Dev.";
  if (ioper == STAT_MINI)  return "Minimum";
  if (ioper == STAT_MAXI)  return "Maximum";
  if (ioper == STAT_SUM)   return "Sum";
  if (ioper == STAT_PROP)  return "Proportion";
  if (ioper == STAT_QUANT) return "Quantile";
  if (ioper == STAT_T)     return "Tonnage";
  if (ioper == STAT_Q)     return "Metal";
  if (ioper == STAT_M)     return "Grade";
  if (ioper == STAT_B)     return "Benefit";
  return "Unknown";
}

VectorInt statsList(const VectorString& opers)
{
  int noper = static_cast<int> (opers.size());
  VectorInt iopers(noper);
  for (int i = 0; i < noper; i++)
  {
    iopers[i] = _statIdentify(opers[i]);
    if (iopers[i] == STAT_UNKNOWN) return VectorInt();
  }
  return iopers;
}

VectorString statsNames(const VectorInt & iopers)
{
  VectorString names;
  for (int i = 0; i < (int) iopers.size(); i++)
  {
    names.push_back(statsName(iopers[i]));
  }
  return names;
}

void dbStatisticsVariables(Db *db,
                           const VectorInt& iatts,
                           const VectorInt& iopers,
                           int iattn,
                           double vmin,
                           double vmax,
                           double proba)
{
  int noper = static_cast<int> (iopers.size());
  if (noper <= 0) return;
  int natt  = static_cast<int> (iatts.size());
  if (natt <= 0) return;

  /* Loop on the samples */

  VectorDouble local(natt);
  for (int iech=0; iech<db->getSampleNumber(); iech++)
  {
    if (! db->isActive(iech)) continue;

    /* Loop on the variables */

    int neff     = 0;
    int nperc    = 0;
    double mean  = 0.;
    double var   = 0.;
    double stdv  = 0.;
    double sum   = 0.;
    double metal = 0.;
    double mini  =  1.e30;
    double maxi  = -1.e30;
    for (int iatt=0; iatt<natt; iatt++)
    {
      int jatt = iatts[iatt];
      double value = db->getArray(iech,jatt);
      if (FFFF(value)) continue;

      local[neff] = value;
      neff ++;
      mean += value;
      sum  += value;
      var  += value * value;
      if (value < mini) mini = value;
      if (value > maxi) maxi = value;
      if (! FFFF(vmin) && value < vmin) continue;
      if (! FFFF(vmax) && value > vmax) continue;
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
        if (iopers[i] == STAT_NUM)
          tab = (double) neff;
        else if (iopers[i] == STAT_MEAN)
          tab = mean;
        else if (iopers[i] == STAT_VAR)
          tab = var;
        else if (iopers[i] == STAT_STDV)
          tab = stdv;
        else if (iopers[i] == STAT_MINI)
          tab = mini;
        else if (iopers[i] == STAT_MAXI)
          tab = maxi;
        else if (iopers[i] == STAT_SUM)
          tab = sum;
        else if (iopers[i] == STAT_PROP)
          tab = (double) nperc / (double) neff;
        else if (iopers[i] == STAT_QUANT)
          tab = _getQuantile(local, neff, proba);
        else if (iopers[i] == STAT_T)
          tab = (double) nperc / (double) neff;
        else if (iopers[i] == STAT_Q)
          tab = metal / (double) neff;
        else if (iopers[i] == STAT_M)
          tab = (nperc > 0) ? metal / (double) nperc : TEST;
        else if (iopers[i] == STAT_B)
          tab = (!FFFF(vmin)) ? (metal - vmin) / (double) neff : TEST;
        else
          return;
      }
      db->setArray(iech, iattn + i, tab);
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
** \param[in]  iopers     List of the operator ranks
** \param[in]  flagIso    Restrain statistics to isotopic samples
** \param[in]  proba      Probability value (between 0 and 1)
** \param[in]  vmin       Minimum threshold
** \param[in]  vmax       Maximum threshold
**
*****************************************************************************/
VectorDouble dbStatisticsMono(Db *db,
                              const VectorInt& iatts,
                              const VectorInt& iopers,
                              bool flagIso,
                              double proba,
                              double vmin,
                              double vmax)
{
  int noper = static_cast<int> (iopers.size());
  int natt  = static_cast<int> (iatts.size());
  int nech  = db->getSampleNumber();

  // Find the Isotopic samples (optional)

  VectorDouble tab;
  VectorDouble local(nech,0.);
  VectorInt accept(nech,0);
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
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

  for (int iatt=0; iatt<natt; iatt++)
  {
    int neff     = 0;
    int nperc    = 0;
    double mean  = 0.;
    double var   = 0.;
    double stdv  = 0.;
    double sum   = 0.;
    double metal = 0.;
    double mini  =  1.e30;
    double maxi  = -1.e30;

    /* Loop on the samples */

    for (int iech=0; iech<db->getSampleNumber(); iech++)
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
      if (iopers[i] == STAT_NUM) tab.push_back((double) neff);
      if (neff > 0)
      {
        if (iopers[i] == STAT_MEAN)  tab.push_back(mean);
        if (iopers[i] == STAT_VAR)   tab.push_back(var);
        if (iopers[i] == STAT_STDV)  tab.push_back(stdv);
        if (iopers[i] == STAT_MINI)  tab.push_back(mini);
        if (iopers[i] == STAT_MAXI)  tab.push_back(maxi);
        if (iopers[i] == STAT_SUM)   tab.push_back(sum);
        if (iopers[i] == STAT_PROP)  tab.push_back((double) nperc / (double) neff);
        if (iopers[i] == STAT_QUANT) tab.push_back(_getQuantile(local, neff, proba));
        if (iopers[i] == STAT_T)     tab.push_back((double) nperc / (double) neff);
        if (iopers[i] == STAT_Q)     tab.push_back(metal / (double) neff);
        if (iopers[i] == STAT_M)
          tab.push_back((nperc > 0) ? metal / (double) nperc : TEST);
        if (iopers[i] == STAT_B)
          tab.push_back((!FFFF(vmin)) ? (metal - vmin) / (double) neff : TEST);
      }
      else
      {
        if (iopers[i] == STAT_MEAN)  tab.push_back(TEST);
        if (iopers[i] == STAT_VAR)   tab.push_back(TEST);
        if (iopers[i] == STAT_STDV)  tab.push_back(TEST);
        if (iopers[i] == STAT_MINI)  tab.push_back(TEST);
        if (iopers[i] == STAT_MAXI)  tab.push_back(TEST);
        if (iopers[i] == STAT_SUM)   tab.push_back(TEST);
        if (iopers[i] == STAT_PROP)  tab.push_back(TEST);
        if (iopers[i] == STAT_QUANT) tab.push_back(TEST);
        if (iopers[i] == STAT_T)     tab.push_back(TEST);
        if (iopers[i] == STAT_Q)     tab.push_back(TEST);
        if (iopers[i] == STAT_M)     tab.push_back(TEST);
        if (iopers[i] == STAT_B)     tab.push_back(TEST);
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
    messerr("This function requires the number of variables (%d) to be equal to 1",
            db->getLocatorNumber(ELoc::Z));
    return props;
  }
  int nech = db->getSampleNumber();
  int nfac = db->getFaciesNumber();

  // Calculate the proportions

  props.resize(nfac,0.);
  int neff = 0;
  for (int iech=0; iech<nech; iech++)
  {
    if (! db->isActiveAndDefined(iech,0)) continue;
    int ifac = (int) db->getVariable(iech,0);
    if (ifac <= 0) continue;
    props[ifac-1] += 1.;
    neff++;
  }

  // Normalization

  if (neff > 0)
  {
    for (int ifac = 0 ; ifac < nfac; ifac++)
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
    messerr("This function requires the number of variables (%d) to be equal to 1",
            db->getLocatorNumber(ELoc::Z));
    return TEST;
  }

  // Calculate the proportions

  double prop = 0.;
  int neff = 0;
  for (int iech=0; iech<db->getSampleNumber(); iech++)
  {
    if (! db->isActiveAndDefined(iech,0)) continue;
    int ifac = (int) db->getVariable(iech,0);
    if (ifac == 1) prop += 1.;
    neff++;
  }

  // Normalization

  if (neff > 0) prop = prop / neff;
  return prop;
}

/****************************************************************************/
/*!
**  Print the correlation matrix for a set of variables of a Db
**
** \return  Error Return code
**
** \param[in]  db          Db structure
** \param[in]  iatts       Vector of attribute ranks
** \param[in]  flagIso    Restrain statistics to isotopic samples
**
*****************************************************************************/
VectorDouble dbStatisticsMulti(Db *db, const VectorInt& iatts, bool flagIso)
{
  int natt = static_cast<int> (iatts.size());

  /* Preliminary checks */

  if (natt <= 1)
    messerr("Correlation matrix will not be printed for a single variable");

  /* Core allocation */

  VectorDouble data(natt,0.);
  VectorDouble mean(natt,0.);
  VectorDouble var(natt,0.);
  VectorDouble num(natt,0.);
  VectorDouble cov(natt * natt,0.);

  /* Loop on the samples */

  int numiso = 0;
  for (int iech=0; iech<db->getSampleNumber(); iech++)
  {
    if (! db->isActive(iech)) continue;

    /* Look for isotopic sample */

    int nundef = 0;
    for (int iatt=0; iatt<natt; iatt++)
    {
      data[iatt] = db->getArray(iech,iatts[iatt]);
      if (FFFF(data[iatt])) nundef++;
    }
    if (flagIso && nundef > 0) continue;

    /* Calculate the 1-point statistics */

    for (int iatt=0; iatt<natt; iatt++)
    {
      if (FFFF(data[iatt])) continue;
      num[iatt]  += 1.;
      mean[iatt] += data[iatt];
      var[iatt]  += data[iatt] * data[iatt];
    }

    if (nundef > 0) continue;
    numiso ++;
    int ijatt = 0;
      for (int iatt=0; iatt<natt; iatt++)
        for (int jatt=0; jatt<natt; jatt++)
        {
          cov[ijatt] += data[iatt] * data[jatt];
          ijatt++;
        }
  }

  /* Normalization */

  for (int iatt=0; iatt<natt; iatt++)
  {
    if (num[iatt] > 0)
    {
      mean[iatt] /= num[iatt];
      var[iatt]  /= num[iatt];
      var[iatt]  -= mean[iatt] * mean[iatt];
      if (var[iatt] <= 0) var[iatt] = 0.;
    }
  }
  if (numiso > 0)
  {
    int ijatt = 0;
    for (int iatt=0; iatt<natt; iatt++)
      for (int jatt=0; jatt<natt; jatt++)
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
** \param[in]  iopers      List of the operator ranks
** \param[in]  varnames    List of variables
** \param[in]  title       Title for the printout (optional)
**
*****************************************************************************/
String statisticsMonoPrint(const VectorDouble& stats,
                           const VectorInt& iopers,
                           const VectorString& varnames,
                           const String& title)
{
  int noper = static_cast<int> (iopers.size());
  int natt  = static_cast<int> (varnames.size());
  std::stringstream sstr;

  // Constitute the vector of row and column names
  VectorString colnames = statsNames(iopers);

  // Printout the matrix
  sstr << toMatrix(title,colnames,varnames,false,noper,natt,stats,true);

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
String statisticsMultiPrint(const VectorDouble& stats,
                            const VectorString& varnames,
                            const String& title)
{
  int natt = static_cast<int> (varnames.size());
  std::stringstream sstr;

  sstr << toMatrix(title, VectorString(), VectorString(), true,
                   natt, natt, stats, true);

  return sstr.str();
}

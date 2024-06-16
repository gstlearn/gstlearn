/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "geoslib_old_f.h"

#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Stats/Classical.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "Basic/AStringable.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Model/Model.hpp"
#include "Space/SpaceRN.hpp"
#include "Space/SpaceTarget.hpp"
#include "Variogram/VarioParam.hpp"
#include "Variogram/Vario.hpp"
#include "Polygon/Polygons.hpp"

#include <math.h>
#include <Matrix/Table.hpp>
#include <string.h>

/****************************************************************************/
/*!
 **  Check the operator name
 **
 ** \return  1 if the operator is valid; 0 otherwise
 **
 ** \param[in]  oper A EStatOption item
 ** \param[in]  flag_multi  1 if multivariate operator is authorized
 ** \param[in]  flag_indic  1 if indicator ("plus","minus","zero") is authorized
 ** \param[in]  flag_sum    1 if sum of variable is authorized
 ** \param[in]  flag_median 1 if median is authorized
 ** \param[in]  flag_qt     1 if QT ("ore","metal") is authorized
 **
 ** \remarks If an error occurred, the message is printed
 **
 *****************************************************************************/
bool _operStatisticsCheck(const EStatOption &oper,
                          int flag_multi,
                          int flag_indic,
                          int flag_sum,
                          int flag_median,
                          int flag_qt)
{
  bool valid;

  /* Initializations */

  valid = false;

  /* Monovariate check */

  if (oper == EStatOption::NUM)  valid = true;
  if (oper == EStatOption::MEAN) valid = true;
  if (oper == EStatOption::VAR)  valid = true;
  if (oper == EStatOption::CORR) valid = true;
  if (oper == EStatOption::STDV) valid = true;
  if (oper == EStatOption::MINI) valid = true;
  if (oper == EStatOption::MAXI) valid = true;
  if (flag_sum)
  {
    if (oper == EStatOption::SUM) valid = true;
  }
  if (flag_median)
  {
    if (oper == EStatOption::MEDIAN) valid = true;
  }

  /* Multivariate check */

  if (flag_multi)
  {
    if (oper == EStatOption::MEAN2) valid = true;
    if (oper == EStatOption::VAR2)  valid = true;
    if (oper == EStatOption::STDV2) valid = true;
    if (flag_sum)
    {
      if (oper == EStatOption::SUM2) valid = true;
    }
  }

  /* Indicator check */

  if (flag_indic)
  {
    if (oper == EStatOption::PLUS)  valid = true;
    if (oper == EStatOption::MOINS) valid = true;
    if (oper == EStatOption::ZERO)  valid = true;
  }

  /* QT variables check */

  if (flag_qt)
  {
    if (oper == EStatOption::ORE)   valid = true;
    if (oper == EStatOption::METAL) valid = true;
  }

  if (!valid) messerr("Invalid operator");

  return (valid);
}

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
  int ifac = (int) dbin->getLocVariable(ELoc::Z,rank, 0);
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
  int ifac1 = (int) dbin->getLocVariable(ELoc::Z,dbin->getGrid().indiceToRank(indg), 0);
  indg[pos] += orient;
  int ifac2 = (int) dbin->getLocVariable(ELoc::Z,dbin->getGrid().indiceToRank(indg), 0);
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
 ** \return  true if the operator is mentioned; false otherwise
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
  if (noper == 0) return true;
  for (int i = 0; i < noper; i++)
  {
    if (opers[i] == refe) return true;
  }
  return false;
}

/****************************************************************************/
/*!
 **  Copy the multivariate into monovariate statistics (before printout)
 **
 ** \param[in]  ncol  Dimension of the (square) matrix
 ** \param[in,out]  tab   Array to be refactored
 **
 *****************************************************************************/
void _refactor(int ncol, VectorDouble& tab)
{
  int ix, iy;

  for (ix = 0; ix < ncol; ix++)
    for (iy = 0; iy < ncol; iy++)
      tab[ix * ncol + iy] = tab[iy * ncol + iy];
}

/****************************************************************************/
/*!
 **  Copy the multivariate or monovariate statistics into the returned array
 **
 ** \param[in]  nx    First dimension of the matrix
 ** \param[in]  ny    First dimension of the matrix
 ** \param[in]  tabin Array to be refactored
 **
 ** \param[out]  tabout Array to be refactored
 **
 *****************************************************************************/
void _copyResults(int nx,
                  int ny,
                  const VectorDouble &tabin,
                  VectorDouble &tabout)
{
  int ix, iy, lec;

  for (ix = lec = 0; ix < nx; ix++)
    for (iy = 0; iy < ny; iy++, lec++)
      tabout[lec] = tabin[lec];
}

/****************************************************************************/
/*!
 **  Calculate the indices of the cell neighboring a target cell
 **
 ** \param[in]  ndim   Space dimension
 ** \param[in]  radius Neighborhood radius
 ** \param[in]  rank0  Rank of the neighbor
 ** \param[in]  indg0  Array of indices of the target cell
 **
 ** \param[out] indg   Array of indices of the neighboring cell
 **
 *****************************************************************************/
void _neighboringCell(int ndim,
                      int radius,
                      int rank0,
                      const VectorInt &indg0,
                      VectorInt &indg)
{
  int nei1d, value, divid, count, reste, ratio, idim;

  /* Initializations */

  nei1d = 2 * radius + 1;
  count = (int) pow(nei1d, (double) ndim);
  value = rank0;

  /* Loop on the space dimensions */

  divid = count;
  for (int jdim = 0; jdim < ndim; jdim++)
  {
    idim = ndim - jdim - 1;
    divid /= nei1d;
    ratio = value / divid;
    reste = value - ratio * divid;
    value = reste;
    indg[idim] = indg0[idim] + ratio - radius;
  }
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

  VH::sortInPlace(tab, true, ntab);
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

VectorString statOptionToName(const std::vector<EStatOption>& opers)
{
  VectorString names;
  for (int i = 0; i < (int) opers.size(); i++)
  {
    EStatOption oper = opers[i];
    names.push_back(oper.getKey());
  }
  return names;
}

std::vector<EStatOption> KeysToStatOptions(const VectorString& opers)
{
  std::vector<EStatOption> options;

  for (int i = 0; i < (int) opers.size(); i++)
  {
    EStatOption opt = EStatOption::fromKey(opers[i]);
    if (opt != EStatOption::UNKNOWN) options.push_back(opt);
  }
  return options;
}

/**
 * \copydoc STATS_2
 *
 * @return Store several statistics calculated on a set of variables of a Db and store them
 * in this same Db in variables already created.
 * These functions should not be used in Target Language.
 */
void dbStatisticsVariables(Db *db,
                           const VectorString &names,
                           const std::vector<EStatOption> &opers,
                           int iptr0,
                           double proba,
                           double vmin,
                           double vmax)
{
  int noper = static_cast<int>(opers.size());
  if (noper <= 0) return;
  if (names.empty()) return;
  VectorInt iuids = db->getUIDs(names);
  int niuid = static_cast<int>(iuids.size());

  /* Loop on the samples */

  VectorDouble local(niuid);
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
    for (int iuid = 0; iuid < niuid; iuid++)
    {
      int juid = iuids[iuid];
      double value = db->getArray(iech, juid);
      if (FFFF(value)) continue;

      local[neff] = value;
      neff++;
      mean += value;
      sum  += value;
      var  += value * value;
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

/**
 * \copydoc STATS_0
 *
 * @param opers      List of the operator ranks
 * @param flagIso    Restrain statistics to isotopic samples
 * @param proba      Probability value (between 0 and 1)
 * @param vmin       Minimum threshold
 * @param vmax       Maximum threshold
 *
 * @return A Table containing the results
 */
Table dbStatisticsMono(Db *db,
                       const VectorString &names,
                       const std::vector<EStatOption> &opers,
                       bool flagIso,
                       double proba,
                       double vmin,
                       double vmax,
                       const String& title)
{
  VectorInt iuids = db->getUIDs(names);
  int niuid = static_cast<int>(iuids.size());
  int noper = static_cast<int>(opers.size());
  int nech = db->getSampleNumber();

  // Find the Isotopic samples (optional)

  VectorDouble tab;
  VectorDouble local(nech, 0.);
  VectorBool accept(nech, true);
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    accept[iech] = false;
    if (! db->isActive(iech)) continue;
    accept[iech] = true;
    int nundef = 0;
    for (int iuid = 0; iuid < niuid; iuid++)
    {
      double value = db->getArray(iech, iuids[iuid]);
      if (FFFF(value)) nundef++;
    }
    if (flagIso && nundef > 0) accept[iech] = false;
  }

  /* Loop on the attributes */

  for (int iuid = 0; iuid < niuid; iuid++)
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
      if (! accept[iech]) continue;
      double value = db->getArray(iech, iuids[iuid]);
      // Skip TEST values (this is necessary when flagIso is set to FALSE)
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

  Table table;
  if (title.empty())
    table.setSkipTitle(true);
  else
    table.setTitle(title);
  table.setSkipDescription(true);
  table.resetFromVD(niuid, noper, tab, false);

  for (int irow=0; irow<niuid; irow++)
    table.setRowName(irow, db->getNameByUID(iuids[irow]));
  for (int icol=0; icol<noper; icol++)
    table.setColumnName(icol, opers[icol].getDescr());

  return table;
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
    int ifac = (int) db->getLocVariable(ELoc::Z,iech, 0);
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
 **
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
    int ifac = (int) db->getLocVariable(ELoc::Z,iech, 0);
    if (ifac == 1) prop += 1.;
    neff++;
  }

  // Normalization

  if (neff > 0) prop = prop / neff;
  return prop;
}

/**
 * \copydoc STATS_0
 *
 * @param flagIso    Restrain statistics to isotopic samples
 *
 * @return A Table containing the correlation matrix
 */
Table dbStatisticsCorrel(Db *db, const VectorString &names, bool flagIso, const String& title)
{
  VectorInt iuids = db->getUIDs(names);
  int niuid = static_cast<int>(iuids.size());

  /* Core allocation */

  VectorDouble data(niuid, 0.);
  VectorDouble mean(niuid, 0.);
  VectorDouble var(niuid, 0.);
  VectorDouble num(niuid, 0.);
  VectorDouble cov(niuid * niuid, 0.);

  /* Loop on the samples */

  int numiso = 0;
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;

    /* Look for isotopic sample */

    int nundef = 0;
    for (int iuid = 0; iuid < niuid; iuid++)
    {
      data[iuid] = db->getArray(iech, iuids[iuid]);
      if (FFFF(data[iuid])) nundef++;
    }
    if (flagIso && nundef > 0) continue;

    /* Calculate the 1-point statistics */

    for (int iuid = 0; iuid < niuid; iuid++)
    {
      if (FFFF(data[iuid])) continue;
      num[iuid] += 1.;
      mean[iuid] += data[iuid];
      var[iuid] += data[iuid] * data[iuid];
    }

    if (nundef > 0) continue;
    numiso++;
    int ijuid = 0;
    for (int iuid = 0; iuid < niuid; iuid++)
      for (int juid = 0; juid < niuid; juid++)
      {
        cov[ijuid] += data[iuid] * data[juid];
        ijuid++;
      }
  }

  /* Normalization */

  for (int iuid = 0; iuid < niuid; iuid++)
  {
    if (num[iuid] > 0)
    {
      mean[iuid] /= num[iuid];
      var[iuid] /= num[iuid];
      var[iuid] -= mean[iuid] * mean[iuid];
      if (var[iuid] <= 0) var[iuid] = 0.;
    }
  }
  if (numiso > 0)
  {
    int ijuid = 0;
    for (int iuid = 0; iuid < niuid; iuid++)
      for (int juid = 0; juid < niuid; juid++)
      {
        cov[ijuid] /= numiso;
        cov[ijuid] -= mean[iuid] * mean[juid];
        cov[ijuid] /= sqrt(var[iuid] * var[juid]);
        ijuid++;
      }
  }

  // Store the results in the symmetric square matrix
  VectorString namloc = db->getNames(names);
  int nvar = (int) namloc.size();

  Table table;
  if (title.empty())
    table.setSkipTitle(true);
  else
    table.setTitle(title);
  table.setSkipDescription(true);
  table.resetFromVD(nvar, nvar, cov, false);
  for (int ivar = 0; ivar < nvar; ivar++)
    table.setColumnName(ivar, db->getNameByUID(iuids[ivar]));
  for (int ivar = 0; ivar < nvar; ivar++)
    table.setRowName(ivar, db->getNameByUID(iuids[ivar]));

  return table;
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

/**
 * \copydoc STATS_0
 *
 * @param opers      List of the operator ranks
 * @param flagIso    Restrain statistics to isotopic samples
 * @param flagCorrel Print the correlation matrix
 * @param radix      Radix given to the printout
 */
void dbStatisticsPrint(const Db *db,
                       const VectorString &names,
                       const std::vector<EStatOption> &opers,
                       bool flagIso,
                       bool flagCorrel,
                       const String &title,
                       const String &radix)
{
  VectorInt iuids = db->getUIDs(names);
  if (iuids.empty()) return;

  char string[50];
  int ncol = static_cast<int>(iuids.size());

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
      data[icol] = db->getArray(iech, iuids[icol]);
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
    ijcol = 0;
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
    _getRowname(radix, ncol, icol, db_name_get_by_att(db, iuids[icol]), string);
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
    _getRowname(radix, ncol, icol, db_name_get_by_att(db, iuids[icol]), string);
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

/**
 * Sphering procedure
 * @param X Input Data vector
 * @return The Sphering matrix (or nullptr if problem)
 *
 * @remark When performing the (forward) sphering, you must perform the following operation
 * @remark        X <- prodMatMat(X, S)
 */
MatrixSquareGeneral* sphering(const AMatrix* X)
{
  if (X->empty()) return nullptr;
  int nech = X->getNRows();
  int nvar = X->getNCols();

  AMatrix* TX = X->transpose();
  AMatrix* prod = MatrixFactory::prodMatMat(TX, X);
  MatrixSquareSymmetric* prodsym = dynamic_cast<MatrixSquareSymmetric*>(prod);
  if (prodsym == nullptr) return nullptr;

  prodsym->prodScalar(1. / (double) nech);
  if (prodsym->computeEigen()) return nullptr;
  VectorDouble eigen_values = prodsym->getEigenValues();
  MatrixSquareGeneral* S = prodsym->getEigenVectors()->clone();

  // Invert the sign of the second Eigen vector (for compatibility with R output)
  for (int ivar = 0; ivar < nvar ; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
    {
      double signe = (jvar < nvar-1) ? 1 : -1;
      S->setValue(ivar, jvar,
                 signe * S->getValue(ivar, jvar) / sqrt(eigen_values[jvar]));
    }

  delete TX;
  delete prod;

  return S;
}

/**
 * \copydoc STATS_1
 *
 * @return Vector of results
 *
 */
VectorDouble dbStatisticsPerCell(Db *db,
                                 DbGrid *dbgrid,
                                 const EStatOption &oper,
                                 const String& name1,
                                 const String& name2,
                                 const VectorDouble &cuts)
{
  VectorDouble result;
  int iuid = db->getUID(name1);
  int juid = -1;
  if (! name2.empty()) juid = db->getUID(name2);
  double z1 = 0.;
  double z2 = 0.;
  int nxyz = dbgrid->getSampleNumber();
  int ncut = (int) cuts.size();
  int ndim = dbgrid->getNDim();
  if (juid < 0) juid = iuid;

  bool flag1 = false;
  bool flag2 = false;
  bool flag_denorm = false;
  bool flag_q = false;
  bool flag_t = false;
  bool flag_s1 = false;
  bool flag_s2 = false;
  bool flag_v1 = false;
  bool flag_v2 = false;
  bool flag_v12 = false;
  bool flag_mini = false;
  bool flag_maxi = false;

  /* Check the operator validity */

  if (! _operStatisticsCheck(oper, 1, 0, 1, 0, 1)) return VectorDouble();

  /* Set the relevant flags */

  if (oper == EStatOption::NUM)
    flag1 = 1;
  else if (oper == EStatOption::MEAN)
    flag1 = flag_s1 = 1;
  else if (oper == EStatOption::SUM)
    flag1 = flag_s1 = flag_denorm = 1;
  else if (oper == EStatOption::STDV)
    flag1 = flag_s1 = flag_v1 = 1;
  else if (oper == EStatOption::VAR)
    flag1 = flag_s1 = flag_v1 = 1;
  else if (oper == EStatOption::MEAN2)
    flag2 = flag_s2 = 1;
  else if (oper == EStatOption::SUM2)
    flag2 = flag_s2 = flag_denorm = 1;
  else if (oper == EStatOption::STDV2)
    flag2 = flag_s2 = flag_v2 = 1;
  else if (oper == EStatOption::VAR2)
    flag2 = flag_s2 = flag_v2 = 1;
  else if (oper == EStatOption::COV)
    flag2 = flag_s1 = flag_s2 = flag_v12 = 1;
  else if (oper == EStatOption::CORR)
    flag2 = flag_s1 = flag_s2 = flag_v1 = flag_v2 = flag_v12 = 1;
  else if (oper == EStatOption::MINI)
    flag1 = flag_mini = 1;
  else if (oper == EStatOption::MAXI)
    flag1 = flag_maxi = 1;
  else if (oper == EStatOption::ORE)
    flag1 = flag_t = 1;
  else if (oper == EStatOption::METAL)
    flag1 = flag_q = 1;
  else
    return (1);

  /* Core allocation */

  VectorDouble coor(ndim);
  VectorInt indg(ndim);
  VectorDouble nn(nxyz, 0.);
  VectorDouble s1;
  VectorDouble s2;
  VectorDouble v1;
  VectorDouble v2;
  VectorDouble v12;
  VectorDouble mini;
  VectorDouble maxi;
  VectorDouble cutval;
  if (flag_s1)
    s1.resize(nxyz, 0.);
  if (flag_s2)
    s2.resize(nxyz, 0.);
  if (flag_v1)
    v1.resize(nxyz, 0.);
  if (flag_v2)
    v2.resize(nxyz, 0.);
  if (flag_v12)
    v12.resize(nxyz, 0.);
  if (flag_mini)
    mini.resize(nxyz, TEST);
  if (flag_maxi)
    maxi.resize(nxyz, TEST);
  if (flag_t || flag_q)
  {
    if (ncut <= 0)
    {
      messerr("For this calculation, 'cuts' must be provided");
      return VectorDouble();
    }
    cutval.resize(nxyz * ncut, 0.);
  }

  /* Loop on the samples */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;

    /* Check the variable(s) */

    if (flag1)
    {
      z1 = db->getArray(iech, iuid);
      if (FFFF(z1)) continue;
    }
    if (flag2)
    {
      z2 = db->getArray(iech, juid);
      if (FFFF(z2)) continue;
    }

    /* Check the location of the data in the grid */

    db->getCoordinatesPerSampleInPlace(iech, coor);
    int iad = dbgrid->getGrid().coordinateToRank(coor);
    if (iad < 0 || iad >= nxyz) continue;
    nn[iad]++;
    if (flag_s1) s1[iad] += z1;
    if (flag_s2) s2[iad] += z2;
    if (flag_v1) v1[iad] += z1 * z1;
    if (flag_v2) v2[iad] += z2 * z2;
    if (flag_v12) v12[iad] += z1 * z2;
    if (flag_mini)
    {
      if (FFFF(mini[iad]) || z1 < mini[iad]) mini[iad] = z1;
    }
    if (flag_maxi)
    {
      if (FFFF(maxi[iad]) || z1 > maxi[iad]) maxi[iad] = z1;
    }
    if (flag_t)
    {
      for (int icut = 0; icut < ncut; icut++)
        if (z1 >= cuts[icut]) cutval[icut + iech * ncut] += 1.;
    }
    if (flag_q)
    {
      for (int icut = 0; icut < ncut; icut++)
        if (z1 >= cuts[icut]) cutval[icut + iech * ncut] += z1;
    }
  }

  /* Normalization */

  for (int i = 0; i < nxyz; i++)
  {
    double ratio = nn[i];
    if (ratio <= 0)
    {
      if (flag_s1) s1[i] = TEST;
      if (flag_s2) s2[i] = TEST;
      if (flag_v1) v1[i] = TEST;
      if (flag_v2) v2[i] = TEST;
      if (flag_v12) v12[i] = TEST;
      if (flag_mini) mini[i] = TEST;
      if (flag_maxi) maxi[i] = TEST;
      if (flag_t || flag_q) for (int icut = 0; icut < ncut; icut++)
        cutval[icut + i * ncut] = TEST;
    }
    else
    {
      if (flag_s1) s1[i] /= ratio;
      if (flag_s2) s2[i] /= ratio;
      if (flag_v1)
      {
        v1[i] = v1[i] / ratio - s1[i] * s1[i];
        v1[i] = (v1[i] < 0.) ? 0. : sqrt(v1[i]);
      }
      if (flag_v2)
      {
        v2[i] = v2[i] / ratio - s2[i] * s2[i];
        v2[i] = (v2[i] < 0.) ? 0. : sqrt(v2[i]);
      }
      if (flag_v12)
      {
        v12[i] = v12[i] / ratio - s1[i] * s2[i];
      }
      if (flag_denorm)
      {
        if (flag_s1) s1[i] *= ratio;
        if (flag_s2) s2[i] *= ratio;
      }
      if (flag_t || flag_q)
      {
        for (int icut = 0; icut < ncut; icut++)
          cutval[icut + i * ncut] /= ratio;
      }
    }
  }

  /* Dispatch according to the type of result expected */

  result.resize(nxyz * MAX(1, ncut), 0.);
  for (int i = 0; i < nxyz; i++)
  {
    if (oper == EStatOption::NUM)
      result[i] = nn[i];
    else if (oper == EStatOption::MEAN)
      result[i] = s1[i];
    else if (oper == EStatOption::SUM)
      result[i] = s1[i];
    else if (oper == EStatOption::STDV)
      result[i] = v1[i];
    else if (oper == EStatOption::VAR)
      result[i] = v1[i] * v1[i];
    else if (oper == EStatOption::MEAN2)
      result[i] = s2[i];
    else if (oper == EStatOption::SUM2)
      result[i] = s2[i];
    else if (oper == EStatOption::STDV2)
      result[i] = v2[i];
    else if (oper == EStatOption::VAR2)
      result[i] = v2[i] * v2[i];
    else if (oper == EStatOption::COV)
      result[i] = v12[i];
    else if (oper == EStatOption::CORR)
    {
      if (v1[i] > 0. && v2[i] > 0.) result[i] = v12[i] / (v1[i] * v2[i]);
    }
    else if (oper == EStatOption::MINI)
      result[i] = mini[i];
    else if (oper == EStatOption::MAXI)
      result[i] = maxi[i];
    else if (oper == EStatOption::ORE)
      for (int icut = 0; icut < ncut; icut++)
        result[i + icut * nxyz] = cutval[icut + i * ncut];
    else if (oper == EStatOption::METAL)
      for (int icut = 0; icut < ncut; icut++)
        result[i + icut * nxyz] = cutval[icut + i * ncut];
    else
      return VectorDouble();
  }

  return result;
}


/**
 * \copydoc STATS_0
 *
 * @param oper       Operator
 * @param flagMono   When True, statistics by variable; otherwise, statistics by pair of variables
 *
 * @return A Table containing the results
 */
Table dbStatisticsMulti(Db *db,
                        const VectorString &names,
                        const EStatOption &oper,
                        bool flagMono,
                        const String& title)
{
  VectorDouble result;

  /* Initializations */

  VectorInt cols = db->getUIDs(names);
  int ncol = (int) cols.size();
  int nech = db->getSampleNumber();
  int ncol2 = ncol * ncol;

  /* Check that all variables are defined */

  for (int icol = 0; icol < ncol; icol++)
  {
    int jcol = cols[icol];
    if (!db->isColIdxValid(jcol))
    {
      messerr("Error: Variable %d is not defined", cols[icol]);
      return Table();
    }
  }

  /* Check the validity of the operator */

  if (! _operStatisticsCheck(oper, 0, 1, 0, 0, 0)) return Table();

  /* Core allocation */

  if (flagMono)
    result.resize(ncol);
  else
    result.resize(ncol2);

  VectorDouble num(ncol2, 0.);
  VectorDouble m1(ncol2, 0.);
  VectorDouble m2(ncol2, 0.);
  VectorDouble v1(ncol2, 0.);
  VectorDouble v2(ncol2, 0.);
  VectorDouble v12(ncol2, 0.);
  VectorDouble mini(ncol2, 0.);
  VectorDouble maxi(ncol2, 0.);
  VectorDouble plus(ncol2, 0.);
  VectorDouble moins(ncol2, 0.);
  VectorDouble zero(ncol2, 0.);

  /* Initializations */

  for (int i = 0; i < ncol * ncol; i++)
  {
    num[i] = m1[i] = m2[i] = v1[i] = v2[i] = v12[i] = 0.;
    plus[i] = moins[i] = zero[i] = 0;
    mini[i] = 1.e30;
    maxi[i] = -1.e30;
  }

  /* Loop on the samples */

  for (int iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    double weight = db->getWeight(iech);

    /* Loop on the first variable */

    for (int icol1 = 0; icol1 < ncol; icol1++)
    {
      int jcol1 = cols[icol1];
      double val1 = db->getArray(iech, jcol1);
      if (FFFF(val1)) continue;

      /* Loop on the second variable */

      for (int icol2 = 0; icol2 < ncol; icol2++)
      {
        int jcol2 = cols[icol2];
        double val2 = db->getArray(iech, jcol2);
        if (FFFF(val2)) continue;

        /* Update statistics */

        int iad = icol1 * ncol + icol2;
        num[iad] += weight;
        m1[iad] += weight * val1;
        m2[iad] += weight * val2;
        v1[iad] += weight * val1 * val1;
        v2[iad] += weight * val2 * val2;
        v12[iad] += weight * val1 * val2;
        if (val1 < mini[iad]) mini[iad] = val1;
        if (val1 > maxi[iad]) maxi[iad] = val1;
        if (val1 > 0) plus[iad] += 1;
        if (val1 < 0) moins[iad] += 1;
        if (val1 == 0) zero[iad] += 1;
      }
    }
  }

  /* Normalization */

  for (int icol1 = 0; icol1 < ncol; icol1++)
    for (int icol2 = 0; icol2 < ncol; icol2++)
    {
      int iad = icol1 * ncol + icol2;
      if (num[iad] <= 0)
      {
        m1[iad] = TEST;
        m2[iad] = TEST;
        v1[iad] = TEST;
        m2[iad] = TEST;
        v12[iad] = TEST;
        mini[iad] = TEST;
        maxi[iad] = TEST;
        plus[iad] = TEST;
        moins[iad] = TEST;
        zero[iad] = TEST;
      }
      else
      {
        m1[iad] /= num[iad];
        m2[iad] /= num[iad];
        v1[iad] = v1[iad] / num[iad] - m1[iad] * m1[iad];
        v2[iad] = v2[iad] / num[iad] - m2[iad] * m2[iad];
        v12[iad] = v12[iad] / num[iad] - m1[iad] * m2[iad];

        if (oper == EStatOption::STDV)
        {
          v12[iad] = (v12[iad] > 0) ? sqrt(v12[iad]) : 0.;
        }
        if (oper == EStatOption::CORR)
        {
          v12[iad] = (v1[iad] > 0 && v2[iad] > 0) ? v12[iad] / sqrt(v1[iad] * v2[iad]) : 0.;
        }
      }
    }

  /* Printout */

  int nx = 0;
  int ny = 0;
  if (flagMono)
  {
    nx = 1;
    ny = ncol;
    _refactor(ncol, num);
    _refactor(ncol, m1);
    _refactor(ncol, v1);
    _refactor(ncol, v12);
    _refactor(ncol, mini);
    _refactor(ncol, maxi);
    _refactor(ncol, plus);
    _refactor(ncol, moins);
    _refactor(ncol, zero);
  }
  else
  {
    nx = ncol;
    ny = ncol;
  }

  /* Set the return array */

  if (oper == EStatOption::NUM)
    _copyResults(nx, ny, num, result);
  else if (oper == EStatOption::MEAN)
    _copyResults(nx, ny, m1, result);
  else if (oper == EStatOption::VAR)
    _copyResults(nx, ny, v12, result);
  else if (oper == EStatOption::CORR)
    _copyResults(nx, ny, v12, result);
  else if (oper == EStatOption::STDV)
    _copyResults(nx, ny, v12, result);
  else if (oper == EStatOption::MINI)
    _copyResults(nx, ny, mini, result);
  else if (oper == EStatOption::MAXI)
    _copyResults(nx, ny, maxi, result);
  else if (oper == EStatOption::PLUS)
    _copyResults(nx, ny, plus, result);
  else if (oper == EStatOption::MOINS)
    _copyResults(nx, ny, moins, result);
  else if (oper == EStatOption::ZERO)
    _copyResults(nx, ny, zero, result);
  else
    messageAbort("This error should never happen");

  // Load the results in a Table

  Table table;
  if (title.empty())
    table.setTitle(oper.getDescr());
  else
    table.setTitle(concatenateStrings(":",title,oper.getDescr()));
  table.setSkipDescription(true);

  if (flagMono)
  {
    table.resetFromVD(ncol, 1, result, false);
  }
  else
  {
    table.resetFromVD(ncol, ncol, result, false);
    for (int icol=0; icol<ncol; icol++)
      table.setColumnName(icol, db->getNameByUID(cols[icol]));
  }
  for (int irow=0; irow<ncol; irow++)
      table.setRowName(irow, db->getNameByUID(cols[irow]));

  return table;
}

/****************************************************************************/
/*!
 **  Calculates the monovariate statistics within cells of a grid
 **
 ** \return  Error return code
 **
 ** \param[in]  db     Db for the points
 ** \param[in]  dbgrid Db for the grid
 ** \param[in]  names  Vector of target variable names
 ** \param[in]  oper   A EStatOption item
 ** \param[in]  radius Neighborhood radius
 ** \param[in]  iptr0  Storage address (first variable)
 **
 *****************************************************************************/
int dbStatisticsInGridTool(Db *db,
                           DbGrid *dbgrid,
                           const VectorString& names,
                           const EStatOption &oper,
                           int radius,
                           int iptr0)
{
  int iptm = -1;
  int iptn = -1;
  int nxyz = dbgrid->getSampleNumber();
  int ndim = dbgrid->getNDim();
  VectorInt iuids = db->getUIDs(names);
  int nuid = (int) iuids.size();
  int count = (int) pow(2. * radius + 1., (double) ndim);

  /* Check the validity of the requested function */

  if (! _operStatisticsCheck(oper, 0, 1, 0, 1, 0)) return 1;

  /* Create the attributes */

  if (oper == EStatOption::MEAN || oper == EStatOption::VAR || oper == EStatOption::STDV)
    iptn = dbgrid->addColumnsByConstant(1, 0.);
  if (oper == EStatOption::VAR || oper == EStatOption::STDV)
    iptm = dbgrid->addColumnsByConstant(1, 0.);

  /* Core allocation */

  VectorDouble coor(ndim);
  VectorInt indg0(ndim);
  VectorInt indg(ndim);
  VectorDouble medtab;
  if (oper == EStatOption::MEDIAN) medtab.resize(ndim);

  /* Loop on the variables */

  for (int iuid = 0; iuid < nuid; iuid++)
  {

    /* Create the output attribute */

    int juid = iuids[iuid];
    int iptr = iptr0 + iuid;

    /* Loop on the samples */

    int nmed = 0;
    for (int iech = 0; iech < db->getSampleNumber(); iech++)
    {

      /* Read a sample */

      if (!db->isActive(iech)) continue;
      db->getCoordinatesPerSampleInPlace(iech, coor);
      if (dbgrid->getGrid().coordinateToIndicesInPlace(coor, indg0, true)) continue;
      double value = db->getArray(iech, juid);
      if (FFFF(value)) continue;

      /* Loop on the neighboring cells */

      for (int ic = 0; ic < count; ic++)
      {
        _neighboringCell(ndim, radius, ic, indg0, indg);
        int iad = dbgrid->getGrid().indiceToRank(indg);
        if (iad < 0 || iad >= nxyz) continue;

        if (oper == EStatOption::NUM)
        {
          dbgrid->updArray(iad, iptr, EOperator::ADD, 1.);
        }
        else if (oper == EStatOption::MEAN)
        {
          dbgrid->updArray(iad, iptn, EOperator::ADD, 1.);
          dbgrid->updArray(iad, iptr, EOperator::ADD, value);
        }
        else if (oper == EStatOption::VAR || oper == EStatOption::STDV)
        {
          dbgrid->updArray(iad, iptn, EOperator::ADD, 1.);
          dbgrid->updArray(iad, iptm, EOperator::ADD, value);
          dbgrid->updArray(iad, iptr, EOperator::ADD, value * value);
        }
        else if (oper == EStatOption::MINI)
        {
          if (value < dbgrid->getArray(iad, iptr))
            dbgrid->setArray(iad, iptr, value);
        }
        else if (oper == EStatOption::MAXI)
        {
          if (value > dbgrid->getArray(iad, iptr))
            dbgrid->setArray(iad, iptr, value);
        }
        else if (oper == EStatOption::MEDIAN)
        {
          medtab[nmed++] = value;
        }
        else if (oper == EStatOption::PLUS)
        {
          if (value > 0.) dbgrid->updArray(iad, iptr, EOperator::ADD, 1.);
        }
        else if (oper == EStatOption::MOINS)
        {
          if (value < 0.) dbgrid->updArray(iad, iptr, EOperator::ADD, 1.);
        }
        else if (oper == EStatOption::ZERO)
        {
          if (value == 0.) dbgrid->updArray(iad, iptr, EOperator::ADD, 1.);
        }
        else
        {
          value = 0.;
        }
      }
    }

    /* Normalization */

    if (oper == EStatOption::MEAN)
    {
      for (int i = 0; i < nxyz; i++)
      {
        double ratio = dbgrid->getArray(i, iptn);
        if (ratio <= 0.)
          dbgrid->setArray(i, iptr, TEST);
        else
          dbgrid->updArray(i, iptr, EOperator::DIVIDE, ratio);
      }
    }
    else if (oper == EStatOption::VAR || oper == EStatOption::STDV)
    {
      for (int i = 0; i < nxyz; i++)
      {
        double ratio = dbgrid->getArray(i, iptn);
        if (ratio <= 0.)
          dbgrid->setArray(i, iptr, TEST);
        else
        {
          double mean = dbgrid->getArray(i, iptm) / ratio;
          double value = dbgrid->getArray(i, iptr) / ratio - mean * mean;
          if (value < 0) value = 0.;
          if (oper == EStatOption::VAR)
            dbgrid->setArray(i, iptr, value);
          else
            dbgrid->setArray(i, iptr, sqrt(value));
        }
      }
    }
    else if (oper == EStatOption::MEDIAN)
    {
      for (int i = 0; i < nxyz; i++)
      {
        double value = (nmed > 0) ? medtab[nmed / 2] : TEST;
        dbgrid->setArray(i, iptr, value);
      }
    }
    else if (oper == EStatOption::MINI)
    {
      for (int i = 0; i < nxyz; i++)
      {
        double value = dbgrid->getArray(i, iptr);
        if (value == 1.e30) dbgrid->setArray(i, iptr, TEST);
      }
    }
    else if (oper == EStatOption::MAXI)
    {
      for (int i = 0; i < nxyz; i++)
      {
        double value = dbgrid->getArray(i, iptr);
        if (value == -1.e30) dbgrid->setArray(i, iptr, TEST);
      }
    }
  }

  /* Delete auxiliary attributes for local calculations */

  if ((oper == EStatOption::MEAN || oper == EStatOption::VAR || oper == EStatOption::STDV)
      && iptn > 0) dbgrid->deleteColumnByUID(iptn);
  if ((oper == EStatOption::VAR || oper == EStatOption::STDV) && iptm > 0)
    dbgrid->deleteColumnByUID(iptm);

  return 0;
}

/****************************************************************************/
/*!
 **  Evaluate the correlation
 **     Correl(Z1(x) , Z2(x))
 **
 ** \return  Array of the indices of pairs of samples (or VectorVectorInt())
 **
 ** \param[in]  db1          Db descriptor (first variable)
 ** \param[in]  db2          Db descriptor (second variable for flag.same=T)
 ** \param[in]  name1        Name of the first variable
 ** \param[in]  name2        Name of the second variable
 ** \param[in]  flagFrom1    Start numbering of indices from 1 if True
 ** \param[in]  verbose      Verbose flag
 **
 ** \remarks The two input Db must match exactly (same number of samples with
 ** \remarks same set of coordinates and same optional selection)
 **
 ** \remarks The returned Vector of Vector of integer 'indices' contain
 ** \remarks the set of indices of the pairs of samples.
 ** \remarks Its contents is i1,j1,i2,j2,...
 ** \remarks The indices are numbered starting from 0
 **
 *****************************************************************************/
VectorVectorInt correlationPairs(Db *db1,
                                 Db *db2,
                                 const String& name1,
                                 const String& name2,
                                 bool flagFrom1,
                                 bool verbose)
{
  VectorVectorInt indices;

  /* Initializations */

  if (db1 == nullptr) return indices;
  if (db2 == nullptr) return indices;
  if (db1->getNDim() != db2->getNDim() ||
      db1->getSampleNumber(true) != db2->getSampleNumber(true))
  {
    messerr("The two input 'db' are not compatible");
    return indices;
  }

  int nech = db1->getSampleNumber();
  int ndim = db1->getNDim();
  int shift = (flagFrom1) ? 1 : 0;
  SpaceRN space(ndim);
  SpaceTarget T1(&space);
  SpaceTarget T2(&space);

  /* Regular correlation */

  indices.resize(2);
  int nb = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    if (!db1->isActive(iech)) continue;
    double val1 = db1->getValue(name1, iech);
    if (FFFF(val1)) continue;
    double val2 = db2->getValue(name2, iech);
    if (FFFF(val2)) continue;

    indices[0].push_back(iech + shift);
    indices[1].push_back(iech + shift);
    nb++;
  }

  /* Messages */

  if (nb <= 0)
  {
    messerr("No sample found where all variables are defined");
    return indices;
  }
  else
  {
    if (verbose)
    {
      message("Total number of samples = %d\n", nech);
      message("Number of samples defined = %d\n", (int) nb);
    }
  }
  return indices;
}

/****************************************************************************/
/*!
 **  Evaluate the shifted correlation calculated as follows:
 **     Correl(Z1(x) , Z2(x+h))
 **
 ** \return  Vector of indices (or VectorVectorInt())
 **
 ** \param[in]  db           Db descriptor
 ** \param[in]  name1        Name of the first variable
 ** \param[in]  name2        Name of the second variable
 ** \param[in]  varioparam   pointer to a VarioParam structure
 ** \param[in]  ipas         Rank of the lag of interest
 ** \param[in]  idir         Rank of the direction of interest (within VarioParam)
 ** \param[in]  verbose      Verbose flag
 **
 ** \remarks The returned Vector of Vector of integer 'indices' contain
 ** \remarks the set of indices of the pairs of samples.
 ** \remarks Its contents is i1,j1,i2,j2,...
 ** \remarks The indices are numbered starting from 1
 **
 *****************************************************************************/
VectorVectorInt hscatterPairs(Db *db,
                              const String& name1,
                              const String& name2,
                              VarioParam *varioparam,
                              int ipas,
                              int idir,
                              bool verbose)
{
  VectorVectorInt indices;
  double dist = 0.;

  // Preliminary checks

  if (db == nullptr) return indices;
  if (varioparam == nullptr) return indices;
  if (idir < 0 || idir >= varioparam->getDirectionNumber()) return indices;

  /* Initializations */

  const DirParam dirparam = varioparam->getDirParam(idir);
  int nech = db->getSampleNumber();
  int ndim = db->getNDim();
  SpaceRN space(ndim);
  SpaceTarget T1(&space);
  SpaceTarget T2(&space);
  indices.resize(2);

  // Creating a local Vario structure from VarioParam (in order to constitute the BiTargetCheck list)

  Vario *vario = Vario::create(*varioparam);
  vario->setDb(db);
  if (vario->prepare()) return 1;

  // Local variables to speed up calculations
  bool hasSel = db->hasLocVariable(ELoc::SEL);

  int nb = 0;
  for (int iech = 0; iech < nech - 1; iech++)
  {
    if (hasSel && !db->isActive(iech)) continue;
    double val1 = db->getValue(name1, iech);
    if (FFFF(val1)) continue;
    db->getSampleAsST(iech, T1);

    for (int jech = iech + 1; jech < nech; jech++)
    {
      if (hasSel && !db->isActive(jech)) continue;
      double val2 = db->getValue(name2, jech);
      if (FFFF(val2)) continue;
      db->getSampleAsST(jech, T2);

      // Reject the point as soon as one BiTargetChecker is not correct
      if (!vario->keepPair(0, T1, T2, &dist)) continue;

      /* Get the rank of the lag */

      int ipasloc = dirparam.getLagRank(dist);
      if (IFFFF(ipasloc)) continue;
      if (ipas != ipasloc) continue;

      /* Point update */

      indices[0].push_back(iech);
      indices[1].push_back(jech);
      nb++;
    }
  }

  /* Messages */

  if (nb <= 0)
  {
    messerr("No sample found where all variables are defined");
  }
  else
  {
    if (verbose)
    {
      message("Total number of samples = %d\n", nech);
      message("Number of pairs used for translated correlation = %d\n", (int) nb);
    }
  }
  return indices;
}

/****************************************************************************/
/*!
 **  Identify samples from scatter plot when included within a polygon
 **
 ** \return  Error return code
 **
 ** \param[in]  db1          Db descriptor (first variable)
 ** \param[in]  db2          Db descriptor (second variable for flag.same=T)
 ** \param[in]  icol1        Rank of the first column
 ** \param[in]  icol2        Rank of the second column
 ** \param[in]  polygon      Polygons structure
 **
 ** \remarks The two input Db must match exactly (same number of samples with
 ** \remarks same set of coordinates and same optional selection)
 **
 *****************************************************************************/
int correlationIdentify(Db *db1,
                        Db *db2,
                        int icol1,
                        int icol2,
                        Polygons *polygon)
{
  if (db1 == nullptr) return (1);
  if (db2 == nullptr) return (1);
  int nech = db1->getSampleNumber();
  int number = 0;

  /* Correlation */

  for (int iech = 0; iech < nech; iech++)
  {
    if (!db1->isActive(iech)) continue;
    double val1 = db1->getArray(iech, icol1);
    if (FFFF(val1)) continue;
    double val2 = db2->getArray(iech, icol2);
    if (FFFF(val2)) continue;

    /* Check of the sample belongs to the polygon */

    VectorDouble coor(3, TEST);
    coor[0] = val1;
    coor[1] = val2;
    if (!polygon->inside(coor, false)) continue;

    /* Print the reference of the sample */

    if (number == 0) mestitle(0, "Samples selected from scatter plot");
    message("Sample #%d - Variable #1=%lf - Variable #2=%lf\n", iech + 1, val1, val2);
    number++;
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Evaluate the experimental conditional expectation
 **
 ** \param[in]  db1           Db descriptor (for target variable)
 ** \param[in]  db2           Db descriptor (for auxiliary variables)
 ** \param[in]  icol1         Rank of the target variable
 ** \param[in]  icol2         Rank of the explanatory variable
 ** \param[in]  mini          Minimum value for the explanaroty variable
 ** \param[in]  maxi          Maximum value for the explanaroty variable
 ** \param[in]  nclass        Number of classes
 ** \param[in]  verbose       Verbose flag
 **
 *****************************************************************************/
VectorVectorDouble condexp(Db *db1,
                           Db *db2,
                           int icol1,
                           int icol2,
                           double mini,
                           double maxi,
                           int nclass,
                           bool verbose)
{
  VectorVectorDouble xycond(2);
  xycond[0].resize(nclass);
  xycond[1].resize(nclass);
  VectorInt ncond(nclass,0);

  /* Loop on the samples */

  for (int iech = 0; iech < db1->getSampleNumber(); iech++)
  {
    if (!db1->isActive(iech)) continue;
    double val1 = db1->getArray(iech, icol1);
    if (FFFF(val1)) continue;
    double val2 = db2->getArray(iech, icol2);
    if (FFFF(val2)) continue;
    if (val2 < mini || val2 > maxi) continue;

    int rank = int((nclass - 1.) * (val2 - mini) / (maxi - mini));

    xycond[0][rank] += val1;
    xycond[1][rank] += val2;
    ncond[rank]++;
  }

  /* Normation */

  for (int iclass = 0; iclass < nclass; iclass++)
  {
    if (ncond[iclass] <= 0)
    {
      xycond[0][iclass] = TEST;
      xycond[1][iclass] = TEST;
    }
    else
    {
      xycond[0][iclass] /= ncond[iclass];
      xycond[1][iclass] /= ncond[iclass];
    }
  }

  /* Optional printout */

  if (verbose)
  {
    message("Experimental Conditional Expectation\n");
    for (int iclass = 0; iclass < nclass; iclass++)
    {
      if (ncond[iclass] > 0)
        message("Class %2d : V1=%lf V2=%lf\n", iclass + 1, xycond[0][iclass],
                xycond[1][iclass]);
    }
  }
  return xycond;
}

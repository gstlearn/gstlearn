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

#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Stats/Classical.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "Basic/AStringable.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Model/Model.hpp"

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

/**
 * Calculate the coefficients of the Deming regression (with 2 variables)
 * @param x Vector for the first variable
 * @param y Vector for the second variable
 * @param delta ratio of error variances (s_y^2 / s_x^2)
 * @return Vector of coefficients for the equation
 * @return y = beta[0] + beta[1] * x
 * @remark Both input vectors are assumed to contain valid values
 * @remark From: https://en.wikipedia.org/wiki/Deming_regression
 */
VectorDouble regrDeming(const VectorDouble &x,
                        const VectorDouble &y,
                        double delta)
{
  VectorDouble beta(2,TEST);
  int nech = (int) x.size();
  if (nech <= 1) return beta;

  double xmean = 0.;
  double ymean = 0.;
  for (int iech = 0; iech < nech; iech++)
  {
    xmean += x[iech];
    ymean += y[iech];
  }
  xmean /= (double) nech;
  ymean /= (double) nech;

  double sxx = 0.;
  double sxy = 0.;
  double syy = 0.;
  for (int iech = 0; iech < nech; iech++)
  {
    double deltax = (x[iech] - xmean);
    double deltay = (y[iech] - ymean);
    sxx += deltax * deltax;
    sxy += deltax * deltay;
    syy += deltay * deltay;
  }
  sxx /= (double) nech;
  sxy /= (double) nech;
  syy /= (double) nech;

  double T = syy - delta * sxx;
  beta[1] = (T + sqrt(T * T + 4. * delta * sxy * sxy)) / (2. * sxy);
  beta[0] = ymean - beta[1] * xmean;
  return beta;
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
  VectorInt iuids = db->getUIDs(names);
  return dbStatisticsMonoByUID(db, iuids, opers, flagIso, proba, vmin, vmax);
}

GSTLEARN_EXPORT Table dbStatisticsMonoT(Db *db,
                                        const VectorString &names,
                                        const std::vector<EStatOption> &opers,
                                        bool flagIso,
                                        double proba,
                                        double vmin,
                                        double vmax)
{
  VectorInt iuids = db->getUIDs(names);
  VectorDouble stats = dbStatisticsMonoByUID(db, iuids, opers, flagIso, proba, vmin, vmax);
  int nrows = (int) iuids.size();
  int ncols = (int) opers.size();
  Table table = Table();
  table.reset(nrows, ncols, stats, false, false);
  for (int irow=0; irow<nrows; irow++)
    table.setRowName(irow, db->getNameByUID(iuids[irow]));

  for (int icol=0; icol<ncols; icol++)
    table.setColumnName(icol, opers[icol].getDescr());

  return table;
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
VectorDouble dbStatisticsMonoByUID(Db *db,
                                   const VectorInt &iatts,
                                   const std::vector<EStatOption> &opers,
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
  VectorInt iuids = db->getUIDs(names);
  return dbStatisticsMultiByUID(db, iuids, flagIso);
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
VectorDouble dbStatisticsMultiByUID(Db *db, const VectorInt &iatts, bool flagIso)
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

  sstr << toMatrix(title, statOptionToName(opers), names, false, noper, natt, stats, true);

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
                       const std::vector<EStatOption> &opers,
                       bool flagIso,
                       bool flagCorrel,
                       const String &title,
                       const String &radix)
{
  VectorInt iuids = db->getUIDs(names);
  if (iuids.size() <= 0) return;
  dbStatisticsPrintByUID(db, iuids, opers, flagIso, flagCorrel, title, radix);
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
void dbStatisticsPrintByUID(const Db *db,
                            const VectorInt &iatts_arg,
                            const std::vector<EStatOption> &opers,
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

/**
 * Sphering procedure
 * @param X Input Data vector
 * @return The Sphering matrix (or nullptr if problem)
 *
 * @remark When performing the (forward) sphering, you must perform the following operation
 * @remark        X <- prodMatrix(X, S)
 */
MatrixRectangular* sphering(const AMatrix* X)
{
  if (X->isEmpty()) return nullptr;
  int nech = X->getNRows();
  int nvar = X->getNCols();

  AMatrix* TX = X->transpose();
  AMatrix* prod = prodMatrix(TX, X);
  prod->prodScalar(1. / (double) nech);

  VectorDouble eigen_values(nvar);
  VectorDouble eigen_vectors(nvar * nvar);
  if (matrix_eigen(prod->getValues().data(), nvar,
                   eigen_values.data(), eigen_vectors.data()))
    return nullptr;

  // Invert the sign of the second Eigen vector (for compatibility with R output)
  MatrixRectangular* S = new MatrixRectangular(nvar, nvar);
  S->setValues(eigen_vectors,true);
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

/****************************************************************************/
/*!
 **  Calculates the statistics of points within cells of a grid
 **
 ** \return  Error return code
 **
 ** \param[in]  db         Db for the points
 ** \param[in]  dbgrid     Db for the grid
 ** \param[in]  oper       A EStatOption item
 ** \param[in]  iatt       Rank of the first attribute
 ** \param[in]  jatt       Rank of the second attribute
 ** \param[in]  cuts       Array of cutoffs (optional)
 **
 *****************************************************************************/
VectorDouble dbStatisticsPerCell(Db *db,
                                 DbGrid *dbgrid,
                                 const EStatOption &oper,
                                 int iatt,
                                 int jatt,
                                 const VectorDouble &cuts)
{
  VectorDouble result;
  double z1 = 0.;
  double z2 = 0.;
  int nxyz = dbgrid->getSampleNumber();
  int ncut = (int) cuts.size();
  int ndim = dbgrid->getNDim();

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
  bool flag_maxi =false;

  /* Check the operator validity */

  if (! _operStatisticsCheck(oper, 1, 0, 1, 0, 1)) return result;

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
    cutval.resize(nxyz * ncut, 0.);

  /* Loop on the samples */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;

    /* Check the variable(s) */

    if (flag1)
    {
      z1 = db->getArray(iech, iatt);
      if (FFFF(z1)) continue;
    }
    if (flag2)
    {
      z2 = db->getArray(iech, jatt);
      if (FFFF(z2)) continue;
    }

    /* Check the location of the data in the grid */

    db->getCoordinatesInPlace(iech, coor);
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
      return result;
  }

  return 0;
}

/****************************************************************************/
/*!
 **  Calculate the multivariate statistics between different variables of a Db
 **
 ** \return  Resulting array (dimension; ncol if flag_mono or ncol*ncol)
 **
 ** \param[in]  db   Db structure
 ** \param[in]  oper A StatOption item
 ** \param[in]  cols Ranks of the variables
 ** \param[in]  flagMono 1 for a monovariate output
 ** \param[in]  verbose  1 for a verbose output
 **
 *****************************************************************************/
VectorDouble dbStatisticsMulti(Db *db,
                               const EStatOption &oper,
                               const VectorInt &cols,
                               bool flagMono,
                               bool verbose)
{
  VectorDouble result;

  /* Initializations */

  int nech = db->getSampleNumber();
  int ncol = static_cast<int>(cols.size());
  int ncol2 = ncol * ncol;

  /* Check that all variables are defined */

  for (int icol = 0; icol < ncol; icol++)
  {
    int jcol = cols[icol];
    if (!db->isColIdxValid(jcol))
    {
      messerr("Error: Variable %d is not defined", cols[icol]);
      return result;
    }
  }

  /* Check the validity of the operator */

  if (! _operStatisticsCheck(oper, 0, 1, 0, 0, 0)) return result;

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

  /* Optional printout */

  if (verbose)
  {
    message("\n");
    if (oper == EStatOption::NUM)
      print_matrix("Matrix of Number of defined samples", 0, 1, nx, ny, NULL,
                   num.data());
    else if (oper == EStatOption::MEAN)
      print_matrix("Matrix of Variable Means", 0, 1, nx, ny, NULL, m1.data());
    else if (oper == EStatOption::VAR)
      print_matrix("Matrix of Variable Variances", 0, 1, nx, ny, NULL, v12.data());
    else if (oper == EStatOption::CORR)
      print_matrix("Matrix of Variable Correlations", 0, 1, nx, ny, NULL, v12.data());
    else if (oper == EStatOption::STDV)
      print_matrix("Matrix of Variable Standard Deviations", 0, 1, nx, ny, NULL,
                   v12.data());
    else if (oper == EStatOption::MINI)
      print_matrix("Matrix of Variable Minima", 0, 1, nx, ny, NULL, mini.data());
    else if (oper == EStatOption::MAXI)
      print_matrix("Matrix of Variable Maxima", 0, 1, nx, ny, NULL, maxi.data());
    else if (oper == EStatOption::PLUS)
      print_matrix("Matrix of Number of positive samples", 0, 1, nx, ny, NULL,
                   plus.data());
    else if (oper == EStatOption::MOINS)
      print_matrix("Matrix of Number of negative samples", 0, 1, nx, ny, NULL,
                   moins.data());
    else if (oper == EStatOption::ZERO)
      print_matrix("Matrix of Number of zero samples", 0, 1, nx, ny, NULL,
                   zero.data());
    else
      messageAbort("This error should never happen");
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

  return result;
}

/****************************************************************************/
/*!
 **  Calculates the monovariate statistics within cells of a grid
 **
 ** \return  Error return code
 **
 ** \param[in]  db     Db for the points
 ** \param[in]  dbgrid Db for the grid
 ** \param[in]  oper   A EStatOption item
 ** \param[in]  cols   Ranks of the variables
 ** \param[in]  radius Neighborhood radius
 ** \param[in]  iptr0  Storage address (first variable)
 **
 *****************************************************************************/
int dbStatisticsInGrid(Db *db,
                       DbGrid *dbgrid,
                       const EStatOption &oper,
                       const VectorInt &cols,
                       int radius,
                       int iptr0)
{
  int iptm = -1;
  int iptn = -1;
  int nxyz = dbgrid->getSampleNumber();
  int ndim = dbgrid->getNDim();
  int ncol = (int) cols.size();
  int count = (int) pow(2. * radius + 1., (double) ndim);

  /* Check that all variables are defined */

  for (int icol = 0; icol < ncol; icol++)
  {
    int jcol = cols[icol];
    if (!db->isColIdxValid(jcol))
    {
      messerr("Error: Variable %d is not defined", cols[icol]);
      return 1;
    }
  }

  /* Check the validity of the requested function */

  if (! _operStatisticsCheck(oper, 0, 1, 0, 1, 0)) return 1;

  /* Create and initialize the new attributes */

  double valdef = 0.;
  if (oper == EStatOption::MINI) valdef = 1.e30;
  if (oper == EStatOption::MAXI) valdef = -1.e30;

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

  for (int icol = 0; icol < ncol; icol++)
  {

    /* Create the output attribute */

    int jcol = cols[icol];
    int iptr = iptr0 + icol;
    if (iptn > 0) db_attribute_init(dbgrid, 1, iptn, 0.);
    if (iptm > 0) db_attribute_init(dbgrid, 1, iptm, 0.);
    db_attribute_init(dbgrid, 1, iptr, valdef);

    /* Loop on the samples */

    int nmed = 0;
    for (int iech = 0; iech < db->getSampleNumber(); iech++)
    {

      /* Read a sample */

      if (!db->isActive(iech)) continue;
      db->getCoordinatesInPlace(iech, coor);
      if (dbgrid->getGrid().coordinateToIndicesInPlace(coor, indg0)) continue;
      double value = db->getArray(iech, jcol);
      if (FFFF(value)) continue;

      /* Loop on the neighboring cells */

      for (int ic = 0; ic < count; ic++)
      {
        _neighboringCell(ndim, radius, ic, indg0, indg);
        int iad = dbgrid->getGrid().indiceToRank(indg);
        if (iad < 0 || iad >= nxyz) continue;

        if (oper == EStatOption::NUM)
        {
          dbgrid->updArray(iad, iptr, 0, 1.);
        }
        else if (oper == EStatOption::MEAN)
        {
          dbgrid->updArray(iad, iptn, 0, 1.);
          dbgrid->updArray(iad, iptr, 0, value);
        }
        else if (oper == EStatOption::VAR || oper == EStatOption::STDV)
        {
          dbgrid->updArray(iad, iptn, 0, 1.);
          dbgrid->updArray(iad, iptm, 0, value);
          dbgrid->updArray(iad, iptr, 0, value * value);
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
          if (value > 0.) dbgrid->updArray(iad, iptr, 0, 1.);
        }
        else if (oper == EStatOption::MOINS)
        {
          if (value < 0.) dbgrid->updArray(iad, iptr, 0, 1.);
        }
        else if (oper == EStatOption::ZERO)
        {
          if (value == 0.) dbgrid->updArray(iad, iptr, 0, 1.);
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
          dbgrid->updArray(i, iptr, 3, ratio);
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


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
#include "Stats/Regression.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"

Regression::Regression()
    : AStringable(),
      _count(0),
      _nvar(0),
      _flagCst(true),
      _coeffs(),
      _variance(0.),
      _varres(0.)
{
}

Regression::Regression(const Regression& r)
    : AStringable(r),
      _count(r._count),
      _nvar(r._nvar),
      _flagCst(r._flagCst),
      _coeffs(r._coeffs),
      _variance(r._variance),
      _varres(r._varres)
{
}

Regression& Regression::operator=(const Regression& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _count = r._count;
    _nvar = r._nvar;
    _flagCst = r._flagCst;
    _coeffs = r._coeffs;
    _variance = r._variance;
    _varres = r._varres;
  }
  return *this;
}

Regression::~Regression()
{
}

String Regression::toString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);

  std::stringstream sstr;

  sstr << toTitle(1, "Linear Regression");
  sstr << "- Calculated on " << _count << " active values" << std::endl;

  int ecr = 0;
  int nvar = _nvar;
  if (_flagCst) nvar--;

  if (_flagCst)
    sstr << "- Constant term           = " << _coeffs[ecr++] << std::endl;
  for (int ivar = 0; ivar < nvar; ivar++)
    sstr << "- Explanatory Variable #" << ivar+1 << " = " << _coeffs[ecr++] << std::endl;

  sstr << "- Initial variance        = " << _variance << std::endl;
  sstr << "- Variance of residuals   = " << _varres << std::endl;

  return sstr.str();
}

bool _regressionCheck(Db *db1,
                      int icol0,
                      const VectorInt &icols,
                      int mode,
                      Db *db2,
                      const Model *model)
{
  int ncol = (int) icols.size();
  int nfex = db2->getLocNumber(ELoc::F);

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
                     int flagCst,
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
      if (flagCst) x[ecr++] = 1.;
      for (int icol = 0; icol < (int) icols.size(); icol++)
        x[ecr++] = db2->getArray(iech, icols[icol]);
      break;

    case 1:
      nfex = db2->getLocNumber(ELoc::F);
      *value = db1->getLocVariable(ELoc::Z,iech, 0);
      if (flagCst) x[ecr++] = 1.;
      for (int i = 0; i < nfex; i++)
        x[ecr++] = db2->getLocVariable(ELoc::F,iech, i);
      break;

    case 2:
      nbfl = model->getDriftNumber();
      *value = db1->getLocVariable(ELoc::Z,iech, 0);
      for (int i = 0; i < nbfl; i++)
         x[ecr++] = model->evalDrift(db2, iech, i, ECalcMember::LHS);
      break;
  }

  bool flagTest = false;
  for (int i = 0; i < (int) x.size() && !flagTest; i++)
    flagTest = FFFF(x[i]);
  return (FFFF(*value) || flagTest);
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
VectorDouble regressionDeming(const VectorDouble &x,
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
 **  Evaluate the regression
 **
 ** \return  Error return code
 **
 ** \param[in,out]  db1        Db descriptor (for target variable)
 ** \param[in]  nameResp       Name of the target variable
 ** \param[in]  nameAux        Vector of names of the explanatory variables
 ** \param[in]  mode           Type of calculation
 ** \li                        0 : standard multivariate case
 ** \li                        1 : using external drifts
 ** \li                        2 : using standard drift functions (in 'model')
 ** \param[in]  flagCst        The constant is added as explanatory variable
 ** \param[in]  db2            Db descriptor (for auxiliary variables)
 ** \param[in]  model          Model (only used for Drift functions if mode==2)
 **
 ** \remark  The flag_mode indicates the type of regression calculation:
 ** \remark  0 : V[icol] as a function of V[icols[i]]
 ** \remark  1 : Z1 as a function of the different Fi's
 **
 *****************************************************************************/
Regression regression(Db *db1,
                      const String &nameResp,
                      const VectorString &nameAux,
                      int mode,
                      bool flagCst,
                      Db *db2,
                      const Model *model)
{
  Regression regr;

  if (db1 == nullptr) return regr;
  if (db2 == nullptr) db2 = db1;

  int icol0 = db1->getUID(nameResp);
  VectorInt icols = db2->getUIDs(nameAux);

  int nfex = db2->getLocNumber(ELoc::F);
  int nech = db1->getSampleNumber();
  int ncol = (int) icols.size();
  int size = 0;
  switch (mode)
  {
    case 0:
      size = ncol;
      if (flagCst) size++;
      break;
    case 1:
      size = nfex;
      if (flagCst) size++;
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

    if (_regressionLoad(db1, db2, iech, icol0, icols, mode, flagCst, model,
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
  regr.setCount(number);
  regr.setNvar(size);
  regr.setFlagCst(flagCst);
  regr.setCoeffs(x);
  regr.setVariance(prod / number - mean * mean);

  /* Calculate the residuals */

  for (int i = 0; i < size; i++)
  {
    prod -= 2. * x[i] * b[i];
    for (int j = 0; j < size; j++)
      prod += x[i] * x[j] * a.getValue(i,j);
  }
  regr.setVarres(prod / number);

  return regr;
}

/****************************************************************************/
/*!
 **  Evaluate the regression
 **
 ** \return  Error return code
 **
 ** \param[in,out]  db1        Db descriptor (for target variable)
 ** \param[in]  iptr0          Storing address (already allocated)
 ** \param[in]  nameResp       Name of the target variable
 ** \param[in]  nameAux        Vector of names of the explanatory variables
 ** \param[in]  mode           Type of calculation
 ** \li                        0 : standard multivariate case
 ** \li                        1 : using external drifts
 ** \li                        2 : using standard drift functions (mode==2)
 ** \param[in]  flagCst        The constant is added as explanatory variable]
 ** \param[in]  db2            Db descriptor (for auxiliary variables)
 ** \param[in]  model          Model descriptor
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
int Regression::apply(Db *db1,
                      int iptr0,
                      const String &nameResp,
                      const VectorString &nameAux,
                      int mode,
                      bool flagCst,
                      Db *db2,
                      const Model *model)
{
  if (db2 == nullptr) db2 = db1;
  int icol0 = db1->getUID(nameResp);
  VectorInt icols;
  if (! nameAux.empty()) icols = db2->getUIDs(nameAux);

  /* Store the regression error at sample points */

  int size = (int) getCoeffs().size();
  double value = 0;
  VectorDouble x(size);

  for (int iech = 0; iech < db1->getSampleNumber(); iech++)
  {
    if (db1->isActive(iech))
    {
      /* Get the information for the current sample */

      if (_regressionLoad(db1, db2, iech, icol0, icols, mode, flagCst, model,
                          &value, x))
      {
        value = TEST;
      }
      else
      {
        for (int i = 0; i < size; i++)
          value -= x[i] * getCoeff(i);
      }
    }
    db1->setArray(iech, iptr0, value);
  }
  return 0;
}

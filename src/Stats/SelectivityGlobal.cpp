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
#include "Db/Db.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Vector.hpp"

#include <math.h>
#include <Stats/SelectivityGlobal.hpp>

SelectivityGlobal::SelectivityGlobal(int nclass)
    : Selectivity(nclass),
      _Test(nclass),
      _Qest(nclass),
      _Best(nclass),
      _Mest(nclass),
      _Tstd(nclass),
      _Qstd(nclass)
{
  ut_vector_fill(_Test, TEST);
  ut_vector_fill(_Qest, TEST);
  ut_vector_fill(_Best, TEST);
  ut_vector_fill(_Mest, TEST);
  ut_vector_fill(_Tstd, TEST);
  ut_vector_fill(_Qstd, TEST);
}

SelectivityGlobal::SelectivityGlobal(const SelectivityGlobal &m)
    : Selectivity(m),
      _Test(m._Test),
      _Qest(m._Qest),
      _Best(m._Best),
      _Mest(m._Mest),
      _Tstd(m._Tstd),
      _Qstd(m._Qstd)
{

}

SelectivityGlobal& SelectivityGlobal::operator=(const SelectivityGlobal &m)
{
  if (this != &m)
  {
    Selectivity::operator=(m);
    _Test = m._Test;
    _Qest = m._Qest;
    _Best = m._Best;
    _Mest = m._Mest;
    _Tstd = m._Tstd;
    _Qstd = m._Qstd;
  }
  return *this;
}

SelectivityGlobal::~SelectivityGlobal()
{

}
void SelectivityGlobal::setBest(int iclass, double best)
{
  if (! _isValid(iclass)) return;
  _Best[iclass] = best;
}

void SelectivityGlobal::setMest(int iclass, double mest)
{
  if (! _isValid(iclass)) return;
  _Mest[iclass] = mest;
}

void SelectivityGlobal::setQest(int iclass, double qest)
{
  if (! _isValid(iclass)) return;
  _Qest[iclass] = qest;
}

void SelectivityGlobal::setQstd(int iclass, double qstd)
{
  if (! _isValid(iclass)) return;
  _Qstd[iclass] = qstd;
}

void SelectivityGlobal::setTest(int iclass, double test)
{
  if (! _isValid(iclass)) return;
  _Test[iclass] = test;
}

void SelectivityGlobal::setTstd(int iclass, double tstd)
{
  if (! _isValid(iclass)) return;
  _Tstd[iclass] = tstd;
}

double SelectivityGlobal::getBest(int iclass) const
{
  if (! _isValid(iclass)) return(TEST);
  return _Best[iclass];
}
double SelectivityGlobal::getMest(int iclass) const
{
  if (! _isValid(iclass)) return(TEST);
  return _Mest[iclass];
}
double SelectivityGlobal::getQest(int iclass) const
{
  if (! _isValid(iclass)) return(TEST);
  return _Qest[iclass];
}
double SelectivityGlobal::getQstd(int iclass) const
{
  if (! _isValid(iclass)) return(TEST);
  return _Qstd[iclass];
}
double SelectivityGlobal::getTest(int iclass) const
{
  if (! _isValid(iclass)) return(TEST);
  return _Test[iclass];
}
double SelectivityGlobal::getTstd(int iclass) const
{
  if (! _isValid(iclass)) return(TEST);
  return _Tstd[iclass];
}

bool SelectivityGlobal::_isValid(int iclass) const
{
  if (iclass < 0 || iclass >= getNCuts())
  {
    mesArg("Selectivity Class", iclass, getNCuts());
    return false;
  }
  return true;
}

SelectivityGlobal* SelectivityGlobal::createFromDb(const VectorDouble& zcuts, const Db* db)
{
  SelectivityGlobal* calcul = nullptr;
  if (zcuts.empty())
  {
    messerr("You must define 'zcuts'");
    return calcul;
  }
  if (db == nullptr)
  {
    messerr("You must provide a valid 'Db'");
    return calcul;
  }
  if (db->getVariableNumber() != 1)
  {
    messerr("The 'Db' must contain a single variable");
    return calcul;
  }

  // Extract the array of data and weights

  VectorDouble tab = db->getColumnByLocator(ELoc::Z, 0, true);
  VectorDouble wtab;
  if (db->hasWeight())
    wtab = db->getColumnByLocator(ELoc::W, 0, true);

  return createFromTab(zcuts, tab, wtab);
}

SelectivityGlobal* SelectivityGlobal::createFromTab(const VectorDouble& zcuts,
                                        const VectorDouble& tab,
                                        const VectorDouble& weights)
{
  SelectivityGlobal* calcul = nullptr;
  if (zcuts.empty())
  {
    messerr("You must define 'cuts'");
    return calcul;
  }
  if (tab.empty())
  {
    messerr("You must provide a valid 'tab'");
    return calcul;
  }
  int nech = (int) tab.size();
  VectorDouble wtab = weights;
  if (wtab.empty())
  {
    wtab.resize(nech,1);
  }
  else
  {
    if (nech != (int) wtab.size())
    {
      messerr("Arguments 'tab' and 'weights' should have same dimension");
      return calcul;
    }
  }

  // Allocate the returned structrue

  int nclass = (int) zcuts.size();
  calcul = new SelectivityGlobal(nclass);

  // Perform calculations

  double tonref = ut_vector_cumul(wtab);
  for (int iclass = 0; iclass < nclass; iclass++)
  {
    double coupure = zcuts[iclass];

    double tonnage = 0.;
    double metal = 0.;
    for (int iech = 0; iech < nech; iech++)
    {
      double x = tab[iech];
      double w = wtab[iech];
      double indic = (x >= coupure);

      tonnage += w * indic;
      metal   += x * w * indic;
    }

    tonnage /= tonref;
    metal   /= tonref;
    double benefit = metal - tonnage * coupure;
    double grade;
    if (tonnage <= 0.)
    {
      grade = TEST;
    }
    else
    {
      grade = metal / tonnage;
    }
    calcul->setZcut(iclass, coupure);
    calcul->setTest(iclass, tonnage);
    calcul->setQest(iclass, metal);
    calcul->setBest(iclass, benefit);
    calcul->setMest(iclass, grade);
  }
  return calcul;
}

/****************************************************************************/
/*!
 **  From the cutoff, tonnage and metal quantity, derive
 **  the conventional benefit and the average recovered grade
 **
 *****************************************************************************/
void SelectivityGlobal::calculateBenefitAndGrade()
{
  int iclass;
  double zval, tval, qval;
  int nclass = getNCuts();

  for (iclass = 0; iclass < nclass; iclass++)
  {
    zval = getZcut(iclass);
    tval = getTest(iclass);
    qval = getQest(iclass);
    setBest(iclass, qval - zval * tval);
    setMest(iclass ,(ABS(tval) < EPSILON6) ? TEST : qval / tval);
  }
}

/****************************************************************************/
/*!
 **  Calculate and print the Gini index
 **
 *****************************************************************************/
void SelectivityGlobal::dumpGini()
{
  int nclass = getNCuts();

  double gini = 1.;
  for (int iclass = 0; iclass < nclass - 1; iclass++)
    gini -= ((getTest(iclass)
        - getTest(iclass + 1))
             * (getQest(iclass) + getQest(iclass + 1)));

  message("Gini calculated on %d classes\n", nclass);
  message("Value of the Gini index = %lf\n", gini);
}

/*****************************************************************************/
/*!
 **  Correct the order relationship for Tonnage
 **
 *****************************************************************************/
void SelectivityGlobal::correctTonnageOrder()
{
  int nclass = getNCuts();
  VectorDouble ta(nclass);
  VectorDouble tb(nclass);

  for (int iclass = nclass - 1; iclass >= 0; iclass--)
  {
    double auxval = getTest(iclass);
    if (iclass < nclass - 1) auxval = MAX(ta[iclass + 1], auxval);
    ta[iclass] = MIN(1., MAX(0., auxval));
  }

  for (int iclass = 0; iclass < nclass; iclass++)
  {
    double auxval = getTest(iclass);
    if (iclass > 0) auxval = MIN(tb[iclass - 1], auxval);
    tb[iclass] = MAX(0., MIN(1., auxval));
  }

  for (int iclass = 0; iclass < nclass; iclass++)
    setTest(iclass, 0.5 * (ta[iclass] + tb[iclass]));
}

/****************************************************************************/
/*!
 **  Interpolate the Grade-Tonnage curves
 **
 ** \param[in] zcuts    Array of cutoffs
 ** \param[in] calest   Input Selectivity
 ** \param[in] verbose  Verbose flag
 **
 *****************************************************************************/
SelectivityGlobal* SelectivityGlobal::createInterpolation(const VectorDouble& zcuts,
                                              const SelectivityGlobal& calest,
                                              bool verbose)
{
  double tval, qval;

  int nclass = calest.getNCuts();
  int ncuts = (int) zcuts.size();

  SelectivityGlobal* calcut = new SelectivityGlobal(ncuts);
  for (int icut = 0; icut < ncuts; icut++)
  {
    double zval = zcuts[icut];
    calcut->setZcut(icut, zval);

    /* Find interval to which cutoffs belongs */

    int iclass = -1;
    for (int jclass = 0; jclass < nclass - 1 && iclass < 0; jclass++)
    {
      double valmin = MIN(calest.getZcut(jclass), calest.getZcut(jclass + 1));
      double valmax = MAX(calest.getZcut(jclass), calest.getZcut(jclass + 1));
      if (zval >= valmin && zval <= valmax) iclass = jclass;
    }

    if (iclass >= 0 && iclass < nclass)
    {

      /* Assuming that cutoffs belongs to the interval the class 'iclass' */

      double zi0 = calest.getZcut(iclass);
      double zi1 = (iclass + 1 > nclass - 1) ? 0. : calest.getZcut(iclass + 1);
      double ti0 = calest.getTest(iclass);
      double ti1 = (iclass + 1 > nclass - 1) ? 0. : calest.getTest(iclass + 1);
      double qi0 = calest.getQest(iclass);
      double qi1 = calest.getQest(iclass + 1);
      calcut->_interpolateInterval(zval, zi0, zi1, ti0, ti1, qi0, qi1, &tval, &qval);
      calcut->setTest(icut, tval);
      calcut->setQest(icut, qval);
    }
    else
    {
      calcut->setTest(icut, 0.);
      calcut->setQest(icut, 0.);
    }
  }

  calcut->calculateBenefitAndGrade();
  if (verbose) calcut->dumpGini();

  return calcut;
}

/*****************************************************************************/
/*!
 **  Interpolate the QT within an interval
 **
 ** \param[in]  zval     Cutoff value
 ** \param[in]  zi0      Lower cutoff of the interval
 ** \param[in]  zi1      Upper cutoff of the interval
 ** \param[in]  ti0      Lower tonnage of the interval
 ** \param[in]  ti1      Upper tonnage of the interval
 ** \param[in]  qi0      Lower metal quantity of the interval
 ** \param[in]  qi1      Upper metal quantity of the interval
 **
 ** \param[out] tval     Tonnage for the current cutoff
 ** \param[out] qval     Metal quantity for the current cutoff
 **
 *****************************************************************************/
void SelectivityGlobal::_interpolateInterval(double zval,
                                       double zi0,
                                       double zi1,
                                       double ti0,
                                       double ti1,
                                       double qi0,
                                       double qi1,
                                       double *tval,
                                       double *qval,
                                       double tol)
{
  double dzi = zi1 - zi0;
  double dti = ti1 - ti0;
  double zmoy = (qi1 - qi0) / (ti1 - ti0);
  double aa0 = (zi1 - zmoy) / (zmoy - zi0);

  if (ABS(zval - zi0) < tol)
  {
    (*tval) = ti0;
    (*qval) = qi0;
    return;
  }

  if (ABS(zval - zi1) < tol)
  {
    (*tval) = ti1;
    (*qval) = qi1;
    return;
  }

  double u = (zval - zi0) / dzi;
  (*tval) = (u <= 0.) ? ti0 : ti0 + dti * pow(u, 1. / aa0);
  (*qval) = (u <= 0.) ? qi0 :
      qi0 + zi0 * ((*tval) - ti0)
      + dzi * dti * pow(u, 1. + 1. / aa0) / (1. + aa0);
}


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
#include "Stats/Selectivity.hpp"
#include "Db/Db.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Vector.hpp"

#include <math.h>

Selectivity::Selectivity(int nclass)
  : _nClass(nclass)
  , _Zcut(nclass)
  , _Test(nclass)
  , _Qest(nclass)
  , _Best(nclass)
  , _Mest(nclass)
  , _Tstd(nclass)
  , _Qstd(nclass)
{
  ut_vector_fill(_Zcut, TEST);
  ut_vector_fill(_Test, TEST);
  ut_vector_fill(_Qest, TEST);
  ut_vector_fill(_Best, TEST);
  ut_vector_fill(_Mest, TEST);
  ut_vector_fill(_Tstd, TEST);
  ut_vector_fill(_Qstd, TEST);
}

Selectivity::Selectivity(const Selectivity &m)
  : _nClass(m._nClass)
  , _Zcut(m._Zcut)
  , _Test(m._Test)
  , _Qest(m._Qest)
  , _Best(m._Best)
  , _Mest(m._Mest)
  , _Tstd(m._Tstd)
  , _Qstd(m._Qstd)
{

}

Selectivity& Selectivity::operator=(const Selectivity &m)
{
  if (this != &m)
  {
    _nClass = m._nClass;;
    _Zcut = m._Zcut;
    _Test = m._Test;
    _Qest = m._Qest;
    _Best = m._Best;
    _Mest = m._Mest;
    _Tstd = m._Tstd;
    _Qstd = m._Qstd;
  }
  return *this;
}

Selectivity::~Selectivity()
{

}
void Selectivity::setBest(int iclass, double best)
{
  if (! _isValid(iclass)) return;
  _Best[iclass] = best;
}

void Selectivity::setMest(int iclass, double mest)
{
  if (! _isValid(iclass)) return;
  _Mest[iclass] = mest;
}

void Selectivity::setQest(int iclass, double qest)
{
  if (! _isValid(iclass)) return;
  _Qest[iclass] = qest;
}

void Selectivity::setQstd(int iclass, double qstd)
{
  if (! _isValid(iclass)) return;
  _Qstd[iclass] = qstd;
}

void Selectivity::setTest(int iclass, double test)
{
  if (! _isValid(iclass)) return;
  _Test[iclass] = test;
}

void Selectivity::setTstd(int iclass, double tstd)
{
  if (! _isValid(iclass)) return;
  _Tstd[iclass] = tstd;
}

void Selectivity::setZcut(int iclass, double zcut)
{
  if (! _isValid(iclass)) return;
  _Zcut[iclass] = zcut;
}

double Selectivity::getBest(int iclass) const
{
  if (! _isValid(iclass)) return(TEST);
  return _Best[iclass];
}
double Selectivity::getMest(int iclass) const
{
  if (! _isValid(iclass)) return(TEST);
  return _Mest[iclass];
}
double Selectivity::getQest(int iclass) const
{
  if (! _isValid(iclass)) return(TEST);
  return _Qest[iclass];
}
double Selectivity::getQstd(int iclass) const
{
  if (! _isValid(iclass)) return(TEST);
  return _Qstd[iclass];
}
double Selectivity::getTest(int iclass) const
{
  if (! _isValid(iclass)) return(TEST);
  return _Test[iclass];
}
double Selectivity::getTstd(int iclass) const
{
  if (! _isValid(iclass)) return(TEST);
  return _Tstd[iclass];
}
double Selectivity::getZcut(int iclass) const
{
  if (! _isValid(iclass)) return(TEST);
  return _Zcut[iclass];
}

bool Selectivity::_isValid(int iclass) const
{
  if (iclass < 0 || iclass >= _nClass)
  {
    mesArg("Selectivity Class", iclass, _nClass);
    return false;
  }
  return true;
}

Selectivity* Selectivity::createFromDb(const VectorDouble& cuts, const Db* db)
{
  Selectivity* calcul = nullptr;
  if (cuts.empty())
  {
    messerr("You must define 'cuts'");
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
    tab = db->getColumnByLocator(ELoc::W, 0, true);

  return createFromTab(cuts, tab, wtab);
}

Selectivity* Selectivity::createFromTab(const VectorDouble& cuts,
                                        const VectorDouble& tab,
                                        const VectorDouble& weights)
{
  Selectivity* calcul = nullptr;
  if (cuts.empty())
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

  int nclass = (int) cuts.size();
  calcul = new Selectivity(nclass);

  // Perform calculations

  double tonref = ut_vector_cumul(wtab);
  for (int iclass = 0; iclass < nclass; iclass++)
  {
    double coupure = cuts[iclass];

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
void Selectivity::calculateBenefitGrade()
{
  int iclass;
  double zval, tval, qval;
  int nclass = getNClass();

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
void Selectivity::dumpGini()
{
  int nclass = getNClass();

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
void Selectivity::correctTonnageOrder()
{
  int nclass = getNClass();
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

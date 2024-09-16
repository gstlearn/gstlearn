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
#include "Model/Constraints.hpp"
#include "Model/ConsItem.hpp"
#include "Basic/Utilities.hpp"

#include <math.h>

Constraints::Constraints(double constantSillValue, const VectorDouble& constantSills)
    : AStringable(),
      _constantSillValue(constantSillValue),
      _constantSills(constantSills),
      _consItems()
{
}

Constraints::Constraints(const Constraints &m)
    : AStringable(m),
      _constantSillValue(m._constantSillValue),
      _constantSills(m._constantSills),
      _consItems()
{
  for (const auto& e: m._consItems)
  {
    _consItems.push_back(e->clone());
  }
}

Constraints& Constraints::operator=(const Constraints &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _constantSillValue = m._constantSillValue;
    _constantSills = m._constantSills;
    for (const auto& e: m._consItems)
    {
      _consItems.push_back(e->clone());
    }
  }
  return *this;
}

Constraints::~Constraints()
{
  for (int i=0; i<(int) _consItems.size(); i++)
    delete _consItems[i];
}

void Constraints::addItem(const ConsItem* item)
{
  _consItems.push_back(item->clone());
}

void Constraints::addItemFromParamId(const EConsElem &elem,
                                     int icov,
                                     int iv1,
                                     int iv2,
                                     const EConsType &type,
                                     double value)
{
  ConsItem* item = ConsItem::createFromParamId(icov, elem, type, value, 0, iv1, iv2);
  _consItems.push_back(item);
}

String Constraints::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  int nitem = static_cast<int> (_consItems.size());

  if (nitem > 0)
    sstr << toTitle(0, "Constraints to be fulfilled in Fitting procedure");

  for (int i = 0; i < nitem; i++)
  {
    sstr << "Constraint #" << i + 1 << std::endl;
    sstr << _consItems[i]->toString(strfmt);
  }

  if (! FFFF(getConstantSillValue()))
    sstr << "- Constraints on the sills =" << getConstantSillValue() << std::endl;

  return sstr.str();
}

int Constraints::isDefinedForSill() const
{
  if (_consItems.size() <= 0) return(0);
  for (int i=0; i<(int) _consItems.size(); i++)
  {
    if (_consItems[i]->getType() == EConsElem::SILL) return(1);
  }
  return(0);
}

void Constraints::modifyConstraintsForSill()
{
  for (int i=0; i<(int) getConsItemNumber(); i++)
  {
    const ConsItem* consitem = getConsItems(i);
    if (consitem->getType() != EConsElem::SILL) continue;
    if (consitem->getValue() > 0) setValue(i,sqrt(consitem->getValue()));
  }
}

void Constraints::setValue(int item, double value)
{
  _consItems[item]->setValue(value);
}

void Constraints::expandConstantSill(int nvar)
{
  _constantSills.resize(nvar,_constantSillValue);
}

bool Constraints::isConstraintSillDefined() const
{
  if (! FFFF(_constantSillValue)) return true;
  if (! _constantSills.empty()) return true;
  return false;
}

/****************************************************************************/
/*!
 **  Print the Auto Fitting Constraints Structure
 **
 ** \param[in]  constraints  Constraints structure
 **
 *****************************************************************************/
void constraints_print(const Constraints& constraints)
{
  constraints.display();
}

/****************************************************************************/
/*!
 **  If a constraint concerns a sill, take its square root
 **  as it corresponds to a constraints on AIC (not on a sill directly)
 **  due to the fact that it will be processed in FOXLEG (not in GOULARD)
 **  This transform only makes sense for MONOVARIATE case (the test should
 **  have been performed beforehand)
 **
 ** \return Error code (if the sill constraint is negative)
 **
 ** \param[in]  constraints  Constraints structure
 **
 *****************************************************************************/
int modify_constraints_on_sill(Constraints& constraints)

{
  int ncons = (int)constraints.getConsItemNumber();
  for (int i = 0; i < ncons; i++)
  {
    const ConsItem* consitem = constraints.getConsItems(i);
    if (consitem->getType() != EConsElem::SILL) continue;
    if (consitem->getValue() < 0) return (1);
    constraints.setValue(i, sqrt(consitem->getValue()));

    // For constraints on the Sill in monovariate case,
    // Add a constraints on AIC for lower bound
    if (consitem->getIV1() == 0 && consitem->getIV2() == 0 &&
        consitem->getIcase() == EConsType::UPPER)
    {
      ConsItem* consjtem = new ConsItem(*consitem);
      consjtem->setValue(-consjtem->getValue());
      consjtem->setIcase(EConsType::LOWER);
      constraints.addItem(consjtem);
    }
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Return the constraint value (if defined) or TEST
 **
 ** \return Returned value or TEST
 **
 ** \param[in,out]  constraints  Constraints structure
 ** \param[in]      icase        Parameter type (EConsType)
 ** \param[in]      igrf         Rank of the Gaussian Random Function
 ** \param[in]      icov         Rank of the structure (starting from 0)
 ** \param[in]      icons        Type of the constraint (EConsElem)
 ** \param[in]      iv1          Rank of the first variable
 ** \param[in]      iv2          Rank of the second variable
 **
 *****************************************************************************/
double constraints_get(const Constraints& constraints,
                       const EConsType& icase,
                       int igrf,
                       int icov,
                       const EConsElem& icons,
                       int iv1,
                       int iv2)
{
  if (!constraints.isDefined()) return (TEST);

  for (int i = 0; i < (int)constraints.getConsItemNumber(); i++)
  {
    const ConsItem* item = constraints.getConsItems(i);
    if (item->getIGrf() != igrf || item->getICov() != icov ||
        item->getType() != icons || item->getIV1() != iv1)
      continue;
    if (icons == EConsElem::SILL && item->getIV2() != iv2) continue;

    if (item->getIcase() == EConsType::EQUAL)
    {
      if (icase == EConsType::LOWER || icase == EConsType::UPPER)
        return (item->getValue());
    }
    else
    {
      if (icase == item->getIcase()) return (item->getValue());
    }
  }
  return (TEST);
}

/****************************************************************************/
/*!
 **  Add constraints to the Option_AutoFit structure
 **
 ** \return Error return code
 **
 ** \param[in]  constraints  Constraints structure
 ** \param[in]  constantSill Constant value for the Sill as a constraint
 **
 *****************************************************************************/
int add_sill_constraints(Constraints& constraints, double constantSill)
{
  constraints.setConstantSillValue(constantSill);

  return (0);
}

/****************************************************************************/
/*!
 **  Add constraints (all equal to 1) to the Option_AutoFit structure
 **
 ** \return Error return code
 **
 ** \param[in]  constraints   Constraints structure
 **
 *****************************************************************************/
int add_unit_sill_constraints(Constraints& constraints)
{
  constraints.setConstantSillValue(1.);
  return (0);
}

/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
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
  for (auto e: m._consItems)
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
    for (auto e: m._consItems)
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

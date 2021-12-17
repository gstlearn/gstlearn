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
#include "Model/Constraints.hpp"
#include "Model/ConsItem.hpp"
#include "geoslib_f.h"

#include <math.h>

Constraints::Constraints()
    : AStringable(),
      _consItems()
{
}

Constraints::Constraints(const Constraints &m)
    : AStringable(m),
      _consItems()
{
  for (auto e: m._consItems)
  {
    _consItems.push_back(dynamic_cast<ConsItem*>(e->clone()));
  }
}

Constraints& Constraints::operator=(const Constraints &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    for (auto e: m._consItems)
    {
      _consItems.push_back(dynamic_cast<ConsItem*>(e->clone()));
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
  _consItems.push_back(dynamic_cast<ConsItem*>(item->clone()));
}

String Constraints::toString(int level) const
{
  std::stringstream sstr;
  int nitem = static_cast<int> (_consItems.size());

  if (nitem > 0)
    sstr << toTitle(0, "Constraints to be fulfilled in Fitting procedure");

  for (int i = 0; i < nitem; i++)
  {
    sstr << "Constraint #" << i + 1 << std::endl;
    sstr << _consItems[i]->toString(level);
  }
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

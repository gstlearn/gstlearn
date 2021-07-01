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

Constraints::Constraints()
    : _consItems()
{
}

Constraints::Constraints(const Constraints &m)
    : _consItems(m._consItems)
{

}

Constraints& Constraints::operator=(const Constraints &m)
{
  if (this != &m)
  {
    _consItems = m._consItems;
  }
  return *this;
}

Constraints::~Constraints()
{
  for (int i=0; i<(int) _consItems.size(); i++)
    delete _consItems[i];
}

ConsItem* Constraints::addItem()
{
  int nitems = _consItems.size();
  _consItems.resize(nitems + 1);
  ConsItem* consitem = new (ConsItem);
  _consItems[nitems] = consitem;
  return consitem;
}

String Constraints::toString(int level) const
{
  std::stringstream sstr;
  int nitem = _consItems.size();

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
    if (_consItems[i]->getType() == CONS_SILL) return(1);
  }
  return(0);
}

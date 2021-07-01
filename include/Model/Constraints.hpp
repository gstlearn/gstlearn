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
#pragma once

#include "Basic/Vector.hpp"
#include "Model/ConsItem.hpp"
#include "Basic/AStringable.hpp"

class Constraints : public AStringable
{
public:
  Constraints();
  Constraints(const Constraints &m);
  Constraints& operator= (const Constraints &m);
  virtual ~Constraints();

  ConsItem* addItem();
  virtual String toString(int level = 0) const override;

  int isDefined() const { return _consItems.size() > 0; }
  int isDefinedForSill() const;
  int getConsItemNumber() const { return _consItems.size(); }

  const std::vector<ConsItem*>& getConsItems() const { return _consItems; }
  ConsItem* getConsItems(int i) const { return _consItems[i]; }

private:
  std::vector<ConsItem *> _consItems;
};

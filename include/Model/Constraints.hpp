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

#include "gstlearn_export.hpp"

// WARNING: Make this include list as small as possible!
#include "Basic/AStringable.hpp"

#include <vector>

class ConsItem;

class GSTLEARN_EXPORT Constraints : public AStringable
{
public:
  Constraints();
  Constraints(const Constraints &m);
  Constraints& operator= (const Constraints &m);
  virtual ~Constraints();

  void addItem(const ConsItem* item);
  virtual String toString(int level = 0) const override;

  int isDefined() const { return _consItems.size() > 0; }
  int isDefinedForSill() const;
  int getConsItemNumber() const { return static_cast<int>(_consItems.size()); }

  const std::vector<ConsItem*>& getConsItems() const { return _consItems; }
  const ConsItem* getConsItems(int i) const { return _consItems[i]; }
  void modifyConstraintsForSill();

  // Pipe to Consitem
  void setValue(int item, double value);

private:
  std::vector<ConsItem *> _consItems;
};

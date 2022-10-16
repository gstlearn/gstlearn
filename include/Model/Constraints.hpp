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

#include "Enum/EConsElem.hpp"
#include "Enum/EConsType.hpp"

#include "Basic/AStringable.hpp"

#include <vector>

class ConsItem;

class GSTLEARN_EXPORT Constraints : public AStringable
{
public:
  Constraints(double constantSillValue = TEST, const VectorDouble& constantSills = VectorDouble());
  Constraints(const Constraints &m);
  Constraints& operator= (const Constraints &m);
  virtual ~Constraints();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  void addItem(const ConsItem* item);
  void addItemFromParamId(const EConsElem &elem = EConsElem::UNKNOWN,
                          int icov = 0,
                          int iv1 = 0,
                          int iv2 = 0,
                          const EConsType &type = EConsType::DEFAULT,
                          double value = 0.);

  int isDefined() const { return _consItems.size() > 0; }
  int isDefinedForSill() const;
  int getConsItemNumber() const { return static_cast<int>(_consItems.size()); }

  const std::vector<ConsItem*>& getConsItems() const { return _consItems; }
  const ConsItem* getConsItems(int i) const { return _consItems[i]; }
  void modifyConstraintsForSill();

  double getConstantSillValue() const { return _constantSillValue; }
  const VectorDouble& getConstantSills() const { return _constantSills; }
  double getConstantSills(int ivar) const { return _constantSills[ivar]; }
  void setConstantSillValue(double value) { _constantSillValue = value; }
  void setConstantSills(const VectorDouble& constantSills) { _constantSills = constantSills; }
  void expandConstantSill(int nvar);
  bool isConstraintSillDefined() const;

  // Pipe to Consitem
  void setValue(int item, double value);

private:
  double _constantSillValue;       /* Constant Sill as a constraint */
  VectorDouble _constantSills;     /* Array of constant Sills (expanded to the variables) */
  std::vector<ConsItem *> _consItems;
};

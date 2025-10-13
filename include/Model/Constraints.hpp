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
#pragma once

#include "gstlearn_export.hpp"

#include "Enum/EConsElem.hpp"
#include "Enum/EConsType.hpp"

#include "Basic/AStringable.hpp"

#include <vector>

namespace gstlrn
{ 
class ConsItem;

class GSTLEARN_EXPORT Constraints : public AStringable
{
public:
  Constraints(double constantSillValue = TEST, const VectorDouble& constantSills = VectorDouble());
  Constraints(const Constraints &m);
  Constraints& operator= (const Constraints &m);
  virtual ~Constraints();

  String toString(const AStringFormat* strfmt = nullptr) const override;

  void addItem(const ConsItem* item);
  void addItemFromParamId(const EConsElem &elem = EConsElem::fromKey("UNKNOWN"),
                          Id icov = 0,
                          Id iv1 = 0,
                          Id iv2 = 0,
                          const EConsType &type = EConsType::fromKey("DEFAULT"),
                          double value = 0.);

  Id isDefined() const { return _consItems.size() > 0; }
  Id isDefinedForSill() const;
  Id getNConsItem() const { return static_cast<Id>(_consItems.size()); }

  const std::vector<ConsItem*>& getConsItems() const { return _consItems; }
  const ConsItem* getConsItems(Id i) const { return _consItems[i]; }
  void modifyConstraintsForSill();

  double getConstantSillValue() const { return _constantSillValue; }
  const VectorDouble& getConstantSills() const { return _constantSills; }
  double getConstantSills(Id ivar) const { return _constantSills[ivar]; }
  void setConstantSillValue(double value) { _constantSillValue = value; }
  void setConstantSills(const VectorDouble& constantSills) { _constantSills = constantSills; }
  void expandConstantSill(Id nvar);
  bool isConstraintSillDefined() const;

  // Pipe to Consitem
  void setValue(Id item, double value);

private:
  double _constantSillValue;       /* Constant Sill as a constraint */
  VectorDouble _constantSills;     /* Vector of constant Sills (expanded to the number of variables) */
  std::vector<ConsItem *> _consItems;
};

GSTLEARN_EXPORT double constraints_get(const Constraints& constraints,
                                       const EConsType& icase,
                                       Id igrf,
                                       Id icov,
                                       const EConsElem& icons,
                                       Id v1,
                                       Id v2);
GSTLEARN_EXPORT void constraints_print(const Constraints& constraints);
GSTLEARN_EXPORT Id modify_constraints_on_sill(Constraints& constraints);
GSTLEARN_EXPORT Id add_unit_sill_constraints(Constraints& constraints);
}
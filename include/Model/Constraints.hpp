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
  void addItemFromParamId(const EConsElem &elem = EConsElem::fromKey("UNKNOWN"),
                          int icov = 0,
                          int iv1 = 0,
                          int iv2 = 0,
                          const EConsType &type = EConsType::fromKey("DEFAULT"),
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
  VectorDouble _constantSills;     /* Vector of constant Sills (expanded to the number of variables) */
  std::vector<ConsItem *> _consItems;
};

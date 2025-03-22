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

#include "Basic/ICloneable.hpp"
#include "Covariances/TabNoStat.hpp"
#include "Enum/EConsElem.hpp"

class GSTLEARN_EXPORT TabNoStatSills : public TabNoStat
{
  public:
  TabNoStatSills();
  TabNoStatSills(const TabNoStatSills &m);
  TabNoStatSills& operator= (const TabNoStatSills &m);
  virtual ~TabNoStatSills();

  IMPLEMENT_CLONING(TabNoStatSills)

  bool isDefinedForVariance() const;
  int getNSills()  const;
 
  String toString(const AStringFormat* strfmt = nullptr) const override;
  String toStringInside(const AStringFormat* strfmt = nullptr,int i = 0) const;
protected:
private:
  bool _isValid(const EConsElem &econs) const override;

};

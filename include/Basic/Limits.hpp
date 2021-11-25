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
#include "Basic/Interval.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/NamingConvention.hpp"

class GSTLEARN_EXPORT Limits : public AStringable
{
public:
  Limits();
  Limits(const VectorDouble& mini,
         const VectorDouble& maxi,
         const VectorBool& incmini = VectorBool(),
         const VectorBool& incmaxi = VectorBool());
  Limits(const VectorDouble& bounds);
  Limits(int nclass);
  Limits(const Limits &m);
  Limits& operator=(const Limits &m);
  virtual ~Limits();

  virtual String toString(int level = 0) const override;

  int getLimitNumber() const { return static_cast<int>(_bounds.size()); }
  std::vector<Interval>& getBounds() { return _bounds; }
  VectorDouble getLowerBounds() const;
  VectorDouble getUpperBounds() const;
  VectorBool   getLowerIncluded() const;
  VectorBool   getUpperIncluded() const;
  bool isInside(double value) const;
  bool empty() const { return _bounds.empty(); }

  int toCategory(Db* db,
                 const String& name = String(),
                 NamingConvention namconv = NamingConvention("Category"));

  int toIndicator(Db* db,
                  const String& name = String(),
                  int OptionIndicator = 1,
                  NamingConvention namconv = NamingConvention("Indicator"));

  int toCategory(Db* db, int iatt, NamingConvention namconv);
  int toIndicator(Db* db,
                  int iatt,
                  int OptionIndicator = 1,
                  NamingConvention namconv = NamingConvention("Indicator"));

private:
  std::vector<Interval> _bounds;
};

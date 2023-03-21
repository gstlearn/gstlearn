/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/VectorT.hpp"
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
  Limits(const VectorDouble& bounds, bool addFromZero = false);
  Limits(int nclass);
  Limits(const Limits &m);
  Limits& operator=(const Limits &m);
  virtual ~Limits();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static Limits* create(const VectorDouble& mini,
                        const VectorDouble& maxi,
                        const VectorBool& incmini = VectorBool(),
                        const VectorBool& incmaxi = VectorBool());
  static Limits* create(const VectorDouble& bounds, bool addFromZero = false);
  static Limits* create(int nclass);

  int getLimitNumber() const { return static_cast<int>(_bounds.size()); }
  std::vector<Interval>& getBounds() { return _bounds; }
  VectorDouble getBound(int iclass, int mode = 0) const;
  VectorDouble getLowerBounds() const;
  VectorDouble getUpperBounds() const;
  VectorBool   getLowerIncluded() const;
  VectorBool   getUpperIncluded() const;
  bool isInside(double value) const;
  bool empty() const { return _bounds.empty(); }

  int toCategory(Db* db,
                 const String& name = String(),
                 const NamingConvention& namconv = NamingConvention("Category")) const;
  int toIndicator(Db* db,
                  const String& name = String(),
                  int OptionIndicator = 1,
                  bool flagBelow = false,
                  bool flagAbove = false,
                  const NamingConvention& namconv = NamingConvention("Indicator")) const;
  VectorDouble statistics(Db *db,
                          const String &name,
                          int optionStat = 1,
                          bool flagBelow = false,
                          bool flagAbove = false);
  int toCategoryByAttribute(Db* db, int iatt, const NamingConvention& namconv) const;
  int toIndicatorByAttribute(Db* db,
                             int iatt,
                             int OptionIndicator = 1,
                             bool flagBelow = false,
                             bool flagAbove = false,
                             const NamingConvention& namconv = NamingConvention("Indicator")) const;

private:
  std::vector<Interval> _bounds;
};

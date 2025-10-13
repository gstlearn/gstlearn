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

#include "Basic/VectorT.hpp"
#include "Basic/Interval.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/NamingConvention.hpp"

namespace gstlrn
{
class GSTLEARN_EXPORT Limits : public AStringable
{
public:
  Limits();
  Limits(const VectorDouble& mini,
         const VectorDouble& maxi,
         const VectorBool& incmini = VectorBool(),
         const VectorBool& incmaxi = VectorBool());
  Limits(const VectorDouble& bounds, bool addFromZero = false);
  Limits(Id nclass);
  Limits(const Limits &m);
  Limits& operator=(const Limits &m);
  virtual ~Limits();

  String toString(const AStringFormat* strfmt = nullptr) const override;

  static Limits* create(const VectorDouble& mini,
                        const VectorDouble& maxi,
                        const VectorBool& incmini = VectorBool(),
                        const VectorBool& incmaxi = VectorBool());
  static Limits* create(const VectorDouble& bounds, bool addFromZero = false);
  static Limits* create(Id nclass);

  Id getNLimit() const { return static_cast<Id>(_bounds.size()); }
  std::vector<Interval>& getBounds() { return _bounds; }
  VectorDouble getBound(Id iclass, Id mode = 0) const;
  VectorDouble getLowerBounds() const;
  VectorDouble getUpperBounds() const;
  VectorBool   getLowerIncluded() const;
  VectorBool   getUpperIncluded() const;
  bool isInside(double value) const;
  bool empty() const { return _bounds.empty(); }

  Id toCategory(Db* db,
                 const String& name = "",
                 const NamingConvention& namconv = NamingConvention("Category")) const;
  Id toIndicator(Db* db,
                  const String& name = "",
                  Id OptionIndicator = 1,
                  bool flagBelow = false,
                  bool flagAbove = false,
                  const NamingConvention& namconv = NamingConvention("Indicator")) const;
  VectorDouble statistics(Db *db,
                          const String &name,
                          Id optionStat = 1,
                          bool flagBelow = false,
                          bool flagAbove = false) const;
  Id toCategoryByAttribute(Db* db, Id iatt, const NamingConvention& namconv) const;
  Id toIndicatorByAttribute(Db* db,
                             Id iatt,
                             Id OptionIndicator = 1,
                             bool flagBelow = false,
                             bool flagAbove = false,
                             const NamingConvention& namconv = NamingConvention("Indicator")) const;
private:
  static Id _computeCategory(Db* db,
                              Id iatt,
                              const VectorDouble& mini,
                              const VectorDouble& maxi,
                              const VectorBool& incmini,
                              const VectorBool& incmaxi,
                              const NamingConvention& namconv);
  static Id _computeIndicator(Db* db,
                               Id iatt,
                               Id flag_indic,
                               const VectorDouble& mini,
                               const VectorDouble& maxi,
                               const VectorBool& incmini,
                               const VectorBool& incmaxi,
                               bool flagBelow,
                               bool flagAbove,
                               const NamingConvention& namconv);
  static VectorDouble _computeLimitStatistics(Db* db,
                                              Id iatt,
                                              const VectorDouble& mini,
                                              const VectorDouble& maxi,
                                              const VectorBool& incmini,
                                              const VectorBool& incmaxi,
                                              Id optionStat,
                                              bool flagBelow,
                                              bool flagAbove);
  static Id _check_bound_consistency(const VectorDouble& mini,
                                      const VectorDouble& maxi,
                                      const VectorBool& incmini,
                                      const VectorBool& incmaxi,
                                      Id* nclass_arg);

private:
  std::vector<Interval> _bounds;
};
}
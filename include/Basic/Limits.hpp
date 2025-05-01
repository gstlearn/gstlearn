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

  int getNLimit() const { return static_cast<int>(_bounds.size()); }
  std::vector<Interval>& getBounds() { return _bounds; }
  VectorDouble getBound(int iclass, int mode = 0) const;
  VectorDouble getLowerBounds() const;
  VectorDouble getUpperBounds() const;
  VectorBool   getLowerIncluded() const;
  VectorBool   getUpperIncluded() const;
  bool isInside(double value) const;
  bool empty() const { return _bounds.empty(); }

  int toCategory(Db* db,
                 const String& name = "",
                 const NamingConvention& namconv = NamingConvention("Category")) const;
  int toIndicator(Db* db,
                  const String& name = "",
                  int OptionIndicator = 1,
                  bool flagBelow = false,
                  bool flagAbove = false,
                  const NamingConvention& namconv = NamingConvention("Indicator")) const;
  VectorDouble statistics(Db *db,
                          const String &name,
                          int optionStat = 1,
                          bool flagBelow = false,
                          bool flagAbove = false) const;
  int toCategoryByAttribute(Db* db, int iatt, const NamingConvention& namconv) const;
  int toIndicatorByAttribute(Db* db,
                             int iatt,
                             int OptionIndicator = 1,
                             bool flagBelow = false,
                             bool flagAbove = false,
                             const NamingConvention& namconv = NamingConvention("Indicator")) const;
private:
  static int _computeCategory(Db* db,
                              int iatt,
                              const VectorDouble& mini,
                              const VectorDouble& maxi,
                              const VectorBool& incmini,
                              const VectorBool& incmaxi,
                              const NamingConvention& namconv);
  static int _computeIndicator(Db* db,
                               int iatt,
                               int flag_indic,
                               const VectorDouble& mini,
                               const VectorDouble& maxi,
                               const VectorBool& incmini,
                               const VectorBool& incmaxi,
                               bool flagBelow,
                               bool flagAbove,
                               const NamingConvention& namconv);
  static VectorDouble _computeLimitStatistics(Db* db,
                                              int iatt,
                                              const VectorDouble& mini,
                                              const VectorDouble& maxi,
                                              const VectorBool& incmini,
                                              const VectorBool& incmaxi,
                                              int optionStat,
                                              bool flagBelow,
                                              bool flagAbove);
  static int _check_bound_consistency(const VectorDouble& mini,
                                      const VectorDouble& maxi,
                                      const VectorBool& incmini,
                                      const VectorBool& incmaxi,
                                      int* nclass_arg);

private:
  std::vector<Interval> _bounds;
};

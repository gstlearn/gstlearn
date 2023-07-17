/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f_private.h"
#include "Db/Db.hpp"
#include "Basic/Limits.hpp"
#include "Basic/Interval.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/NamingConvention.hpp"

Limits::Limits()
    : AStringable(),
      _bounds()
{
}

Limits::Limits(const Limits &m)
    : AStringable(m),
      _bounds(m._bounds)
{

}

Limits::Limits(const VectorDouble& mini,
               const VectorDouble& maxi,
               const VectorBool& incmini,
               const VectorBool& incmaxi)
{
  if (mini.size() != maxi.size())
    my_throw("Arguments 'mini' and 'maxi' should have the same dimension");
  int nclass = static_cast<int> (mini.size());
  if (nclass <= 0)
    my_throw("You must define at least one item in 'mini' and 'maxi'");
  if (incmini.size() != 0 && (int) incmini.size() != nclass)
    my_throw("Arguments 'incmini' and 'mini' should have the same dimension");
  if (incmaxi.size() != 0 && (int) incmaxi.size() != nclass)
    my_throw("Arguments 'incmaxi' and 'maxi' should have the same dimension");

  _bounds.clear();
  for (int i = 0; i < nclass; i++)
  {
    bool incmini_loc = (incmini.empty()) ? true : incmini[i];
    bool incmaxi_loc = (incmaxi.empty()) ? false : incmaxi[i];
    Interval bd = Interval(mini[i],maxi[i],incmini_loc,incmaxi_loc);
    _bounds.push_back(bd);
  }
}

/**
 * Create the limits from a list of bounds. Intervals are delimited between two given bounds: [z_i, z_(i+1)[
 * @param bounds    List of cutoffs used to create the limits. The number of limits is equal to the number of elements in the 'bounds' vector minus one.
 * @param addFromZero When TRUE, add a class from 0 to bounds[0]
 */
Limits::Limits(const VectorDouble& bounds, bool addFromZero)
{
  if (bounds.size() == 1)
  {
    // Same as next constructor
    int nclass = static_cast<int> (bounds[0]);
    _bounds.clear();
    for (int i = 0; i < nclass; i++)
    {
      Interval bd = Interval(i + 0.5, i + 1.5);
      _bounds.push_back(bd);
    }
  }
  else
  {
    int nclass = static_cast<int> (bounds.size()) - 1;
    if (nclass <= 0)
      my_throw("The argument 'bounds' should have at least 2 items");

    _bounds.clear();

    // Add the first class from 0 to bounds[0] (optional)

    if (addFromZero && bounds[0] > 0)
    {
      Interval bd = Interval(0., bounds[0]);
      _bounds.push_back(bd);
    }

    // Store the remaining classes
    for (int i = 0; i < nclass; i++)
    {
      Interval bd;
      if (bounds[i] == bounds[i + 1])
        bd = Interval(bounds[i], bounds[i + 1], true, true);
      else
        bd = Interval(bounds[i], bounds[i + 1]);
      _bounds.push_back(bd);
    }
  }
}

Limits::Limits(int nclass)
{
  if (nclass <= 0)
    my_throw("The argument 'nclass' should be strictly positive");

  _bounds.clear();
  for (int i = 0; i < nclass; i++)
  {
    Interval bd = Interval(i + 0.5, i + 1.5);
    _bounds.push_back(bd);
  }
}

Limits& Limits::operator=(const Limits &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _bounds = m._bounds;
  }
  return *this;
}

Limits::~Limits()
{

}

Limits* Limits::create(const VectorDouble& mini,
                       const VectorDouble& maxi,
                       const VectorBool& incmini,
                       const VectorBool& incmaxi)
{
  Limits* limits = new Limits(mini, maxi, incmini, incmaxi);
  return limits;
}

Limits* Limits::create(const VectorDouble& bounds, bool addFromZero)
{
  Limits* limits = new Limits(bounds, addFromZero);
  return limits;
}

Limits* Limits::create(int nclass)
{
  Limits* limits = new Limits(nclass);
  return limits;
}

String Limits::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  for (int i = 0; i < (int) _bounds.size(); i++)
    sstr << "Bound( " << i+1 << " ) : " << _bounds[i].toString() << std::endl;

  return sstr.str();
}

/**
 * Retrieve the set of bounds or one bound
 * @param iclass Rank of the class
 * @param mode   0 for both bounds; 1 for lower bound; 2 for upper bound
 * @return The vector of bound values
 */
VectorDouble Limits::getBound(int iclass, int mode) const
{
  VectorDouble bounds;
  if (iclass < 0 || iclass >= getLimitNumber()) return bounds;

  if (mode == 0 || mode == 1)
    bounds.push_back(_bounds[iclass].getVmin());
  if (mode == 0 || mode == 2)
    bounds.push_back(_bounds[iclass].getVmax());
  return bounds;
}

VectorDouble Limits::getLowerBounds() const
{
  int nclass = getLimitNumber();
  VectorDouble lower(nclass);
  for (int i = 0; i < nclass; i++)
    lower[i] = _bounds[i].getVmin();
  return lower;
}

VectorDouble Limits::getUpperBounds() const
{
  int nclass = getLimitNumber();
  VectorDouble upper(nclass);
  for (int i = 0; i < nclass; i++)
    upper[i] = _bounds[i].getVmax();
  return upper;
}

VectorBool Limits::getLowerIncluded() const
{
  int nclass = getLimitNumber();
  VectorBool mininc(nclass);
  for (int i = 0; i < nclass; i++)
    mininc[i] = _bounds[i].getMinIncluded();
  return mininc;
}

VectorBool Limits::getUpperIncluded() const
{
  int nclass = getLimitNumber();
  VectorBool maxinc(nclass);
  for (int i = 0; i < nclass; i++)
    maxinc[i] = _bounds[i].getMaxIncluded();
  return maxinc;
}

bool Limits::isInside(double value) const
{
  for (int i = 0; i < getLimitNumber(); i++)
  {
    if (! _bounds[i].isInside(value)) return false;
  }
  return true;
}

int Limits::toCategoryByAttribute(Db* db,
                                  int iatt,
                                  const NamingConvention& namconv) const
{
  return _db_category(db, iatt, getLowerBounds(), getUpperBounds(),
                      getLowerIncluded(), getUpperIncluded(), namconv);
}

int Limits::toCategory(Db* db,
                       const String& name,
                       const NamingConvention& namconv) const
{
  int iatt = db->getUID(name);
  if (iatt < 0) return 1;
  return toCategoryByAttribute(db, iatt, namconv);
}

/**
 * Create indicators variables on the intervals defined by the limits for a given variable in a Db.  
 * Note:
 *
 * - If OptionIndicator is 1, the Db-class will contain the new indicator variables.
 * There are as many new variables as they are classes.
 * Each sample of the indicator variable for class 'iclass' is set to 1 if the sample belongs to this class
 * or 0 otherwise.
 *
 * - If OptionIndicator is 0, the Db-class will contain one variable such that each sample contains
 * the average of the variable calculated over the samples whose value belong to this class.

 * @param db                 Db containing the variable to be discretized (from which the indicators are computed)
 * @param name               Name of the variable in the Db to be discretized.
 * @param OptionIndicator    When 1, the function assignes the indicator variables.
 *							             When 0, the function assignes the average of the class.
 * @param flagBelow          When True, consider samples below lowest bound
 * @param flagAbove          When True, consider samples above highest bound
 * @param namconv            Naming convention
 *
 * @return
 */
int Limits::toIndicator(Db* db,
                        const String& name,
                        int OptionIndicator,
                        bool flagBelow,
                        bool flagAbove,
                        const NamingConvention& namconv) const
{
  int iatt = db->getUID(name);
  if (iatt < 0) return 1;
  return toIndicatorByAttribute(db, iatt, OptionIndicator, flagBelow, flagAbove, namconv);
}

int Limits::toIndicatorByAttribute(Db *db,
                                   int iatt,
                                   int OptionIndicator,
                                   bool flagBelow,
                                   bool flagAbove,
                                   const NamingConvention &namconv) const
{
  return _db_indicator(db, iatt, OptionIndicator,
                       getLowerBounds(), getUpperBounds(),
                       getLowerIncluded(), getUpperIncluded(),
                       flagBelow, flagAbove, namconv);
}

/**
 * Calculate the statistics per Class
 * @param db         Target Db
 * @param name       Name of the Target Variable
 * @param optionStat 1 for Mean; 2 for Proportions
 * @param flagBelow  When TRUE, add a class for samples below lowest bound
 * @param flagAbove  When TRUE, add a class for samples above highest bound
 * @return
 */
VectorDouble Limits::statistics(Db* db, const String& name,
                                int optionStat, bool flagBelow, bool flagAbove)
{
  int iatt = db->getUID(name);
  if (iatt < 0) return 1;
  return _db_limits_statistics(db, iatt,
                               getLowerBounds(), getUpperBounds(),
                               getLowerIncluded(), getUpperIncluded(),
                               optionStat, flagBelow, flagAbove);
}

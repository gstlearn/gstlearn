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
#include "Basic/Limits.hpp"
#include "Basic/Interval.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/Utilities.hpp"
#include "Db/Db.hpp"
#include "geoslib_define.h"

namespace gstlrn
{
Limits::Limits()
  : AStringable()
  , _bounds()
{
}

Limits::Limits(const Limits& m)
  : AStringable(m)
  , _bounds(m._bounds)
{
}

Limits::Limits(const VectorDouble& mini,
               const VectorDouble& maxi,
               const VectorBool& incmini,
               const VectorBool& incmaxi)
{
  _bounds.clear();
  if (mini.size() != maxi.size())
  {
    messerr("Arguments 'mini' and 'maxi' should have the same dimension. Limits empty");
    return;
  }
  Id nclass = static_cast<Id>(mini.size());
  if (nclass <= 0)
  {
    messerr("You must define at least one item in 'mini' and 'maxi'. Limits empty");
    return;
  }
  if (incmini.size() != 0 && static_cast<Id>(incmini.size()) != nclass)
  {
    messerr("Arguments 'incmini' and 'mini' should have the same dimension. Limits empty");
    return;
  }
  if (incmaxi.size() != 0 && static_cast<Id>(incmaxi.size()) != nclass)
  {
    messerr("Arguments 'incmaxi' and 'maxi' should have the same dimension. Limits empty");
    return;
  }

  for (Id i = 0; i < nclass; i++)
  {
    bool incmini_loc = (incmini.empty()) ? true : incmini[i];
    bool incmaxi_loc = (incmaxi.empty()) ? false : incmaxi[i];
    Interval bd(mini[i], maxi[i], incmini_loc, incmaxi_loc);
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
    Id nclass = static_cast<Id>(bounds[0]);
    _bounds.clear();
    for (Id i = 0; i < nclass; i++)
    {
      Interval bd(i + 0.5, i + 1.5);
      _bounds.push_back(bd);
    }
  }
  else
  {
    Id nclass = static_cast<Id>(bounds.size()) - 1;
    if (nclass <= 0)
      my_throw("The argument 'bounds' should have at least 2 items");

    _bounds.clear();

    // Add the first class from 0 to bounds[0] (optional)

    if (addFromZero && bounds[0] > 0)
    {
      Interval bd(0., bounds[0]);
      _bounds.push_back(bd);
    }

    // Store the remaining classes
    for (Id i = 0; i < nclass; i++)
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

Limits::Limits(Id nclass)
{
  if (nclass <= 0)
    my_throw("The argument 'nclass' should be strictly positive");

  _bounds.clear();
  for (Id i = 0; i < nclass; i++)
  {
    Interval bd(i + 0.5, i + 1.5);
    _bounds.push_back(bd);
  }
}

Limits& Limits::operator=(const Limits& m)
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
  auto* limits = new Limits(mini, maxi, incmini, incmaxi);
  return limits;
}

Limits* Limits::create(const VectorDouble& bounds, bool addFromZero)
{
  auto* limits = new Limits(bounds, addFromZero);
  return limits;
}

Limits* Limits::create(Id nclass)
{
  auto* limits = new Limits(nclass);
  return limits;
}

String Limits::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  for (Id i = 0; i < static_cast<Id>(_bounds.size()); i++)
    sstr << "Bound( " << i + 1 << " ) : " << _bounds[i].toString() << std::endl;

  return sstr.str();
}

/**
 * Retrieve the set of bounds or one bound
 * @param iclass Rank of the class
 * @param mode   0 for both bounds; 1 for lower bound; 2 for upper bound
 * @return The vector of bound values
 */
VectorDouble Limits::getBound(Id iclass, Id mode) const
{
  VectorDouble bounds;
  if (iclass < 0 || iclass >= getNLimit()) return bounds;

  if (mode == 0 || mode == 1)
    bounds.push_back(_bounds[iclass].getVmin());
  if (mode == 0 || mode == 2)
    bounds.push_back(_bounds[iclass].getVmax());
  return bounds;
}

VectorDouble Limits::getLowerBounds() const
{
  auto nclass = getNLimit();
  VectorDouble lower(nclass);
  for (Id i = 0; i < nclass; i++)
    lower[i] = _bounds[i].getVmin();
  return lower;
}

VectorDouble Limits::getUpperBounds() const
{
  auto nclass = getNLimit();
  VectorDouble upper(nclass);
  for (Id i = 0; i < nclass; i++)
    upper[i] = _bounds[i].getVmax();
  return upper;
}

VectorBool Limits::getLowerIncluded() const
{
  auto nclass = getNLimit();
  VectorBool mininc(nclass);
  for (Id i = 0; i < nclass; i++)
    mininc[i] = _bounds[i].getMinIncluded();
  return mininc;
}

VectorBool Limits::getUpperIncluded() const
{
  auto nclass = getNLimit();
  VectorBool maxinc(nclass);
  for (Id i = 0; i < nclass; i++)
    maxinc[i] = _bounds[i].getMaxIncluded();
  return maxinc;
}

bool Limits::isInside(double value) const
{
  for (Id i = 0; i < getNLimit(); i++)
  {
    if (!_bounds[i].isInside(value)) return false;
  }
  return true;
}

Id Limits::toCategoryByAttribute(Db* db,
                                 Id iatt,
                                 const NamingConvention& namconv) const
{
  return _computeCategory(db, iatt, getLowerBounds(), getUpperBounds(),
                          getLowerIncluded(), getUpperIncluded(), namconv);
}

Id Limits::toCategory(Db* db,
                      const String& name,
                      const NamingConvention& namconv) const
{
  Id iatt = db->getUID(name);
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
Id Limits::toIndicator(Db* db,
                       const String& name,
                       Id OptionIndicator,
                       bool flagBelow,
                       bool flagAbove,
                       const NamingConvention& namconv) const
{
  Id iatt = db->getUID(name);
  if (iatt < 0) return 1;
  return toIndicatorByAttribute(db, iatt, OptionIndicator, flagBelow, flagAbove, namconv);
}

Id Limits::toIndicatorByAttribute(Db* db,
                                  Id iatt,
                                  Id OptionIndicator,
                                  bool flagBelow,
                                  bool flagAbove,
                                  const NamingConvention& namconv) const
{
  return _computeIndicator(db, iatt, OptionIndicator, getLowerBounds(),
                           getUpperBounds(), getLowerIncluded(),
                           getUpperIncluded(), flagBelow, flagAbove, namconv);
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
VectorDouble Limits::statistics(Db* db,
                                const String& name,
                                Id optionStat,
                                bool flagBelow,
                                bool flagAbove) const
{
  Id iatt = db->getUID(name);
  if (iatt < 0) return 1;
  return _computeLimitStatistics(db, iatt, getLowerBounds(), getUpperBounds(),
                                 getLowerIncluded(), getUpperIncluded(),
                                 optionStat, flagBelow, flagAbove);
}

/****************************************************************************/
/*!
 **  Convert the contents of the ivar-th continuous variable
 **  into a categorical array
 **
 ** \return  Error return code
 **
 ** \param[in]  db      Db structure
 ** \param[in]  iatt    Rank of the attribute
 ** \param[in]  mini    Array containing the minima per class
 ** \param[in]  maxi    Array containing the maxima per class
 ** \param[in]  incmini Array containing the inclusive flag for minima per class
 ** \param[in]  incmaxi Array containing the inclusive flag for maxima per class
 ** \param[in]  namconv Naming convention
 **
 ** \remark When the array mini or maxi are not provided, then
 ** \remark mini[iclass] = Id(iclass-1)
 ** \remark maxi[iclass[ = Id(iclass)
 **
 *****************************************************************************/
Id Limits::_computeCategory(Db* db,
                            Id iatt,
                            const VectorDouble& mini,
                            const VectorDouble& maxi,
                            const VectorBool& incmini,
                            const VectorBool& incmaxi,
                            const NamingConvention& namconv)
{
  // Determination of the number of classes

  Id nclass;
  if (_check_bound_consistency(mini, maxi, incmini, incmaxi, &nclass))
    return 1;

  /* Create the variable */

  Id iptr = db->addColumnsByConstant(1, TEST);
  if (iptr < 0) return (1);

  /* Loop on the samples */

  for (Id iech = 0; iech < db->getNSample(); iech++)
  {
    if (!db->isActive(iech)) continue;
    double value = db->getArray(iech, iatt);
    if (FFFF(value)) continue;

    /* Loop on the limits classes */

    Id ival = 0;
    for (Id iclass = 0; iclass < nclass; iclass++)
    {
      double minival = (mini.empty()) ? iclass + 0.5 : mini[iclass];
      double maxival = (maxi.empty()) ? iclass + 1.5 : maxi[iclass];
      if (!FFFF(minival))
      {
        Id flag = (incmini.empty()) ? 1 : static_cast<Id>(incmini[iclass]);
        if ((flag == 0 && value <= minival) || (flag == 1 && value < minival))
          continue;
      }
      if (!FFFF(maxival))
      {
        Id flag = (incmaxi.empty()) ? 0 : static_cast<Id>(incmaxi[iclass]);
        if ((flag == 0 && value >= maxival) || (flag == 1 && value > maxival))
          continue;
      }
      ival = iclass + 1;
    }

    /* Set the returning value */

    db->setArray(iech, iptr, static_cast<double>(ival));
  }

  namconv.setNamesAndLocators(db, iatt, db, iptr);

  return (0);
}

/****************************************************************************/
/*!
 **  Create indicator variables
 **
 ** \return  Error returned code
 **
 ** \param[in]  db      Db structure
 ** \param[in]  iatt    Rank of the target variable
 ** \param[in]  flag_indic Type of variable(s) to be stored:
 **                     1 the indicator variable
 **                     0 the mean variable per class
 ** \param[in]  mini    Array containing the minima per class
 ** \param[in]  maxi    Array containing the maxima per class
 ** \param[in]  incmini Array containing the inclusive flag for minima per class
 ** \param[in]  incmaxi Array containing the inclusive flag for maxima per class
 ** \param[in]  flagBelow If True, consider the values below the lower bound
 ** \param[in]  flagAbove If True, consider the values above the upper bound
 ** \param[in]  namconv Naming convention
 **
 ** \remark When both arrays mini and maxi are not provided, then:
 ** \remark mini[iclass] = iclass-1
 ** \remark maxi[iclass[ = iclass
 **
 *****************************************************************************/
Id Limits::_computeIndicator(Db* db,
                             Id iatt,
                             Id flag_indic,
                             const VectorDouble& mini,
                             const VectorDouble& maxi,
                             const VectorBool& incmini,
                             const VectorBool& incmaxi,
                             bool flagBelow,
                             bool flagAbove,
                             const NamingConvention& namconv)
{
  // Determination of the number of classes

  Id nclass;
  if (_check_bound_consistency(mini, maxi, incmini, incmaxi, &nclass))
    return 1;

  /* Core allocation */

  VectorInt count(nclass);
  VectorInt flagmin(nclass);
  VectorInt flagmax(nclass);
  VectorDouble mean(nclass);
  VectorDouble minival(nclass);
  VectorDouble maxival(nclass);
  for (Id iclass = 0; iclass < nclass; iclass++)
  {
    count[iclass]   = 0;
    mean[iclass]    = 0.;
    minival[iclass] = (mini.empty()) ? iclass + 0.5 : mini[iclass];
    maxival[iclass] = (maxi.empty()) ? iclass + 1.5 : maxi[iclass];
    flagmin[iclass] = (incmini.empty()) ? 1 : static_cast<Id>(incmini[iclass]);
    flagmax[iclass] = (incmaxi.empty()) ? 0 : static_cast<Id>(incmaxi[iclass]);
  }
  Id nbelow     = 0;
  Id nabove     = 0;
  double mbelow = 0.;
  double mabove = 0.;

  /* Find extrema of all classes to sort too small or too large samples */

  double zmini = MAXIMUM_BIG;
  double zmaxi = MINIMUM_BIG;
  for (Id iclass = 0; iclass < nclass; iclass++)
  {
    if (!FFFF(minival[iclass])) zmini = MIN(zmini, minival[iclass]);
    if (!FFFF(maxival[iclass])) zmaxi = MAX(zmaxi, maxival[iclass]);
  }
  if (FFFF(zmini)) zmini = MINIMUM_BIG;
  if (FFFF(zmaxi)) zmaxi = MAXIMUM_BIG;

  /* Create the variables */

  Id iptr_indic = -1;
  Id iptr_below = -1;
  Id iptr_above = -1;
  Id iptr_mean  = -1;
  if (flag_indic)
  {
    if (flagBelow)
    {
      iptr_below = db->addColumnsByConstant(1, 0.);
      if (iptr_below < 0) return 1;
    }
    iptr_indic = db->addColumnsByConstant(nclass, 0.);
    if (iptr_indic < 0) return 1;
    if (flagAbove)
    {
      iptr_above = db->addColumnsByConstant(1, 0.);
      if (iptr_above < 0) return 1;
    }
  }
  else
  {
    iptr_mean = db->addColumnsByConstant(1, TEST);
    if (iptr_mean < 0) return 1;
  }

  /* Loop on the samples */

  for (Id iech = 0; iech < db->getNSample(); iech++)
  {
    if (!db->isActive(iech)) continue;
    double value = db->getArray(iech, iatt);
    if (FFFF(value)) continue;

    /* Loop on the limit classes */

    Id found = -1;
    for (Id iclass = 0; iclass < nclass; iclass++)
    {
      Id belong = 1;
      if (!FFFF(minival[iclass]))
      {
        if ((flagmin[iclass] == 0 && value <= minival[iclass]) ||
            (flagmin[iclass] == 1 && value < minival[iclass]))
          belong = 0;
      }
      if (!FFFF(maxival[iclass]))
      {
        if ((flagmax[iclass] == 0 && value >= maxival[iclass]) ||
            (flagmax[iclass] == 1 && value > maxival[iclass]))
          belong = 0;
      }

      /* Store the indicator (if required) */

      if (flag_indic)
        db->setArray(iech, iptr_indic + iclass, belong);

      if (belong)
      {
        mean[iclass] += value;
        count[iclass]++;
        found = iclass;
      }
    }

    // Store the Below of Above indicator

    if (flag_indic)
    {
      if (flagBelow)
        db->setArray(iech, iptr_below, value < zmini);

      if (flagAbove)
        db->setArray(iech, iptr_above, value > zmaxi);
    }

    /* Update the statistics for classes below and above */

    Id iclass;
    if (found < 0)
    {
      if (value < zmini)
      {
        nbelow += 1;
        iclass = -1;
        mbelow += value;
      }
      else
      {
        nabove += 1;
        iclass = nclass;
        mabove += value;
      }
    }
    else
    {
      iclass = found;
    }
    if (!flag_indic)
      db->setArray(iech, iptr_mean, static_cast<double>(iclass));
  }

  /* Calculate the statistics */

  for (Id iclass = 0; iclass < nclass; iclass++)
  {
    if (count[iclass] <= 0)
      mean[iclass] = TEST;
    else
      mean[iclass] /= static_cast<double>(count[iclass]);
  }
  if (nbelow > 0) mbelow /= static_cast<double>(nbelow);
  if (nabove > 0) mabove /= static_cast<double>(nabove);
  if (!flagBelow) mbelow = TEST;
  if (!flagAbove) mabove = TEST;

  /* Assign the mean variable per class */

  if (!flag_indic)
  {
    for (Id iech = 0; iech < db->getNSample(); iech++)
    {
      if (!db->isActive(iech)) continue;
      double value = db->getArray(iech, iptr_mean);
      if (FFFF(value)) continue;
      Id iclass = static_cast<Id>(value);
      if (iclass < 0)
        db->setArray(iech, iptr_mean, mbelow);
      else if (iclass < nclass)
        db->setArray(iech, iptr_mean, mean[iclass]);
      else
        db->setArray(iech, iptr_mean, mabove);
    }
  }

  // Naming convention
  if (flag_indic == 1)
  {
    if (flagBelow)
      namconv.setNamesAndLocators(db, iatt, db, iptr_below, "Below", 1);
    namconv.setNamesAndLocators(db, iatt, db, iptr_indic, "Class", nclass);
    if (flagAbove)
      namconv.setNamesAndLocators(db, iatt, db, iptr_above, "Above", 1);
  }
  else
    namconv.setNamesAndLocators(db, iatt, db, iptr_mean, "Mean", 1);

  return (0);
}

/****************************************************************************/
/*!
 **  Calculate statistics per class
 **
 ** \return  Array of statistics
 **
 ** \param[in]  db      Db structure
 ** \param[in]  iatt    Rank of the target variable
 ** \param[in]  mini    Array containing the minima per class
 ** \param[in]  maxi    Array containing the maxima per class
 ** \param[in]  incmini Array containing the inclusive flag for minima per class
 ** \param[in]  incmaxi Array containing the inclusive flag for maxima per class
 ** \param[in]  optionStat 1 for Proportions; 2 for Mean
 ** \param[in]  flagBelow If True, consider the values below the lower bound
 ** \param[in]  flagAbove If True, consider the values above the upper bound
 **
 *****************************************************************************/
VectorDouble Limits::_computeLimitStatistics(Db* db,
                                             Id iatt,
                                             const VectorDouble& mini,
                                             const VectorDouble& maxi,
                                             const VectorBool& incmini,
                                             const VectorBool& incmaxi,
                                             Id optionStat,
                                             bool flagBelow,
                                             bool flagAbove)
{
  // Determination of the number of classes

  Id nclass;
  if (_check_bound_consistency(mini, maxi, incmini, incmaxi, &nclass))
    return 1;

  /* Core allocation */

  VectorInt count(nclass);
  VectorInt flagmin(nclass);
  VectorInt flagmax(nclass);
  VectorDouble mean(nclass);
  VectorDouble minival(nclass);
  VectorDouble maxival(nclass);
  for (Id iclass = 0; iclass < nclass; iclass++)
  {
    count[iclass]   = 0;
    mean[iclass]    = 0.;
    minival[iclass] = (mini.empty()) ? iclass + 0.5 : mini[iclass];
    maxival[iclass] = (maxi.empty()) ? iclass + 1.5 : maxi[iclass];
    flagmin[iclass] = (incmini.empty()) ? 1 : static_cast<Id>(incmini[iclass]);
    flagmax[iclass] = (incmaxi.empty()) ? 0 : static_cast<Id>(incmaxi[iclass]);
  }
  Id nbelow     = 0;
  Id nabove     = 0;
  double mbelow = 0.;
  double mabove = 0.;

  /* Find extrema of all classes to sort too small or too large samples */

  double zmini = MAXIMUM_BIG;
  double zmaxi = MINIMUM_BIG;
  for (Id iclass = 0; iclass < nclass; iclass++)
  {
    if (!FFFF(minival[iclass])) zmini = MIN(zmini, minival[iclass]);
    if (!FFFF(maxival[iclass])) zmaxi = MAX(zmaxi, maxival[iclass]);
  }
  if (FFFF(zmini)) zmini = MINIMUM_BIG;
  if (FFFF(zmaxi)) zmaxi = MAXIMUM_BIG;

  /* Loop on the samples */

  Id nactive = 0;
  for (Id iech = 0; iech < db->getNSample(); iech++)
  {
    if (!db->isActive(iech)) continue;
    double value = db->getArray(iech, iatt);
    if (FFFF(value)) continue;

    /* Loop on the limit classes */

    Id found = -1;
    for (Id iclass = 0; iclass < nclass; iclass++)
    {
      Id belong = 1;
      if (!FFFF(minival[iclass]))
      {
        if ((flagmin[iclass] == 0 && value <= minival[iclass]) ||
            (flagmin[iclass] == 1 && value < minival[iclass]))
          belong = 0;
      }
      if (!FFFF(maxival[iclass]))
      {
        if ((flagmax[iclass] == 0 && value >= maxival[iclass]) ||
            (flagmax[iclass] == 1 && value > maxival[iclass]))
          belong = 0;
      }

      if (belong)
      {
        mean[iclass] += value;
        count[iclass]++;
        found = iclass;
      }
    }

    /* Update the statistics for classes below and above */

    if (found < 0)
    {
      if (value < zmini)
      {
        nbelow += 1;
        mbelow += value;
      }
      else
      {
        nabove += 1;
        mabove += value;
      }
    }
    nactive++;
  }

  /* Calculate the statistics */

  for (Id iclass = 0; iclass < nclass; iclass++)
  {
    if (count[iclass] <= 0)
      mean[iclass] = TEST;
    else
      mean[iclass] /= static_cast<double>(count[iclass]);
  }
  if (nbelow > 0) mbelow /= static_cast<double>(nbelow);
  if (nabove > 0) mabove /= static_cast<double>(nabove);

  /* Returning the results */

  VectorDouble stats;

  if (optionStat == 1)
  {
    if (flagBelow) stats.push_back(static_cast<double>(nbelow) / static_cast<double>(nactive));
    for (Id iclass = 0; iclass < nclass; iclass++)
      stats.push_back(static_cast<double>(count[iclass]) / static_cast<double>(nactive));
    if (flagAbove) stats.push_back(static_cast<double>(nabove) / static_cast<double>(nactive));
  }
  else
  {
    if (flagBelow) stats.push_back(mbelow);
    for (Id iclass = 0; iclass < nclass; iclass++)
      stats.push_back(mean[iclass]);
    if (flagAbove) stats.push_back(mabove);
  }
  return stats;
}

/****************************************************************************/
/*!
 **  Check consistency between the different bounds vectors
 **  and returns the number of classes
 **
 ** \return  Error return code
 **
 ** \param[in]  mini    Array containing the minima per class
 ** \param[in]  maxi    Array containing the maxima per class
 ** \param[in]  incmini Array containing the inclusive flag for minima per class
 ** \param[in]  incmaxi Array containing the inclusive flag for maxima per class
 **
 ** \param[out] nclass_arg Number of classes
 **
 *****************************************************************************/
Id Limits::_check_bound_consistency(const VectorDouble& mini,
                                    const VectorDouble& maxi,
                                    const VectorBool& incmini,
                                    const VectorBool& incmaxi,
                                    Id* nclass_arg)
{
  Id nclass = 0;
  if (!mini.empty())
  {
    if (nclass > 0 && nclass != static_cast<Id>(mini.size()))
    {
      messerr("Wrong dimension of 'mini'(%d). It should be %d",
              static_cast<Id>(mini.size()), nclass);
      return 1;
    }
    nclass = static_cast<Id>(mini.size());
  }
  if (!maxi.empty())
  {
    if (nclass > 0 && nclass != static_cast<Id>(maxi.size()))
    {
      messerr("Wrong dimension of 'maxi'(%d). It should be %d",
              static_cast<Id>(maxi.size()), nclass);
      return 1;
    }
    nclass = static_cast<Id>(maxi.size());
  }
  if (!incmini.empty())
  {
    if (nclass > 0 && nclass != static_cast<Id>(incmini.size()))
    {
      messerr("Wrong dimension of 'incmini'(%d). It should be %d",
              incmini.size(), nclass);
      return 1;
    }
    nclass = static_cast<Id>(incmini.size());
  }
  if (!incmaxi.empty())
  {
    if (nclass > 0 && nclass != static_cast<Id>(incmaxi.size()))
    {
      messerr("Wrong dimension of 'incmaxi'(%d). It should be %d",
              incmaxi.size(), nclass);
      return 1;
    }
    nclass = static_cast<Id>(incmaxi.size());
  }
  if (nclass <= 0)
  {
    messerr("You must define at least one valid limit");
    return 1;
  }
  *nclass_arg = nclass;
  return 0;
}

} // namespace gstlrn
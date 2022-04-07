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
#include "geoslib_f.h"
#include "geoslib_f_private.h"
#include "Db/Db.hpp"
#include "Basic/Limits.hpp"
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
    throw("Arguments 'mini' and 'maxi' should have the same dimension");
  int nclass = static_cast<int> (mini.size());
  if (nclass <= 0)
    throw("You must define at least one item in 'mini' and 'maxi'");
  if (incmini.size() != 0 && (int) incmini.size() != nclass)
    throw("Arguments 'incmini' and 'mini' should have the same dimension");
  if (incmaxi.size() != 0 && (int) incmaxi.size() != nclass)
    throw("Arguments 'incmaxi' and 'maxi' should have the same dimension");

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
 * Create the limits from a list of bounds. Intervals are demilited between two given bounds: [z_i, z_(i+1)[
 * @param bounds    list of cutoffs used to create the limits. The number of limits is equal to the number of elements in the 'bounds' vector minus one. 
 * @return
 */
Limits::Limits(const VectorDouble& bounds)
{
  int nclass = static_cast<int> (bounds.size()) - 1;
  if (nclass <= 0)
    throw("The argument 'bounds' should have at least 2 items");

  _bounds.clear();
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

Limits::Limits(int nclass)
{
  if (nclass <= 0)
    throw("The argument 'nclass' should be strictly positive");

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

String Limits::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  for (int i = 0; i < (int) _bounds.size(); i++)
    sstr << "Bound( " << i+1 << " ) : " << _bounds[i].toString() << std::endl;

  return sstr.str();
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

int Limits::toCategoryByAttribute(Db* db, int iatt, const NamingConvention& namconv)
{
  return _db_category(db, iatt, getLowerBounds(), getUpperBounds(),
                      getLowerIncluded(), getUpperIncluded(), namconv);
}

int Limits::toCategory(Db* db, const String& name,
                       const NamingConvention& namconv)
{
  int iatt = db->getUID(name);
  if (iatt < 0) return 1;
  return toCategoryByAttribute(db, iatt, namconv);
}

/**
 * Create indicators variables on the intervals defined by the limits for a given variable in a Db.  
 * Note:  
 * If OptionIndicator is 1, the Db-class will contain the new indicator variables. There are as many new variables as they are classes. Each sample of the indicator variable for class 'iclass' is
 * set to 1 if the sample belongs to this class or 0 otherwise.  
 * If OptionIndicator is 0, the Db-class will contain a variable such that each sample contains the variable average calculated over the samples whose value belong to this class.

 * @param db                 Db containing the variable to be discretized (from which the indicators are computed)
 * @param name               Name of the variable in the Db to be discretized.
 * @param OptionIndicator    When 1, the function calculates the indicator variables.
 *							 When 0, it calculates the discretized variable: each sample contains the average value of the variable within the class to which it belongs.
 * @return
 */
int Limits::toIndicator(Db* db,
                        const String& name,
                        int OptionIndicator,
                        const NamingConvention& namconv)
{
  int iatt = db->getUID(name);
  if (iatt < 0) return 1;
  return toIndicatorByAttribute(db, iatt, OptionIndicator, namconv);
}

int Limits::toIndicatorByAttribute(Db* db,
                        int iatt,
                        int OptionIndicator,
                        const NamingConvention& namconv)
{
  return _db_indicator(db, iatt, OptionIndicator, getLowerBounds(),
                       getUpperBounds(), getLowerIncluded(),
                       getUpperIncluded(), namconv);
}

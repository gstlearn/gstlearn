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
#include "Db/Db.hpp"
#include "Basic/Limits.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/NamingConvention.hpp"
#include "geoslib_f.h"
#include "geoslib_f_private.h"

Limits::Limits()
    : _bounds()
{
}

Limits::Limits(const Limits &m)
    : _bounds(m._bounds)
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
    _bounds = m._bounds;
  }
  return *this;
}

Limits::~Limits()
{

}

String Limits::toString(int level) const
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

int Limits::toCategory(Db* db, int iatt, NamingConvention namconv)
{
  return _db_category(db, iatt, getLowerBounds(), getUpperBounds(),
                     getLowerIncluded(), getUpperIncluded(), namconv);
}

int Limits::toCategory(Db* db, const String& name, NamingConvention namconv)
{
  VectorInt iatts = db->ids(name, true);
  if (iatts.empty()) return 1;
  return toCategory(db, iatts[0], namconv);
}

int Limits::toIndicator(Db* db,
                        const String& name,
                        int OptionIndicator,
                        NamingConvention namconv)
{
  VectorInt iatts = db->ids(name, true);
  if (iatts.empty()) return 1;
  return toIndicator(db, iatts[0], OptionIndicator, namconv);
}

int Limits::toIndicator(Db* db,
                        int iatt,
                        int OptionIndicator,
                        NamingConvention namconv)
{
  return _db_indicator(db, iatt, OptionIndicator, getLowerBounds(),
                      getUpperBounds(), getLowerIncluded(),
                      getUpperIncluded(), namconv);
}

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
#include "Fractures/Environ.hpp"
#include "Basic/AStringable.hpp"

Environ::Environ(double xmax,
                 double ymax,
                 double deltax,
                 double deltay,
                 double xextend,
                 double mean,
                 double stdev)
  : AStringable(),
    _xmax(xmax),
    _ymax(ymax),
    _deltax(deltax),
    _deltay(deltay),
    _xextend(xextend),
    _mean(mean),
    _stdev(stdev),
    _families(),
    _faults()
{
}

Environ::Environ(const Environ& r)
    : AStringable(r),
      _xmax(r._xmax),
      _ymax(r._ymax),
      _deltax(r._deltax),
      _deltay(r._deltay),
      _xextend(r._xextend),
      _mean(r._mean),
      _stdev(r._stdev),
      _families(r._families),
      _faults(r._faults)
{
}

Environ& Environ::operator=(const Environ& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _xmax = r._xmax;
    _ymax = r._ymax;
    _deltax = r._deltax;
    _deltay = r._deltay;
    _xextend = r._xextend;
    _mean = r._mean;
    _stdev = r._stdev;
    _families = r._families;
    _faults = r._faults;
  }
  return *this;
}

Environ::~Environ()
{
}

Environ* Environ::create(double xmax,
                         double ymax,
                         double deltax,
                         double deltay,
                         double xextend,
                         double mean,
                         double stdev)
{
  return new Environ(xmax, ymax, deltax, deltay, xextend, mean, stdev);
}

String Environ::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

   /* General characteristics */

   sstr << toTitle(0, "Geometry");
   sstr << "Field extension (horizontal)    = " << _xmax << std::endl;
   sstr << "Field extension (vertical)      = " << _ymax << std::endl;
   sstr << "Field dilation (horizontal)     = " << _deltax << std::endl;
   sstr << "Field dilation (vertical)       = " << _deltay << std::endl;
   sstr << "Mean of thickness law           = " << _mean << std::endl;
   sstr << "St. dev. of thickness law       = " << _stdev << std::endl;
   sstr << "Number of families              = " << getNFamilies() << std::endl;
   sstr << "Number of faults                = " << getNFaults() << std::endl;

   /* Loop on the families */

  for (int j = 0; j < getNFamilies(); j++)
  {
    sstr << toTitle(2, "Family #%d/%d", j + 1, getNFamilies());
    sstr << _families[j].toString(strfmt);
  }

  /* Loop on the faults */

  for (int i = 0; i < getNFaults(); i++)
  {
    mestitle(1, "Fault #%d/%d", i + 1, getNFaults());
    sstr << _faults[i].toString(strfmt);
  }

  return sstr.str();
}

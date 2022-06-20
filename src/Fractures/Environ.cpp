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
#include "Basic/ASerializable.hpp"

Environ::Environ(double xmax,
                 double ymax,
                 double deltax,
                 double deltay,
                 double xextend,
                 double mean,
                 double stdev)
  : AStringable(),
    ASerializable(),
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
      ASerializable(r),
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
    ASerializable::operator=(r);
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

int Environ::dumpToNF(const String& neutralFilename, bool verbose) const
{
  std::ofstream os;
  int ret = 1;
  if (_fileOpenWrite(neutralFilename, "Fracture_Environ", os, verbose))
  {
    ret = _serialize(os, verbose);
    if (ret && verbose) messerr("Problem writing in the Neutral File.");
    os.close();
  }
  return ret;
}

/**
 * Create a Environ by loading the contents of a Neutral File
 *
 * @param neutralFilename Name of the Neutral File
 * @param verbose         Verbose
 */
Environ* Environ::createFromNF(const String& neutralFilename, bool verbose)
{
  Environ* environ = nullptr;
  std::ifstream is;
  if (_fileOpenRead(neutralFilename, "Fracture_Environ", is, verbose))
  {
    environ = new Environ;
    if (environ->_deserialize(is, verbose))
    {
      if (verbose) messerr("Problem reading the Neutral File.");
      delete environ;
      environ = nullptr;
    }
    is.close();
  }
  return environ;
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

int Environ::_deserialize(std::istream& is, bool verbose)
{
  int nfamilies, nfaults;
  bool ret = true;

  ret = ret && _recordRead<int>(is, "Number of families", nfamilies);
  ret = ret && _recordRead<int>(is, "Number of main faults", nfaults);
  ret = ret && _recordRead<double>(is, "Maximum horizontal distance", _xmax);
  ret = ret && _recordRead<double>(is, "Maximum vertical distance", _ymax);
  ret = ret && _recordRead<double>(is, "Dilation along the horizontal axis", _deltax);
  ret = ret && _recordRead<double>(is, "Dilation along the vertical axis", _deltay);
  ret = ret && _recordRead<double>(is, "Mean of thickness distribution", _mean);
  ret = ret && _recordRead<double>(is, "Stdev of thickness distribution", _stdev);
  if (ret) return 1;

  for (int ifam = 0; ifam < nfamilies; ifam++)
  {
    Family family = Family();
    ret = ret && family._deserialize(is, verbose);
    addFamily(family);
  }

  for (int ifault = 0; ifault < nfaults; ifault++)
  {
    Fault fault = Fault();
    ret = ret && fault._deserialize(is, verbose);
    addFault(fault);
  }

  return 0;
}

int Environ::_serialize(std::ostream& os, bool verbose) const
{
  bool ret = true;
  ret = ret && _recordWrite<int>(os, "Number of families", getNFamilies());
  ret = ret && _recordWrite<int>(os, "Number of main faults", getNFaults());
  ret = ret && _recordWrite<double>(os, "Maximum horizontal distance", _xmax);
  ret = ret && _recordWrite<double>(os, "Maximum vertical distance", _ymax);
  ret = ret && _recordWrite<double>(os, "Dilation along the horizontal axis", _deltax);
  ret = ret && _recordWrite<double>(os, "Dilation along the vertical axis", _deltay);
  ret = ret && _recordWrite<double>(os, "Mean of thickness distribution", _mean);
  ret = ret && _recordWrite<double>(os, "Stdev of thickness distribution", _stdev);

  for (int ifam = 0; ifam < getNFamilies(); ifam++)
  {
    ret = ret && _commentWrite(os, "Characteristics of family");
    const Family& family = getFamily(ifam);
    ret = ret && family._serialize(os, verbose);
  }

  /* Loop on the main faults */

  for (int ifault = 0; ifault < getNFaults(); ifault++)
  {
    ret = ret && _commentWrite(os, "Characteristics of main fault");
    const Fault& fault = getFault(ifault);
    ret = ret && fault._serialize(os, verbose);
  }

  return 0;
}

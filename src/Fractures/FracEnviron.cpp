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
#include "Fractures/FracEnviron.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"

FracEnviron::FracEnviron(double xmax,
                         double ymax,
                         double deltax,
                         double deltay,
                         double mean,
                         double stdev)
  : AStringable(),
    ASerializable(),
    _xmax(xmax),
    _ymax(ymax),
    _deltax(deltax),
    _deltay(deltay),
    _mean(mean),
    _stdev(stdev),
    _families(),
    _faults()
{
}

FracEnviron::FracEnviron(const FracEnviron& r)
    : AStringable(r),
      ASerializable(r),
      _xmax(r._xmax),
      _ymax(r._ymax),
      _deltax(r._deltax),
      _deltay(r._deltay),
      _mean(r._mean),
      _stdev(r._stdev),
      _families(r._families),
      _faults(r._faults)
{
}

FracEnviron& FracEnviron::operator=(const FracEnviron& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    ASerializable::operator=(r);
    _xmax = r._xmax;
    _ymax = r._ymax;
    _deltax = r._deltax;
    _deltay = r._deltay;
    _mean = r._mean;
    _stdev = r._stdev;
    _families = r._families;
    _faults = r._faults;
  }
  return *this;
}

FracEnviron::~FracEnviron()
{
}

/**
 * Create a Environ by loading the contents of a Neutral File
 *
 * @param neutralFilename Name of the Neutral File
 * @param verbose         Verbose
 */
FracEnviron* FracEnviron::createFromNF(const String& neutralFilename, bool verbose)
{
  FracEnviron* envir = nullptr;
  std::ifstream is;
  envir = new FracEnviron();
  bool success = false;
  if (envir->_fileOpenRead(neutralFilename, is, verbose))
  {
    success =  envir->deserialize(is, verbose);
  }
  if (! success)
  {
    delete envir;
    envir = nullptr;
  }
  return envir;
}

FracEnviron* FracEnviron::create(double xmax,
                                 double ymax,
                                 double deltax,
                                 double deltay,
                                 double mean,
                                 double stdev)
{
  return new FracEnviron(xmax, ymax, deltax, deltay, mean, stdev);
}

String FracEnviron::toString(const AStringFormat* strfmt) const
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

  for (int i = 0; i < getNFamilies(); i++)
  {
    sstr << toTitle(2, "Family #%d/%d", i + 1, getNFamilies());
    sstr << _families[i].toString(strfmt);
  }

  /* Loop on the faults */

  for (int i = 0; i < getNFaults(); i++)
  {
    sstr << toTitle(2, "Fault #%d/%d", i + 1, getNFaults());
    sstr << _faults[i].toString(strfmt);
  }

  return sstr.str();
}

bool FracEnviron::_deserialize(std::istream& is, bool verbose)
{
  int nfamilies = 0;
  int nfaults = 0;
  bool ret = true;
  ret = ret && _recordRead<int>(is, "Number of families", nfamilies);
  ret = ret && _recordRead<int>(is, "Number of main faults", nfaults);
  ret = ret && _recordRead<double>(is, "Maximum horizontal distance", _xmax);
  ret = ret && _recordRead<double>(is, "Maximum vertical distance", _ymax);
  ret = ret && _recordRead<double>(is, "Dilation along the horizontal axis", _deltax);
  ret = ret && _recordRead<double>(is, "Dilation along the vertical axis", _deltay);
  ret = ret && _recordRead<double>(is, "Mean of thickness distribution", _mean);
  ret = ret && _recordRead<double>(is, "Stdev of thickness distribution", _stdev);
  if (! ret) return ret;

  for (int ifam = 0; ret && ifam < nfamilies; ifam++)
  {
    FracFamily family;
    ret = ret && family.deserialize(is, verbose);
    if (ret) addFamily(family);
  }

  for (int ifault = 0; ret && ifault < nfaults; ifault++)
  {
    FracFault fault;
    ret = ret && fault.deserialize(is, verbose);
    if (ret) addFault(fault);
  }
  return ret;
}

bool FracEnviron::_serialize(std::ostream& os, bool verbose) const
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

  for (int ifam = 0; ret && ifam < getNFamilies(); ifam++)
  {
    ret = ret && _commentWrite(os, "Characteristics of family");
    const FracFamily& family = getFamily(ifam);
    ret = ret && family.serialize(os, verbose);
  }

  /* Loop on the main faults */

  for (int ifault = 0; ret && ifault < getNFaults(); ifault++)
  {
    ret = ret && _commentWrite(os, "Characteristics of main fault");
    const FracFault& fault = getFault(ifault);
    ret = ret && fault.serialize(os, verbose);
  }
  return ret;
}

double FracEnviron::getXextend() const
{
  return _xmax + 2. * _deltax;
}

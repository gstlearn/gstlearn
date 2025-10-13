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
#include "Fractures/FracEnviron.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/SerializeHDF5.hpp"

namespace gstlrn
{
FracEnviron::FracEnviron(double xmax,
                         double ymax,
                         double deltax,
                         double deltay,
                         double mean,
                         double stdev)
  : AStringable()
  , ASerializable()
  , _xmax(xmax)
  , _ymax(ymax)
  , _deltax(deltax)
  , _deltay(deltay)
  , _mean(mean)
  , _stdev(stdev)
  , _families()
  , _faults()
{
}

FracEnviron::FracEnviron(const FracEnviron& r)
  : AStringable(r)
  , ASerializable(r)
  , _xmax(r._xmax)
  , _ymax(r._ymax)
  , _deltax(r._deltax)
  , _deltay(r._deltay)
  , _mean(r._mean)
  , _stdev(r._stdev)
  , _families(r._families)
  , _faults(r._faults)
{
}

FracEnviron& FracEnviron::operator=(const FracEnviron& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    ASerializable::operator=(r);
    _xmax     = r._xmax;
    _ymax     = r._ymax;
    _deltax   = r._deltax;
    _deltay   = r._deltay;
    _mean     = r._mean;
    _stdev    = r._stdev;
    _families = r._families;
    _faults   = r._faults;
  }
  return *this;
}

FracEnviron::~FracEnviron()
{
}

/**
 * Create a Environ by loading the contents of a Neutral File
 *
 * @param NFFilename Name of the Neutral File
 * @param verbose    Verbose
 */
FracEnviron* FracEnviron::createFromNF(const String& NFFilename, bool verbose)
{
  auto* envir = new FracEnviron();
  if (envir->_fileOpenAndDeserialize(NFFilename, verbose)) return envir;
  delete envir;
  return nullptr;
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

  for (Id i = 0; i < getNFamilies(); i++)
  {
    sstr << toTitle(2, "Family #%d/%d", i + 1, getNFamilies());
    sstr << _families[i].toString(strfmt);
  }

  /* Loop on the faults */

  for (Id i = 0; i < getNFaults(); i++)
  {
    sstr << toTitle(2, "Fault #%d/%d", i + 1, getNFaults());
    sstr << _faults[i].toString(strfmt);
  }

  return sstr.str();
}

bool FracEnviron::_deserializeAscii(std::istream& is, bool verbose)
{
  Id nfamilies = 0;
  Id nfaults   = 0;
  bool ret      = true;
  ret           = ret && _recordRead<Id>(is, "Number of families", nfamilies);
  ret           = ret && _recordRead<Id>(is, "Number of main faults", nfaults);
  ret           = ret && _recordRead<double>(is, "Maximum horizontal distance", _xmax);
  ret           = ret && _recordRead<double>(is, "Maximum vertical distance", _ymax);
  ret           = ret && _recordRead<double>(is, "Dilation along the horizontal axis", _deltax);
  ret           = ret && _recordRead<double>(is, "Dilation along the vertical axis", _deltay);
  ret           = ret && _recordRead<double>(is, "Mean of thickness distribution", _mean);
  ret           = ret && _recordRead<double>(is, "Stdev of thickness distribution", _stdev);
  if (!ret) return ret;

  for (Id ifam = 0; ret && ifam < nfamilies; ifam++)
  {
    FracFamily family;
    ret = ret && family._deserializeAscii(is, verbose);
    if (ret) addFamily(family);
  }

  for (Id ifault = 0; ret && ifault < nfaults; ifault++)
  {
    FracFault fault;
    ret = ret && fault._deserializeAscii(is, verbose);
    if (ret) addFault(fault);
  }
  return ret;
}

bool FracEnviron::_serializeAscii(std::ostream& os, bool verbose) const
{
  bool ret = true;
  ret      = ret && _recordWrite<Id>(os, "Number of families", getNFamilies());
  ret      = ret && _recordWrite<Id>(os, "Number of main faults", getNFaults());
  ret      = ret && _recordWrite<double>(os, "Maximum horizontal distance", _xmax);
  ret      = ret && _recordWrite<double>(os, "Maximum vertical distance", _ymax);
  ret      = ret && _recordWrite<double>(os, "Dilation along the horizontal axis", _deltax);
  ret      = ret && _recordWrite<double>(os, "Dilation along the vertical axis", _deltay);
  ret      = ret && _recordWrite<double>(os, "Mean of thickness distribution", _mean);
  ret      = ret && _recordWrite<double>(os, "Stdev of thickness distribution", _stdev);

  for (Id ifam = 0; ret && ifam < getNFamilies(); ifam++)
  {
    ret                      = ret && _commentWrite(os, "Characteristics of family");
    const FracFamily& family = getFamily(ifam);
    ret                      = ret && family._serializeAscii(os, verbose);
  }

  /* Loop on the main faults */

  for (Id ifault = 0; ret && ifault < getNFaults(); ifault++)
  {
    ret                    = ret && _commentWrite(os, "Characteristics of main fault");
    const FracFault& fault = getFault(ifault);
    ret                    = ret && fault._serializeAscii(os, verbose);
  }
  return ret;
}

double FracEnviron::getXextend() const
{
  return _xmax + 2. * _deltax;
}
#ifdef HDF5
bool FracEnviron::_deserializeH5(H5::Group& grp, [[maybe_unused]] bool verbose)
{
  auto fracG = SerializeHDF5::getGroup(grp, "FracEnviron");
  if (!fracG) return false;

  /* Read the grid characteristics */
  bool ret      = true;
  Id nfamilies = 0;
  Id nfaults   = 0;

  ret = ret && SerializeHDF5::readValue(*fracG, "NFamilies", nfamilies);
  ret = ret && SerializeHDF5::readValue(*fracG, "NFaults", nfaults);
  ret = ret && SerializeHDF5::readValue(*fracG, "Xmax", _xmax);
  ret = ret && SerializeHDF5::readValue(*fracG, "Ymax", _ymax);
  ret = ret && SerializeHDF5::readValue(*fracG, "Deltax", _deltax);
  ret = ret && SerializeHDF5::readValue(*fracG, "Deltay", _deltay);
  ret = ret && SerializeHDF5::readValue(*fracG, "Mean", _mean);
  ret = ret && SerializeHDF5::readValue(*fracG, "Stdev", _stdev);

  // Loop on the families
  for (Id ifam = 0; ret && ifam < getNFamilies(); ifam++)
  {
    String locName = "Family" + std::to_string(ifam);
    auto famG      = SerializeHDF5::getGroup(*fracG, locName);
    if (!famG) return false;

    FracFamily family;
    ret = ret && family._deserializeH5(*famG, verbose);
    if (ret) addFamily(family);
  }

  // Loop on the main faults
  for (Id ifault = 0; ret && ifault < getNFaults(); ifault++)
  {
    String locName = "Fault" + std::to_string(ifault);
    auto faultG    = SerializeHDF5::getGroup(*fracG, locName);
    if (!faultG) return false;

    FracFault fault;
    ret = ret && fault._deserializeH5(*faultG, verbose);
    if (ret) addFault(fault);
  }

  return ret;
}

bool FracEnviron::_serializeH5(H5::Group& grp, [[maybe_unused]] bool verbose) const
{
  auto fracG = grp.createGroup("FracEnviron");

  bool ret = true;

  ret = ret && SerializeHDF5::writeValue(fracG, "NFamilies", getNFamilies());
  ret = ret && SerializeHDF5::writeValue(fracG, "NFaults", getNFaults());
  ret = ret && SerializeHDF5::writeValue(fracG, "Xmax", _xmax);
  ret = ret && SerializeHDF5::writeValue(fracG, "Ymax", _ymax);
  ret = ret && SerializeHDF5::writeValue(fracG, "Deltax", _deltax);
  ret = ret && SerializeHDF5::writeValue(fracG, "Deltay", _deltay);
  ret = ret && SerializeHDF5::writeValue(fracG, "Mean", _mean);
  ret = ret && SerializeHDF5::writeValue(fracG, "Stdev", _stdev);

  // Loop on the families
  auto famsG = fracG.createGroup("Families");
  for (Id ifam = 0, nfam = getNFamilies(); ret && ifam < nfam; ifam++)
  {
    const FracFamily& family = getFamily(ifam);
    String locName           = "Family" + std::to_string(ifam);
    auto famG                = famsG.createGroup(locName);
    ret                      = ret && family._serializeH5(famsG, verbose);
  }

  // Loop on the main faults
  auto faultsG = fracG.createGroup("Fauts");
  for (Id ifault = 0, nfault = getNFaults(); ret && ifault < nfault; ifault++)
  {
    const FracFault& fault = getFault(ifault);
    String locName         = "Fault" + std::to_string(ifault);
    auto faultG            = faultsG.createGroup(locName);
    ret                    = ret && fault._serializeH5(faultsG, verbose);
  }

  return ret;
}
#endif
}
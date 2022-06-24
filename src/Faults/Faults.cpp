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
#include "Faults/Faults.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/Utilities.hpp"
#include "geoslib_f.h"

Faults::Faults()
  : AStringable(),
    ASerializable(),
    _faults()
{
}

Faults::Faults(const Faults& r)
    : AStringable(r),
      ASerializable(r),
      _faults(r._faults)
{
}

Faults& Faults::operator=(const Faults& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    ASerializable::operator=(r);
    _faults = r._faults;
  }
  return *this;
}

Faults::~Faults()
{
}

String Faults::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  int nfaults = getNFaults();
  sstr << "Number of Faults = " << nfaults << std::endl;

  for (int i = 0; i < nfaults; i++)
  {
    sstr << "Fault #" << i+1 << std::endl;
    sstr << _faults[i].toString(strfmt);
  }
  return sstr.str();
}

bool Faults::_serialize(std::ostream& os, bool verbose) const
{
  bool ret = true;
  ret = ret && _recordWrite<int>(os, "Number of Faults", getNFaults());
  for (int i = 0; ret && i < getNFaults(); i++)
    ret = ret && _faults[i].serialize(os, verbose);
  return ret;
}

bool Faults::_deserialize(std::istream& is, bool verbose)
{
  int nfaults;
  bool ret = true;
  ret = ret && _recordRead<int>(is, "Number of Faults", nfaults);

  for (int i = 0; ret && i < nfaults; i++)
  {
    Line2D fault;
    ret = ret && fault.deserialize(is, verbose);
    addFault(fault);
  }
  return ret;
}

Faults* Faults::createFromNF(const String& neutralFilename, bool verbose)
{
  Faults* faults = nullptr;
  std::ifstream is;
  faults = new Faults();
  bool success = false;
  if (faults->_fileOpenRead(neutralFilename, is, verbose))
  {
    success =  faults->deserialize(is, verbose);
  }
  if (! success)
  {
    delete faults;
    faults = nullptr;
  }
  return faults;
}

void Faults::addFault(const Line2D& fault)
{
  _faults.push_back(fault);
}


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
#include "Faults/Faults.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Geometry/GeometryHelper.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/SpacePoint.hpp"

namespace gstlrn
{
Faults::Faults()
  : AStringable()
  , ASerializable()
  , _faults()
{
}

Faults::Faults(const Faults& r)
  : AStringable(r)
  , ASerializable(r)
  , _faults(r._faults)
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
  auto nfaults = getNFaults();
  if (nfaults <= 0) return sstr.str();

  sstr << "Number of Faults = " << nfaults << std::endl;

  for (Id i = 0; i < nfaults; i++)
  {
    sstr << "Fault #" << i + 1 << std::endl;
    sstr << _faults[i].toString(strfmt);
  }
  return sstr.str();
}

bool Faults::_serializeAscii(std::ostream& os, bool verbose) const
{
  bool ret = true;
  ret      = ret && _recordWrite<Id>(os, "Number of Faults", getNFaults());
  for (Id i = 0; ret && i < getNFaults(); i++)
    ret = ret && _faults[i]._serializeAscii(os, verbose);
  return ret;
}

bool Faults::_deserializeAscii(std::istream& is, bool verbose)
{
  Id nfaults = 0;
  bool ret    = true;
  ret         = ret && _recordRead<Id>(is, "Number of Faults", nfaults);

  for (Id i = 0; ret && i < nfaults; i++)
  {
    PolyLine2D fault;
    ret = ret && fault._deserializeAscii(is, verbose);
    addFault(fault);
  }
  return ret;
}

Faults* Faults::createFromNF(const String& NFFilename, bool verbose)
{
  auto* faults = new Faults();
  if (faults->_fileOpenAndDeserialize(NFFilename, verbose)) return faults;
  delete faults;
  return nullptr;
}

void Faults::addFault(const PolyLine2D& fault)
{
  _faults.push_back(fault);
}

bool Faults::isSplitByFaultSP(const SpacePoint& P1, const SpacePoint& P2) const
{
  // This is limited to 2D case in RN

  if (getDefaultSpaceType() != ESpaceType::RN || P1.getNDim() != 2)
  {
    messerr("This is limited to 2-D case in RN");
    return false;
  }

  double xt1 = P1.getCoord(0);
  double yt1 = P1.getCoord(1);
  double xt2 = P2.getCoord(0);
  double yt2 = P2.getCoord(1);
  return isSplitByFault(xt1, yt1, xt2, yt2);
}

bool Faults::isSplitByFault(double xt1, double yt1, double xt2, double yt2) const
{
  double xint, yint;

  // Segment bounding box

  double xtmin = MIN(xt1, xt2);
  double xtmax = MAX(xt1, xt2);
  double ytmin = MIN(yt1, yt2);
  double ytmax = MAX(yt1, yt2);

  // Loop on the Fault polylines

  for (Id ifault = 0, nfault = getNFaults(); ifault < nfault; ifault++)
  {
    const PolyLine2D& fault = getFault(ifault);

    const VectorDouble& x = fault.getX();
    const VectorDouble& y = fault.getY();

    // Check if bounding boxes overlap

    if (VH::maximum(x) < xtmin) continue;
    if (xtmax < VH::minimum(x)) continue;
    if (VH::maximum(y) < ytmin) continue;
    if (ytmax < VH::minimum(y)) continue;

    // Loop on the segments of the polyline

    double x1 = x[0];
    double y1 = y[0];
    for (Id ip = 1, np = fault.getNPoints(); ip < np; ip++)
    {
      const double x2 = x[ip];
      const double y2 = y[ip];
      if (GH::segmentIntersect(x1, y1, x2, y2, xt1, yt1, xt2, yt2, &xint, &yint)) return true;
      x1 = x2;
      y1 = y2;
    }
  }
  return false;
}

#ifdef HDF5
bool Faults::_deserializeH5(H5::Group& grp, [[maybe_unused]] bool verbose)
{
  auto faultG = SerializeHDF5::getGroup(grp, "Faults");
  if (!faultG) return false;

  /* Read the grid characteristics */
  bool ret    = true;
  Id nfaults = 0;

  ret = ret && SerializeHDF5::readValue(*faultG, "NFaults", nfaults);

  auto faultsG = SerializeHDF5::getGroup(*faultG, "Lines");
  if (!faultsG) return false;
  for (Id i = 0; ret && i < nfaults; i++)
  {
    String locName = "Line" + std::to_string(i);
    auto lineG     = SerializeHDF5::getGroup(*faultsG, locName);
    if (!lineG) return false;

    PolyLine2D fault;
    ret = ret && fault._deserializeH5(*lineG, verbose);
    addFault(fault);
  }
  return ret;
}

bool Faults::_serializeH5(H5::Group& grp, [[maybe_unused]] bool verbose) const
{
  auto faultG = grp.createGroup("Faults");

  bool ret = true;

  ret = ret && SerializeHDF5::writeValue(faultG, "NFaults", getNFaults());

  auto faultsG = faultG.createGroup("Lines");
  for (Id ifault = 0, nfaults = getNFaults(); ret && ifault < nfaults; ifault++)
  {
    String locName = "Line" + std::to_string(ifault);
    auto lineG     = faultsG.createGroup(locName);

    ret = ret && _faults[ifault]._serializeH5(lineG, verbose);
  }

  return ret;
}
#endif
}
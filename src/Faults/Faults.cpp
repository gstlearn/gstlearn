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

#include "Geometry/GeometryHelper.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/SpacePoint.hpp"

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
  if (nfaults <= 0) return sstr.str();

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
  int nfaults = 0;
  bool ret = true;
  ret = ret && _recordRead<int>(is, "Number of Faults", nfaults);

  for (int i = 0; ret && i < nfaults; i++)
  {
    PolyLine2D fault;
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

bool Faults::isSplitByFault(double xt1,double yt1, double xt2, double yt2) const
{
  double xint, yint;

  // Segment bounding box

  double xtmin = MIN(xt1, xt2);
  double xtmax = MAX(xt1, xt2);
  double ytmin = MIN(yt1, yt2);
  double ytmax = MAX(yt1, yt2);

  // Loop on the Fault polylines

  for (int ifault = 0, nfault = getNFaults(); ifault < nfault; ifault++)
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
    for (int ip = 1, np = fault.getNPoints(); ip < np; ip++)
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

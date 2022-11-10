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
#include "Basic/AStringable.hpp"
#include "Simulation/SimuPartitionParam.hpp"

#include <math.h>

SimuPartitionParam::SimuPartitionParam(int nbtuba,
                                       double intensity,
                                       const VectorDouble& dilate)
    : AStringable(),
      _nbtuba(nbtuba),
      _intensity(intensity),
      _dilate(dilate)
{
}

SimuPartitionParam::SimuPartitionParam(const SimuPartitionParam &r)
    : AStringable(r),
      _nbtuba(r._nbtuba),
      _intensity(r._intensity),
      _dilate(r._dilate)
{
}

SimuPartitionParam& SimuPartitionParam::operator=(const SimuPartitionParam &r)
{
  if (this != &r)
  {
    AStringable::operator =(r);
    _nbtuba = r._nbtuba;
    _intensity = r._intensity;
    _dilate = r._dilate;
  }
  return *this;
}

SimuPartitionParam::~SimuPartitionParam()
{
}

String SimuPartitionParam::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << "Intensity of Poisson Law = " << _intensity << std::endl;
  sstr << "Number of Bands used for valuation simulation = " << _nbtuba << std::endl;
  if (! _dilate.empty())
    sstr << toVector("Dilation (used for Poisson)",_dilate);

  return sstr.str();
}

double SimuPartitionParam::getDilate(int idim) const
{
  if (_dilate.empty()) return 0.;
  return _dilate[idim];
}


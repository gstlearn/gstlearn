/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_old_f.h"

#include "Neigh/ANeighParam.hpp"
#include "Morpho/Morpho.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "Db/Db.hpp"

ANeighParam::ANeighParam(bool flag_xvalid, const ASpace* space)
    : ASpaceObject(space),
      ASerializable()
{
}

ANeighParam& ANeighParam::operator=(const ANeighParam& r)
{
  if (this != &r)
  {
    ASpaceObject::operator=(r);
    ASerializable::operator=(r);
   }
  return *this;
}

ANeighParam::~ANeighParam()
{
}

ANeighParam::ANeighParam(const ANeighParam& r)
    : ASpaceObject(r),
      ASerializable(r)
{
}

bool ANeighParam::isConsistent(const ASpace* /*space*/) const
{
  return true;
}

String ANeighParam::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << "Space dimension = " << getNDim() << std::endl;

  return sstr.str();
}

bool ANeighParam::_deserialize(std::istream& is, bool /*verbose*/)
{
  int ndim = 0;

  bool ret = true;
  ret = ret && _recordRead<int>(is, "Space Dimension", ndim);
  if (ret) setNDim(ndim);
  return ret;
}

bool ANeighParam::_serialize(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;
  ret = ret && _recordWrite<int>(os, "Space Dimension", getNDim());

  return ret;
}

bool ANeighParam::_isDimensionValid(int idim) const
{
  if (idim < 0 || idim >= (int) getNDim())
  {
    messerr("Error in 'idim'(%d). It should lie within [0,%d[",idim,getNDim());
    return false;
  }
  return true;
}

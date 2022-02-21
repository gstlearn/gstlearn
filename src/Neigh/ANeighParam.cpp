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
#include "Neigh/ANeighParam.hpp"
#include "Morpho/Morpho.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Db/Db.hpp"
#include "geoslib_f.h"
#include "geoslib_old_f.h"

ANeighParam::ANeighParam(int ndim, bool flag_xvalid)
    : AStringable(),
      ASerializable(),
      _nDim(ndim),
      _flagXvalid(flag_xvalid),
      _flagContinuous(0),
      _distCont(0.)
{
}

ANeighParam& ANeighParam::operator=(const ANeighParam& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    ASerializable::operator=(r);
    _nDim = r._nDim;
    _flagXvalid = r._flagXvalid;
    _flagContinuous = r._flagContinuous;
    _distCont = r._distCont;
   }
  return *this;
}

ANeighParam::~ANeighParam()
{
}

ANeighParam::ANeighParam(const ANeighParam& r)
    : AStringable(r),
      ASerializable(r),
      _nDim(r._nDim),
      _flagXvalid(r._flagXvalid),
      _flagContinuous(r._flagContinuous),
      _distCont(r._distCont)
{
}

String ANeighParam::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  /* NeighUniqueborhood options */

  sstr << toTitle(0,"Neighborhood characteristics");

  sstr << "Space dimension = " << getNDim() << std::endl;
  if (getFlagXvalid() != 0)
    sstr << "The Cross-Validation Option is switched ON" << std::endl;

  return sstr.str();
}

int ANeighParam::_deserialize(std::istream& is, bool /*verbose*/)
{
  int ndim, flag_xvalid;

  bool ret = _recordRead<int>(is, "Space Dimension", ndim);
  ret = ret && _recordRead<int>(is, "Cross-validation flag", flag_xvalid);
  if (! ret) return 1;

  setNDim(ndim);
  setFlagXvalid(flag_xvalid);

  return 0;
}

int ANeighParam::_serialize(std::ostream& os, bool /*verbose*/) const
{
  bool ret = _recordWrite<int>(os, "Space Dimension", getNDim());
  ret = ret && _recordWrite<int>(os, "Cross-Validation flag", getFlagXvalid());

  return ret ? 0 : 1;
}

bool ANeighParam::_isDimensionValid(int idim) const
{
  if (idim < 0 || idim >= _nDim)
  {
    messerr("Error in 'idim'(%d). It should lie within [0,%d[",idim,_nDim);
    return false;
  }
  return true;
}

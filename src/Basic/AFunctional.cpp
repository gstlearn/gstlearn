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
#include "Basic/AFunctional.hpp"
#include "Basic/VectorNumT.hpp"
#include "Db/Db.hpp"

AFunctional::AFunctional(int ndim)
    : _ndim(ndim)
{
}

AFunctional::AFunctional(const AFunctional &m)
    : _ndim(m._ndim)
{
}

AFunctional& AFunctional::operator=(const AFunctional &m)
{
  if (this != &m)
  {
    _ndim = m._ndim;
  }
  return *this;
}

AFunctional::~AFunctional()
{
}

VectorDouble AFunctional::getFunctionValues(const Db *db, bool useSel) const
{
  if (db == nullptr) return VectorDouble();
  if (_ndim != db->getNDim())
  {
    messerr("You cannot evaluate the function on input Db: they do not have the same Space Dimension");
    return VectorDouble();
  }
  VectorDouble coor(_ndim);
  VectorDouble vec;

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (useSel && ! db->isActive(iech)) continue;

    for (int idim = 0; idim < _ndim; idim++)
      coor[idim] = db->getCoordinate(iech,idim);

    double value = getFunctionValue(coor);
    vec.push_back(value);
  }
  return vec;
}

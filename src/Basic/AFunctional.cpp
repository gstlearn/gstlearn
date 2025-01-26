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

  for (int iech = 0; iech < db->getNSample(); iech++)
  {
    if (useSel && ! db->isActive(iech)) continue;

    for (int idim = 0; idim < _ndim; idim++)
      coor[idim] = db->getCoordinate(iech,idim);

    double value = getFunctionValue(coor);
    vec.push_back(value);
  }
  return vec;
}

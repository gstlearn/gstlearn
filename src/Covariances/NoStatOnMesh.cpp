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
#include "Covariances/NoStatOnMesh.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "geoslib_define.h"

namespace gstlrn
{
NoStatOnMesh::NoStatOnMesh(std::shared_ptr<const Db> dbref,
                           const String& colname)
  : NoStatArray(dbref, colname)
{
}

void NoStatOnMesh::informMeshByMesh(const AMesh* amesh, bool verbose)
{
  DECLARE_UNUSED(verbose);
  _informFieldByMesh(amesh, _tabmesh, true);
}

void NoStatOnMesh::informMeshByApex(const AMesh* amesh, bool verbose)
{
  DECLARE_UNUSED(verbose);
  _informFieldByMesh(amesh, _tabvertices, false);
}

void NoStatOnMesh::_informFieldByMesh(const AMesh* amesh, VectorDouble& tab, bool onMeshes)
{
  Id iatt = _dbNoStat->getUID(_colName);
  if (iatt < 0)
  {
    messerr("The Non-stationary attribute  is not defined in _dbNoStat anymore");
    return;
  }

  /* Loop on the point samples */

  Id nech = onMeshes ? amesh->getNMeshes() : amesh->getNApices();

  tab.resize(nech);

  for (Id iech = 0; iech < nech; iech++)
  {
    tab[iech] = _dbNoStat->getArray(iech, iatt);
  }
}

} // namespace gstlrn

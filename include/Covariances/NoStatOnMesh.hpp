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
#pragma once

#include "Covariances/NoStatArray.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
  class GSTLEARN_EXPORT NoStatOnMesh: public NoStatArray
  {
  public:
    NoStatOnMesh() {};
    NoStatOnMesh(std::shared_ptr<const Db> dbref, const String& colname);

    NoStatOnMesh(const NoStatOnMesh& m) = delete;
    NoStatOnMesh& operator=(const NoStatOnMesh& m) = delete;

    void informMeshByMesh(const AMesh* amesh, bool verbose = false) override;
    void informMeshByApex(const AMesh* amesh, bool verbose = false) override;

    virtual ~NoStatOnMesh() {};

  private:
    void
      _informFieldByMesh(const AMesh* amesh, VectorDouble& tab, bool onMeshes);
  };
} // namespace gstlrn

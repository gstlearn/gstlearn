/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Mesh/SphTriangle.hpp"
#include "Mesh/MeshSpherical.hpp"

/**
 * Meshing defined in the Euclidean space
 */
class GSTLEARN_EXPORT MeshSphericalExt: public MeshSpherical
{
public:
  MeshSphericalExt();
  MeshSphericalExt(const MeshSphericalExt &m);
  MeshSphericalExt& operator=(const MeshSphericalExt &m);
  virtual ~MeshSphericalExt();

  int resetFromDb(Db *dbin,
                  Db *dbout,
                  const String &triswitch = "nqQ",
                  bool verbose = false);
  AMesh* spde_mesh_load(Db *dbin,
                        Db *dbout = nullptr,
                        const VectorDouble &gext = VectorDouble(),
                        const String &triswitch = "-r2",
                        bool verbose = false);

private:
  AMesh* _load2DSph(bool verbose, Db *dbin, Db *dbout, const String &triswitch);
  void _meshesSphLoadVertices(SphTriangle *t);
};

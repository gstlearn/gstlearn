/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Mesh/MeshSpherical.hpp"
#include "ExternalTools/sphtriangle.hpp"

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

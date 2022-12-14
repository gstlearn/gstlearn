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

#include "Mesh/MeshEStandard.hpp"

/**
 * Meshing defined in the Euclidean space
 */
class GSTLEARN_EXPORT MeshEStandardExt: public MeshEStandard
{
public:
  MeshEStandardExt();
  MeshEStandardExt(const MeshEStandardExt &m);
  MeshEStandardExt& operator=(const MeshEStandardExt &m);
  virtual ~MeshEStandardExt();

  int resetFromDb(Db *dbin,
                  Db *dbout = nullptr,
                  const VectorDouble &dilate = VectorDouble(),
                  const String &triswitch = "nqQ",
                  bool verbose = false);
  AMesh* spde_mesh_load(Db *dbin,
                        Db *dbout = nullptr,
                        const VectorDouble &gext = VectorDouble(),
                        const String &triswitch = "nqQ",
                        bool verbose = false);

private:
  int _create1D(int verbose, Db *dbin, Db *dbout, const VectorDouble &dilate);
  int _create2D(bool verbose,
                Db *dbin,
                Db *dbout,
                const VectorDouble &dilate,
                const String& triswitch);
  int _create3D(bool verbose,
                Db *dbin,
                Db *dbout,
                const VectorDouble &dilate,
                const String& triswitch);
  AMesh* _load1D(bool verbose, Db *dbin, Db *dbout, const VectorDouble &gext);
  AMesh* _load2D(bool verbose,
                 Db *dbin,
                 Db *dbout,
                 const VectorDouble &gext,
                 const String& triswitch);
  AMesh* _load3D(bool verbose,
                 Db *dbin,
                 Db *dbout,
                 const VectorDouble &gext,
                 const String &triswitch);
};

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

#include "gstlearn_export.hpp"

#include "Db/DbGrid.hpp"
#include "Basic/ICloneable.hpp"
#include "Mesh/MeshETurbo.hpp"

namespace gstlrn
{
/**
 * \brief
 * Class containing the Data Information organized as a Turbo Meshing
 *
 * This class is derived from the DbGrid class, with a specific decoration: samples
 * or 'nodes' are connected via oriented the 'meshing' information which is
 * stored in this class.
 *
 * Note that in this particular case of Turbo Meshing, the Db must be organized
 * as a grid and the information regarding the meshing is minimum.
 */
class GSTLEARN_EXPORT DbMeshTurbo: public DbGrid
{
public:
  DbMeshTurbo(const VectorInt& nx              = VectorInt(),
              const VectorDouble& dx           = VectorDouble(),
              const VectorDouble& x0           = VectorDouble(),
              const VectorDouble& angles       = VectorDouble(),
              const ELoadBy& order             = ELoadBy::fromKey("SAMPLE"),
              const VectorDouble& tab          = VectorDouble(),
              const VectorString& names        = VectorString(),
              const VectorString& locatorNames = VectorString(),
              bool flag_polarized              = false,
              bool verbose                     = false,
              Id mode                         = 1);
  DbMeshTurbo(const DbMeshTurbo& r);
  DbMeshTurbo& operator=(const DbMeshTurbo& r);
  virtual ~DbMeshTurbo();

public:
  /// ICloneable interface
  IMPLEMENT_CLONING(DbMeshTurbo)

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Db Interface
  bool isMesh() const override { return true; }
  bool mayChangeSampleNumber() const override { return false; }
  bool isConsistent() const override;

  static DbMeshTurbo* create(const VectorInt& nx,
                             const VectorDouble& dx           = VectorDouble(),
                             const VectorDouble& x0           = VectorDouble(),
                             const VectorDouble& angles       = VectorDouble(),
                             const ELoadBy& order             = ELoadBy::fromKey("SAMPLE"),
                             const VectorDouble& tab          = VectorDouble(),
                             const VectorString& names        = VectorString(),
                             const VectorString& locatorNames = VectorString(),
                             bool flag_polarized              = false,
                             bool verbose                     = false);
  static DbMeshTurbo* createFromNF(const String& NFFilename, bool verbose = true);

  Id getNApices() const { return _mesh.getNApices(); }
  Id getNMeshes() const { return _mesh.getNMeshes(); }
  Id getApex(Id imesh, Id rank) const { return _mesh.getApex(imesh, rank); }
  double getCoor(Id imesh, Id rank, Id idim) const
  {
    return _mesh.getCoor(imesh, rank, idim);
  }
  void getCoordinatesPerMeshInPlace(Id imesh, Id rank, VectorDouble& coords) const
  {
    _mesh.getCoordinatesPerMeshInPlace(imesh, rank, coords);
  }
  double getApexCoor(Id i, Id idim) const
  {
    return _mesh.getApexCoor(i, idim);
  }
  void getApexCoordinatesInPlace(Id i, VectorDouble& coords) const
  {
    _mesh.getApexCoordinatesInPlace(i, coords);
  }
  VectorDouble getCoordinatesPerMesh(Id imesh, Id idim, bool flagClose = false) const
  {
    return _mesh.getCoordinatesPerMesh(imesh, idim, flagClose);
  }
  
protected:
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "DbMeshTurbo"; }

private:
  MeshETurbo _mesh;
};
}
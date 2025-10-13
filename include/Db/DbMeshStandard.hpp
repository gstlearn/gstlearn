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

#include "Db/Db.hpp"
#include "Basic/ICloneable.hpp"
#include "Mesh/MeshEStandard.hpp"

namespace gstlrn
{
/**
 * \brief
 * Class containing the Data Information organized as a General Meshing
 *
 * This class is derived from the Db class, with a specific decoration: samples
 * or 'nodes' are connected via oriented the 'meshing' information which is
 * stored in this class.
 */
class GSTLEARN_EXPORT DbMeshStandard: public Db
{
public:
  DbMeshStandard(Id ndim                         = 0,
                 Id napexpermesh                 = 1,
                 const VectorDouble& apices       = VectorDouble(),
                 const VectorInt& meshes          = VectorInt(),
                 const ELoadBy& order             = ELoadBy::fromKey("SAMPLE"),
                 const VectorDouble& tab          = VectorDouble(),
                 const VectorString& names        = VectorString(),
                 const VectorString& locatorNames = VectorString(),
                 bool verbose                     = false);
  DbMeshStandard(const DbMeshStandard& r);
  DbMeshStandard& operator=(const DbMeshStandard& r);
  virtual ~DbMeshStandard();

public:
  /// ICloneable interface
  IMPLEMENT_CLONING(DbMeshStandard)

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Db Interface
  bool isMesh() const override { return true; }
  bool mayChangeSampleNumber() const override { return false; }
  bool isConsistent() const override;

  static DbMeshStandard*
  create(Id ndim,
         Id napexpermesh,
         const VectorDouble& apices,
         const VectorInt& meshes,
         const ELoadBy& order             = ELoadBy::fromKey("SAMPLE"),
         const VectorDouble& tab          = VectorDouble(),
         const VectorString& names        = VectorString(),
         const VectorString& locatorNames = VectorString(),
         bool verbose                     = false);
  static DbMeshStandard* createFromNF(const String& NFFilename, bool verbose = true);
  static DbMeshStandard*
  createFromExternal(const MatrixDense& apices,
                     const MatrixInt& meshes,
                     const ELoadBy& order      = ELoadBy::fromKey("SAMPLE"),
                     const VectorDouble& tab   = VectorDouble(),
                     const VectorString& names = VectorString(),
                     const VectorString& locatorNames = VectorString(),
                     bool verbose                     = false);

  Id getNApices() const { return _mesh.getNApices(); }
  Id getNMeshes() const { return _mesh.getNMeshes(); }
  Id getApex(Id imesh, Id rank) const { return _mesh.getApex(imesh, rank); }
  double getCoor(Id imesh, Id rank, Id idim) const;
  void getCoordinatesPerMeshInPlace(Id imesh, Id rank, VectorDouble& coords) const;
  double getApexCoor(Id i, Id idim) const;
  void getApexCoordinatesInPlace(Id i, VectorDouble& coords) const;
  VectorDouble getCoordinatesPerMesh(Id imesh, Id idim, bool flagClose = false) const;
  
protected:
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "DbMeshStandard"; }

private:
  MeshEStandard _mesh;
};
}
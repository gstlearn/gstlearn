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
  DbMeshStandard(int ndim                         = 0,
                 int napexpermesh                 = 1,
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
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Db Interface
  bool isMesh() const override { return true; }
  bool mayChangeSampleNumber() const override { return false; }
  bool isConsistent() const override;

  static DbMeshStandard*
  create(int ndim,
         int napexpermesh,
         const VectorDouble& apices,
         const VectorInt& meshes,
         const ELoadBy& order             = ELoadBy::fromKey("SAMPLE"),
         const VectorDouble& tab          = VectorDouble(),
         const VectorString& names        = VectorString(),
         const VectorString& locatorNames = VectorString(),
         bool verbose                     = false);
  static DbMeshStandard* createFromNF(const String& neutralFilename,
                                      bool verbose = true);
  static DbMeshStandard*
  createFromExternal(const MatrixRectangular& apices,
                     const MatrixInt& meshes,
                     const ELoadBy& order      = ELoadBy::fromKey("SAMPLE"),
                     const VectorDouble& tab   = VectorDouble(),
                     const VectorString& names = VectorString(),
                     const VectorString& locatorNames = VectorString(),
                     bool verbose                     = false);

  int getNApices() const { return _mesh.getNApices(); }
  int getNMeshes() const { return _mesh.getNMeshes(); }
  int getApex(int imesh, int rank) const { return _mesh.getApex(imesh, rank); }
  double getCoor(int imesh, int rank, int idim) const;
  void getCoordinatesInPlace(int imesh, int rank, VectorDouble& coords) const;
  double getApexCoor(int i, int idim) const;
  void getApexCoordinatesInPlace(int i, VectorDouble& coords) const;
  VectorDouble getCoordinatesPerMesh(int imesh, int idim, bool flagClose = false) const;
  
protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os,
                          bool verbose = false) const override;
  String _getNFName() const override { return "DbMeshStandard"; }

private:
  MeshEStandard _mesh;
};

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

#include "Mesh/AMesh.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixInt.hpp"

namespace gstlrn
{

/**
 * Meshing defined in the Spherical Space
 */
class GSTLEARN_EXPORT MeshSpherical : public AMesh
{
public:
  MeshSpherical(const MatrixDense& apices = MatrixDense(),
                const MatrixInt& meshes = MatrixInt());
  MeshSpherical(const MeshSpherical &m);
  MeshSpherical& operator= (const MeshSpherical &m);
  virtual ~MeshSpherical();

  /// Interface to AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface to AMesh
  Id getNApices() const override;
  Id getNMeshes() const override;
  double getMeshSize(Id imesh) const override;
  Id getApex(Id imesh, Id rank) const override;
  double getCoor(Id imesh, Id rank, Id idim) const override;
  double getApexCoor(Id i, Id idim) const override;
  Id getEmbeddedNDim() const override { return 3; }
  void getEmbeddedCoorPerMesh(Id imesh, Id ic, VectorDouble& coords) const override;
  void getEmbeddedCoorPerApex(Id iapex, VectorDouble& coords) const override;
  void getBarycenterInPlace(Id imesh, vect coord) const override;

  static MeshSpherical* createFromNF(const String& NFFilename, bool verbose = true);
  static MeshSpherical* create(const MatrixDense& apices = MatrixDense(),
                               const MatrixInt& meshes         = MatrixInt());

  Id reset(Id ndim,
           Id napexpermesh,
           const VectorDouble& apices,
           const VectorInt& meshes,
           bool byCol,
           bool verbose = false);
  Id getVariety() const override { return 1; }

  const MatrixDense& getApices() const { return _apices; }
  const MatrixInt& getMeshes() const { return _meshes; }
  VectorVectorInt getMeshesAsVVI() const {return _meshes.getMatrix();}

protected:
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os,bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "MeshSpherical"; }

private:
  void _defineBoundingBox();
  VectorDouble _defineUnits() const;
  Id _recopy(const MeshSpherical& m);
  static double _closestValue(double ref, double coor, double period);
  void _checkConsistency() const;
  bool _weightsInMesh(const VectorDouble& coor,
                      const VectorVectorDouble& corners,
                      double meshsize,
                      VectorDouble& weights,
                      double eps = EPSILON5) const override;
  static void _getCoordOnSphere(double longitude,
                                double latitude,
                                VectorDouble& coords);

private:
  MatrixDense _apices; // Dimension: NRow=napices; Ncol=Ndim(=2)
  MatrixInt         _meshes; // Dimension: Nrow=Nmesh; Ncol=NApexPerMesh
};
}

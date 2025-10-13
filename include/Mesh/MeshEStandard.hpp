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

#include "Basic/VectorNumT.hpp"
#include "Mesh/AMesh.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixInt.hpp"

namespace gstlrn
{
class MeshETurbo;

/**
 * Standard Meshing defined in the Euclidean space
 */
class GSTLEARN_EXPORT MeshEStandard: public AMesh
{
public:
  MeshEStandard();
  MeshEStandard(const MeshEStandard &m);
  MeshEStandard& operator=(const MeshEStandard &m);
  virtual ~MeshEStandard();

  /// Interface to AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface for AMesh
  Id     getNApices() const override;
  Id     getNMeshes() const override;
  Id     getApex(Id imesh, Id rank) const override;
  double  getCoor(Id imesh, Id rank, Id idim) const override;
  double  getApexCoor(Id i, Id idim) const override;
  double  getMeshSize(Id imesh) const override;
  static MeshEStandard* createFromNF(const String& NFFilename, bool verbose = true);
  static MeshEStandard* createFromExternal(const MatrixDense& apices,
                                           const MatrixInt& meshes,
                                           bool verbose = false);

  VectorInt    getMeshList() const { return _meshes.getValues(); }
  VectorDouble getPointList(bool byCol = true) const;
  Id reset(const MatrixDense& apices,
            const MatrixInt& meshes,
            bool verbose = false);
  Id reset(Id ndim,
            Id napexpermesh,
            const VectorDouble &apices,
            const VectorInt &meshes,
            bool byCol = true,
            bool verbose = false);
  Id resetFromTurbo(const MeshETurbo &turbo, bool verbose = false);

  const MatrixDense& getApices() const { return _apices; }
  const MatrixInt& getMeshes() const { return _meshes; }

protected:
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os,bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "MeshEStandard"; }
  void _defineBoundingBox(void);

private:
  bool _coorInMesh(const VectorDouble& coor,
                   Id imesh,
                   double meshsize,
                   VectorDouble& weights) const;
  void _deallocate();
  Id  _recopy(const MeshEStandard &m);
  void _checkConsistency() const;
  void _setApex(Id imesh, Id rank, Id value);
  void _validate();

private:
  MatrixDense _apices; // Dimension: NRow=napices; Ncol=Ndim
  MatrixInt   _meshes; // Dimension: Nrow=Nmesh; Ncol=NApexPerMesh

  friend class DbMeshStandard;
};
}
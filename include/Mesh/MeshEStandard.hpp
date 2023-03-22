/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
#include "Basic/ASerializable.hpp"
#include "Mesh/AMesh.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixInt.hpp"

class MeshETurbo;

/**
 * Meshing defined in the Euclidean space
 */
class GSTLEARN_EXPORT MeshEStandard: public AMesh
{
public:
  MeshEStandard();
  MeshEStandard(const MeshEStandard &m);
  MeshEStandard& operator=(const MeshEStandard &m);
  virtual ~MeshEStandard();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface for AMesh
  int     getNApices() const override;
  int     getNMeshes() const override;
  int     getApex(int imesh, int rank) const override;
  double  getCoor(int imesh, int rank, int idim) const override;
  double  getApexCoor(int i, int idim) const override;
  double  getMeshSize(int imesh) const override;
  cs*     getMeshToDb(const Db *db, bool verbose = false) const override;

  static MeshEStandard* createFromNF(const String& neutralFilename,
                                     bool verbose = true);
  static MeshEStandard* createFromExternal(const MatrixRectangular& apices,
                                           const MatrixInt& meshes,
                                           bool verbose = false);

  VectorInt    getMeshList() const { return _meshes.getValues(); }
  VectorDouble getPointList(bool byCol = true) const;
  int reset(const MatrixRectangular& apices,
            const MatrixInt& meshes,
            bool verbose = false);
  int reset(int ndim,
            int napexpermesh,
            const VectorDouble &apices,
            const VectorInt &meshes,
            bool byCol = true,
            bool verbose = false);
  int reset(int ndim,
            int napexpermesh,
            int npoints,
            int nmeshes,
            const double *apices,
            const int *meshes,
            bool byCol = true,
            bool verbose = false);
  int resetFromTurbo(const MeshETurbo &turbo, bool verbose = false);

  const MatrixRectangular& getApices() const { return _apices; }
  const MatrixInt& getMeshes() const { return _meshes; }

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os,bool verbose = false) const override;
  String _getNFName() const override { return "MeshEStandard"; }
  void _defineBoundingBox(void);

private:
  VectorDouble _defineUnits() const;
  VectorDouble _defineContainers() const;
  bool _coorInMeshContainer(const VectorDouble& coor,
                            int imesh,
                            const VectorDouble& container) const;
  bool _coorInMesh(const VectorDouble& coor,
                   int imesh,
                   double meshsize,
                   VectorDouble& weights) const;
  void _setContainer(VectorDouble &container,
                     int imesh,
                     int idim,
                     double vmin,
                     double vmax) const;
  void _getContainer(const VectorDouble& container,
                     int imesh,
                     int idim,
                     double* vmin,
                     double* vmax) const;
  void _printContainers(const VectorDouble& container) const;
  void _deallocate();
  int  _recopy(const MeshEStandard &m);
  void _checkConsistency() const;
  void _setApex(int imesh, int rank, int value);
  void _validate();

private:
  MatrixRectangular _apices; // Dimension: NRow=napices; Ncol=Ndim
  MatrixInt         _meshes; // Dimension: Nrow=Nmesh; Ncol=NApexPerMesh
};

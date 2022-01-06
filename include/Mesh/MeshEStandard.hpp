/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* Created on: 9 avr. 2019 by N. Desassis                                     */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/Vector.hpp"
#include "Mesh/AMesh.hpp"
#include "Matrix/MatrixRectangular.hpp"

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

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int    getNApices() const override;
  int    getNMeshes() const override;
  int    getApex(int imesh, int rank) const override;
  double getCoor(int imesh, int rank, int idim) const override;
  double getApexCoor(int i, int idim) const override;
  double getMeshSize(int imesh) const override;
  VectorInt    getMeshList() const { return _meshes; }
  VectorDouble getPointList(bool byCol = true) const;

  void   getDuplicates(int verbose,
                       Db *dbin,
                       Db *dbout,
                       int *nbdupl,
                       int **dupl1,
                       int **dupl2) const;
  int resetFromDb(Db *dbin,
                  Db *dbout,
                  const VectorDouble& dilate = VectorDouble(),
                  const String& triswitch = "Q",
                  bool verbose = false);
  int reset(const MatrixRectangular& apices,
            const VectorInt& meshes,
            bool verbose = false);
  int resetOldStyle(int ndim,
                    const VectorDouble& apices,
                    const VectorInt& meshes,
                    bool verbose = false);
  cs*        getMeshToDb(const Db *db, int verbose = 0) const override;
  double*    interpolateMeshToDb(Db *db, double* mtab) const override;
  int        convertFromOldMesh(SPDE_Mesh* s_mesh, int verbose);

private:
  int _create1D(int ndim_ref,
                int verbose,
                Db *dbin,
                Db *dbout,
                const VectorDouble& dilate);
  int _create2D(int ndim_ref,
                int verbose,
                Db *dbin,
                Db *dbout,
                const VectorDouble& dilate,
                const char *triswitch);
  int _create3D(int ndim_ref,
                int verbose,
                Db *dbin,
                Db *dbout,
                const VectorDouble& dilate,
                const char *triswitch);
  void    _defineBoundingBox();
  void    _defineUnits();
  double* _defineContainers() const;
  bool    _coorInMeshContainer(double* coor, int imesh, double* container) const;
  bool    _coorInMesh(double* coor,
                   int imesh,
                   double meshsize,
                   double* weights) const;
  void    _setContainer(double* container,
                     int imesh,
                     int idim,
                     double vmin,
                     double vmax) const;
  void _getContainer(double* container,
                     int imesh,
                     int idim,
                     double* vmin,
                     double* vmax) const;
  void _printContainers(double* container) const;
  void _deallocate();
  int  _recopy(const MeshEStandard &m);
  void _checkConsistency() const;

private:
  MatrixRectangular _apices; // Dimension: NRow=napices; Ncol=Ndim
  VectorInt         _meshes; // TODO MatrixRectangular of Int. Dimension: Nrow=Nmesh; Ncol=NApexPerMesh
  VectorDouble      _units;
};

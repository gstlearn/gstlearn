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
#include "Basic/Vector.hpp"
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
  int     getApex(int imesh, int rank, bool inAbsolute = true) const override;
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
  void   getDuplicates(int verbose,
                       Db *dbin,
                       Db *dbout,
                       int *nbdupl,
                       int **dupl1,
                       int **dupl2) const;
  int resetFromDb(Db *dbin,
                  Db *dbout = nullptr,
                  const VectorDouble& dilate = VectorDouble(),
                  const String& triswitch = "Q",
                  bool verbose = false);
  int reset(const MatrixRectangular& apices,
            const MatrixInt& meshes,
            bool verbose = false);
  int resetOldStyle(int ndim,
                    int napexpermesh,
                    const VectorDouble& apices,
                    const VectorInt& meshes,
                    bool verbose = false);
  int        resetFromTurbo(const MeshETurbo& turbo, bool verbose = false);
  int        convertFromOldMesh(SPDE_Mesh* s_mesh, int verbose);

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os,bool verbose = false) const override;
  String _getNFName() const override { return "MeshEStandard"; }

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

private:
  MatrixRectangular _apices; // Dimension: NRow=napices; Ncol=Ndim
  MatrixInt         _meshes; // Dimension: Nrow=Nmesh; Ncol=NApexPerMesh
};

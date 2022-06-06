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
#include "Mesh/AMesh.hpp"
#include "Matrix/MatrixRectangular.hpp"

/**
 * Meshing defined in 3D for a 2D manifold
 */
class GSTLEARN_EXPORT MeshManifold : public AMesh {

public:
  MeshManifold();
  MeshManifold(const MeshManifold &m);
  MeshManifold& operator= (const MeshManifold &m);
  virtual ~MeshManifold();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface to AMesh
  int     getNApices() const override;
  int     getNMeshes() const override;
  int     getApex(int imesh, int rank) const override;
  double  getCoor(int imesh, int rank, int idim) const override;
  double  getApexCoor(int i, int idim) const override;
  int     getEmbeddedNDim() const override { return 3; }
  void    getEmbeddedCoorPerMesh(int imesh, int ic, VectorDouble& coords) const override;
  void    getEmbeddedCoorPerApex(int iapex, VectorDouble& coords) const override;

  int dumpToNF(const String& neutralFilename, bool verbose = false) const override;
  static MeshManifold* createFromNF(const String& neutralFilename,
                                     bool verbose = false);
  void    getDuplicates(Db *dbin, Db *dbout,
                        int *nbdupl,int **dupl1,int **dupl2, int verbose=0) const;
  cs* getMeshToDb(const Db *db,
                  bool fatal = false,
                  bool verbose = false) const override;
  int getVariety() const { return 2; }

  VectorInt getMeshes() const {return _meshes;}

protected:
  virtual int _deserialize(std::istream& is, bool verbose = false) override;
  virtual int _serialize(std::ostream& os,bool verbose = false) const override;

private:
  void    _defineBoundingBox();
  int     _recopy(const MeshManifold &m);

private:
  MatrixRectangular _apices;
  VectorInt         _meshes; // TODO Transform it into MatrixRectangular of Int
  VectorDouble      _units;
};

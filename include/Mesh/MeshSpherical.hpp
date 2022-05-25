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
 * Meshing defined in the Spherical Space
 */
class GSTLEARN_EXPORT MeshSpherical : public AMesh {

public:
	MeshSpherical();
  MeshSpherical(const MeshSpherical &m);
  MeshSpherical& operator= (const MeshSpherical &m);
	virtual ~MeshSpherical();

	/// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface to AMesh
  int     getNApices() const override;
  int     getNMeshes() const override;
  double  getMeshSize(int imesh) const override;
  int     getApex(int imesh, int rank) const override;
  double  getCoor(int imesh, int rank, int idim) const override;
  double  getApexCoor(int i, int idim) const override;
  int     getEmbeddedDim() const override { return 3; }
  void    getEmbeddedCoor(int imesh, int ic, VectorDouble& coords) const override;

  void    getDuplicates(Db *dbin, Db *dbout,
                        int *nbdupl,int **dupl1,int **dupl2, int verbose=0) const;
  cs* getMeshToDb(const Db *db,
                  bool fatal = false,
                  bool verbose = false) const override;
  double* interpolateMeshToDb(Db *db, double* mtab) const override;
  VectorInt getMeshes() const {return _meshes;}
  int     reset(Db* dbin,Db *dbout,const String& triswitch, int verbose);

private:
  void    _defineBoundingBox();
  void    _defineUnits();
  bool    _coorInMesh(double* coor,int imesh,double* weights) const;
  int     _recopy(const MeshSpherical &m);

private:
  MatrixRectangular _apices;
  VectorInt         _meshes; // TODO Transform it into MatrixRectangular of Int
  VectorDouble      _units;
};

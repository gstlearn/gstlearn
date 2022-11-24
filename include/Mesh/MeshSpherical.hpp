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
#include "Matrix/MatrixInt.hpp"

/**
 * Meshing defined in the Spherical Space
 */
class GSTLEARN_EXPORT MeshSpherical : public AMesh {

public:
	MeshSpherical(const MatrixRectangular& apices = MatrixRectangular(),
	              const MatrixInt& meshes = MatrixInt());
  MeshSpherical(const MeshSpherical &m);
  MeshSpherical& operator= (const MeshSpherical &m);
	virtual ~MeshSpherical();

	/// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface to AMesh
  int     getNApices() const override;
  int     getNMeshes() const override;
  double  getMeshSize(int imesh) const override;
  int     getApex(int imesh, int rank, bool inAbsolute = true) const override;
  double  getCoor(int imesh, int rank, int idim) const override;
  double  getApexCoor(int i, int idim) const override;
  int     getEmbeddedNDim() const override { return 3; }
  void    getEmbeddedCoorPerMesh(int imesh, int ic, VectorDouble& coords) const override;
  void    getEmbeddedCoorPerApex(int iapex, VectorDouble& coords) const override;

  static MeshSpherical* createFromNF(const String &neutralFilename,
                                     bool verbose = true);
  static MeshSpherical* create(const MatrixRectangular &apices = MatrixRectangular(),
                               const MatrixInt &meshes = MatrixInt());

  cs*     getMeshToDb(const Db *db, bool verbose = false) const override;
  int     getVariety() const { return 1; }

  VectorVectorInt getMeshes() const {return _meshes.getMatrix();}
  int     resetFromDb(Db* dbin,Db *dbout,const String& triswitch, int verbose);

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os,bool verbose = false) const override;
  String _getNFName() const override { return "MeshSpherical"; }

private:
  void _defineBoundingBox();
  VectorDouble _defineUnits() const;
  bool _coorInMesh(const VectorDouble& coor,
                   int imesh,
                   double meshsize,
                   VectorDouble& weights,
                   bool flag_approx = true) const;
  int _recopy(const MeshSpherical &m);
  double _closestValue(double ref, double coor, double period) const;

private:
  MatrixRectangular _apices; // Dimension: NRow=napices; Ncol=Ndim(=2)
  MatrixInt         _meshes; // Dimension: Nrow=Nmesh; Ncol=NApexPerMesh
};

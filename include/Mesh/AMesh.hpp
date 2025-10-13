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

#include "Basic/ASerializable.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"

namespace gstlrn
{
class MatrixDense;
class ProjMatrix;
class MatrixInt;
class Db;
class GSTLEARN_EXPORT AMesh: public AStringable, public ASerializable
{
public:
  AMesh();
  AMesh(const AMesh& m);
  AMesh& operator=(const AMesh& m);
  virtual ~AMesh();

  /// Interface to AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface for AMesh
  /*! Returns the number of apex per mesh */
  virtual Id getNApexPerMesh() const { return _nDim + 1; }
  /*! Returns the number of apices */
  virtual Id getNApices() const = 0;
  /*! Returns the number of meshes */
  virtual Id getNMeshes() const = 0;
  /*! Returns the rank of apex 'rank' for mesh 'imesh' */
  virtual Id getApex(Id imesh, Id rank) const = 0;
  /*! Returns coordinate 'idim' of apex 'rank' of mesh 'imesh' */
  virtual double getCoor(Id imesh, Id rank, Id idim) const = 0;
  /*! Returns coordinate 'idim' of apex 'rank' of mesh 'imesh' */
  virtual void getCoordinatesPerMeshInPlace(Id imesh, Id rank, VectorDouble& coords) const;
  /*! Returns coordinate 'idim' of apex 'i' */
  virtual double getApexCoor(Id i, Id idim) const = 0;
  /*! Returns coordinates of apex 'i' */
  virtual void getApexCoordinatesInPlace(Id i, VectorDouble& coords) const;
  /*! Returns the mesh size */
  virtual double getMeshSize(Id imesh) const = 0;
  /*! Initialize the Sparse Matrix for projecting the Db on a Mesh */
  virtual void resetProjFromDb(ProjMatrix* m, const Db* db, Id rankZ = -1, bool verbose = false) const;

  /*! Returns the space variety */
  virtual Id getVariety() const { return 0; }
  virtual Id getEmbeddedNDim() const { return _nDim; }
  virtual void getEmbeddedCoorPerMesh(Id imesh, Id ic, VectorDouble& coords) const;
  virtual void getEmbeddedCoorPerApex(Id iapex, VectorDouble& coords) const;
  virtual void getBarycenterInPlace(Id imesh, vect coord) const;

  /*! Returns the Sparse Matrix for projecting the Mesh to a Db */
  ProjMatrix* createProjMatrix(const Db* db, Id rankZ = -1, bool verbose = false) const;

  /*! Returns the space dimension */
  Id getNDim() const { return _nDim; }
  /*! Returns the minimum of the Bounding box for a given space dimension */
  double getExtendMin(Id idim) const { return _extendMin[idim]; }
  /*! Returns the maximum of the Bounding box for a given space dimension */
  double getExtendMax(Id idim) const { return _extendMax[idim]; }
  /*! Returns the Vector of Extrema of the Bounding Box */
  VectorDouble getExtrema(Id idim) const;
  /*! Returns the list of apexes and meshes */
  void getElements(MatrixDense& apices, MatrixInt& meshes) const;

  Id isCompatibleDb(const Db* db) const;
  VectorDouble getMeshSizes() const;

  /*! Print the list of meshes and apices */
  void printMesh(Id imesh0 = -1) const;
  void printMeshes(Id level = 0, Id nline_max = -1) const;
  /*! Returns Vector of Apex coordinates for space index */
  VectorDouble getCoordinatesPerApex(Id idim) const;
  /*! Returns the list of indices of Meshes sharing the same Apex */
  VectorInt getMeshByApexPair(Id apex1, Id apex2) const;
  /*! Returns the vector of coordinates for a mesh */
  VectorDouble getCoordinatesPerMesh(Id imesh, Id idim, bool flagClose = false) const;
  /*! Returns the coordinates of an Apex */
  VectorDouble getApexCoordinates(Id iapex) const;

  VectorVectorDouble getCoordinatesPerMesh(Id imesh) const;
  VectorVectorDouble getEmbeddedCoordinatesPerMesh(Id imesh = 0) const;
  void getEmbeddedCoordinatesPerMeshInPlace(Id imesh, VectorVectorDouble& vec) const;
  VectorVectorDouble getEmbeddedCoordinatesPerApex() const;

  VectorDouble getDistances(Id iapex0, const VectorInt& japices = VectorInt()) const;

  VectorVectorDouble getAllCoordinates() const;
  MatrixDense getAllApices() const;
  MatrixInt getAllMeshes() const;

  double getCenterCoordinate(Id imesh, Id idim) const;
  VectorVectorDouble getAllCenterCoordinates() const;

  VectorVectorInt getNeighborhoodPerMesh() const;
  VectorVectorInt getNeighborhoodPerApex() const;
  static void dumpNeighborhood(std::vector<VectorInt>& Vmesh, Id nline_max = 1);

protected:
  void _setNDim(Id ndim) { _nDim = ndim; }
  Id _setExtend(const VectorDouble& extendmin, const VectorDouble& extendmax);
  virtual bool _weightsInMesh(const VectorDouble& coor,
                              const VectorVectorDouble& corners,
                              double meshsize,
                              VectorDouble& weights,
                              double eps = EPSILON5) const;
  double _getMeshUnit(const VectorVectorDouble& corners) const;

protected:
  void _recopy(const AMesh& m);
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "AMesh"; }

private:
  bool _isSpaceDimensionValid(Id idim) const;
  void _printMeshListByIndices(Id nline_max = -1) const;
  void _printMeshListByCoordinates(Id nline_max = -1) const;
  Id _findBarycenter(const VectorDouble& target,
                     const VectorDouble& units,
                     Id nb_neigh,
                     VectorInt& neighs,
                     VectorDouble& weight) const;
  VectorDouble _defineUnits(void) const;

private:
  Id _nDim;
  VectorDouble _extendMin;
  VectorDouble _extendMax;
};

typedef std::vector<const AMesh*> VectorMeshes;

GSTLEARN_EXPORT void dumpMeshes(const VectorMeshes& meshes);

} // namespace gstlrn

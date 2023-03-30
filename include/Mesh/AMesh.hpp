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
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"

class MatrixRectangular;
class MatrixInt;
class Db;
class cs; /// TODO : Dependency to csparse to be removed

class GSTLEARN_EXPORT AMesh : public AStringable, public ASerializable
{

public:
	AMesh();
  AMesh(const AMesh &m);
  AMesh& operator= (const AMesh &m);
	virtual ~AMesh();

	/// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface for AMesh
  /*! Returns the number of apex per mesh */
  virtual int getNApexPerMesh() const { return _nDim + 1; }
  /*! Returns the number of apices */
  virtual int getNApices() const = 0;
  /*! Returns the number of meshes */
  virtual int getNMeshes() const = 0;
  /*! Returns the rank of apex 'rank' for mesh 'imesh' */
  virtual int getApex(int imesh, int rank) const = 0;
  /*! Returns coordinate 'idim' of apex 'rank' of mesh 'imesh' */
  virtual double getCoor(int imesh, int rank, int idim) const = 0;
  /*! Returns coordinate 'idim' of apex 'rank' of mesh 'imesh' */
  virtual void getCoordinatesInPlace(int imesh, int rank, VectorDouble& coords) const;
  /*! Returns coordinate 'idim' of apex 'i' */
  virtual double getApexCoor(int i, int idim) const = 0;
  /*! Returns coordinates of apex 'i' */
  virtual void getApexCoordinatesInPlace(int i, VectorDouble& coords) const;
  /*! Returns the mesh size */
  virtual double getMeshSize(int imesh) const = 0;
#ifndef SWIG
  /*! Returns the Sparse Matrix for projecting a Mesh to a Db */
  virtual cs* getMeshToDb(const Db *db,
                          bool verbose = false) const = 0;
#endif
  /*! Returns the space variety */
  virtual int  getVariety() const { return 0; }
  virtual int  getEmbeddedNDim() const { return _nDim; }
  virtual void getEmbeddedCoorPerMesh(int imesh, int ic, VectorDouble& coords) const;
  virtual void getEmbeddedCoorPerApex(int iapex, VectorDouble& coords) const;

  /*! Returns the space dimension */
  int getNDim() const { return _nDim; }
  /*! Returns the minimum of the Bounding box for a given space dimension */
  double getExtendMin(int idim) const { return _extendMin[idim]; }
  /*! Returns the maximum of the Bounding box for a given space dimension */
  double getExtendMax(int idim) const { return _extendMax[idim]; }
  /*! Returns the Vector of Extrema of the Bounding Box */
  VectorDouble getExtrema(int idim) const;
  /*! Returns the list of apexes and meshes */
  void getElements(MatrixRectangular& apices, MatrixInt& meshes) const;

  int  isCompatibleDb(const Db *db) const;
  VectorDouble getMeshSizes() const;

  /*! Print the list of meshes and apices */
  void printMesh(int imesh0) const;
  void printMeshes(int level=0, int nline_max=-1) const;
  /*! Returns Vector of Apex coordinates for space index */
  VectorDouble getCoordinates(int idim) const;
  /*! Returns the list of indices of Meshes sharing the same Apex */
  VectorInt getMeshByApexPair(int apex1, int apex2) const;
  /*! Returns the vector of coordinates for a mesh */
  VectorDouble getCoordinatesPerMesh(int imesh, int idim, bool flagClose=false) const;
  /*! Returns the coordinates of an Apex */
  VectorDouble getApexCoordinates(int iapex) const;

  VectorVectorDouble getCoordinatesPerMesh(int imesh) const;
  VectorVectorDouble getEmbeddedCoordinatesPerMesh(int imesh = 0) const;
  void getEmbeddedCoordinatesPerMesh(int imesh, VectorVectorDouble& coors) const;
  VectorVectorDouble getEmbeddedApexCoordinates() const;

  VectorDouble getDistances(int iapex0, const VectorInt& japices = VectorInt());

  VectorVectorDouble getAllCoordinates() const;

  /// TODO : replace by VectorVectorInt ?
  std::vector<VectorInt> getNeighborhoodPerMesh() const;
  std::vector<VectorInt> getNeighborhoodPerApex() const;
  void dumpNeighborhood(std::vector<VectorInt>& Vmesh);

protected:
  void _setNDim(int ndim) { _nDim = ndim; }
  int  _setExtend(const VectorDouble extendmin, const VectorDouble extendmax);
  bool _weightsInMesh(const VectorDouble& coor,
                      const VectorVectorDouble& corners,
                      double meshsize,
                      VectorDouble& weights,
                      double eps = EPSILON5) const;;
  double _getMeshUnit(const VectorVectorDouble& corners) const;

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os,bool verbose = false) const override;
  String _getNFName() const override { return "AMesh"; }

private:
  void _recopy(const AMesh &m);
  bool _isSpaceDimensionValid(int idim) const;
  void _printMeshListByIndices(int nline_max = -1) const;
  void _printMeshListByCoordinates(int nline_max = -1) const;

private:
  int          _nDim;
  VectorDouble _extendMin;
  VectorDouble _extendMax;
};

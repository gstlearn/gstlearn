/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"
#include "csparse_d.h"

class MatrixRectangular;
class Db;
class SPDE_Mesh;

class GSTLEARN_EXPORT AMesh : public AStringable
{

public:
	AMesh();
  AMesh(const AMesh &m);
  AMesh& operator= (const AMesh &m);
	virtual ~AMesh();

	/// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /*! Returns the space variety */
  int getVariety() const { return _variety; }
  /*! Returns the space dimension */
  int getNDim() const { return _nDim; }
  /*! Set the Variety */
  void setVariety(int variety) { _variety = variety; }
  /*! Set the Space dimension */
  void setNDim(int ndim) { _nDim = ndim; }
  /*! Returns the minimum of the Bounding box for a given space dimension */
  double getExtendMin(int idim) const { return _extendMin[idim]; }
  /*! Returns the maximum of the Bounding box for a given space dimension */
  double getExtendMax(int idim) const { return _extendMax[idim]; }
  /*! Returns the Vector of Extrema of the Bounding Box */
  VectorDouble getExtrema(int idim) const;
  /*! Returns the set of apexes and meshes */
  void getElements(MatrixRectangular& apices, VectorInt& meshes) const;

  int  setExtend(const VectorDouble extendmin, const VectorDouble extendmax);
  void getDuplicates(int verbose, Db *dbin, Db *dbout,
                     int *nbdupl,int **dupl1,int **dupl2) const;
  int  isCompatibleDb(const Db *db) const;
  VectorDouble getMeshSizes() const;

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
  /*! Returns coordinate 'idim' of apex 'i' */
  virtual double getApexCoor(int i, int idim) const = 0;
  /*! Returns the mesh size */
  virtual double getMeshSize(int imesh) const = 0;
  /*! Returns the Sparse Matrix for projecting a Mesh to a Db */
  virtual cs* getMeshToDb(const Db *db, int verbose = 0) const = 0;
  /*! Interpolates an array from Mesh to Db */
  virtual double* interpolateMeshToDb(Db *db, double* mtab) const = 0;
  /*! Print the list of meshes and apices */
  void printMeshes(int imesh) const;

  /*! Convert from New Mesh into Old Mesh */
  SPDE_Mesh* _convertToOldMesh(AMesh* a_mesh) const;
  /*! Returns Vector of Apex coordinates for space index */
  VectorDouble getCoordinates(int idim) const;
  /*! Returns the list of indices of Meshes sharing the same Apex */
  VectorInt getMeshByApexPair(int apex1, int apex2) const;
  /*! Returns the vector of coordinates for an apex */
  VectorDouble getCoordinatesPerMesh(int imesh, int idim, bool flagClose=false) const;

  std::vector<VectorInt> getNeighborhoodPerMesh() const;
  std::vector<VectorInt> getNeighborhoodPerApex() const;
  void dumpNeighborhood(std::vector<VectorInt>& Vmesh);

private:
  void _recopy(const AMesh &m);
  bool _isSpaceDimensionValid(int idim) const;

private:
  int          _variety;
  int          _nDim;
  VectorDouble _extendMin;
  VectorDouble _extendMax;
};

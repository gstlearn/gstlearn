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

#include "../Basic/Grid.hpp"

class MatrixRectangular;
class Db;
class CovAniso;

/**
 * Meshing defined as a Turbo based on a Regular Grid
 * It actually avoids storing all the meshing information
 * and produces quicker methods
 */
class GSTLEARN_EXPORT MeshETurbo: public AMesh
{
public:
  MeshETurbo();
  MeshETurbo(const VectorInt& nx,
             const VectorDouble& dx = VectorDouble(),
             const VectorDouble& x0 = VectorDouble(),
             const VectorDouble& rotmat = VectorDouble(),
             bool flag_polarized = true,
             int verbose = 0);
  MeshETurbo(const Db* db, int verbose = 0);
  MeshETurbo(const MeshETurbo &m);
  MeshETurbo& operator=(const MeshETurbo &r);
  virtual ~MeshETurbo();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int    getNApices() const override;
  int    getNMeshes() const override;
  int    getApex(int imesh, int rank) const override;
  double getCoor(int imesh, int rank, int idim) const override;
  double getApexCoor(int i, int idim) const override;
  double getMeshSize(int imesh) const override;
  void   setPolarized(bool flag) { _isPolarized = flag; }
  void   setMaskArrayFromInt(int* array);
  void   setMaskArrayFromDouble(double* array);
  int initFromExtend(const VectorDouble& extendmin,
                     const VectorDouble& extendmax,
                     const VectorDouble& cellsize,
                     const VectorDouble& rotmat = VectorDouble(),
                     bool flag_polarized = true,
                     int verbose = 0);
  int initFromGrid(const VectorInt& nx,
                   const VectorDouble& dx = VectorDouble(),
                   const VectorDouble& x0 = VectorDouble(),
                   const VectorDouble& rotmat = VectorDouble(),
                   bool flag_polarized = true,
                   int verbose = 0);
  int initFromCova(const CovAniso& cova,
                   const Db* field,
                   double ratio,
                   int nbExt = 0,
                   bool useSel = true,
                   int verbose = 0);
  bool isNodeMasked(int iabs) const;
  cs*  getMeshToDb(const Db *db, int verbose = 0) const override;
  double* interpolateMeshToDb(Db *db, double* mtab) const override;

  const Grid& getGrid() const
  {
    return _grid;
  }

private:
  int  _defineGrid(const VectorDouble& cellsize);
  void _setNumberElementPerCell();
  int  _getPolarized(VectorInt indg) const;
  int  _addWeights(const int verbose,
                  const int icas,
                  const VectorInt& indg0,
                  VectorInt& indgg,
                  const VectorDouble& coor,
                  VectorInt& indices,
                  double *rhs,
                  double *lambda) const;

private:
  void _deallocate();
  void _fromMeshToIndex(int imesh, int *node, int *icas) const;

private:
  Grid _grid;
  int   _nPerCell;
  bool  _isPolarized;
  bool  _isMaskDefined;
  bool* _maskGrid;
};

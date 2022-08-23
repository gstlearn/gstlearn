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

#include "Basic/Grid.hpp"

class MatrixRectangular;
class DbGrid;
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
  MeshETurbo(const DbGrid* dbgrid, int verbose = 0);
  MeshETurbo(const MeshETurbo &m);
  MeshETurbo& operator=(const MeshETurbo &r);
  virtual ~MeshETurbo();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface to AMesh
  int     getNApices() const override;
  int     getNMeshes() const override;
  int     getApex(int imesh, int rank) const override;
  int     getAbsoluteMeshRank(int imesh) const override;
  double  getCoor(int imesh, int rank, int idim) const override;
  double  getApexCoor(int i, int idim) const override;
  double  getMeshSize(int imesh) const override;
  cs*     getMeshToDb(const Db *db, bool fatal = false,
                      bool verbose = false) const override;

  void   setPolarized(bool flag) { _isPolarized = flag; }

  static MeshETurbo* createFromNF(const String& neutralFilename,
                                     bool verbose = true);
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
                   const VectorDouble& sel = VectorDouble(),
                   bool flag_polarized = true,
                   int verbose = 0);
  int initFromCova(const CovAniso& cova,
                   const DbGrid* field,
                   double ratio,
                   int nbExt = 0,
                   bool useSel = true,
                   int verbose = 0);
  const Grid& getGrid() const { return _grid; }

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
  void _deallocate();
  void _fromMeshToIndex(int imesh, int *node, int *icas) const;
  void _fromSelToMeshingMask(const VectorDouble& sel);
  int  _nmeshInCompleteGrid() const;
  bool _isMaskDefined() const { return (! _rankActiveMesh.empty()); }
  int  _fromActiveToAbsolute(int imesh) const;

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os,bool verbose = false) const override;
  String _getNFName() const override { return "MeshETurbo"; }

private:
  Grid  _grid;
  int   _nPerCell;
  bool  _isPolarized;
  VectorInt _rankActiveMesh;
};

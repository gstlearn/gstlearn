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

#include "Basic/VectorNumT.hpp"
#include "Basic/Indirection.hpp"
#include "Basic/Grid.hpp"
#include "Mesh/AMesh.hpp"

class MatrixRectangular;
class DbGrid;
class CovAniso;
class cs;

/**
 * Meshing defined as a Turbo based on a Regular Grid
 * It actually avoids storing all the meshing information
 * and produces faster methods
 */
class GSTLEARN_EXPORT MeshETurbo: public AMesh
{
public:
  MeshETurbo(int mode = 1);
  MeshETurbo(const VectorInt& nx,
             const VectorDouble& dx = VectorDouble(),
             const VectorDouble& x0 = VectorDouble(),
             const VectorDouble& angles = VectorDouble(),
             bool flag_polarized = false,
             bool verbose = false,
             int mode = 1);
  MeshETurbo(const DbGrid *dbgrid,
             bool flag_polarized = false,
             bool verbose = false,
             int mode = 1);
  MeshETurbo(const MeshETurbo &r);
  MeshETurbo& operator=(const MeshETurbo &r);
  virtual ~MeshETurbo();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface to AMesh
  int     getNApices() const override;
  int     getNMeshes() const override;
  int     getApex(int imesh, int rank) const override;
  double  getCoor(int imesh, int rank, int idim) const override;
  void    getCoordinatesPerMeshInPlace(int imesh, int rank, VectorDouble& coords) const override;
  double  getApexCoor(int i, int idim) const override;
  void    getApexCoordinatesInPlace(int i, VectorDouble& coords) const override;
  double  getMeshSize(int imesh) const override;
  void    resetProjMatrix(ProjMatrix* m, const Db *db, int rankZ = -1, bool verbose = false) const override;
  void    setPolarized(bool flag) { _isPolarized = flag; }

  static MeshETurbo* create(const VectorInt& nx,
                            const VectorDouble& dx = VectorDouble(),
                            const VectorDouble& x0 = VectorDouble(),
                            const VectorDouble& angles = VectorDouble(),
                            bool flag_polarized = false,
                            bool verbose = false);
  static MeshETurbo* createFromNF(const String &neutralFilename,
                                  bool verbose = true);
  static MeshETurbo* createFromGrid(const DbGrid *dbgrid,
                                    bool flag_polarized = false,
                                    bool verbose = false,
                                    int mode = 1);
  static MeshETurbo* createFromGridInfo(const Grid *grid,
                                        bool flag_polarized = false,
                                        bool verbose = false,
                                        int mode = 1);
  static MeshETurbo* createFromCova(const CovAniso& cova,
                                    const Db* field,
                                    double ratio,
                                    int nbExt = 0,
                                    bool isPolarized = false,
                                    bool useSel = true,
                                    bool flagNoStatRot = false,
                                    int nxmax = 300,
                                    bool verbose = false);

  int initFromExtend(const VectorDouble& extendmin,
                     const VectorDouble& extendmax,
                     const VectorDouble& cellsize,
                     const VectorDouble& rotmat = VectorDouble(),
                     bool flag_polarized = false,
                     bool verbose = false);
  int initFromGridByMatrix(const VectorInt& nx,
                           const VectorDouble& dx     = VectorDouble(),
                           const VectorDouble& x0     = VectorDouble(),
                           const VectorDouble& rotmat = VectorDouble(),
                           const VectorDouble& sel    = VectorDouble(),
                           bool flag_polarized        = false,
                           bool verbose               = false);
  int initFromGridByAngles(const VectorInt& nx,
                           const VectorDouble& dx     = VectorDouble(),
                           const VectorDouble& x0     = VectorDouble(),
                           const VectorDouble& angles = VectorDouble(),
                           const VectorDouble& sel    = VectorDouble(),
                           bool flag_polarized        = false,
                           bool verbose               = false);
  int initFromCova(const CovAniso& cova,
                   const Db* field,
                   double ratio,
                   int nbExt = 0,
                   bool isPolarized = false,
                   bool useSel = true,
                   bool flagNoStatRot = false,
                   int nxmax = 300,
                   bool verbose = false);
  const Grid& getGrid() const { return _grid; }

  const Indirection& getGridIndirect() const { return _gridIndirect; }
  const Indirection& getMeshIndirect() const { return _meshIndirect; }
  void getApexIndicesInPlace(int i, VectorInt& indg) const;

private:
  int _defineGrid(const VectorDouble& cellsize);
  void _setNElementPerCell();
  int _getPolarized(const constvectint indg) const;
  int _addWeights(int icas,
                  const constvectint indg0,
                  const constvect coor,
                  const vectint indices,
                  const vect lambda,
                  bool verbose = false) const;
  void _deallocate();
  void _getGridFromMesh(int imesh, int *node, int *icas) const;
  void _buildMaskInMeshing(const VectorDouble& sel);
  int  _nmeshInCompleteGrid() const;
  bool _addElementToTriplet(NF_Triplet& NF_T,
                            int iech,
                            const VectorDouble& coor,
                            const VectorInt& indg0,
                            bool verbose) const;
  int _initFromGridInternal(const VectorDouble& sel,
                            bool flag_polarized,
                            bool verbose);

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os,bool verbose = false) const override;
  String _getNFName() const override { return "MeshETurbo"; }

private:
  Grid  _grid;
  int   _nPerCell;
  bool  _isPolarized;
  Indirection _meshIndirect;
  Indirection _gridIndirect;
};

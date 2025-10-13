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

#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
#include "Basic/Indirection.hpp"
#include "Basic/Grid.hpp"
#include "Mesh/AMesh.hpp"



namespace gstlrn
{
class DbGrid;
class CovAniso;
/**
 * Meshing defined as a Turbo based on a Regular Grid
 * It actually avoids storing all the meshing information
 * and produces faster methods
 */
class GSTLEARN_EXPORT MeshETurbo: public AMesh
{
public:
  MeshETurbo(Id mode = 1);
  MeshETurbo(const VectorInt& nx,
             const VectorDouble& dx = VectorDouble(),
             const VectorDouble& x0 = VectorDouble(),
             const VectorDouble& angles = VectorDouble(),
             bool flag_polarized = false,
             bool verbose = false,
             Id mode = 1);
  MeshETurbo(const DbGrid *dbgrid,
             bool flag_polarized = false,
             bool verbose = false,
             Id mode = 1);
  MeshETurbo(const MeshETurbo &r);
  MeshETurbo& operator=(const MeshETurbo &r);
  virtual ~MeshETurbo();

  /// Interface to AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface to AMesh
  Id     getNApices() const override;
  Id     getNMeshes() const override;
  Id     getApex(Id imesh, Id rank) const override;
  double  getCoor(Id imesh, Id rank, Id idim) const override;
  void    getCoordinatesPerMeshInPlace(Id imesh, Id rank, VectorDouble& coords) const override;
  double  getApexCoor(Id i, Id idim) const override;
  void    getApexCoordinatesInPlace(Id i, VectorDouble& coords) const override;
  double  getMeshSize(Id imesh) const override;
  void    resetProjFromDb(ProjMatrix* m, const Db *db, Id rankZ = -1, bool verbose = false) const override;
  void    setPolarized(bool flag) { _isPolarized = flag; }

  static MeshETurbo* create(const VectorInt& nx,
                            const VectorDouble& dx = VectorDouble(),
                            const VectorDouble& x0 = VectorDouble(),
                            const VectorDouble& angles = VectorDouble(),
                            bool flag_polarized = false,
                            bool verbose = false);
  static MeshETurbo* createFromNF(const String &NFFilename, bool verbose = true);
  static MeshETurbo* createFromGrid(const DbGrid *dbgrid,
                                    bool flag_polarized = false,
                                    bool verbose = false,
                                    Id mode = 1);
  static MeshETurbo* createFromGridInfo(const Grid *grid,
                                        bool flag_polarized = false,
                                        bool verbose = false,
                                        Id mode = 1);
  static MeshETurbo* createFromCova(const CovAniso& cova,
                                    const Db* field,
                                    double ratio,
                                    Id nbExt = 0,
                                    bool isPolarized = false,
                                    bool useSel = true,
                                    Id nxmax = 300,
                                    bool verbose = false);

  Id initFromExtend(const VectorDouble& extendmin,
                     const VectorDouble& extendmax,
                     const VectorDouble& cellsize,
                     const VectorDouble& rotmat = VectorDouble(),
                     bool flag_polarized = false,
                     bool verbose = false);
  Id initFromGridByMatrix(const VectorInt& nx,
                           const VectorDouble& dx     = VectorDouble(),
                           const VectorDouble& x0     = VectorDouble(),
                           const VectorDouble& rotmat = VectorDouble(),
                           const VectorDouble& sel    = VectorDouble(),
                           bool flag_polarized        = false,
                           bool verbose               = false);
  Id initFromGridByAngles(const VectorInt& nx,
                           const VectorDouble& dx     = VectorDouble(),
                           const VectorDouble& x0     = VectorDouble(),
                           const VectorDouble& angles = VectorDouble(),
                           const VectorDouble& sel    = VectorDouble(),
                           bool flag_polarized        = false,
                           bool verbose               = false);
  Id initFromCova(const CovAniso& cova,
                   const Db* field,
                   double ratio,
                   Id nbExt = 0,
                   bool isPolarized = false,
                   bool useSel = true,
                   Id nxmax = 300,
                   bool verbose = false);
  const Grid& getGrid() const { return _grid; }

  const Indirection& getGridIndirect() const { return _gridIndirect; }
  const Indirection& getMeshIndirect() const { return _meshIndirect; }
  void getApexIndicesInPlace(Id i, VectorInt& indg) const;
  Id getMeshFromCoordinates(const VectorDouble& coor,
                             VectorInt& indices,
                             VectorDouble& lambdas) const;

private:
  Id _defineGrid(const VectorDouble& cellsize);
  void _setNElementPerCell();
  Id _getPolarized(const constvectint indg) const;
  Id _addWeights(Id icas,
                  const constvectint indg0,
                  const constvect coor,
                  const vectint indices,
                  const vect lambda,
                  bool verbose = false) const;
  void _deallocate();
  void _getGridFromMesh(Id imesh, Id* node, Id* icas) const;
  void _buildMaskInMeshing(const VectorDouble& sel);
  Id  _nmeshInCompleteGrid() const;
  bool _addElementToTriplet(NF_Triplet& NF_T,
                            Id iech,
                            const VectorDouble& coor,
                            const VectorInt& indg0,
                            bool verbose) const;
  Id _initFromGridInternal(const VectorDouble& sel,
                            bool flag_polarized,
                            bool verbose);

protected:
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os,bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "MeshETurbo"; }

private:
  Grid  _grid;
  Id   _nPerCell;
  bool  _isPolarized;
  Indirection _meshIndirect;
  Indirection _gridIndirect;

  /// factor allocations
  mutable std::vector<Id>    _indg;
  mutable std::vector<Id>    _indices;
  mutable std::vector<double> _lambdas;
  mutable std::vector<double> _rhs;
  mutable std::vector<Id>    _indgg;

  friend class DbMeshTurbo;
};

GSTLEARN_EXPORT bool isTurbo(const VectorMeshes& meshes);
}

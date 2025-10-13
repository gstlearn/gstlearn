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
#include "Geometry/Rotation.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
class MatrixSquare;

class GSTLEARN_EXPORT Grid: public AStringable, public ASerializable
{

public:
  Grid(Id ndim                = 0,
       const VectorInt& nx    = VectorInt(),
       const VectorDouble& x0 = VectorDouble(),
       const VectorDouble& dx = VectorDouble());
  Grid(const Grid& r);
  Grid& operator=(const Grid& r);
  virtual ~Grid();

public:
  void initThread() const;
  static VectorInt gridIndices(const VectorInt& nx,
                               const String& string,
                               bool startFromZero = true,
                               bool invert        = true,
                               bool verbose       = false);
  static Id generateMirrorIndex(Id nx, Id ix);

  void resetFromSpaceDimension(Id ndim);
  void resetFromGrid(Grid* grid);
  Id resetFromVector(const VectorInt& nx        = VectorInt(),
                     const VectorDouble& dx     = VectorDouble(),
                     const VectorDouble& x0     = VectorDouble(),
                     const VectorDouble& angles = VectorDouble());
  void setX0(Id idim, double value);
  void setDX(Id idim, double value);
  void setNX(Id idim, Id value);
  void setRotationByMatrix(const MatrixSquare& rotmat);
  void setRotationByVector(const VectorDouble& rotmat);
  void setRotationByAngles(const VectorDouble& angles);
  void setRotationByAngle(double angle);

  Id getNDim() const { return _nDim; }
  double getX0(Id idim) const;
  double getDX(Id idim) const;
  Id getNX(Id idim) const;
  Id getNTotal() const;
  double getCellSize() const;
  double getExtend(Id idim, bool flagCell = false) const;
  double getVolume(bool flagCell = false) const;
  VectorDouble getExtends(bool flagCell = false) const;

  /// Interface to AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  void copyParams(Id mode, const Grid& gridaux);
  double getCoordinate(Id rank, Id idim, bool flag_rotate = true) const;
  VectorDouble getCoordinatesByRank(Id rank, bool flag_rotate = true) const;
  void getCoordinatesByRankInPlace(VectorDouble& coor, Id rank, bool flag_rotate) const;
  VectorDouble getCoordinatesByIndice(const VectorInt& indice,
                                      bool flag_rotate               = true,
                                      const VectorInt& shift         = VectorInt(),
                                      const VectorDouble& dxsPerCell = VectorDouble()) const;
  void getCoordinatesByIndiceInPlace(VectorDouble& coor,
                                     const VectorInt& indice,
                                     bool flag_rotate               = true,
                                     const VectorInt& shift         = VectorInt(),
                                     const VectorDouble& dxsPerCell = VectorDouble()) const;
  VectorDouble getCoordinatesByCorner(const VectorInt& icorner) const;
  void getCoordinatesByCornerInPlace(VectorDouble& coor, const VectorInt& icorner) const;
  VectorDouble getCellCoordinatesByCorner(Id node,
                                          const VectorInt& shift         = VectorInt(),
                                          const VectorDouble& dxsPerCell = VectorDouble()) const;
  void getCellCoordinatesByCornerInPlace(VectorDouble& coor,
                                         Id node,
                                         const VectorInt& shift         = VectorInt(),
                                         const VectorDouble& dxsPerCell = VectorDouble()) const;
#ifndef SWIG
  double indiceToCoordinate(Id idim0,
                            const constvectint indice,
                            const constvect percent = {},
                            bool flag_rotate        = true) const;
#endif // SWIG
  VectorDouble indicesToCoordinate(const VectorInt& indice,
                                   const VectorDouble& percent = VectorDouble()) const;
#ifndef SWIG
  void indicesToCoordinateInPlace(const constvectint indice,
                                  const vect coor,
                                  const constvect percent = {},
                                  bool flag_rotate        = true) const;
#endif // SWIG
  double rankToCoordinate(Id idim0,
                          Id rank,
                          const VectorDouble& percent = VectorDouble()) const;
  VectorDouble rankToCoordinates(Id rank,
                                 const VectorDouble& percent = VectorDouble()) const;
  void rankToCoordinatesInPlace(Id rank,
                                VectorDouble& coor,
                                const VectorDouble& percent = VectorDouble()) const;
#ifndef SWIG
  Id indiceToRank(const constvectint indice) const;
  void rankToIndice(Id rank, vectint indices, bool minusOne = false) const;
#endif // SWIG
  VectorInt coordinateToIndices(const VectorDouble& coor,
                                bool centered = false,
                                double eps    = EPSILON6) const;
  void coordinateToIndicesInPlace(VectorInt& indices,
                                  const VectorDouble& coor,
                                  bool centered,
                                  double eps) const;
  Id coordinateToIndicesInPlace(const VectorDouble& coor,
                                VectorInt& indice,
                                bool centered = false,
                                double eps    = EPSILON6) const;
  Id coordinateToRank(const VectorDouble& coor,
                      bool centered = false,
                      double eps    = EPSILON6) const;
  VectorInt getCenterIndices(bool flagSup = false) const;
  void getCenterIndicesInPlace(VectorInt& indices, bool flagSup) const;
  VectorInt generateGridIndices(const String& string,
                                bool startFromZero = true,
                                bool invert        = true,
                                bool verbose       = false) const;
  bool sampleBelongsToCell(constvect coor,
                           constvect center,
                           const VectorDouble& dxsPerCell) const;
  bool sampleBelongsToCell(const VectorDouble& coor,
                           const VectorDouble& center,
                           const VectorDouble& dxsPerCell = VectorDouble()) const;
  bool sampleBelongsToCell(const VectorDouble& coor,
                           Id rank,
                           const VectorDouble& dxsPerCell = VectorDouble()) const;
  VectorDouble getRotAngles() const { return _rotation.getAngles(); }
  VectorDouble getRotMat() const { return _rotation.getMatrixDirect().getValues(); }
  double getRotAngle(Id idim) const { return _rotation.getAngle(idim); }
  VectorInt getNXs() const { return _nx; }
  VectorDouble getX0s() const { return _x0; }
  VectorDouble getDXs() const { return _dx; }
  const Rotation& getRotation() const { return _rotation; }
  void setRotation(const Rotation& rotation) { _rotation = rotation; }
  bool isSame(const Grid& grid) const;
  bool isSameMesh(const Grid& grid) const;
  bool isRotated() const { return _rotation.isRotated(); }
  bool isSameRotation(const Grid& grid) const { return _rotation.isSame(grid.getRotation()); }
  VectorDouble getAxis(Id idim) const;

  void iteratorInit(const VectorInt& order = VectorInt());
  VectorInt iteratorNext(void);
  void iteratorNext(std::vector<Id>& indices);
  bool empty() const;

  void dilate(Id mode,
              const VectorInt& nshift,
              VectorInt& nx,
              VectorDouble& dx,
              VectorDouble& x0) const;
  void multiple(const VectorInt& nmult,
                bool flagCell,
                VectorInt& nx,
                VectorDouble& dx,
                VectorDouble& x0) const;
  void divider(const VectorInt& nmult,
               bool flagCell,
               VectorInt& nx,
               VectorDouble& dx,
               VectorDouble& x0) const;
  Id getMirrorIndex(Id idim, Id ix) const;
  bool isInside(const VectorInt& indices) const;

protected:
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "Grid"; }

private:
  const MatrixSquare& _getRotMat() const { return _rotation.getMatrixDirect(); }
  const MatrixSquare& _getRotInv() const { return _rotation.getMatrixInverse(); }
  void _allocate();
  void _recopy(const Grid& r);
  bool _isSpaceDimensionValid(Id idim) const;

private:
  Id _nDim;
  VectorInt _nx;
  VectorDouble _x0;
  VectorDouble _dx;
  Rotation _rotation;

  // Iterator

  Id _iter;
  Id _nprod;
  std::vector<Id> _counts;
  VectorInt _order;
  VectorInt _indices;
  mutable VectorInt _dummy;
  // Some working vectors, defined in order to avoid too many allocations

  friend class DbGrid;
};
} // namespace gstlrn
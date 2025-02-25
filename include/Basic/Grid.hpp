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
#include "geoslib_define.h"
#include "Geometry/Rotation.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"

class GridOld;
class MatrixSquareGeneral;

class GSTLEARN_EXPORT Grid: public AStringable, public ASerializable
{

public:
  Grid(int ndim = 0,
       const VectorInt &nx = VectorInt(),
       const VectorDouble &x0 = VectorDouble(),
       const VectorDouble &dx = VectorDouble());
  Grid(const Grid &r);
  Grid& operator=(const Grid &r);
  virtual ~Grid();

public:
  static VectorInt gridIndices(const VectorInt &nx,
                               const String &string,
                               bool startFromZero = true,
                               bool invert = true,
                               bool verbose = false);
  static int generateMirrorIndex(int nx, int ix);

  void resetFromSpaceDimension(int ndim);
  void resetFromGrid(Grid* grid);
  int  resetFromVector(const VectorInt& nx = VectorInt(),
                       const VectorDouble& dx = VectorDouble(),
                       const VectorDouble& x0 = VectorDouble(),
                       const VectorDouble& angles = VectorDouble());
  void    setX0(int idim,double value);
  void    setDX(int idim,double value);
  void    setNX(int idim,int    value);
  void    setRotationByMatrix(const MatrixSquareGeneral& rotmat);
  void    setRotationByVector(const VectorDouble& rotmat);
  void    setRotationByAngles(const VectorDouble& angles);
  void    setRotationByAngle(double angle);

  int     getNDim() const { return _nDim; }
  double  getX0(int idim) const;
  double  getDX(int idim) const;
  int     getNX(int idim) const;
  int     getNTotal() const;
  double  getCellSize() const;
  double  getExtend(int idim, bool flagCell = false) const;
  double  getVolume(bool flagCell = false) const;
  VectorDouble getExtends(bool flagCell = false) const;

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  void    copyParams(int mode, const Grid& gridaux);
  double  getCoordinate(int rank, int idim, bool flag_rotate=true) const;
  VectorDouble getCoordinatesByRank(int rank, bool flag_rotate=true) const;
  VectorDouble getCoordinatesByIndice(const VectorInt &indice,
                                      bool flag_rotate = true,
                                      const VectorInt& shift = VectorInt(),
                                      const VectorDouble& dxsPerCell = VectorDouble()) const;
  VectorDouble getCoordinatesByCorner(const VectorInt& icorner) const;
  VectorDouble getCellCoordinatesByCorner(int node,
                                          const VectorInt& shift = VectorInt(),
                                          const VectorDouble& dxsPerCell = VectorDouble()) const;
#ifndef SWIG
  double indiceToCoordinate(int idim0,
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
  double rankToCoordinate(int idim0,
                          int rank,
                          const VectorDouble& percent = VectorDouble()) const;
  VectorDouble rankToCoordinates(int rank,
                                 const VectorDouble& percent = VectorDouble()) const;
  void rankToCoordinatesInPlace(int rank,
                                VectorDouble& coor,
                                const VectorDouble& percent = VectorDouble()) const;
#ifndef SWIG
  int     indiceToRank(const constvectint indice) const;
  void    rankToIndice(int rank, vectint indices, bool minusOne = false) const;
#endif // SWIG
  VectorInt coordinateToIndices(const VectorDouble &coor,
                                bool centered = false,
                                double eps = EPSILON6) const;
  int coordinateToIndicesInPlace(const VectorDouble &coor,
                                 VectorInt &indice,
                                 bool centered = false,
                                 double eps = EPSILON6) const;
  int coordinateToRank(const VectorDouble &coor,
                       bool centered = false,
                       double eps = EPSILON6) const;
  VectorInt getCenterIndices(bool flagSup = false) const;
  VectorInt generateGridIndices(const String &string,
                                bool startFromZero = true,
                                bool invert = true,
                                bool verbose = false) const;
  bool sampleBelongsToCell(constvect coor,
                           constvect center,
                           const VectorDouble& dxsPerCell) const;
  bool sampleBelongsToCell(const VectorDouble& coor,
                           const VectorDouble& center,
                           const VectorDouble &dxsPerCell = VectorDouble()) const;
  bool sampleBelongsToCell(const VectorDouble& coor,
                           int rank,
                           const VectorDouble &dxsPerCell = VectorDouble()) const;
  VectorDouble    getRotAngles() const { return _rotation.getAngles(); }
  VectorDouble    getRotMat() const { return _rotation.getMatrixDirect().getValues(); }
  double getRotAngle(int idim) const { return _rotation.getAngle(idim); }
  VectorInt       getNXs() const { return _nx; }
  VectorDouble    getX0s() const { return _x0; }
  VectorDouble    getDXs() const { return _dx; }
  const Rotation&       getRotation() const { return _rotation; }
  bool  isSame(const Grid& grid) const;
  bool  isSameMesh(const Grid& grid) const;
  bool  isRotated() const { return _rotation.isRotated(); }
  bool  isSameRotation(const Grid& grid) const { return _rotation.isSame(grid.getRotation()); }
  VectorDouble getAxis(int idim) const;

  void iteratorInit(const VectorInt& order = VectorInt());
  VectorInt iteratorNext(void);
  void iteratorNext(std::vector<int>& indices);
  bool empty() const;

  void dilate(int mode,
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
  int getMirrorIndex(int idim, int ix) const;
  bool isInside(const VectorInt& indices) const;

  /// Interface for ASerializable
  bool _deserialize(std::istream& is, bool verbose = false) override;
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serialize(std::ostream& os, bool verbose = false) const override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
  String _getNFName() const override { return "Grid"; }

private:
  const MatrixSquareGeneral& _getRotMat() const { return _rotation.getMatrixDirect(); }
  const MatrixSquareGeneral& _getRotInv() const { return _rotation.getMatrixInverse(); }
  void _allocate();
  void _recopy(const Grid &r);
  bool _isSpaceDimensionValid(int idim) const;

private:
  int          _nDim;
  VectorInt    _nx;
  VectorDouble _x0;
  VectorDouble _dx;
  Rotation     _rotation;

  // Iterator
  int       _iter;
  int       _nprod;
  std::vector<int> _counts;
  VectorInt _order;
  VectorInt _indices;

  // Some working vectors, defined in order to avoid too many allocations
  mutable VectorInt    _iwork0;
  mutable VectorDouble _work1;
  mutable std::vector<double> _work2;
};

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
#include "geoslib_define.h"
#include "Basic/Rotation.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Basic/AStringable.hpp"

class GridOld;

class GSTLEARN_EXPORT GridC : public AStringable
{

public:
  GridC(int ndim = 0,
        const VectorInt& nx    = VectorInt(),
        const VectorDouble& x0 = VectorDouble(),
        const VectorDouble& dx = VectorDouble());
  GridC(const GridC &r);
  GridC& operator= (const GridC &r);
	virtual ~GridC();

public:
  void    init(int ndim);
  void    init(GridC* grid);
  int     init(const VectorInt&    nx = VectorInt(),
               const VectorDouble& dx = VectorDouble(),
               const VectorDouble& x0 = VectorDouble(),
               const VectorDouble& angles = VectorDouble());
  void    setX0(int idim,double value);
  void    setDX(int idim,double value);
  void    setNX(int idim,int    value);
  void    setRotationFromMatrix(const MatrixSquareGeneral& rotmat);
  void    setRotationFromMatrix(const VectorDouble& rotmat);
  void    setRotationFromAngles(VectorDouble angles);
  void    setRotationFromAngles(double angle);

  int     getNDim() const { return _nDim; }
  double  getX0(int idim) const;
  double  getDX(int idim) const;
  int     getNX(int idim) const;
  int     getNTotal() const;
  double  getCellSize() const;

  virtual String toString(int level = 0) const override;

  void    copyParams(int mode, const GridC& gridaux);
  double  getCoordinate(int rank, int idim, bool flag_rotate=true) const;
  VectorDouble getCoordinate(int rank, bool flag_rotate=true) const;
  VectorDouble getCoordinate(const VectorInt& indice, bool flag_rotate=true) const;
  VectorDouble getCoordinateFromCorner(const VectorInt& icorner) const;
  double indiceToCoordinate(int idim0,
                            const VectorInt& indice,
                            const VectorDouble& percent = VectorDouble()) const;
  void indiceToCoordinate(const VectorInt& indice,
                          VectorDouble& coor,
                          const VectorDouble& percent) const;
  VectorDouble indiceToCoordinate(const VectorInt& indice,
                                  const VectorDouble& percent = VectorDouble()) const;
  double rankToCoordinate(int idim0,
                          int rank,
                          const VectorDouble& percent = VectorDouble()) const;
  VectorDouble rankToCoordinate(int rank,
                                const VectorDouble& percent = VectorDouble()) const;
  void rankToCoordinate(int rank,
                        VectorDouble& coor,
                        const VectorDouble& percent = VectorDouble()) const;
  int     indiceToRank(const VectorInt& indice) const;
  void    rankToIndice(int node,VectorInt& indice, bool minusOne = false) const;
  int     coordinateToIndice(const VectorDouble& coor,
                             VectorInt& indice,
                             double eps = EPSILON6) const;
  int     coordinateToRank(const VectorDouble& coor, double eps = EPSILON6) const;

  const VectorDouble    getRotAngles() const { return _rotation.getAngles(); }
  const VectorDouble    getRotMat() const { return _rotation.getMatrixDirect().getValues(); }
  double getRotAngles(int idim) const { return _rotation.getAngles(idim); }
  const VectorInt       getNX() const { return _nx; }
  const VectorDouble    getX0() const { return _x0; }
  const VectorDouble    getDX() const { return _dx; }
  const Rotation&       getRotation() const { return _rotation; }
  bool  isSame(const GridC& grid) const;
  bool  isSameMesh(const GridC& grid) const;
  bool  isRotated() const { return _rotation.isRotated(); }
  bool  isSameRotation(const GridC& grid) const { return _rotation.isSame(grid.getRotation()); }
  VectorDouble getAxis(int idim) const;

  void iteratorInit(const VectorInt& order = VectorInt());
  VectorInt iteratorNext(void);
  bool empty() const;

  void dilate(int mode,
              const VectorInt& nshift,
              VectorInt& nx,
              VectorDouble& dx,
              VectorDouble& x0) const;
  void multiple(const VectorInt& nmult,
                int flag_cell,
                VectorInt& nx,
                VectorDouble& dx,
                VectorDouble& x0) const;
  void divider(const VectorInt& nmult,
               int flag_cell,
               VectorInt& nx,
               VectorDouble& dx,
               VectorDouble& x0) const;

private:
  const MatrixSquareGeneral _getRotMat() const { return _rotation.getMatrixDirect(); }
  const MatrixSquareGeneral _getRotInv() const { return _rotation.getMatrixInverse(); }
  void _allocate();
  void _recopy(const GridC &r);
  bool _isValid(int idim) const;

private:
  int          _nDim;
  VectorInt    _nx;
  VectorDouble _x0;
  VectorDouble _dx;
  Rotation     _rotation;

  // Iterator
  int       _iter;
  int       _nprod;
  VectorInt _counts;
  VectorInt _order;
  VectorInt _indices;
};

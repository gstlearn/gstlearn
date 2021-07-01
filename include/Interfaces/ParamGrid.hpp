#ifndef PARAM_GRID_HPP
#define PARAM_GRID_HPP

#include "Interfaces/interface_d.hpp"
#include "Basic/GridC.hpp"

class ParamGrid {
public:
  ParamGrid();
  ParamGrid(VectorInt nx, VectorDouble x0, VectorDouble dx, VectorDouble Rotation, int cell_order);
  ~ParamGrid();

  VectorInt     getNx() const;
  VectorDouble  getX0() const;
  VectorDouble  getDx() const;
  VectorDouble  getRotation() const;
  int           getCellOrder() const;
  VectorDouble  getValues(int i) const;
  void          fromGeoslib(GridC grid);
  void          reset();

private:
  VectorInt     _Nx;
  VectorDouble  _X0;
  VectorDouble  _Dx;
  VectorDouble  _Rotation;
  int           _CellOrder;
};

#endif

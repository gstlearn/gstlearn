#pragma once

#include "gstlearn_export.hpp"
#include "Interfaces/interface_d.hpp"

#include "Basic/Grid.hpp"
#include "Db/ELoadBy.hpp"

class GSTLEARN_EXPORT ParamGrid {
public:
  ParamGrid();
  ParamGrid(VectorInt nx,
            VectorDouble x0,
            VectorDouble dx,
            VectorDouble Rotation,
            const ELoadBy& cell_order);
  ~ParamGrid();

  VectorInt      getNx() const;
  VectorDouble   getX0() const;
  VectorDouble   getDx() const;
  VectorDouble   getRotation() const;
  const ELoadBy& getCellOrder() const;
  VectorDouble   getValues(int i) const;
  void           fromGeoslib(Grid grid);
  void           reset();

private:
  VectorInt     _Nx;
  VectorDouble  _X0;
  VectorDouble  _Dx;
  VectorDouble  _Rotation;
  ELoadBy       _CellOrder;
};


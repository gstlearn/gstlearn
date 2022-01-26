#include "Interfaces/ParamGrid.hpp"
#include <iostream>
ParamGrid::ParamGrid()
:_Nx()
, _X0()
, _Dx()
, _Rotation()
, _CellOrder()
{
}

ParamGrid::ParamGrid(VectorInt nx,
                     VectorDouble x0,
                     VectorDouble dx,
                     VectorDouble Rotation,
                     const ELoadBy& cell_order)
:_Nx(nx)
, _X0(x0)
,_Dx(dx)
, _Rotation(Rotation)
,_CellOrder(cell_order)
{
}

ParamGrid::~ParamGrid()
{
}


VectorInt  ParamGrid::getNx() const
{
  return(_Nx);
}

VectorDouble  ParamGrid::getX0() const
{
  return(_X0);
}

VectorDouble  ParamGrid::getDx() const
{
  return(_Dx);
}

VectorDouble  ParamGrid::getRotation() const
{
  return(_Rotation);
}

const ELoadBy& ParamGrid::getCellOrder() const
{
  return(_CellOrder);
}

/****************************************************************************/
/*
**  Function returning a Vector of double containing the coordinate of the 
**  desired dimension
**  
**  \param[in] i    dimension desired.
**
*****************************************************************************/
VectorDouble ParamGrid::getValues(int i) const
{
  int    nx = getNx()[i];
  double x0 = getX0()[i];
  double dx = getDx()[i];
  VectorDouble res(nx);

  for (int i = 0; i < nx; i++)
  {
    res[i] = x0 + dx * i;
  }
  return(res);
}

/****************************************************************************/
/*
**  Clear all the attribute of ParamGrid
**
*****************************************************************************/
void ParamGrid::reset()
{
  _Nx.clear();
  _X0.clear();
  _Dx.clear();
  _Rotation.clear();
  _CellOrder = ELoadBy::COLUMN;
}

/****************************************************************************/
/*
**  fill a p_grid from a geoslib Grid object
**
*****************************************************************************/
void ParamGrid::fromGeoslib(Grid grid)
{
  int i = 0;
  reset();
  _CellOrder = ELoadBy::COLUMN;
  while ( i < grid.getNDim())
  {
    _Nx.push_back(grid.getNX(i));
    _X0.push_back(grid.getX0(i));
    _Dx.push_back(grid.getDX(i));
    if (grid.isRotated())
      _Rotation.push_back(grid.getRotAngle(i));
    i++;
  }
}

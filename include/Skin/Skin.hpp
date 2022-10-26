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
#include "Skin/ISkinFunctions.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/Vector.hpp"

class GSTLEARN_EXPORT Skin
{
public:
  Skin(const ISkinFunctions* skf, DbGrid* dbgrid = nullptr);
  Skin(const Skin& r);
  Skin& operator=(const Skin& r);
  virtual ~Skin();

  int gridShift(int lec, int dir);
  int init(bool verbose = false);
  int remains(bool verbose = false);
  void getNext(int *rank, int *ipos);
  int unstack(int rank0, int ipos0);
  void skinPrint();

private:
  double _getWeight(int ipos, int idir);
  int    _gridShift(const VectorInt& indg0, int dir);
  void   _cellDelete(int rank);
  int    _cellAlreadyFilled(int ipos);
  void   _cellModify(int rank, double energy);
  int    _cellAdd(int ipos, double energy);
  int    _getNDim() const;

private:
  const ISkinFunctions* _skf;
  DbGrid* _dbgrid;
  int _nxyz;
  int _nval;
  int _date;
  int _nvalMax;
  double  _total;
  double  _totalMax;
  VectorInt _address;
  VectorDouble _energy;
};

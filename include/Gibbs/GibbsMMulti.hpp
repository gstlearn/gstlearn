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
#include "Gibbs/GibbsMulti.hpp"
#include "Basic/Vector.hpp"

class Db;
class Model;

class GSTLEARN_EXPORT GibbsMMulti: public GibbsMulti
{
  typedef struct
  {
    int _size;
    int _nbgh;
    int _pivot;
    VectorInt _mvRanks;
    VectorVectorDouble _weights;
  } WgtVect;

public:
  GibbsMMulti();
  GibbsMMulti(Db* db, Model* model);
  GibbsMMulti(const GibbsMMulti &r);
  GibbsMMulti& operator=(const GibbsMMulti &r);
  virtual ~GibbsMMulti();

  void update(VectorVectorDouble& y,
              int isimu,
              int ipgs,
              int iter) override;
  int covmatAlloc(bool verbose) override;

  void setEps(double eps) { _eps = eps; }
  void setStoreTables(bool storeTables) { _storeTables = storeTables; }

private:
  void _display(int iact) const;
  void _display() const;
  int  _getVariableNumber() const;
  void _tableStore(int mode, const cs* Cmat);
  void _getWeights(int iech, WgtVect& area) const;
  int  _calculateWeights(int iech, WgtVect& area, double tol = EPSILON3) const;
  bool _checkForInternalStorage(bool verbose = false);
  void _storeAllWeights(bool verbose = false);
  int  _getSizeOfArea(const WgtVect& area) const;

private:
  cs*       _Ln;
  VectorInt _Pn;
  double    _eps;
  bool      _storeTables;

  // Mutable arrays (declared to speed up the process)
  mutable VectorDouble _b;
  mutable VectorDouble _x;
  mutable std::vector<WgtVect> _areas;
};

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
#include "Basic/HDF5format.hpp"

#include <vector>

class Db;
class Model;

class GSTLEARN_EXPORT GibbsMMulti: public GibbsMulti
{
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
  void cleanup() override;

  bool getFlagStoreInternal() const { return _flagStoreInternal; }
  void setFlagStoreInternal(bool flagStoreInternal) { _flagStoreInternal = flagStoreInternal; }

private:
  int  _getVariableNumber() const;
  void _tableStore(int mode, const cs* Cmat);
  void _getWeights(int iact0) const;
  int  _calculateWeights(int iact0, double tol = EPSILON3) const;
  int  _storeAllWeights(bool verbose = false);
  int  _getSizeOfWeights(const VectorDouble& weights) const;
  void _getEstimate(int ipgs0,
                    int ivar0,
                    int iact0,
                    int icase,
                    VectorVectorDouble& y,
                    double *yk,
                    double *vark) const;

private:
  cs*        _Ln;
  VectorInt  _Pn;
  double     _eps;
  bool       _storeTables;
  HDF5format _hdf5;
  bool       _flagStoreInternal;

  // Mutable arrays (declared to speed up the process)
  mutable VectorDouble _b;
  mutable VectorDouble _x;
  mutable VectorDouble _weights;
  mutable VectorVectorDouble _areas;
};

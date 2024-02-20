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

#include "Gibbs/GibbsMulti.hpp"
#include "Basic/HDF5format.hpp"

#include <vector>

class Db;
class Model;
class cs; // TODO to be removed

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
  int covmatAlloc(bool verbose, bool verboseTimer = false) override;

  void setEps(double eps) { _eps = eps; }
  void cleanup() override;

  bool getFlagStoreInternal() const { return _flagStoreInternal; }
  void setFlagStoreInternal(bool flagStoreInternal) ;

private:
  int  _getVariableNumber() const;
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
  MatrixSparse* _Cmat;
  double     _eps;
  HDF5format _hdf5;
  bool       _flagStoreInternal;

  // Mutable arrays (declared to speed up the process)
  mutable VectorDouble _b;
  mutable VectorDouble _x;
  mutable VectorDouble _weights;
  mutable VectorVectorDouble _areas;
};

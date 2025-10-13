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
#include "Basic/VectorNumT.hpp"

#include <string>

#define TO_ndim 2
#define TO_ncorner 3
#define TO_npercell 2

namespace gstlrn
{ 
typedef struct
{
  VectorInt rows;
  VectorInt cols;
  VectorDouble values;
} TripletND;

/**
 * \brief Turbo Optimizer for a specific 2-D environment,
 * \brief with an isotropic Matérn Model
 */
class GSTLEARN_EXPORT TurboOptimizer
{
private:
  bool _isCalculated;
  Id _nx;
  Id _ny;
  double _dx;
  double _dy;
  double _x0;
  double _y0;
  double _scale;
  double _sill;
  Id _param;
  Id _poncif;
  Id _center;
  Id _nxred;
  Id _half;
  Id _flagOne;
  VectorDouble _Blin;
  VectorDouble _TildeC_T;
  VectorDouble _Lambda_T;
  VectorDouble _S_T;
  VectorDouble _Q_T;

public:
  TurboOptimizer(Id nx = 2,
                 Id ny = 2,
                 double dx = 1.,
                 double dy = 1.,
                 double x0 = 0.,
                 double y0 = 0.,
                 double scale = 1.,
                 double sill = 1.,
                 Id param = 1,
                 Id flagOne = 1);
  TurboOptimizer(const TurboOptimizer &tbo);
  TurboOptimizer& operator=(const TurboOptimizer &tbo);
  virtual ~TurboOptimizer();

  void setGrid(Id nx = 2,
               Id ny = 2,
               double dx = 1.,
               double dy = 1.,
               double x0 = 0.,
               double y0 = 0.);
  void setModelByRange(double range = 1., double sill = 1., Id param = 1);
  void setModelByScale(double scale = 1., double sill = 1., Id param = 1);
  void setEnviron(Id flagOne = 1);
  void run(bool verbose = false);
  VectorDouble getBlin() const;
  VectorDouble getTildeC() const;
  VectorDouble getLambda() const;
  TripletND getS() const;
  TripletND getQ() const;
  TripletND interpolate(const VectorDouble& x,
                        const VectorDouble& y) const;

  VectorInt interpolate_rows(const VectorDouble& x,
                             const VectorDouble& y) const
  {
    return interpolate(x, y).rows;
  }
  VectorInt interpolate_cols(const VectorDouble& x,
                             const VectorDouble& y) const
  {
    return interpolate(x, y).cols;
  }
  VectorDouble interpolate_values(const VectorDouble& x,
                                  const VectorDouble& y) const
  {
    return interpolate(x, y).values;
  }

  VectorInt getQ_rows() const
  {
    return getQ().rows;
  }
  VectorInt getQ_cols() const
  {
    return getQ().cols;
  }
  VectorDouble getQ_values() const
  {
    return getQ().values;
  }

  void printClass() const;
  void printMeshes() const;
  void printS(Id nper_batch = 5,
              Id row_begin = 0,
              Id row_end = 0,
              Id col_begin = 0,
              Id col_end = 0) const;
  void printQ(Id nper_batch = 5,
              Id row_begin = 0,
              Id row_end = 0,
              Id col_begin = 0,
              Id col_end = 0) const;

  Id getHalf() const
  {
    return _half;
  }
  Id getCenter() const
  {
    return _center;
  }
  Id getNxred() const
  {
    return _nxred;
  }
  Id getPoncif() const
  {
    return _poncif;
  }
  const VectorDouble& getTildeCT() const
  {
    return _TildeC_T;
  }
  const VectorDouble& getLambdaT() const
  {
    return _Lambda_T;
  }
  const VectorDouble& getQT() const
  {
    return _Q_T;
  }
  const VectorDouble& getST() const
  {
    return _S_T;
  }

private:
  Id _getNMeshes() const
  {
    return ((_nx - 1) * (_ny - 1) * TO_npercell);
  }
  Id _getNVertices() const
  {
    return _nx * _ny;
  }
  Id _getNVertices_red() const
  {
    return _nxred * _nxred;
  }
  double _getMeshSize() const
  {
    return (_dx * _dy / TO_npercell);
  }
  Id _getVertex(Id imesh, Id rank) const;
  double _getCoor(Id node, Id idim0) const;
  double _getCoorByMesh(Id imesh, Id rank, Id idim0) const;
  void _fromMeshToIndex(Id imesh, Id *node, Id *icas) const;
  void _rankToIndice(Id rank, VectorInt& indice, bool minusOne) const;
  static Id _MSS(Id icas, Id icorn, Id idim0);
  Id _indiceToRank(VectorInt& indice, bool flag_complete = true) const;
  void _loadHH(VectorDouble& hh) const;
  double _rangeToScale(double range) const;
  Id _coordinateToIndice(double x, double y, VectorInt& indice) const;
  double _indiceToCoordinate(Id idim0, const VectorInt& indice) const;
  static void _printVector(const std::string& title,
                           VectorDouble& uu,
                           I32 width = 10,
                           I32 ndec  = 3);
  static void _printMatrix(const std::string& title,
                           Id nrow,
                           Id ncol,
                           VectorDouble& uu,
                           Id nper_batch,
                           Id row_shift = 0,
                           Id col_shift = 0,
                           I32 width    = 10,
                           I32 ndec     = 6);
  static void _invert_3x3(VectorDouble& uu, VectorDouble& vv, double tol = 1.e-6);
  static void _prodMatrix(Id size,
                          const VectorDouble& aa,
                          const VectorDouble& bb,
                          VectorDouble& cc);
  static void _prodMatrixVector(Id size,
                                const VectorDouble& aa,
                                const VectorDouble& bb,
                                VectorDouble& cc);

  void _updateMargin(Id idim0, VectorInt& indice) const;
  void _getRankInTemplate(VectorInt& indice1,
                          VectorInt& indice2) const;
  Id _determineInternalGrid(bool verbose);
  VectorDouble _buildTildeC() const;
  VectorDouble _buildLambda(const VectorDouble& TildeC) const;
  VectorDouble _buildS(const VectorDouble& TildeC) const;
  VectorDouble _buildBlin() const;
  VectorDouble _buildQ(const VectorDouble& ss,
                       const VectorDouble& blin,
                       const VectorDouble& lambda) const;
  VectorDouble _getVectorFromTemplate(const VectorDouble& vecin) const;
  TripletND _getMatrixFromTemplate(const VectorDouble& matin,
                                 Id nperline) const;
  Id _addWeights(Id icas,
                  double x,
                  double y,
                  const VectorInt& indg0,
                  VectorInt& indices,
                  VectorDouble& lambda) const;
  VectorDouble _expandTripletToMatrix(Id row_begin,
                                             Id row_end,
                                             Id col_begin,
                                             Id col_end,
                                             const TripletND& triplet) const;
};
} // namespace gstlrn

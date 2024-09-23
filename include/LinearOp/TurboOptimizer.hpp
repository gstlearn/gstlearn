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
  int _nx;
  int _ny;
  double _dx;
  double _dy;
  double _x0;
  double _y0;
  double _scale;
  double _sill;
  int _param;
  int _poncif;
  int _center;
  int _nxred;
  int _half;
  int _flagOne;
  VectorDouble _Blin;
  VectorDouble _TildeC_T;
  VectorDouble _Lambda_T;
  VectorDouble _S_T;
  VectorDouble _Q_T;

public:
  TurboOptimizer(int nx = 2,
                 int ny = 2,
                 double dx = 1.,
                 double dy = 1.,
                 double x0 = 0.,
                 double y0 = 0.,
                 double scale = 1.,
                 double sill = 1.,
                 int param = 1,
                 int flagOne = 1);
  TurboOptimizer(const TurboOptimizer &tbo);
  TurboOptimizer& operator=(const TurboOptimizer &tbo);
  virtual ~TurboOptimizer();

  void setGrid(int nx = 2,
               int ny = 2,
               double dx = 1.,
               double dy = 1.,
               double x0 = 0.,
               double y0 = 0.);
  void setModelByRange(double range = 1., double sill = 1., int param = 1);
  void setModelByScale(double scale = 1., double sill = 1., int param = 1);
  void setEnviron(int flagOne = 1);
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
  void printS(int nper_batch = 5,
              int row_begin = 0,
              int row_end = 0,
              int col_begin = 0,
              int col_end = 0) const;
  void printQ(int nper_batch = 5,
              int row_begin = 0,
              int row_end = 0,
              int col_begin = 0,
              int col_end = 0) const;

  int getHalf() const
  {
    return _half;
  }
  int getCenter() const
  {
    return _center;
  }
  int getNxred() const
  {
    return _nxred;
  }
  int getPoncif() const
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
  int _getNMeshes() const
  {
    return ((_nx - 1) * (_ny - 1) * TO_npercell);
  }
  int _getNVertices() const
  {
    return _nx * _ny;
  }
  int _getNVertices_red() const
  {
    return _nxred * _nxred;
  }
  double _getMeshSize() const
  {
    return (_dx * _dy / TO_npercell);
  }
  int _getVertex(int imesh, int rank) const;
  double _getCoor(int node, int idim0) const;
  double _getCoorByMesh(int imesh, int rank, int idim0) const;
  void _fromMeshToIndex(int imesh, int *node, int *icas) const;
  void _rankToIndice(int rank, VectorInt& indice, bool minusOne) const;
  static int _MSS(int icas, int icorn, int idim0);
  int _indiceToRank(VectorInt& indice, bool flag_complete = true) const;
  void _loadHH(VectorDouble& hh) const;
  double _rangeToScale(double range) const;
  int _coordinateToIndice(double x, double y, VectorInt& indice) const;
  double _indiceToCoordinate(int idim0, const VectorInt& indice) const;
  static void _printVector(const std::string& title,
                           VectorDouble& uu,
                           int width = 10,
                           int ndec  = 3);
  static void _printMatrix(const std::string& title,
                           int nrow,
                           int ncol,
                           VectorDouble& uu,
                           int nper_batch,
                           int row_shift = 0,
                           int col_shift = 0,
                           int width     = 10,
                           int ndec      = 6);
  static void _invert_3x3(VectorDouble& uu, VectorDouble& vv, double tol = 1.e-6);
  static void _prodMatrix(int size,
                          const VectorDouble& aa,
                          const VectorDouble& bb,
                          VectorDouble& cc);
  static void _prodMatrixVector(int size,
                                const VectorDouble& aa,
                                const VectorDouble& bb,
                                VectorDouble& cc);

  void _updateMargin(int idim0, VectorInt& indice) const;
  void _getRankInTemplate(VectorInt& indice1,
                          VectorInt& indice2) const;
  int _determineInternalGrid(bool verbose);
  VectorDouble _buildTildeC() const;
  VectorDouble _buildLambda(const VectorDouble& TildeC) const;
  VectorDouble _buildS(const VectorDouble& TildeC) const;
  VectorDouble _buildBlin() const;
  VectorDouble _buildQ(const VectorDouble& ss,
                       const VectorDouble& blin,
                       const VectorDouble& lambda) const;
  VectorDouble _getVectorFromTemplate(const VectorDouble& vecin) const;
  TripletND _getMatrixFromTemplate(const VectorDouble& matin,
                                 int nperline) const;
  int _addWeights(int icas,
                  double x,
                  double y,
                  const VectorInt& indg0,
                  VectorInt& indices,
                  VectorDouble& lambda) const;
  VectorDouble _expandTripletToMatrix(int row_begin,
                                             int row_end,
                                             int col_begin,
                                             int col_end,
                                             const TripletND& triplet) const;
};

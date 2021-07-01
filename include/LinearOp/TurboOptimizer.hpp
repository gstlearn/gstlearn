/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/******************************************************************************/
#pragma once
#include <string>
#include <vector>

#define TO_ndim 2
#define TO_ncorner 3
#define TO_npercell 2

typedef struct
{
  std::vector<int> rows;
  std::vector<int> cols;
  std::vector<double> values;
} TripletND;

/**
 * \brief Turbo Optimizer for a specific 2-D environment,
 * \brief with an isotropic Matérn Model
 */
class TurboOptimizer
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
  std::vector<double> _Blin;
  std::vector<double> _TildeC_T;
  std::vector<double> _Lambda_T;
  std::vector<double> _S_T;
  std::vector<double> _Q_T;

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
  std::vector<double> getBlin() const;
  std::vector<double> getTildeC() const;
  std::vector<double> getLambda() const;
  TripletND getS() const;
  TripletND getQ() const;
  TripletND interpolate(const std::vector<double>& x,
                      const std::vector<double>& y) const;

  std::vector<int> interpolate_rows(const std::vector<double>& x,
                                    const std::vector<double>& y) const
  {
    return interpolate(x, y).rows;
  }
  std::vector<int> interpolate_cols(const std::vector<double>& x,
                                    const std::vector<double>& y) const
  {
    return interpolate(x, y).cols;
  }
  std::vector<double> interpolate_values(const std::vector<double>& x,
                                         const std::vector<double>& y) const
  {
    return interpolate(x, y).values;
  }

  std::vector<int> getQ_rows() const
  {
    return getQ().rows;
  }
  std::vector<int> getQ_cols() const
  {
    return getQ().cols;
  }
  std::vector<double> getQ_values() const
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
  const std::vector<double>& getTildeCT() const
  {
    return _TildeC_T;
  }
  const std::vector<double>& getLambdaT() const
  {
    return _Lambda_T;
  }
  const std::vector<double>& getQT() const
  {
    return _Q_T;
  }
  const std::vector<double>& getST() const
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
  void _rankToIndice(int rank, std::vector<int>& indice, bool minusOne) const;
  int _MSS(int icas, int icorn, int idim0) const;
  int _indiceToRank(std::vector<int>& indice, bool flag_complete = true) const;
  void _loadHH(std::vector<double>& hh) const;
  double _rangeToScale(double range) const;
  int _coordinateToIndice(double x, double y, std::vector<int>& indice) const;
  double _indiceToCoordinate(int idim0, const std::vector<int> indice) const;
  void _printVector(const std::string& title,
                    std::vector<double>& uu,
                    int width = 10,
                    int ndec = 3) const;
  void _printMatrix(const std::string& title,
                    int nrow,
                    int ncol,
                    std::vector<double>& uu,
                    int nper_batch,
                    int row_shift = 0,
                    int col_shift = 0,
                    int width = 10,
                    int ndec = 6) const;
  void _invert_3x3(std::vector<double>& uu,
                   std::vector<double>& vv,
                   double tol = 1.e-6) const;
  void _prodMatrix(int size,
                   const std::vector<double>& aa,
                   const std::vector<double>& bb,
                   std::vector<double> & cc) const;
  void _prodMatVect(int size,
                    const std::vector<double>& aa,
                    const std::vector<double>& bb,
                    std::vector<double> & cc) const;

  void _updateMargin(int idim0, std::vector<int>& indice) const;
  void _getRankInTemplate(std::vector<int>& indice1,
                          std::vector<int>& indice2) const;
  int _determineInternalGrid(bool verbose);
  std::vector<double> _buildTildeC() const;
  std::vector<double> _buildLambda(const std::vector<double> TildeC) const;
  std::vector<double> _buildS(const std::vector<double>& TildeC) const;
  std::vector<double> _buildBlin() const;
  std::vector<double> _buildQ(const std::vector<double>& ss,
                              const std::vector<double>& blin,
                              const std::vector<double>& lambda) const;
  std::vector<double> _getVectorFromTemplate(const std::vector<double>& vecin) const;
  TripletND _getMatrixFromTemplate(const std::vector<double>& matin,
                                 int nperline) const;
  int _addWeights(int icas,
                  double x,
                  double y,
                  const std::vector<int>& indg0,
                  std::vector<int>& indices,
                  std::vector<double>& lambda) const;
  std::vector<double> _expandTripletToMatrix(int row_begin,
                                             int row_end,
                                             int col_begin,
                                             int col_end,
                                             const TripletND& triplet) const;
};

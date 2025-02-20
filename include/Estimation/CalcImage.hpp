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
#include "geoslib_define.h"

#include "Enum/EMorpho.hpp"

#include "Calculators/ACalcInterpolator.hpp"

class DbGrid;
class ANeigh;
class NeighImage;

class GSTLEARN_EXPORT CalcImage: public ACalcInterpolator {
public:
  CalcImage();
  CalcImage(const CalcImage& r)            = delete;
  CalcImage& operator=(const CalcImage& r) = delete;
  virtual ~CalcImage();

  void setFlagFilter(bool flagFilter) { _flagFilter = flagFilter; }
  void setFlagFFT(bool flagFFT) { _flagFFT = flagFFT; }
  void setSeed(int seed) { _seed = seed; }
  void setFlagMorpho(bool flagMorpho) { _flagMorpho = flagMorpho; }
  void setOper(const EMorpho& oper) { _oper = oper; }
  void setOption(int option) { _option = option; }
  void setRadius(const VectorInt& radius) { _radius = radius; }
  void setVmin(double vmin) { _vmin = vmin; }
  void setVmax(double vmax) { _vmax = vmax; }
  void setVerbose(bool verbose) { _verbose = verbose; }
  void setDistErode(bool distErode) { _distErode = distErode; }
  void setNvarMorpho(int nvarMorpho) { _nvarMorpho = nvarMorpho; }
  void setFlagSmooth(bool flagSmooth) { _flagSmooth = flagSmooth; }
  void setSmoothRange(double smoothRange) { _smoothRange = smoothRange; }
  void setSmoothType(int smoothType) { _smoothType = smoothType; }

private:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;

  bool _filterImage(DbGrid* dbgrid, const ModelGeneric* modelgeneric);
  static DbGrid* _buildMarpat(const NeighImage* neigh,
                              const VectorVectorInt& ranks,
                              const MatrixRectangular& wgt);
  static VectorVectorInt _getActiveRanks(const DbGrid* dblocal);

private:
  int _iattOut;

  bool _flagFilter;
  bool _flagFFT;
  int _seed;

  bool _flagMorpho;
  int _nvarMorpho;
  EMorpho _oper;
  double _vmin;
  double _vmax;
  int _option;
  VectorInt _radius;
  bool _distErode;
  bool _verbose;

  bool _flagSmooth;
  int _smoothType;
  double _smoothRange;
};

GSTLEARN_EXPORT int krimage(DbGrid* dbgrid,
                            Model* model,
                            ANeigh* neigh,
                            bool flagFFT                    = false,
                            int seed                        = 13431,
                            const NamingConvention& namconv = NamingConvention("Filtering"));
GSTLEARN_EXPORT int dbMorpho(DbGrid* dbgrid,
                             const EMorpho& oper,
                             double vmin                     = 0.,
                             double vmax                     = 1.5,
                             int option                      = 0,
                             const VectorInt& radius         = VectorInt(),
                             bool flagDistErode              = false,
                             bool verbose                    = false,
                             const NamingConvention& namconv = NamingConvention("Morpho"));
GSTLEARN_EXPORT int dbSmoother(DbGrid* dbgrid,
                               ANeigh* neigh,
                               int type                        = 1,
                               double range                    = 1.,
                               const NamingConvention& namconv = NamingConvention("Smooth"));

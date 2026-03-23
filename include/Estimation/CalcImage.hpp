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

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Enum/EMorpho.hpp"

#include "Calculators/ACalcInterpolator.hpp"

namespace gstlrn
{

  class DbGrid;
  class ANeigh;
  class NeighImage;
  class ModelCovList;
  class Model;

  class GSTLEARN_EXPORT CalcImage: public ACalcInterpolator
  {
  public:
    CalcImage();
    CalcImage(const CalcImage& r) = delete;
    CalcImage& operator=(const CalcImage& r) = delete;
    virtual ~CalcImage();

    void setFlagFilter(bool flagFilter) { _flagFilter = flagFilter; }

    void setFlagFFT(bool flagFFT) { _flagFFT = flagFFT; }

    void setSeed(Id seed) { _seed = seed; }

    void setFlagMorpho(bool flagMorpho) { _flagMorpho = flagMorpho; }

    void setOper(const EMorpho& oper) { _oper = oper; }

    void setOption(Id option) { _option = option; }

    void setRadius(const VectorInt& radius) { _radius = radius; }

    void setVmin(double vmin) { _vmin = vmin; }

    void setVmax(double vmax) { _vmax = vmax; }

    void setVerbose(bool verbose) { _verbose = verbose; }

    void setDistErode(bool distErode) { _distErode = distErode; }

    void setNvarMorpho(Id nvarMorpho) { _nvarMorpho = nvarMorpho; }

    void setFlagSmooth(bool flagSmooth) { _flagSmooth = flagSmooth; }

    void setSmoothRange(double smoothRange) { _smoothRange = smoothRange; }

    void setSmoothType(Id smoothType) { _smoothType = smoothType; }

  private:
    bool _check() override;
    bool _preprocess() override;
    bool _run() override;
    bool _postprocess() override;
    void _rollback() override;

    bool _filterImage(DbGrid* dbgrid, const ModelGeneric* model);
    static void _image_smoother(DbGrid* dbgrid,
                                const NeighImage* neigh,
                                Id type,
                                double range,
                                Id iptr0);
    static DbGrid* _buildMarpat(const NeighImage* neigh,
                                const VectorVectorInt& ranks,
                                const MatrixDense& wgt,
                                Id optionVerbose = 0);
    static VectorVectorInt _getActiveRanks(const DbGrid* dblocal);

  private:
    Id _iattOut;

    bool _flagFilter;
    bool _flagFFT;
    Id _seed;

    bool _flagMorpho;
    Id _nvarMorpho;
    EMorpho _oper;
    double _vmin;
    double _vmax;
    Id _option;
    VectorInt _radius;
    bool _distErode;
    bool _verbose;

    bool _flagSmooth;
    Id _smoothType;
    double _smoothRange;
  };

} // namespace gstlrn

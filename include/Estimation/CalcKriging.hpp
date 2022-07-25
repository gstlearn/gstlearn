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

#include "geoslib_define.h"

#include "Calculators/ACalcInterpolator.hpp"

class Db;
class DbGrid;

class GSTLEARN_EXPORT CalcKriging: public ACalcInterpolator
{
public:
  CalcKriging(bool flag_est = true, bool flag_std = true, bool flag_varZ = false);
  CalcKriging(const CalcKriging &r) = delete;
  CalcKriging& operator=(const CalcKriging &r) = delete;
  virtual ~CalcKriging();

  const EKrigOpt& getCalcul() const
  {
    return _calcul;
  }

  void setCalcul(const EKrigOpt &calcul)
  {
    _calcul = calcul;
  }

  void setMatCl(const VectorVectorDouble &matCl)
  {
    _matCL = matCl;
  }

  const VectorInt& getNdisc() const
  {
    return _ndisc;
  }

  void setNdisc(const VectorInt &ndisc)
  {
    _ndisc = ndisc;
  }

  const VectorInt& getRankColCok() const
  {
    return _rankColCok;
  }

  void setRankColCok(const VectorInt &rankColCok)
  {
    _rankColCok = rankColCok;
  }

private:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;
  int _getNVar() const override;

private:
  bool _flagEst;
  bool _flagStd;
  bool _flagVarZ;

  EKrigOpt  _calcul;
  VectorInt _ndisc;
  VectorInt _rankColCok;
  VectorVectorDouble _matCL;

  int _iptrEst;
  int _iptrStd;
  int _iptrVarZ;
};

GSTLEARN_EXPORT int kriging(Db *dbin,
                            Db *dbout,
                            Model *model,
                            ANeighParam *neighparam,
                            const EKrigOpt &calcul = EKrigOpt::PONCTUAL,
                            bool flag_est = true,
                            bool flag_std = true,
                            bool flag_varz = false,
                            VectorInt ndisc = VectorInt(),
                            VectorInt rank_colcok = VectorInt(),
                            VectorVectorDouble matCL = VectorVectorDouble(),
                            const NamingConvention& namconv = NamingConvention("Kriging"));

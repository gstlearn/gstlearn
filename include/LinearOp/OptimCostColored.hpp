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
#include "OptimCostBinary.hpp"

namespace gstlrn
{
class GSTLEARN_EXPORT OptimCostColored : public OptimCostBinary
{
public:
  OptimCostColored();
  OptimCostColored(Id nprop,
                   PrecisionOp* pmat,
                   const ProjMatrix* projdata,
                   const ProjMatrix* projseis = nullptr,
                   const VectorDouble& propseis = VectorDouble(),
                   const VectorDouble& varseis = VectorDouble());
  OptimCostColored(const OptimCostColored &m);
  OptimCostColored& operator = (const OptimCostColored &m);
	virtual ~OptimCostColored();

  void reset(Id nprop,
             PrecisionOp* pmat,
             const ProjMatrix* projdata,
             const ProjMatrix* projseis = nullptr,
             const VectorDouble& propseis = VectorDouble(),
             const VectorDouble& varseis = VectorDouble());

  VectorVectorDouble minimize(const VectorDouble& facies,
                              const VectorVectorInt& splits = VectorVectorInt(),
                              const VectorDouble& meanprops = VectorDouble(),
                              bool verbose = false,
                              Id maxiter = 100,
                              double eps = 5.e-4);

  VectorVectorInt initSplit(Id nfacies, bool verbose = false) const;
  void   printSplits(const VectorVectorInt& splits = VectorVectorInt()) const;

  void setMeanProps(const VectorDouble& meanProps) { _meanProps = meanProps; }
  void setSplits(const VectorVectorInt& splits) { _splits = splits; }

private:
  void   _getFaciesToIndic(const VectorDouble& facies,
                           const VectorInt&    split,
                           VectorDouble&       indic) const;
  double _getFaciesToProportion(const VectorInt& split) const;
  Id    _checkFacies(const VectorDouble& facies) const;
  Id    _checkSplits(const VectorVectorInt& splits);
  Id    _checkMeanProportions(const VectorDouble& meanprops);
  void   _copyMultProportions(Id level,
                              Id ip,
                              const VectorDouble& propfac,
                              VectorVectorDouble& propfacs);
                                
private:
  Id             _nprop;
  VectorVectorInt _splits;
  VectorDouble    _meanProps;
};
}
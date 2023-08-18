/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "OptimCostBinary.hpp"
#include "Basic/AStringable.hpp"

class GSTLEARN_EXPORT OptimCostColored : public OptimCostBinary
{
public:
  OptimCostColored();
  OptimCostColored(int nprop,
                   PrecisionOp* pmat,
                   const ProjMatrix* projdata,
                   const ProjMatrix* projseis = nullptr,
                   const VectorDouble& propseis = VectorDouble(),
                   const VectorDouble& varseis = VectorDouble());
  OptimCostColored(const OptimCostColored &m);
  OptimCostColored& operator = (const OptimCostColored &m);
	virtual ~OptimCostColored();

  void reset(int nprop,
             PrecisionOp* pmat,
             const ProjMatrix* projdata,
             const ProjMatrix* projseis = nullptr,
             const VectorDouble& propseis = VectorDouble(),
             const VectorDouble& varseis = VectorDouble());

  VectorVectorDouble minimize(const VectorDouble& facies,
                              const VectorVectorInt& splits = VectorVectorInt(),
                              const VectorDouble& meanprops = VectorDouble(),
                              bool verbose = false,
                              int maxiter = 100,
                              double eps = 5.e-4);

  VectorVectorInt initSplit(int nfacies, bool verbose = false) const;
  void   printSplits(const VectorVectorInt& splits = VectorVectorInt()) const;

  void setMeanProps(const VectorDouble& meanProps) { _meanProps = meanProps; }
  void setSplits(const VectorVectorInt& splits) { _splits = splits; }

private:
  void   _getFaciesToIndic(const VectorDouble& facies,
                           const VectorInt&    split,
                           VectorDouble&       indic) const;
  double _getFaciesToProportion(const VectorInt& split) const;
  int    _checkFacies(const VectorDouble& facies) const;
  int    _checkSplits(const VectorVectorInt& splits);
  int    _checkMeanProportions(const VectorDouble& meanprops);
  void   _copyMultProportions(int level,
                              int ip,
                              const VectorDouble& propfac,
                              VectorVectorDouble& propfacs);
                                
private:
  int             _nprop;
  VectorVectorInt _splits;
  VectorDouble    _meanProps;
};

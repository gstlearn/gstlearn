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

  void    init(int                 nprop,
               PrecisionOp*  	     pmat,
               const ProjMatrix*   projdata,
               const ProjMatrix*   projseis = nullptr,
               const VectorDouble&  propseis = VectorDouble(),
               const VectorDouble&  varseis  = VectorDouble());

  VectorVectorDouble minimize(const VectorDouble& facies,
                              const VectorVectorInt& splits = VectorVectorInt(),
                              const VectorDouble& meanprops = VectorDouble(),
                              bool verbose = false,
                              int maxiter = 100,
                              double eps = 5.e-4);

  VectorVectorInt createSplit(int nfacies, bool verbose = false) const;
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

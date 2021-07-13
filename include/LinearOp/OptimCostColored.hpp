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

#include "OptimCostBinary.hpp"

class OptimCostColored : public OptimCostBinary {

public:
	OptimCostColored();
	virtual ~OptimCostColored();

  void    init(int                 nprop,
               PrecisionOp*  	     pmat,
               const ProjMatrix*   projdata,
               const ProjMatrix*   projseis = nullptr,
               const VectorDouble& propseis = VectorDouble(),
               const VectorDouble& varseis  = VectorDouble());

  int     minimize(VectorDouble& facies,
                   std::vector<VectorDouble>& propfacs,
                   std::vector<VectorInt>&    splits,
                   VectorDouble& meanprops,
                   bool          verbose = false,
                   int           maxiter = 100,
                   double        eps = 5.e-4);

private:

  void   _getFaciesToIndic(const VectorDouble& facies,
                           const VectorInt&    split,
                           VectorDouble&       indic) const;
  double _getFaciesToProportion(const VectorInt& split) const;
  int    _checkFacies(const VectorDouble& facies) const;
  int    _checkSplits(const std::vector<VectorInt>& splits);
  int    _checkMeanProportions(const VectorDouble& meanprops);
  void   _copyMultProportions(int level,
                              int ip,
                              const VectorDouble& propfac,
                              std::vector<VectorDouble>& propfacs);
  void   _printSplits() const;
                                
private:
  int                    _nprop;
  std::vector<VectorInt> _splits;
  VectorDouble           _meanProps;
};

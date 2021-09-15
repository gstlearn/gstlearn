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

#include "Gibbs/AGibbs.hpp"
#include "Basic/Vector.hpp"

class Db;
class Model;

class GibbsStandard : public AGibbs
{
public:
  GibbsStandard();
  GibbsStandard(const GibbsStandard &r);
  GibbsStandard& operator=(const GibbsStandard &r);
  virtual ~GibbsStandard();

  int calculInitialize(Db *dbin,
                       Model *model,
                       int isimu,
                       int igrf,
                       int ipgs,
                       bool verbose);
  int calculIteration(Db *dbin,
                      Model *model,
                      int isimu,
                      int ipgs,
                      int igrf,
                      int verbose);
  int calculatePropagation(Db *db,
                           Model *model,
                           int isimu,
                           bool verbose);
  int covmatAlloc(Db *dbin, Model *model, bool verbose);

private:
  VectorDouble _covmat;
};

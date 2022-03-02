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
#include "Gibbs/GibbsMulti.hpp"
#include "Basic/Vector.hpp"

class Db;
class Model;

class GSTLEARN_EXPORT GibbsUMulti : public GibbsMulti
{
public:
  GibbsUMulti();
  GibbsUMulti(Db* db, Model* model);
  GibbsUMulti(const GibbsUMulti &r);
  GibbsUMulti& operator=(const GibbsUMulti &r);
  virtual ~GibbsUMulti();

  void update(VectorVectorDouble& y,
              int isimu,
              int ipgs,
              int iter) override;
  int covmatAlloc(bool verbose, bool verboseTimer = false) override;

private:
  VectorDouble _covmat;
};

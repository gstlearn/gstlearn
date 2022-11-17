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

#include "Gibbs/AGibbs.hpp"

class Db;
class Model;

class GSTLEARN_EXPORT GibbsFactory
{
public:
  GibbsFactory();
  virtual ~GibbsFactory();

  static AGibbs *createGibbs(Db* db,
                             Model* model,
                             bool flagMoving);
  static AGibbs *createGibbs(Db* db,
                             std::vector<Model *> models,
                             double rho,
                             bool flag_propagation);
};

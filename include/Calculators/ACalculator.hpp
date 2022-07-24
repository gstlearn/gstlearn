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

class GSTLEARN_EXPORT ACalculator
{
public:
  ACalculator();
  ACalculator(const ACalculator &r) = delete;
  ACalculator& operator=(const ACalculator &r) = delete;
  virtual ~ACalculator();

  bool run();

protected:
  virtual bool _run() = 0;

  virtual bool _check() const { return true; }
  virtual bool _preprocess()  { return true; }
  virtual bool _postprocess() { return true; }
  virtual void _rollback() { }
};

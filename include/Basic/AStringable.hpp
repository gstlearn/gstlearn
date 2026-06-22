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

#include "geoslib_define.h"

#include "Basic/AStringFormat.hpp"

// Even if not used directly here, this include is used in all subsequent functions
// composing the stream (as required in toString() method)
#include <sstream>

namespace gstlrn
{
  class AMatrix;

  class GSTLEARN_EXPORT AStringable
  {
  public:
    AStringable();
    AStringable(const AStringable& r);
    AStringable& operator=(const AStringable& r);
    virtual ~AStringable();

    virtual String toString(const AStringFormat* strfmt = nullptr) const;

    virtual void display(const AStringFormat* strfmt = nullptr) const final;
#ifndef SWIG // TODO : overload not available in customized SWIG 4.2.3 and more
    virtual void display(Id level) const final;
#endif

    void printConcreteClassName() const;
  };

} // namespace gstlrn

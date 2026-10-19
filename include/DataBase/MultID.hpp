/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2025) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Basic/VectorNumT.hpp"
#include "Enum/ERole.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
  class GSTLEARN_EXPORT MultID
  {
  public:
    MultID(const String& name);
    MultID(const VectorString& names);
    MultID(const VectorInt& icols);
    MultID(const ERole& role);
    MultID(Id icol);

#ifndef SWIG
    template<int N>
    MultID(const char (&name)[N])
      : MultID(String(name))
    {
    }
#endif

    const VectorString& getNames() const { return _names; }

    const VectorInt& getICols() const { return _icols; }

    const ERole& getRole() const { return _role; }

    String getDescr() const;

    void setNames(const VectorString& names) { _names = names; }

    void setICols(const VectorInt& icols) { _icols = icols; }

    void setERole(const ERole& role) { _role = role; }

    bool hasRoleDefined() const;

    /**
     * Factory helpers mainly intended for language bindings (SWIG).
     */
#ifndef SWIG
    static MultID create(const String& name);
    static MultID create(const VectorString& names);
    static MultID create(const VectorInt& icols);
    static MultID create(const ERole& role);
    static MultID create(const MultID& MultID);
#endif

  private:
    VectorInt _icols{};
    VectorString _names{};
    ERole _role{ERole::UNDEFINED};
  };

} // namespace gstlrn

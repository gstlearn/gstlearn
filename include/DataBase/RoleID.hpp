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

#include "gstlearn_export.hpp"

#include "Enum/ERole.hpp"
#include "geoslib_define.h"

#include <cctype>
#include <optional>

namespace gstlrn
{
  /**
   * @brief This class describes the Role that a Column carries in a DbData.
   * It is characterized by:
   * - its 'Role' (from the Enum ERole.hpp)
   * - its 'index' within this Role (warning: this is 0 based)
   *
   * As an example, one can think of a Column where:
   * - the 'Role' is set to ERole::X
   * - the 'index' is set to 2
   * This stands for a Column containing the third coordinate os samples
   *
   * Important remark: 'index' is a à-based number. However, the printout of this Role1D
   * is performed using a 1-based number (i.e. index+1) to be more user friendly.
   * The previous example will be printed as "x3".
   */
  class GSTLEARN_EXPORT RoleID
  {
  public:
    RoleID(const ERole& role = ERole::UNDEFINED, Id index = 0);

    const ERole& getRole() const { return _role; }

    Id getIndex() const { return _index; }

    void setRole(const ERole& role) { _role = role; }

    void setIndex(Id index) { _index = index; }

    bool isEqual(const RoleID& RoleID, bool checkIndex = true) const;

    bool isUnique() const;

    String getName() const;

    static std::optional<RoleID> createFromName(String string);

    void removeRole();

    bool isDefined() const { return _role != ERole::UNDEFINED; }

  private:
    ERole _role{ERole::UNDEFINED};
    Id _index{0};
  };

} // namespace gstlrn

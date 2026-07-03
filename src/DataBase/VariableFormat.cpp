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
#include "DataBase/VariableFormat.hpp"
#include "Enum/EVariableAlign.hpp"
#include "Enum/EVariableNumeric.hpp"
#include "Enum/EVariableTruncate.hpp"
#include "geoslib_define.h"
#include <iomanip>

namespace gstlrn
{
  VariableFormat::VariableFormat(
    Id layout,
    Id precision,
    const EVariableAlign& align,
    const EVariableTruncate& truncate,
    const EVariableNumeric& mode)
    : _layout(layout)
    , _precision(precision)
    , _align(align)
    , _truncate(truncate)
    , _mode(mode)
  {
  }

  std::string VariableFormat::encode() const
  {
    std::string s = "w";

    if (_align == EVariableAlign::LEFT)
      s += "<";
    else
      s += ">";

    s += std::to_string(_layout);

    if (_truncate == EVariableTruncate::LEFT)
      s += "l";
    else
      s += "r";

    if (_mode != EVariableNumeric::NONE)
    {
      s += ":";

      switch (_mode.getValue())
      {
        case EVariableNumeric::E_FIXED: s += "f"; break;
        case EVariableNumeric::E_SCIENTIFIC: s += "e"; break;
        case EVariableNumeric::E_AUTO: s += "g"; break;
        default: break;
      }

      if (_precision >= 0) s += std::to_string(_precision);
    }

    return s;
  }

  std::string VariableFormat::apply(UChar value, bool isUndefined) const
  {
    DECLARE_UNUSED(isUndefined);
    if (_layout <= 0) return "";

    std::string out = (value ? "True" : "False");
    return out;
  }

  void VariableFormat::applyNumericFormat(std::ostringstream& os) const
  {
    if (_mode != EVariableNumeric::NONE)
    {
      switch (_mode.getValue())
      {
        case EVariableNumeric::E_FIXED:
          os << std::fixed;
          if (_precision >= 0) os << std::setprecision(_precision);
          break;

        case EVariableNumeric::E_SCIENTIFIC: os << std::scientific; break;

        case EVariableNumeric::E_AUTO:

        default: break;
      }
    }
  }

} // namespace gstlrn

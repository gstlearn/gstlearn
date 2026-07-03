#pragma once

#include "Enum/EVariableAlign.hpp"
#include "Enum/EVariableNumeric.hpp"
#include "Enum/EVariableTruncate.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
  class GSTLEARN_EXPORT VariableFormat
  {
  public:
    VariableFormat(
      Id layout = 10,
      Id precision = 3,
      const EVariableAlign& align = EVariableAlign::RIGHT,
      const EVariableTruncate& truncate = EVariableTruncate::RIGHT,
      const EVariableNumeric& mode = EVariableNumeric::AUTO);
    VariableFormat(const VariableFormat& r) = default;
    VariableFormat& operator=(const VariableFormat& r) = default;
    ~VariableFormat() = default;

    Id getLayout() const { return _layout; }

    Id getPrecision() const { return _precision; }

    EVariableAlign getAlign() const { return _align; }

    EVariableTruncate getTruncate() const { return _truncate; }

    EVariableNumeric getMode() const { return _mode; }

    std::string encode() const;

    void applyNumericFormat(std::ostringstream& os) const;

    template<class T>
    std::string apply(const T& value, bool isUndefined = true) const
    {
      if (_layout <= 0) return "";

      String out;
      if (isUndefined)
      {
        out = "NA";
      }
      else
      {
        std::ostringstream tmp;

        if constexpr (std::is_arithmetic_v<T>) applyNumericFormat(tmp);

        tmp << value;

        out = tmp.str();
      }

      Id len = static_cast<Id>(out.size());

      if (len > _layout) out.resize(_layout);

      if (len < _layout)
      {
        size_t nfill = _layout - len;

        if (_align == EVariableAlign::LEFT)
          out.append(nfill, ' ');
        else
          out.insert(0, nfill, ' ');
      }

      return out;
    }

    std::string apply(UChar value, bool isUndefined = true) const;

  private:
    Id _layout;
    Id _precision;
    EVariableAlign _align;
    EVariableTruncate _truncate;
    EVariableNumeric _mode;
  };

} // namespace gstlrn

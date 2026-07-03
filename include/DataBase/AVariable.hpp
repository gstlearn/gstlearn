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

#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/Undefined.hpp"
#include "Basic/VectorNumT.hpp"
#include "DataBase/VariableFormat.hpp"

#include <memory>

namespace gstlrn
{
  class VariableDouble;
  class VariableBool;
  class VariableInt;
  class VariableString;

  class GSTLEARN_EXPORT AVariable: public AStringable, public ICloneable
  {
  public:
    AVariable() = default;

    AVariable(
      const String& name,
      size_t ncols = 1,
      const String& unit = "",
      bool naallowed = true,
      const VariableFormat& format = VariableFormat())
      : AStringable()
      , _name(name)
      , _nCols(ncols)
      , _unit(unit)
      , _NAAllowed(naallowed)
      , _format(format)
    {
    }

    ~AVariable() = default;

    /**
     * @brief Returns the count of samples in the variable.
     *
     * @return size_t
     */
    virtual size_t getNSamples() const = 0;
    /**
     * @brief Returns the count of values in the variable (for all Cols) which are not undefined.
     *
     * @return size_t
     */
    virtual size_t getNAllDefined() const = 0;
    /**
     * @brief Get the set of values in the variable for a given column as a VectorDouble object
     *
     * @return VectorDouble
     */
    virtual VectorDouble getValues(Id icol = 0) const = 0;
    /**
     * @brief Get the All Values object (for all columns) as a VectorDouble object
     *
     * @return VectorDouble
     */
    virtual String getType() const = 0;

    static std::unique_ptr<AVariable>
      create(const String& name, const String& type = "double");

    void setName(const String& name) { _name = name; }

    void setUnit(const String& unit) { _unit = unit; }

    void setNCols(size_t ncols) { _nCols = ncols; }

    void setNAAllowed(bool allowed) { _NAAllowed = allowed; }

    void setFormat(const VariableFormat& format) { _format = format; }

    const String& getName() const { return _name; }

    const String& getUnit() const { return _unit; }

    size_t getNCols() const { return _nCols; }

    bool isNAAllowed() const { return _NAAllowed; }

    /**
     * @brief Returns the count of values in the variable (e.g.: number of samples * number of columns)
     *
     * @return size_t
     */
    size_t getNAllValues() const { return getNSamples() * getNCols(); };

  protected:
    static bool _getFlagVariableCheck();
    static void _setFlagVariableCheck(bool flag);

  protected:
    String _name{UNDEF_STRING};
    size_t _nCols{1};
    String _unit{" "};
    bool _NAAllowed{true};
    VariableFormat _format;
  };

} // namespace gstlrn

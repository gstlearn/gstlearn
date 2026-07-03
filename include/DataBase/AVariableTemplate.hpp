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

#include "Basic/Message.hpp"
#include "DataBase/AVariable.hpp"
#include "gstlearn_export.hpp"

#include <iostream>
#include <vector>

namespace gstlrn
{

  template<class T>
  class GSTLEARN_EXPORT AVariableTemplate: public AVariable
  {
  public:
    AVariableTemplate()
      : AVariable()
      , _values()
    {
    }

    AVariableTemplate(
      const String& name,
      size_t ncols = 1,
      const String& unit = "",
      bool naallowed = true,
      const VariableFormat& format = VariableFormat())
      : AVariable(name, ncols, unit, naallowed, format)
      , _values()
    {
    }

    AVariableTemplate(const String& name, const std::vector<T>& values)
      : AVariable(name, 1, "")
      , _values(values)
    {
    }

    size_t getNSamples() const override { return _values.size() / getNCols(); }

    size_t getNDefined(size_t icol) const
    {
      const size_t nsamples = getNSamples();
      const size_t offset = icol * nsamples;

      size_t ndef = 0;
      for (size_t isample = 0; isample < nsamples; ++isample)
      {
        if (!isUndefined(offset + isample)) ++ndef;
      }
      return ndef;
    }

    size_t getNAllDefined() const override
    {
      size_t ndef = 0;
      for (size_t i = 0; i < getNAllValues(); ++i)
      {
        if (!isUndefined(i)) ++ndef;
      }
      return ndef;
    }

    virtual bool isUndefined(size_t i) const = 0;

    /**
     * @brief Resize the variable to contain n samples, filling it with val.
     *
     * @param nsample Total number of samples to resize the variable to
     * @param val Default value to fill the variable with
     *
     * @remarks The resize method will overwrite the current samples of the variable,
     *          and fill it with the specified value. It will not modify the number of columns.
     *          If the new size is smaller than the current size, the variable will be truncated.
     */
    virtual void resize(const size_t nsample, const T& val)
    {
      _values.resize(nsample * getNCols(), val);
    }

    void fill(const T& val) { std::fill(_values.begin(), _values.end(), val); }

    virtual void
      setValue(const T& val, const size_t isample, const size_t icol = 0)
    {

      if (_getFlagVariableCheck())
      {
        if (!checkArg("AVariableTemplate::setValue", isample, getNSamples()))
          return;
        if (!checkArg("AVariableTemplate::setValue", icol, getNCols())) return;
      }

      const size_t nsamples = getNSamples();
      const size_t offset = icol * nsamples;
      size_t index = offset + isample;
      if (isNAAllowed() && isNA(val))
      {
        messerr("Cannot set value to NA, NA values are not allowed");
        return;
      }
      _values.at(index) = val;
    }

    virtual void setValues(const std::vector<T>& vals, const size_t icol = 0)
    {
      if (_getFlagVariableCheck())
      {
        if (!checkArg("AVariableTemplate::setValues", icol, getNCols())) return;
      }
      if (!isNAAllowed())
      {
        // Check that no NA value is provided in the input vector
        for (size_t isample = 0; isample < vals.size(); ++isample)
        {
          if (isNA(vals[isample]))
          {
            messerr(
              "Cannot set value for item %d to NA, NA values are not allowed",
              isample);
            return;
          }
        }
      }
      if (getNSamples() <= 0)
      {
        resize(vals.size(), getNA<T>());
      }

      const size_t nsamples = getNSamples();
      const size_t offset = icol * nsamples;
      for (size_t isample = 0; isample < nsamples; ++isample)
      {
        _values.at(offset + isample) = vals[isample];
      }
    }

    T getValueAsType(const size_t isample, const size_t icol = 0) const
    {
      if (_getFlagVariableCheck())
      {
        if (!checkArg(
              "AVariableTemplate::getValueAsType", isample, getNSamples()))
          return getNA<T>();
        if (!checkArg("AVariableTemplate::getValueAsType", icol, getNCols()))
          return getNA<T>();
      }
      const size_t nsamples = getNSamples();
      const size_t offset = icol * nsamples;
      return _values.at(offset + isample);
    }

    virtual const std::vector<T>& getAllValuesAsType() const { return _values; }

    virtual void display_old() const
    {
      std::cout << "Variable " << _name << " :\n";
      size_t ncols = getNCols();
      size_t i = 0;
      for (size_t icol = 0; icol < ncols; ++icol)
      {
        for (size_t isample = 0; isample < getNSamples(); ++isample, i++)
        {
          if (ncols > 1)
            std::cout << "    value[" << icol << "][" << isample << "] = ";
          else
            std::cout << "    value[" << isample << "] = ";
          std::cout << _format.apply(_values[i], isUndefined(i)) << " "
                    << getUnit() << std::endl;
        }
      }
    }

  protected:
    std::vector<T> _values;
  };

} // namespace gstlrn

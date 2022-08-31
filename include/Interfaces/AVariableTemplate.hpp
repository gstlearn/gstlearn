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
#include "Interfaces/interface_d.hpp"
#include "Interfaces/AVariable.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <sstream>

template<class T>
class GSTLEARN_EXPORT AVariableTemplate : public AVariable
{
  public:

    AVariableTemplate() : AVariable(), _values() {}
    AVariableTemplate(const String& name) : AVariable(name), _values() {}
    AVariableTemplate(const String& name, std::vector<T> values) : AVariable(name), _values(values) {}
    AVariableTemplate(const AVariableTemplate& ref) : AVariable(ref), _values(ref._values) {}
    AVariableTemplate& operator=(const AVariableTemplate& ref)
    {
      if (this != &ref)
      {
        _name = ref._name;
        _values = ref._values;
        
      }
      return(*this);
    }

    virtual ~AVariableTemplate()
    {
    }

    unsigned int getNValues() const override { return (unsigned int) _values.size(); }

    //bool isUndefined(int i) const { return getUndef<T>() == getValueAsType(i); }
    virtual void resize(int n, const T& val) { _values.resize(n, val); }

    virtual void setValue(int i, const T& val) { _checkIdx(i); _values[i] = val; }
    
    virtual void setValues(const std::vector<T>& values) { _values = values; }

    T getValueAsType(int i) const { _checkIdx(i); return _values[i]; }

    virtual const std::vector<T>& getValuesAsType() const { return _values; }

    virtual void display_old() const
    {
      VectorDouble vals = getValues();
      std::cout << "Variable " << _name << " :" << std::endl;

      for (int i=0, n=(int) _values.size(); i<n; i++)
      {
        std::stringstream sstr;
        sstr << "  var[" << i << "]=" << _values[i] << " (" << vals[i] << ")";
        std::cout << sstr.str() << std::endl;
      }
    }

  private:
    void _checkIdx(int i) const
    {
      if (i < 0 || i >= (int)_values.size())
        throw("Index out of range");
    }

  protected:
    std::vector<T> _values;
};


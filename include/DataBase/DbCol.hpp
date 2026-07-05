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
#include "Basic/VectorT.hpp"
#include "gstlearn_export.hpp"

#include <functional>
#include <iostream>
#include <optional>
#include <type_traits>
#include <variant>

namespace gstlrn
{

  /**
   * @brief Generic implementation of a Column inside a Db
   */
  class GSTLEARN_EXPORT ADbCol
  {
  public:
    ADbCol(String&& name)
      : _name{std::move(name)}
      , _data{VectorDouble{}}
    {
    }

    template<typename VectorType>
    ADbCol(String&& name, VectorType&& array)
      : _name{std::move(name)}
      , _data{array}
    {
    }

    template<typename VectorType>
    ADbCol(const String& name, VectorType&& array)
      : _name{name}
      , _data{array}
    {
    }

    const String& getName() const { return this->_name; }

    template<typename T, typename VectorType>
    static constexpr bool isConsistent()
    {
      return std::is_convertible_v<T, typename VectorType::value_type>;
    }

    template<typename T>
    std::optional<T> getValue(const Id i) const
    {
      std::optional<T> v{};
      std::visit(
        [i, &v](auto&& arg)
        {
          using VectorType = std::decay_t<decltype(arg)>;
          if constexpr (isConsistent<T, VectorType>())
          {
            v = arg[i];
          }
        },
        this->_data);
      return v;
    }

    template<typename T>
    bool setValue(const Id i, const T& v)
    {
      bool res = false;
      std::visit(
        [i, &v, &res](auto&& arg)
        {
          using VectorType = std::decay_t<decltype(arg)>;
          if constexpr (isConsistent<T, VectorType>())
          {
            arg[i] = v;
            res = true;
          }
        },
        this->_data);
      return res;
    }

    // template<typename VectorType>
    // std::optional<std::reference_wrapper<const VectorType>> getArray() const
    // {
    //   if (auto* p = std::get_if<VectorType>(&_data)) return *p;
    //   return {};
    // }

#ifndef SWIG
    template<typename VectorType>
    const VectorType* getArray() const
    {
      return std::get_if<VectorType>(&_data);
    }
#endif

  private:
    String _name;
    std::variant<
      VectorDouble,
      VectorFloat,
      VectorInt,
      VectorUChar,
      VectorString,
      VectorBool>
      _data;
  };

  /**
   * @brief Similar to vtkFieldData
   */
  class DbData
  {
  public:
    void mytest() const
    {
      std::cout << "Activation du test" << std::endl;

      std::cout << "DbData has " << this->_cols.size() << " columns\n";
      for (const auto& col: this->_cols)
      {
        std::cout << "Column name: " << col.getName() << "\n";
      }
    }

#ifndef SWIG
    template<typename VectorType>
    void addArray(VectorType&& array, String&& name)
    {
      this->_cols.emplace_back(
        std::move(name), std::forward<VectorType>(array));
    }

    template<typename VectorType>
    void addArray(VectorType&& array, const String& name)
    {
      this->_cols.emplace_back(name, std::forward<VectorType>(array));
    }
#endif

    void removeArray(const String& name)
    {
      const auto col = this->getArray(name);
      const auto id = col ? col->second : -1;
      this->removeArray(id);
    }

    void removeArray(const Id icol)
    {
      this->_cols.erase(this->_cols.begin() + icol);
    }

    std::optional<std::reference_wrapper<ADbCol>> getArray(const Id icol)
    {
      if (icol < 0 || icol >= static_cast<Id>(this->_cols.size()))
      {
        return {};
      }
      return {this->_cols[icol]};
    }

#ifndef SWIG
    std::optional<std::reference_wrapper<const ADbCol>>
      getArray(const Id icol) const
    {
      if (icol < 0 || icol >= static_cast<Id>(this->_cols.size()))
      {
        return {};
      }
      return {this->_cols[icol]};
    }
#endif

    std::optional<std::pair<std::reference_wrapper<ADbCol>, Id>>
      getArray(const String& name)
    {
      Id icol{};
      for (auto& col: this->_cols)
      {
        if (col.getName() == name)
        {
          return {{col, icol}};
        }
        icol++;
      }

      return {};
    }

#ifndef SWIG
    std::optional<String> getArrayName(const Id icol)
    {
      auto arr = this->getArray(icol);
      return arr ? std::optional<String>{arr->get().getName()} : std::nullopt;
    }
#endif

    bool hasArray(const String& name)
    {
      auto arr = this->getArray(name);
      return arr.has_value();
    }

#ifndef SWIG
    template<typename T>
    std::optional<T> getValue(const Id iech, const Id icol) const
    {
      const auto array = this->getArray(icol);
      if (!array) return {};
      const auto val = array->get().getValue<T>(iech);
      return val;
    }
#endif

#ifndef SWIG
    template<typename T>
    bool setValue(const Id iech, const Id icol, const T value)
    {
      const auto array = this->getArray(icol);
      if (!array) return false;
      return array->get().setValue<T>(iech, value);
    }
#endif

  private:
    std::vector<ADbCol> _cols;
  };

} // namespace gstlrn

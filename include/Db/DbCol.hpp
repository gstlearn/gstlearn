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

#include <algorithm>
#include <functional>
#include <map>
#include <optional>
#include <string_view>
#include <type_traits>
#include <variant>

namespace gstlrn
{

  class GSTLEARN_EXPORT Dictionary
  {
  public:
    using Category = std::pair<Id, std::string_view>;

    Dictionary() = default;

    Dictionary(std::map<Id, String>&& data)
      : _data{std::move(data)}
    {
    }

    bool addCategory(const Id key, const String& val)
    {
      for (const auto& p: this->_data)
      {
        if (p.first == key || p.second == val)
        {
          return false;
        }
      }
      this->_data[key] = val;
      return true;
    }

    bool hasCategory(const Category& cat) const
    {
      const auto it = this->_data.find(cat.first);
      if (it == this->_data.end()) return false;
      if (it->second != cat.second) return false;
      return true;
    }

    std::optional<Category> operator[](const Id key) const
    {
      const auto it = this->_data.find(key);
      if (it == this->_data.end()) return {};
      return *it;
    }

  private:
    std::map<Id, String> _data;
  };

  class GSTLEARN_EXPORT VectorCategory
  {
  public:
    using Category = Dictionary::Category;
    using value_type = Category;

    VectorCategory(const size_t n, const Dictionary& dict)
      : _data(n)
      , _dict{dict}
    {
    }

    void resize(const size_t count) { this->_data.resize(count); }

    size_t size() const { return this->_data.size(); }

    std::optional<Category> getCategory(const size_t pos) const
    {
      if (pos >= this->_data.size()) return {};
      return this->_data[pos];
    }

    bool setCategory(const size_t pos, const Category& cat)
    {
      if (pos >= this->_data.size()) return false;
      if (!this->_dict.hasCategory(cat)) return false;
      this->_data[pos] = cat;
      return true;
    }

    std::optional<Category>& operator[](const size_t o)
    {
      return this->_data[o];
    }

    const std::optional<Category>& operator[](const size_t o) const
    {
      return this->_data[o];
    }

  private:
    std::vector<std::optional<Category>> _data;
    Dictionary _dict;
  };

  template<typename VectorType>
  class Array2D
  {
  public:
    using vector_type = VectorType;
    using value_type = typename VectorType::value_type;

    Array2D() = default;

    Array2D(VectorType&& vector)
      : _inner{static_cast<Id>(vector.size())}
      , _outer{1}
      , _buf{std::forward<VectorType>(vector)}
    {
    }

    Array2D(const Id inner, const Id outer = 1, const value_type val = {})
      : _inner{inner}
      , _outer{outer}
      , _buf(outer * inner, val)
    {
    }

    void resize(const Id inner, const Id outer = 1, const value_type val = {})
    {
      this->_inner = inner;
      this->_outer = outer;
      this->_buf.resize(this->_outer * this->_inner, val);
    }

    void addSamples(const Id nnewsamp, const value_type val)
    {
      VectorType newbuf(this->_outer * (this->_inner + nnewsamp), val);
      for (Id o = 0; o < this->_outer; ++o)
      {
        for (Id i = 0; i < this->_inner; ++i)
        {
          newbuf[(o * this->_inner) + i] = this->_buf[(o * this->_inner) + i];
        }
      }
      this->_inner += nnewsamp;
      std::swap(this->_buf, newbuf);
    }

    void deleteSample(const Id isamp)
    {
      VectorType newbuf(this->_buf);
      newbuf.resize(this->_outer * (this->_inner - 1));

      for (Id o = 0; o < this->_outer; ++o)
      {
        Id a{};
        for (Id i = 0; i < this->_inner; ++i)
        {
          if (i == isamp)
          {
            continue;
          }
          newbuf[(o * (this->_inner - 1)) + a] =
            this->_buf[(o * this->_inner) + i];
          a++;
        }
      }

      this->_inner -= 1;
      std::swap(this->_buf, newbuf);
    }

    std::optional<value_type> getValue(const size_t o, const size_t i) const
    {
      if constexpr (std::is_same_v<VectorType, VectorCategory>)
      {
        return this->_buf.getCategory((o * this->_inner) + i);
      }
      else
      {
        return this->_buf[(o * this->_inner) + i];
      }
    }

    bool setValue(const size_t o, const size_t i, const value_type& val)
    {
      if constexpr (std::is_same_v<VectorType, VectorCategory>)
      {
        return this->_buf.setCategory((o * this->_inner) + i, val);
      }
      else
      {
        this->_buf[(o * this->_inner) + i] = val;
        return true;
      }
    }

    Id inner() const { return this->_inner; }

    Id outer() const { return this->_outer; }

    value_type* data() { return this->_buf.data(); }

    const value_type* data() const { return this->_buf.data(); }

    VectorType& getBuffer() { return this->_buf; }

  private:
    Id _inner{};
    Id _outer{};
    VectorType _buf;
    // VectorString _names; // TODO one name per sub-column
  };

  /**
   * @brief Generic implementation of a Column inside a Db
   */
  class GSTLEARN_EXPORT ADbCol
  {
  public:
    ADbCol(String&& name, const Id i = 0, const Id o = 1)
      : _name{std::move(name)}
      , _data{Array2D<VectorDouble>{i, o}}
    {
    }

    template<typename VectorType>
    ADbCol(String&& name, VectorType&& array)
      : _name{std::move(name)}
      , _data{std::forward<VectorType>(array)}
    {
    }

    template<typename VectorType>
    ADbCol(const String& name, VectorType&& array)
      : _name{name}
      , _data{std::forward<VectorType>(array)}
    {
    }

    const String& GetName() const { return this->_name; }

    void SetName(const String& newName) { this->_name = newName; }

    template<typename T, typename VectorType>
    static constexpr bool isConsistent()
    {
      return std::is_same_v<T, typename VectorType::value_type>
          || (std::is_same_v<T, bool> && std::is_same_v<VectorType, VectorBool>)
          || std::is_convertible_v<T, typename VectorType::value_type>;
    }

    template<typename T>
    std::optional<T> getValue(const Id i, const Id o) const
    {
      std::optional<T> v{};
      std::visit(
        [i, o, &v](auto&& arg)
        {
          if (o < 0 || o >= arg.outer() || i < 0 || i >= arg.inner())
          {
            return;
          }
          using VectorType = std::decay_t<decltype(arg)>::vector_type;
          if constexpr (isConsistent<T, VectorType>())
          {
            auto val = arg.getValue(o, i);
            if (val) v = val.value();
          }
        },
        this->_data);
      return v;
    }

    template<typename T>
    bool setValue(const Id i, const Id o, const T& v)
    {
      bool res = false;
      std::visit(
        [i, o, &v, &res](auto&& arg)
        {
          if (o < 0 || o >= arg.outer() || i < 0 || i >= arg.inner())
          {
            return;
          }
          using VectorType = std::decay_t<decltype(arg)>::vector_type;
          if constexpr (isConsistent<T, VectorType>())
          {
            res = arg.setValue(o, i, v);
          }
        },
        this->_data);
      return res;
    }

    template<typename T>
    void addSamples(const Id nnewsamp, const T val)
    {
      std::visit(
        [nnewsamp, val](auto&& arg)
        {
          using VectorType = std::decay_t<decltype(arg)>::vector_type;
          if constexpr (isConsistent<T, VectorType>())
          {
            arg.addSamples(nnewsamp, val);
          }
        },
        this->_data);
    }

    void deleteSample(const Id isamp)
    {
      std::visit([isamp](auto&& arg) { arg.deleteSample(isamp); }, this->_data);
    }

    template<typename VectorType>
    std::optional<std::reference_wrapper<VectorType>> getVector()
    {
      auto* pv = std::get_if<Array2D<VectorType>>(&this->_data);
      if (pv == nullptr) return {};
      return pv->getBuffer();
    }

  private:
    String _name;
    std::variant<
      Array2D<VectorDouble>,
      Array2D<VectorFloat>,
      Array2D<VectorInt>,
      Array2D<VectorUChar>,
      Array2D<VectorString>,
      Array2D<VectorBool>,
      Array2D<VectorCategory>>
      _data;
  };

  /**
   * @brief Similar to vtkFieldData
   */
  class DbData
  {
  public:
    template<typename VectorType>
    void AddArray(VectorType&& array, String&& name)
    {
      Id index{};
      for (auto& col: this->_cols)
      {
        if (col.GetName() == name)
        {
          this->RemoveArray(index);
        }
        index++;
      }
      this->_cols.emplace_back(
        std::move(name), std::forward<VectorType>(array));
    }

    template<typename VectorType>
    void AddArray(VectorType&& array, const String& name)
    {
      Id index{};
      for (auto& col: this->_cols)
      {
        if (col.GetName() == name)
        {
          this->RemoveArray(index);
        }
        index++;
      }
      this->_cols.emplace_back(name, std::forward<VectorType>(array));
    }

    bool RenameArray(const Id index, const String& newName)
    {
      // check if newName already in use
      for (const auto& col: this->_cols)
      {
        if (col.GetName() == newName)
        {
          return false;
        }
      }

      this->_cols[index].SetName(newName);
      return true;
    }

    void RemoveArray(const String& name)
    {
      const auto col = this->GetArray(name);
      const auto id = col ? col->second : -1;
      this->RemoveArray(id);
    }

    void RemoveArray(const Id index)
    {
      if (index < 0 || index >= static_cast<Id>(this->_cols.size()))
      {
        return;
      }

      this->_cols.erase(this->_cols.begin() + index);
    }

    std::optional<std::reference_wrapper<ADbCol>> GetArray(const Id index)
    {
      if (index < 0 || index >= static_cast<Id>(this->_cols.size()))
      {
        return {};
      }
      return {this->_cols[index]};
    }

    std::optional<std::reference_wrapper<const ADbCol>>
      GetArray(const Id index) const
    {
      if (index < 0 || index >= static_cast<Id>(this->_cols.size()))
      {
        return {};
      }
      return {this->_cols[index]};
    }

    std::optional<std::pair<std::reference_wrapper<ADbCol>, Id>>
      GetArray(const String& name)
    {
      Id index{};
      for (auto& col: this->_cols)
      {
        if (col.GetName() == name)
        {
          return {{col, index}};
        }
        index++;
      }

      return {};
    }

    std::optional<String> GetArrayName(const Id index)
    {
      auto arr = this->GetArray(index);
      return arr ? std::optional<String>{arr->get().GetName()} : std::nullopt;
    }

    bool HasArray(const String& name)
    {
      auto arr = this->GetArray(name);
      return arr.has_value();
    }

    template<typename T>
    std::optional<T>
      getValue(const Id iech, const Id icol, const Id isubcol = 0) const
    {
      const auto array = this->GetArray(icol);
      if (!array) return {};
      const auto val = array->get().getValue<T>(iech, isubcol);
      return val;
    }

    template<typename T>
    bool setValue(const Id iech, const Id icol, const Id subcol, const T value)
    {
      const auto array = this->GetArray(icol);
      if (!array) return false;
      return array->get().setValue<T>(iech, subcol, value);
    }

  private:
    std::vector<ADbCol> _cols;
  };

} // namespace gstlrn

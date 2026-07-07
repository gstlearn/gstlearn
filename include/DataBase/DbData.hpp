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

#include "DataBase/DbCol.hpp"

#include <functional>
#include <iostream>
#include <optional>

namespace gstlrn
{
  /**
   * @brief Similar to vtkFieldData
   */
  class DbData
  {
  public:
    /**
     * @brief Produces a summary of the DbData content
     */
    void printContents() const
    {
      Id ncols = this->_cols.size();

      std::cout << "DbData contains " << ncols << " columns\n";
      Id icol = 0;
      for (const auto& col: this->_cols)
      {
        std::cout << "Column " << icol << "/" << ncols << " : "
                  << col.getColName() << std::endl;
        ++icol;
      }
    }

    template<typename VectorType>
    void addArray(VectorType&& array, String&& name)
    {
      Id icol{};
      for (auto& col: this->_cols)
      {
        if (col.getColName() == name)
        {
          this->removeArray(icol);
        }
        icol++;
      }
      this->_cols.emplace_back(
        std::move(name), std::forward<VectorType>(array));

      _buildNameToICol();
    }

    template<typename VectorType>
    void addArray(VectorType&& array, const String& name)
    {
      Id icol{};
      for (auto& col: this->_cols)
      {
        if (col.getColName() == name)
        {
          this->removeArray(icol);
        }
        icol++;
      }
      this->_cols.emplace_back(name, std::forward<VectorType>(array));

      _buildNameToICol();
    }

    void removeArray(const String& name)
    {
      const auto icol = _identifyCol(name);
      if (icol >= 0) this->removeArray(icol);

      _buildNameToICol();
    }

    void removeArray(const Id icol)
    {
      this->_cols.erase(this->_cols.begin() + icol);

      _buildNameToICol();
    }

    std::optional<std::reference_wrapper<DbCol>> getArray(const Id icol)
    {
      if (!_isValidCol(icol)) return std::nullopt;
      return {this->_cols[icol]};
    }

    std::optional<std::reference_wrapper<const DbCol>>
      getArray(const Id icol) const
    {
      if (!_isValidCol(icol)) return std::nullopt;
      return {this->_cols[icol]};
    }

    std::optional<std::reference_wrapper<const DbCol>>
      getArray(const String& name) const
    {
      Id icol = _identifyCol(name);
      if (icol < 0) return std::nullopt;
      return {this->_cols[icol]};
    }

    std::optional<String> getName(const Id icol)
    {
      auto arr = this->getArray(icol);
      return arr ? std::optional<String>{arr->get().getColName()}
                 : std::nullopt;
    }

    bool hasArray(const String& name)
    {
      auto arr = this->getArray(name);
      return arr.has_value();
    }

    template<typename T>
    std::optional<T>
      getValue(const Id isample, const Id icol, const Id iversion = 0) const
    {
      const auto array = this->getArray(icol);
      if (!array) return std::nullopt;
      const auto val = array->get().getValue<T>(isample, iversion);
      return val;
    }

    // template<typename T>
    // std::optional<T>
    //   getValue(const Id isample, const String& name, const Id iversion = 0)
    //     const
    // {
    //   const auto array = this->getArray(name);
    //   if (!array) return std::nullopt;
    //   const auto val = array->get().getValue<T>(isample, iversion);
    //   return val;
    // }

    template<typename T>
    bool setValue(
      const T value,
      const Id isample,
      const Id icol,
      const Id iversion = 0)
    {
      const auto array = this->getArray(icol);
      if (!array) return false;
      return array->get().setValue<T>(value, isample, iversion);
    }

    // template<typename T>
    // bool setValue(
    //   const T value,
    //   const Id isample,
    //   const String& name,
    //   const Id iversion = 0)
    // {
    //   const auto array = this->getArray(name);
    //   if (!array) return false;
    //   return array->get().setValue<T>(value, isample, iversion);
    // }

  private:
    bool _isValidCol(const Id icol) const
    {
      return (icol >= 0 && icol < static_cast<Id>(this->_cols.size()));
    }

    /**
     * @brief Build the internal map from column name to column index
     */
    void _buildNameToICol()
    {
      _nameToICol.clear();

      for (Id icol = 0; icol < static_cast<Id>(_cols.size()); ++icol)
      {
        _nameToICol[_cols[icol].getColName()] = icol;
      }
    }

    Id _identifyCol(const String& name) const
    {
      auto it = _nameToICol.find(name);
      if (it != _nameToICol.end())
      {
        return it->second;
      }
      return -1; // Not found
    }

  private:
    std::vector<DbCol> _cols;
    std::unordered_map<String, Id> _nameToICol;
  };

} // namespace gstlrn

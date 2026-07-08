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
#include "Enum/ERole.hpp"

#include <functional>
#include <iostream>
#include <optional>

namespace gstlrn
{
  struct DbColInfo
  {
    DbCol col;
    ERole role = ERole::UNDEFINED;
    Id rank = 0;

    template<typename VectorType>
    DbColInfo(
      String name,
      VectorType&& array,
      const ERole& role = ERole::UNDEFINED,
      Id rank = 0)
      : col(std::move(name), std::forward<VectorType>(array))
      , role(role)
      , rank(rank)
    {
    }
  };

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
      Id ncols = this->_colInfos.size();

      std::cout << "DbData contains " << ncols << " columns\n";
      Id icol = 0;
      for (const auto& c: this->_colInfos)
      {
        std::cout << "Column " << icol << "/" << ncols << " : "
                  << c.col.getColName();
        if (c.role != ERole::UNDEFINED)
        {
          std::cout << " (Role: " << c.role.getKey() << ", Rank: " << c.rank
                    << ")";
        }
        std::cout << std::endl;
        ++icol;
      }
    }

    template<typename VectorType>
    void addArray(
      VectorType&& array,
      const String& name,
      const ERole& role = ERole::UNDEFINED,
      Id rank = 0)
    {
      Id icol{};
      for (auto& c: this->_colInfos)
      {
        if (c.col.getColName() == name)
        {
          removeArray(icol);
          break;
        }
        icol++;
      }
      _colInfos.emplace_back(name, std::forward<VectorType>(array), role, rank);

      // Update the map as the list of Columns may have been modified
      _buildNameToICol();
    }

    void removeArray(const String& name)
    {
      const auto icol = _identifyCol(name);
      if (icol >= 0) this->removeArray(icol);
    }

    void removeArray(const Id icol)
    {
      this->_colInfos.erase(this->_colInfos.begin() + icol);

      // Update the map as the list of Columns may have been modified
      _buildNameToICol();
    }

    std::optional<std::reference_wrapper<DbCol>> getArray(const Id icol)
    {
      if (!_isValidCol(icol)) return std::nullopt;
      return {this->_colInfos[icol].col};
    }

    std::optional<std::reference_wrapper<const DbCol>>
      getArray(const Id icol) const
    {
      if (!_isValidCol(icol)) return std::nullopt;
      return {this->_colInfos[icol].col};
    }

    std::optional<std::reference_wrapper<const DbCol>>
      getArray(const String& name) const
    {
      Id icol = _identifyCol(name);
      if (icol < 0) return std::nullopt;
      return {this->_colInfos[icol].col};
    }

    std::optional<String> getName(const Id icol) const
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

    std::optional<ERole> getRole(Id icol) const
    {
      if (!_isValidCol(icol)) return std::nullopt;
      return _colInfos[icol].role;
    }

    std::optional<Id> getRank(Id icol) const
    {
      if (!_isValidCol(icol)) return std::nullopt;
      return _colInfos[icol].rank;
    }

    bool setRole(Id icol, const ERole& role, Id rank = 0)
    {
      if (!_isValidCol(icol)) return false;
      _colInfos[icol].role = role;
      _colInfos[icol].rank = rank;
      return true;
    }

  private:
    bool _isValidCol(const Id icol) const
    {
      return (icol >= 0 && icol < static_cast<Id>(this->_colInfos.size()));
    }

    /**
     * @brief Build the internal map from column name to column index
     */
    void _buildNameToICol()
    {
      _nameToICol.clear();

      Id icol = 0;
      for (const auto& info: _colInfos)
      {
        _nameToICol[info.col.getColName()] = icol;
        ++icol;
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
    std::vector<DbColInfo> _colInfos;
    std::unordered_map<String, Id> _nameToICol;
  };

} // namespace gstlrn

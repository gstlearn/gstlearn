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

#include "Basic/Utilities.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorT.hpp"
#include "gstlearn_export.hpp"

#include <optional>
#include <type_traits>
#include <variant>

namespace gstlrn
{

  /**
   * @brief Generic implementation of a Column inside a Db
   */
  class GSTLEARN_EXPORT DbCol
  {
  public:
    DbCol(String&& name)
      : _name{std::move(name)}
      , _nSample{0}
      , _nVersion{1}
      , _data{VectorDouble{}}
    {
    }

    template<typename VectorType>
    DbCol(String&& name, VectorType&& array)
      : _name{std::move(name)}
      , _nSample{static_cast<Id>(array.size())}
      , _nVersion{1}
      , _data{array}
    {
    }

    template<typename VectorType>
    DbCol(const String& name, VectorType&& array)
      : _name{name}
      , _nSample{static_cast<Id>(array.size())}
      , _nVersion{1}
      , _data{array}
    {
    }

    /**
     * @brief Get the Name of the Column
     *
     * @return const String&
     */
    const String& getColName() const { return this->_name; }

    /**
     * @brief Get the number of samples contained in the Column
     *
     * @return const Id
     */
    Id getNSample() const { return this->_nSample; }

    /**
     * @brief Get the number of versions contained in the Column
     *
     * @return const Id
     */
    Id getNVersion() const { return this->_nVersion; }

    bool setNVersion(const Id nVersion)
    {
      if (nVersion <= 0)
      {
        return false;
      }
      if (!isMultiple(getNSample(), nVersion))
      {
        messerr(
          "Dimension of '_data' (%d) should be a multiple of the number of "
          "versions (%d).",
          getNSample(), nVersion);
        return false;
      }
      this->_nVersion = nVersion;
      this->_nSample = getNSample() / nVersion;
      return true;
    }

    /**
     * @brief Check the consistency of input with the type of the array
     *
     * @tparam T
     * @tparam VectorType
     */
    template<typename T, typename VectorType>
    static constexpr bool isConsistent()
    {
      return std::is_convertible_v<T, typename VectorType::value_type>;
    }

    /**
     * @brief Get the Value of an element of the Column
     *
     * @tparam T
     * @param isample Rank of the sample
     * @param iversion Rank of the version
     * @return std::optional<T>, e.g. the value of the sample if it exists, std::nullopt otherwise
     */
    template<typename T>
    std::optional<T> getValue(const Id isample, const Id iversion = 0) const
    {
      if (!this->_isValidIDs(isample, iversion)) return std::nullopt;

      std::optional<T> v{};
      std::visit(
        [isample, iversion, &v](auto&& arg)
        {
          using VectorType = std::decay_t<decltype(arg)>;
          if constexpr (isConsistent<T, VectorType>())
          {
            v = arg[isample];
          }
        },
        this->_data);
      return v;
    }

    /**
     * @brief Set the Value of an element of the Column
     *
     * @tparam T
     * @param v Value to be written
     * @param isample Rank of the sample
     * @param iversion Rank of the version
     */
    template<typename T>
    bool setValue(const T& v, const Id isample, const Id iversion = 0)
    {
      if (!this->_isValidIDs(isample, iversion)) return false;

      bool res = false;
      std::visit(
        [isample, iversion, &v, &res](auto&& arg)
        {
          using VectorType = std::decay_t<decltype(arg)>;
          if constexpr (isConsistent<T, VectorType>())
          {
            arg[isample] = v;
            res = true;
          }
        },
        this->_data);
      return res;
    }

    /**
     * @brief Get the array of a specific type
     *
     * @tparam VectorType The type of the array to get
     * @return const VectorType* Pointer to the array if it exists, nullptr otherwise
     *
     * @remark This method has been added to allow getArray to be used from DbData.
     */
    template<typename VectorType>
    const VectorType* getArray() const
    {
      return std::get_if<VectorType>(&_data);
    }

  private:
    bool _isValidIDs(Id isample, Id iversion = 0) const
    {
      return isample >= 0 && isample < _nSample && iversion >= 0
          && iversion < _nVersion;
    }

  private:
    String _name;
    Id _nSample{0};
    Id _nVersion{0};
    std::variant<
      VectorDouble,
      VectorFloat,
      VectorInt,
      VectorUChar,
      VectorString,
      VectorBool>
      _data;
  };

} // namespace gstlrn

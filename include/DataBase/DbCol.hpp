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
#include "DataBase/Array2D.hpp"
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
    DbCol(String&& name, const Id nsample = 0, const Id nversion = 1)
      : _name{std::move(name)}
      , _data{Array2D<VectorDouble>{nsample, nversion}}
    {
    }

    template<typename VectorType>
    DbCol(String&& name, VectorType&& array)
      : _name{std::move(name)}
      , _data{std::forward<VectorType>(array)}
    {
    }

    template<typename VectorType>
    DbCol(const String& name, VectorType&& array)
      : _name{name}
      , _data{std::forward<VectorType>(array)}
    {
    }

#ifndef SWIG
    DbCol(const DbCol&) = default;
    DbCol(DbCol&&) = default;
    DbCol& operator=(const DbCol&) = default;
    DbCol& operator=(DbCol&&) = default;
    ~DbCol() = default;
#endif

    /**
     * @brief Get the Name of the Column
     *
     * @return const String&
     */
    const String& getColName() const { return this->_name; }

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
      std::optional<T> v{};
      std::visit(
        [isample, iversion, &v](auto&& arg)
        {
          if (iversion < 0 || iversion >= arg.outer() || isample < 0
              || isample >= arg.inner())
          {
            return;
          }
          using VectorType = std::decay_t<decltype(arg)>::vector_type;
          if constexpr (isConsistent<T, VectorType>())
          {
            v = arg[iversion][isample];
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
      bool res = false;
      std::visit(
        [isample, iversion, &v, &res](auto&& arg)
        {
          if (iversion < 0 || iversion >= arg.outer() || isample < 0
              || isample >= arg.inner())
          {
            return;
          }
          using VectorType = std::decay_t<decltype(arg)>::vector_type;
          if constexpr (isConsistent<T, VectorType>())
          {
            arg[iversion][isample] = v;
            res = true;
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

    void deleteSample(const Id isample)
    {
      std::visit(
        [isample](auto&& arg) { arg.deleteSample(isample); }, this->_data);
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
      Array2D<VectorBool>>
      _data;
  };

} // namespace gstlrn

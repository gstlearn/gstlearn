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

#include "gstlearn_export.hpp"

#include "DataBase/Dictionary.hpp"
#include <cstddef>
#include <optional>
#include <vector>

namespace gstlrn
{

  class VectorCategory
  // class GSTLEARN_EXPORT VectorCategory
  {
  public:
    using Category = Dictionary::Category;
    using value_type = Category;

    VectorCategory(
      const size_t nsample = 0,
      const Dictionary& dict = Dictionary())
      : _data(nsample)
      , _dict{dict}
    {
    }

    inline void resize(const size_t count) { this->_data.resize(count); }

    inline size_t size() const { return this->_data.size(); }

    std::optional<Category> getCategory(const size_t isample) const;

    bool setCategory(const size_t isample, const Category& cat);

#ifndef SWIG
    /**
     * @brief Provides access to a sample category.
     *
     * @param isample Sample index.
     *
     * @return Reference to the optional category associated with the sample.
     *
     * @warning No bounds checking is performed.
     */
    std::optional<Category>& operator[](const size_t isample)
    {
      return this->_data[isample];
    }

    /**
     * @brief Provides read-only access to a sample category.
     *
     * @param isample Sample index.
     *
     * @return Constant reference to the optional category associated with
     *         the sample.
     *
     * @warning No bounds checking is performed.
     */
    const std::optional<Category>& operator[](const size_t isample) const
    {
      return this->_data[isample];
    }
#endif

  private:
    std::vector<std::optional<Category>> _data;
    Dictionary _dict;
  };

} // namespace gstlrn

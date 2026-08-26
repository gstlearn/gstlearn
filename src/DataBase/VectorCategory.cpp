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
#include "DataBase/VectorCategory.hpp"

namespace gstlrn
{
  /**
   * @brief Returns the category associated with a sample.
   *
   * @param isample Sample index.
   *
   * @return The category associated with the sample, or an empty optional
   *         if the sample has no category or if the index is out of range.
   */
  std::optional<Dictionary::Category>
    VectorCategory::getCategory(const size_t isample) const
  {
    if (isample >= this->_data.size()) return {};
    return this->_data[isample];
  }

  /**
   * @brief Sets the category associated with a sample.
   *
   * The operation fails if the sample index is out of range or if the
   * specified category is not defined in the associated dictionary.
   *
   * @param isample Sample index.
   * @param cat Category to assign to the sample.
   *
   * @return @c true if the category was successfully assigned,
   *         @c false otherwise.
   */
  bool VectorCategory::setCategory(const size_t isample, const Category& cat)
  {
    if (isample >= this->_data.size()) return false;
    if (!this->_dict.hasCategory(cat)) return false;
    this->_data[isample] = cat;
    return true;
  }
} // namespace gstlrn

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
#include <optional>
#include <type_traits>
#include <variant>

namespace gstlrn
{

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
      VectorType newbuf(this->_outer * (this->_inner - 1));

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

#ifndef SWIG
    std::span<value_type> operator[](const size_t o)
    {
      return {&this->_buf[o * this->_inner], static_cast<size_t>(this->_inner)};
    }

    std::span<const value_type> operator[](const size_t o) const
    {
      return {&this->_buf[o * this->_inner], static_cast<size_t>(this->_inner)};
    }
#endif

    Id inner() const { return this->_inner; }

    Id outer() const { return this->_outer; }

    value_type* data() { return this->_buf.data(); }

    const value_type* data() const { return this->_buf.data(); }

    VectorType& getBuffer() { return this->_buf; }

  private:
    Id _inner{};
    Id _outer{};
    VectorType _buf;
    VectorString _names;
  };

} // namespace gstlrn

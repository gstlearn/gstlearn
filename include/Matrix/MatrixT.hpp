#pragma once

#ifdef USE_BOOST_SPAN
#include <boost/core/span.hpp>
#else
#include <span>
#endif

#include <algorithm>
#include <vector>

namespace gstlrn
{

  template<typename T>
  class MatrixT
  {
  public:
    MatrixT(size_t nrow = 0, size_t ncol = 0, T defaultValue = T())
      : _nrow(nrow)
      , _ncol(ncol)
      , _data(_nrow * _ncol, defaultValue)
    {
    }

    T getValue(size_t row, size_t col) const
    {
      return _data[(row * _ncol) + col];
    }

    void setValue(size_t row, size_t col, const T& value)
    {
      _data[(row * _ncol) + col] = value;
    }

    T* getRowPtr(size_t row) { return &_data[row * _ncol]; }

    const T* getRowPtr(size_t row) const { return &_data[row * _ncol]; }

    void clear()
    {
      _nrow = 0;
      _ncol = 0;
      _data.clear();
    }

#ifdef USE_BOOST_SPAN
    using span = boost::span<T>;
    using constspan = boost::span<const T>;
#else
    using span = std::span<T>;
    using constspan = std::span<const T>;
#endif

    span getRow(size_t row) { return {getRowPtr(row), getRowPtr(row) + _ncol}; }

    constspan getRow(size_t row) const
    {
      return {getRowPtr(row), getRowPtr(row) + _ncol};
    }

    T& operator()(size_t row, size_t col) { return _data[(row * _ncol) + col]; }

    const T& operator()(size_t row, size_t col) const
    {
      return _data[(row * _ncol) + col];
    }

    void resize(size_t nrow, size_t ncol, T defaultValue = T())
    {
      _nrow = nrow;
      _ncol = ncol;
      _data.resize(_nrow * _ncol);
      fill(defaultValue);
    }

    void fill(const T& value) { std::fill(_data.begin(), _data.end(), value); }

    size_t getSize() const { return _nrow * _ncol; }

    bool empty() const { return getSize() == 0; }

    size_t getNRows() const { return _nrow; }

    size_t getNCols() const { return _ncol; }

    const std::vector<T>& getData() const { return _data; }

  private:
    size_t _nrow;
    size_t _ncol;
    std::vector<T> _data; // Stockage contigu en 1D
  };
} // namespace gstlrn

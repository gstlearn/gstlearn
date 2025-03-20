#pragma once
#include <cstddef>
#include <vector>

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

  MatrixT(const MatrixT& other)
    : _nrow(other._nrow)
    , _ncol(other._ncol)
    , _data(other._data)
  {
  }

  MatrixT& operator=(const MatrixT& other)
  {
    if (this != &other)
    {
      _nrow = other._nrow;
      _ncol = other._ncol;
      _data = other._data;
    }
    return *this;
  }

  T* getRowPtr(size_t row)
  {
    return &_data[row * _ncol];
  }
  T& operator()(size_t row, size_t col)
  {
    return _data[row * _ncol + col];
  }
  const T& operator()(size_t row, size_t col) const
  {
    return _data[row * _ncol + col];
  }
  void resize(size_t nrow, size_t ncol, T defaultValue = T())
  {
    _nrow = nrow;
    _ncol = ncol;
    _data.assign(_nrow * _ncol, defaultValue);
  }
  void fill(const T& value)
  {
    std::fill(_data.begin(), _data.end(), value);
  }
  size_t getSize() const { return _nrow * _ncol; }
  size_t getNRows() const { return _nrow; }
  size_t getNCols() const { return _ncol; }

private:
  size_t _nrow; 
  size_t _ncol;
  std::vector<T> _data; // Stockage contigu en 1D
};

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
#pragma once

#include <array>
#ifdef USE_BOOST_SPAN
#  include <boost/core/span.hpp>
#else
#  include <span>
#endif

namespace gstlrn
{
// Template générique pour tableaux 3D à tailles fixes
template<typename T, size_t NX, size_t NY, size_t NZ>
class Array3D
{
private:
  std::array<T, NX * NY * NZ> data {};

public:
  // Accès classique (i,j,k)
  constexpr T& operator()(size_t i, size_t j, size_t k) noexcept
  {
    return data[i * NY * NZ + j * NZ + k];
  }
  constexpr const T& operator()(size_t i, size_t j, size_t k) const noexcept
  {
    return data[i * NY * NZ + j * NZ + k];
  }

  // Accès type arr[i][j][k]
  struct LayerYConst
  {
    const T* base;
    constexpr const T* operator[](size_t j) const noexcept { return base + j * NZ; }
  };

  struct LayerY
  {
    T* base;
    constexpr T* operator[](size_t j) noexcept { return base + j * NZ; }
  };

  constexpr LayerY operator[](size_t i) noexcept
  {
    return LayerY {data.data() + i * NY * NZ};
  }

  constexpr LayerYConst operator[](size_t i) const noexcept
  {
    return LayerYConst {data.data() + i * NY * NZ};
  }

#ifndef SWIG
#ifdef USE_BOOST_SPAN
  using span = boost::span<T>;
  using constspan = boost::span<const T>;
#else
  using span = std::span<T>;
  using constspan = std::span<const T>;
#endif
  // Vue sur la "colonne k" (dernier indice) pour i,j fixes
  span row(size_t i, size_t j) noexcept
  {
    return span(data.data() + i * NY * NZ + j * NZ, NZ);
  }
  constspan row(size_t i, size_t j) const noexcept
  {
    return constspan(data.data() + i * NY * NZ + j * NZ, NZ);
  }
#endif
};

} // namespace gstlrn
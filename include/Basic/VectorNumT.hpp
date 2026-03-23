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

#include "Basic/Undefined.hpp"
#include "Basic/VectorT.hpp"
#include "geoslib_define.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <type_traits>
#include <vector>

namespace gstlrn
{

  // Forward declaration nécessaire pour pouvoir spécialiser le trait final_scalar_type
  template<typename>
  class VectorNumT;

  // Trait récursif (au niveau namespace) pour récupérer le type scalaire final
  template<typename U>
  struct final_scalar_type
  {
    using type = U;
  };

  template<typename U>
  struct final_scalar_type<VectorNumT<U>>
  {
    using type = typename final_scalar_type<U>::type;
  };

  template<typename U>
  using final_scalar_type_t = typename final_scalar_type<U>::type;

  /***************************************************************************
   **
   ** Vector of T values (numerical type).
   ** T type must define copy constructor and assignment operator
   ** T type must override numerical operators (+, -, *, /)
   **
   ***************************************************************************/
  template<typename T>
  class VectorNumT: public VectorT<T>
  {
  public:
    typedef VectorT<T> Parent;
    typedef std::vector<T> Vector;
    typedef typename Vector::value_type value_type;
    typedef typename Vector::size_type size_type;
    typedef typename Vector::iterator iterator;
    typedef typename Vector::const_iterator const_iterator;
    typedef typename Vector::reverse_iterator reverse_iterator;
    typedef typename Vector::const_reverse_iterator const_reverse_iterator;
    using scalar_type = final_scalar_type_t<T>; // type scalaire final concret

  public:
    inline VectorNumT()
      : Parent()
    {
    }

    inline VectorNumT(const Vector& vec)
      : Parent(vec)
    {
    }

    inline VectorNumT(size_type count, const T& value = {})
      : Parent(count, value)
    {
    }

    template<class InputIt>
    inline VectorNumT(InputIt first, InputIt last)
      : Parent(first, last)
    {
    }

    inline VectorNumT(const VectorNumT& other) = default;
#ifndef SWIG
    inline VectorNumT(std::initializer_list<T> init)
      : Parent(init)
    {
    }
#endif
    inline ~VectorNumT() = default;

  private:
    // === Fonction utilitaire : nom du type final ===
    template<typename U>
    static std::string get_final_type_name()
    {
      using BaseT = final_scalar_type_t<U>;

      if constexpr (std::is_same_v<BaseT, int>)
        return "int";
      else if constexpr (std::is_same_v<BaseT, long>)
        return "long";
      else if constexpr (std::is_same_v<BaseT, long long>)
        return "Id";
      else if constexpr (std::is_same_v<BaseT, float>)
        return "float";
      else if constexpr (std::is_same_v<BaseT, double>)
        return "double";
      else if constexpr (std::is_same_v<BaseT, unsigned int>)
        return "unsigned int";
      else
        return "unknown";
    }

    // Helper pour extraire le type scalaire final
    template<typename U>
    struct BaseType
    {
      using type = U;
    };

    template<typename U>
    struct BaseType<std::vector<U>>
    {
      using type = typename BaseType<U>::type;
    };

    template<typename U>
    static bool isScalarNA(const U& val)
    {
      if constexpr (std::is_arithmetic_v<U>)
        return isNA(val);
      else
        return false;
    }

  public:
    inline bool isEqual(const VectorNumT& other, double eps = 1.e-10) const;
    inline bool isConstant();
    template<typename U>
    bool sameDimension(const U& other) const;

    inline double sum() const;
    inline double prod() const;
    inline double minimum(bool flagAbs = false) const;
    inline double maximum(bool flagAbs = false) const;
    inline double mean() const;
    inline double median() const;
    inline double variance(bool scaleByN = false) const;
    inline double stdv(bool scaleByN = false) const;
    inline double norm(Id normType = 2) const;
    inline double norm2() const;

    inline Id count(Id flagDef = 0) const;
    inline void identify() const;
    inline void
      dump(const String& title = String(), bool newLineAfterTitle = true) const;

    inline double innerProduct(const VectorNumT<T>& v, Id size = 0) const;

    inline double normTo(const VectorNumT<T>& other) const;
    inline double correlation(const VectorNumT<T>& other) const;
    inline void normalizeInPlace(Id normType = 2);
  };

  template<typename T>
  bool VectorNumT<T>::isEqual(const VectorNumT& other, double eps) const
  {
    if (other.size() != this->size()) return false;

    for (size_type i = 0, n = this->size(); i < n; i++)
    {
      const auto& a = this->at(i);
      const auto& b = other.at(i);

      if constexpr (std::is_arithmetic_v<T>)
      {
        // Cas scalaire (int, double, etc.)
        if (ABS(a - b) > eps) return false;
      }
      else
      {
        // Cas récursif : T est lui-même un VectorNumT<U>
        if (!a.isEqual(b, eps)) return false;
      }
    }
    return true;
  }

  template<typename T>
  bool VectorNumT<T>::isConstant()
  {
    if (VectorNumT::_v.empty()) return false;
    T refval = VectorNumT::_v.at(0);
    for (auto v: VectorNumT::_v)
    {
      if (v != refval) return false;
    }
    return true;
  }

  template<typename T>
  template<typename U>
  bool VectorNumT<T>::sameDimension(const U& other) const
  {
    // Vérifie que 'other' est un VectorNumT avec le même type d'élément
    if constexpr (!std::is_same_v<U, VectorNumT<T>>)
    {
      return false; // types incompatibles
    }

    if (this->size() != other.size()) return false;

    // Si les éléments sont eux-mêmes des vecteurs, récursion
    if constexpr (!std::is_arithmetic_v<T>)
    {
      for (size_t i = 0; i < this->size(); ++i)
      {
        if (!this->_v[i].sameDimension(other._v[i])) return false;
      }
    }

    return true;
  }

  template<typename T>
  double VectorNumT<T>::maximum(bool flagAbs) const
  {
    double mymax = std::numeric_limits<double>::lowest();
    if (VectorNumT::_v.empty()) return mymax;

    for (const auto& v: VectorNumT::_v)
    {
      if constexpr (std::is_arithmetic_v<T>)
      {
        if (isScalarNA(v)) continue;
        auto myv = v;
        if (flagAbs) myv = ABS(myv);
        if (myv > mymax) mymax = myv;
      }
      else
      {
        auto subv = v.maximum(flagAbs);
        if (subv > mymax) mymax = subv;
      }
    }
    return mymax;
  }

  template<typename T>
  double VectorNumT<T>::sum() const
  {
    double result = 0;
    if (this->_v.empty()) return result;

    for (const auto& v: this->_v)
    {
      if constexpr (std::is_arithmetic_v<T>)
      {
        if (isScalarNA(v)) continue;
        result += static_cast<double>(v);
      }
      else
      {
        result += v.sum();
      }
    }
    return result;
  }

  template<typename T>
  double VectorNumT<T>::minimum(bool flagAbs) const
  {
    double mymin = std::numeric_limits<double>::max();
    if (this->_v.empty()) return mymin;

    for (const auto& v: this->_v)
    {
      if constexpr (std::is_arithmetic_v<T>)
      {
        if (isScalarNA(v)) continue;

        auto val = static_cast<double>(v);
        if (flagAbs) val = std::fabs(val);

        if (val < mymin) mymin = val;
      }
      else
      {
        double subv = v.minimum(flagAbs);
        if (subv < mymin) mymin = subv;
      }
    }

    return mymin;
  }

  template<typename T>
  double VectorNumT<T>::mean() const
  {
    if (this->_v.empty()) return getNA<double>();

    double sum = 0.;
    double count = 0.;

    for (const auto& v: this->_v)
    {
      if constexpr (std::is_arithmetic_v<T>)
      {
        if (isScalarNA(v)) continue;
        sum += v;
        count += 1.;
      }
      else
      {
        sum += v.mean();
        count += 1.;
      }
    }

    if (count > 0.) return sum / count;
    return TEST;
  }

  template<typename T>
  double VectorNumT<T>::median() const
  {
    if (this->_v.empty()) return getNA<double>();

    VectorNumT<double> medians;

    for (const auto& v: this->_v)
    {
      if constexpr (std::is_arithmetic_v<T>)
      {
        if (isScalarNA(v)) continue;
        medians.push_back(static_cast<double>(v));
      }
      else
      {
        medians.push_back(v.median());
      }
    }

    if (medians.empty()) return getNA<double>();

    std::sort(medians.begin(), medians.end());
    size_t n = medians.size();

    if (n % 2 == 1) return medians[n / 2];
    return 0.5 * (medians[n / 2 - 1] + medians[n / 2]);
  }

  template<typename T>
  double VectorNumT<T>::variance(bool scaleByN) const
  {
    if (this->size() <= 0) return static_cast<double>(NAN);

    double sum = 0.;
    double sumsq = 0.;
    double count = 0.;

    for (const auto& v: this->_v)
    {
      if constexpr (std::is_arithmetic_v<T>)
      {
        if (isScalarNA(v)) continue;
        sum += v;
        sumsq += v * v;
        count += 1.;
      }
      else
      {
        double subvar = v.variance(scaleByN);
        double submean = v.mean();
        count += 1.;
        sum += submean;
        sumsq += subvar + submean * submean; // Var(X) + mean(X)^2
      }
    }

    if (count == 0.) return TEST;

    double mean = sum / count;
    double var;

    if (scaleByN)
      var = sumsq / count - mean * mean;
    else
      var = (sumsq - count * mean * mean) / (count - 1.);

    return var;
  }

  template<typename T>
  double VectorNumT<T>::stdv(bool scaleByN) const
  {
    double var = this->variance(scaleByN);
    if (isScalarNA(var)) return getNA<double>();
    return std::sqrt(var);
  }

  template<typename T>
  double VectorNumT<T>::norm2() const
  {
    double value = norm(2);
    return value * value;
  }

  template<typename T>
  double VectorNumT<T>::norm(Id normType) const
  {
    if (this->_v.empty()) return 0.;

    double result = 0.;

    if (normType == 0)
    {
      for (const auto& v: this->_v)
      {
        double val;
        if constexpr (std::is_arithmetic_v<T>)
        {
          if (isScalarNA(v)) continue;
          val = ABS(v);
        }
        else
        {
          val = v.norm(normType);
        }
        if (val > result) result = val;
      }
    }
    else if (normType == 1)
    {
      for (const auto& v: this->_v)
      {
        if constexpr (std::is_arithmetic_v<T>)
        {
          if (isScalarNA(v)) continue;
          result += ABS(v);
        }
        else
        {
          result += v.norm(normType);
        }
      }
    }
    else if (normType == 2)
    {
      for (const auto& v: this->_v)
      {
        double val;
        if constexpr (std::is_arithmetic_v<T>)
        {
          if (isScalarNA(v)) continue;
          val = v;
        }
        else
        {
          val = v.norm(normType);
        }
        result += val * val;
      }
      result = std::sqrt(result);
    }

    return result;
  }

  template<typename T>
  double VectorNumT<T>::normTo(const VectorNumT<T>& other) const
  {
    if (!this->sameDimension(other)) return getNA<double>();

    double normValue = 0.;

    for (size_type i = 0, n = this->size(); i < n; ++i)
    {
      const auto& a = this->at(i);
      const auto& b = other.at(i);

      if constexpr (std::is_arithmetic_v<T>)
      {
        if (isScalarNA(a) || isScalarNA(b)) continue;
        double diff = static_cast<double>(a) - static_cast<double>(b);
        normValue += diff * diff;
      }
      else
      {
        // Récursion pour vecteurs imbriqués
        double subNorm = a.normTo(b);
        normValue += subNorm * subNorm;
      }
    }
    normValue = std::sqrt(normValue);
    return normValue;
  }

  template<typename T>
  double VectorNumT<T>::correlation(const VectorNumT<T>& other) const
  {
    if (!this->sameDimension(other) || this->size() <= 0)
      return getNA<double>();

    double sum1 = 0.;
    double sum2 = 0.;
    double sum11 = 0.;
    double sum22 = 0.;
    double sum12 = 0.;
    Id count = 0;

    for (size_type i = 0, n = this->size(); i < n; ++i)
    {
      const auto& a = this->at(i);
      const auto& b = other.at(i);

      if constexpr (std::is_arithmetic_v<T>)
      {
        if (isScalarNA(a) || isScalarNA(b)) continue;

        sum1 += a;
        sum2 += b;
        sum11 += a * a;
        sum22 += b * b;
        sum12 += a * b;
        count++;
      }
      else
      {
        double subCorr = a.correlation(b);
        if (isNA(subCorr)) continue;

        // pour vecteurs imbriqués, on peut approximer par pondération simple
        // ou utiliser une somme des corrélations pondérée par taille
        sum1 += a.mean();
        sum2 += b.mean();
        sum11 += a.variance(true) + a.mean() * a.mean();
        sum22 += b.variance(true) + b.mean() * b.mean();
        sum12 += a.innerProduct(b) / static_cast<double>(a.size());
        count++;
      }
    }

    if (count <= 0) return TEST;

    double m1 = sum1 / count;
    double m2 = sum2 / count;
    double v11 = sum11 / count - m1 * m1;
    double v22 = sum22 / count - m2 * m2;
    double v12 = sum12 / count - m1 * m2;

    if (v11 <= 0. || v22 <= 0.) return TEST;

    return v12 / std::sqrt(v11 * v22);
  }

  template<typename T>
  double VectorNumT<T>::prod() const
  {
    double result = 1;
    if (this->_v.empty()) return result;

    for (const auto& v: this->_v)
    {
      if constexpr (std::is_arithmetic_v<T>)
      {
        if (isScalarNA(v)) continue;
        result *= static_cast<double>(v);
      }
      else
      {
        result *= v.prod();
      }
    }
    return result;
  }

  template<typename T>
  Id VectorNumT<T>::count(Id flagDef) const
  {
    if (this->_v.empty()) return 0;

    Id cnt = 0;

    for (const auto& v: this->_v)
    {
      if constexpr (std::is_arithmetic_v<T>)
      {
        if ((flagDef == 0) || (flagDef == +1 && !isScalarNA(v))
            || (flagDef == -1 && isScalarNA(v)))
          cnt++;
      }
      else
      {
        cnt += v.count(flagDef);
      }
    }

    return cnt;
  }

  template<typename T>
  double VectorNumT<T>::innerProduct(const VectorNumT<T>& v, Id size) const
  {
    Id n =
      (size > 0) ? size : static_cast<Id>(std::min(this->size(), v.size()));
    if (!this->sameDimension(v)) return 0.;

    double prod = 0.;

    for (Id i = 0; i < n; ++i)
    {
      const auto& a = this->at(i);
      const auto& b = v.at(i);

      if constexpr (std::is_arithmetic_v<T>)
      {
        if (isScalarNA(a) || isScalarNA(b)) continue;
        prod += static_cast<double>(a) * static_cast<double>(b);
      }
      else
      {
        prod += a.innerProduct(b);
      }
    }

    return prod;
  }

  template<typename T>
  void VectorNumT<T>::normalizeInPlace(Id normType)
  {
    double normValue = this->norm(normType);
    if (normValue == 0.) return;

    for (size_type i = 0, n = this->size(); i < n; i++)
    {
      auto& v = this->operator[](i);

      if constexpr (std::is_arithmetic_v<T>)
      {
        if (isScalarNA(v)) continue;
        v = static_cast<T>(v / normValue);
      }
      else
      {
        v.normalizeInPlace(normType);
      }
    }
  }

  /**
   * @brief Identify the VectorNumT
   *
   * @tparam T Input argument (template)
   */
  template<typename T>
  void VectorNumT<T>::identify() const
  {
    if (VectorNumT::size() <= 0) return;

    std::cout << "Vector of";
    if constexpr (!std::is_arithmetic_v<T>) std::cout << " Vector of";
    std::cout << " < " << get_final_type_name<T>() << " > " << std::endl;
  }

  template<typename T>
  void VectorNumT<T>::dump(const String& title, bool newLineAfterTitle) const
  {
    if (VectorNumT::size() <= 0) return;

    if (!title.empty())
    {
      if (newLineAfterTitle)
        std::cout << title << std::endl;
      else
        std::cout << title << " : ";
    }
    // --- Cas 1 : T est un type arithmétique → on affiche directement le vecteur courant
    if constexpr (std::is_arithmetic_v<T>)
    {
      this->display();
    }
    // --- Cas 2 : T est un VectorNumT<U> → on regarde le type U
    else
    {
      using Elem = typename T::value_type;

      // Cas 2a : le sous-type est arithmétique → on appelle display() sur chaque sous-vecteur
      if constexpr (std::is_arithmetic_v<Elem>)
      {
        for (const auto& v: this->_v) v.display();
      }
      // Cas 2b : le sous-type n’est pas arithmétique → on descend récursivement
      else
      {
        for (const auto& v: this->_v) v.dump();
      }
    }
  }

#ifndef SWIG
  template<typename T>
  std::ostream&
    operator<<(std::ostream& os, const VectorNumT<VectorNumT<T>>& vec)
  {
    os << "[";
    for (Id i = 0, n = static_cast<Id>(vec.size()); i < n; i++)
    {
      os << vec.at(i).toString();
      if (i != n - 1) os << " ";
    }
    os << "]";
    return os;
  }
#endif

  typedef VectorNumT<Id> VectorInt;
  typedef VectorNumT<double> VectorDouble;
  typedef VectorNumT<float> VectorFloat;
  typedef VectorNumT<UChar>
    VectorUChar; // Use typedef because swig doesn't like 'unsigned char' in two words
  typedef VectorNumT<VectorInt> VectorVectorInt;
  typedef VectorNumT<VectorDouble> VectorVectorDouble;
  typedef VectorNumT<VectorFloat> VectorVectorFloat;

} // namespace gstlrn

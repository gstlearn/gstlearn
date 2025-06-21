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
#include <cstddef>
#include <utility> // std::pair
#include <tuple>   // structured bindings

class LowerTriangularRange {
public:
    class Iterator {
    public:
        using value_type = std::pair<size_t, size_t>;
        using reference = const value_type;
        using pointer = void;
        using iterator_category = std::input_iterator_tag;
        using difference_type = std::ptrdiff_t;

        Iterator(size_t n = 0, size_t i = 0, size_t j = 0)
            : n_(n), i_(i), j_(j) {}

        reference operator*() const {
            return {i_, j_};
        }

        Iterator& operator++() {
            if (j_ < i_) {
                ++j_;
            } else {
                ++i_;
                j_ = 0;
            }
            return *this;
        }

        bool operator==(const Iterator& other) const {
            return i_ == other.i_;
        }

        bool operator!=(const Iterator& other) const {
            return !(*this == other);
        }

    private:
        size_t n_;
        size_t i_;
        size_t j_;
    };

    explicit LowerTriangularRange(size_t n) : n_(n) {}

    Iterator begin() const { return Iterator(n_, 0, 0); }
    Iterator end() const { return Iterator(n_, n_, 0); }

private:
    size_t n_;
};

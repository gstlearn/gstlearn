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
#include "Covariances/ACov.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
class ACor;

class GSTLEARN_EXPORT CovBase: public ACov
{
public:
    CovBase(const MatrixSquareSymmetric &sills = MatrixSquareSymmetric(), ACor* cor = nullptr);
    CovBase(const CovBase &r);
    CovBase& operator=(const CovBase &r);
    virtual ~CovBase();
private:
    MatrixSquareSymmetric _sills;
    ACor* _cor;
};

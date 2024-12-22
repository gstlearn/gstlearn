/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c3) MINES Paris / ARMINES                                        */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/

#include "Covariances/CovBase.hpp"
#include "Covariances/ACor.hpp"
#include "Covariances/ACov.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"

CovBase::CovBase(const MatrixSquareSymmetric &sills,
                 ACor* cor)
: ACov(cor == nullptr? nullptr : cor->getSpace())
, _sills(sills)
, _cor(cor)
{

}

CovBase::CovBase(const CovBase &r)
: ACov(r)
, _sills(r._sills)
, _cor(r._cor)
{

}

CovBase& CovBase::operator=(const CovBase &r)
{
    if (this != &r)
    {
        ACov::operator=(r);
        _sills = r._sills;
        _cor = r._cor;
    }
    return *this;
}

CovBase::~CovBase()
{

}
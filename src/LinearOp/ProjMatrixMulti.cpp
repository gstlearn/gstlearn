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

#include "LinearOp/ProjMatrixMulti.hpp"
#include "LinearOp/ProjMatrix.hpp"

ProjMatrixMulti::ProjMatrixMulti(const std::vector<ProjMatrix*> &proj, int nvar)
: _projs(proj)
, _apicesNumber(0)
, _pointsNumber(0)
, _nvar(nvar)
{
    for(int istruct = 0; istruct < (int)_projs.size();istruct++)
    {
        _apicesNumber += _projs[istruct]->getApexNumber();
        _pointsNumber += _projs[istruct]->getPointNumber();
    } 
    _apicesNumber *= _nvar;
    _pointsNumber *= _nvar;
}

int  ProjMatrixMulti::_point2mesh(const Eigen::VectorXd& inv,
                                        Eigen::VectorXd& outv) const
{
    return 0;
}
int  ProjMatrixMulti::_mesh2point(const Eigen::VectorXd& inv,
                                        Eigen::VectorXd& outv) const
{
    return 0;
}

int ProjMatrixMulti::getApexNumber() const 
{
    return _apicesNumber;
}

int ProjMatrixMulti::getPointNumber() const
{
    return _pointsNumber;
}
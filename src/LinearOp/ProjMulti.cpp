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

#include "LinearOp/ProjMulti.hpp"
#include "Basic/AStringable.hpp"
#include "LinearOp/IProjMatrix.hpp"
#include "Matrix/VectorEigen.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <algorithm>

int ProjMulti::findFirstNoNullOnRow(int j) const
{
    int i = 0;

    while (i < (int)_projs[j].size() && _projs[j][i]==nullptr)
    {
        i++;
    }
    if (i == (int)_projs[j].size())
    {
        i = -1;
        messerr("All the projectors of row %d are nullptr",j);
    }
    return i;
}

int ProjMulti::findFirstNoNullOnCol(int j) const
{
    int i = 0;
    while (i < (int)_projs.size() && _projs[i][j]==nullptr)
    {
        i++;
    }
    if (i == (int)_projs.size())
    {
        i = -1;
        messerr("All the projectors of column %d are nullptr.",j);
    }
    return i;
}


ProjMulti::ProjMulti(const std::vector<std::vector<const IProjMatrix*>> &proj)
: _projs(proj)
, _pointNumber(0)
, _apexNumber(0)
, _nlatent((int)proj[0].size())
, _nvariable((int)proj.size())
{   
    if (proj.empty())
    {
        messerr("Proj is empty.");
        _makeEmpty();
        return;
    }

    if (_nlatent == 0)
    {
        messerr("There is no projection in line 0.");
        _makeEmpty();
        return;
    }

    for (int i = 1; i < (int)proj.size(); i++)
    {
        if ((int)proj[i].size() != _nlatent)
        {
            messerr("All the elements of proj have to share the same size.");
            messerr("Line %d has %d elements instead of %d.",i,proj[i].size(),_nlatent);
            _makeEmpty();
            return;
        }
    } 

    for (int i = 0; i < (int)_projs.size(); i++)
    {
        int fcol = findFirstNoNullOnRow(i);
        if (fcol == -1)
        {   _makeEmpty();
            return;
        }
        auto npoints = _projs[i][fcol]->getPointNumber();
        for (int j = fcol + 1; j < (int)_projs[0].size(); j++)
        {   if(_projs[i][j] != nullptr)
            {
                if (_projs[i][j]->getPointNumber() != npoints)
                {
                    messerr("Inconsistency between the IProjMatrix Point Numbers.");
                    messerr("Element [%d,%d] should have Point Number = %d  instead of %d.",
                            i,j,npoints,_projs[i][j]->getPointNumber());
                    _makeEmpty();
                    return;
                }
            }
        }
        _pointNumbers.push_back(npoints);
        _pointNumber += npoints;
    }

    for (int j = 0; j < (int)_projs[0].size(); j++)
    {   
        int frow = findFirstNoNullOnCol(j);
        if (frow == -1)
        {   
            _makeEmpty();
            return;
        }
        auto nvertex = _projs[frow][j]->getApexNumber();
        for (int i = frow + 1; i < (int)_projs.size(); i++)
        {   
            if (_projs[i][j] != nullptr)
            {
                if (_projs[i][j]->getApexNumber() != nvertex)
                {
                    messerr("Inconsistency between the IProjMatrix Apex Numbers.");
                    messerr("Element [%d,%d] should have Apex Number = %d  instead of %d.",
                            i,j,nvertex,_projs[i][j]->getApexNumber());
                    _makeEmpty();
                    return;
                }
            }
        }
        _apexNumbers.push_back(nvertex);
        _apexNumber += nvertex;
    }
}

void ProjMulti::_makeEmpty()
{
    _pointNumbers.resize(0);
    _apexNumbers.resize(0);
    _nvariable = 0;
    _nlatent = 0;
    _work.resize(0);
    _workmesh.resize(0);
    _projs.resize(0);
}

int  ProjMulti::_addPoint2mesh(const Eigen::VectorXd& inv,
                                        Eigen::VectorXd& outv) const
{
    int iadvar = 0;
    for (int i = 0; i < _nlatent; i++)
    {   
        int iad = 0;
        int nvertex = _apexNumbers[i];
        _workmesh.resize(nvertex);
        VectorEigen::fill(_workmesh,0.);
        for(int j = 0; j < _nvariable; j++)
        {
            if (_projs[j][i] != nullptr)
            {
                Eigen::Map<const Eigen::VectorXd> view(inv.data()+iad,_pointNumbers[j]);
                _projs[j][i]->addPoint2mesh(view,_workmesh);
            }
            iad += _pointNumbers[j];
        }
        VectorEigen::addInPlace(_workmesh,outv,iadvar);
        iadvar += _apexNumbers[i];
    }
    return 0;
}
int  ProjMulti::_addMesh2point(const Eigen::VectorXd& inv,
                                     Eigen::VectorXd& outv) const
{
    int iadvar = 0;
    for (int i = 0; i < _nvariable; i++)
    {   
        int iad = 0;
        int npoint = _pointNumbers[i];
        _work.resize(npoint);
        VectorEigen::fill(_work,0.);
        for(int j = 0; j < _nlatent; j++)
        {
            if (_projs[i][j] != nullptr)
            {
                Eigen::Map<const Eigen::VectorXd> view(inv.data()+iad,_apexNumbers[j]);
                _projs[i][j]->addMesh2point(view,_work);
            }
            iad += _apexNumbers[j];
        }
        VectorEigen::addInPlace(_work,outv,iadvar);
        iadvar += _pointNumbers[i];
    }
    return 0;
}

int ProjMulti::getApexNumber() const 
{
    return _apexNumber;
}

int ProjMulti::getPointNumber() const
{
    return _pointNumber;
}
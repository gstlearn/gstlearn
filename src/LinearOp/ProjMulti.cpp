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
#include "Basic/VectorHelper.hpp"
#include "LinearOp/IProj.hpp"

int ProjMulti::findFirstNoNullOnRow(int j) const
{
    int i = 0;

    while (i < (int)_projs[j].size() && _projs[j][i]==nullptr)
    {
        i++;
    }
    if (i == (int)_projs[j].size())
    {
        messerr("All the projectors of row %d are nullptr",j);
        return -1;    
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

bool ProjMulti::_checkArg(const std::vector<std::vector<const IProj*>> &projs) const
{
    if (projs.empty())
    {   
        if (!_silent)
            messerr("projs is empty.");
        return true;
    }

    int nvariable = (int)projs.size();
    int nlatent = (int)projs[0].size();
    
    if (nlatent == 0)
    {
        messerr("There is no projection in line 0.");
        return true;
    }

    for (int i = 1; i < nvariable; i++)
    {
        if ((int)projs[i].size() != nlatent)
        {
            messerr("All the elements of proj have to share the same size.");
            messerr("Line %d has %d elements instead of %d.",i,projs[i].size(),nlatent);
            return true;
        }
    } 

    for (int i = 0; i < nvariable; i++)
    {
        int fcol = findFirstNoNullOnRow(i);
        if (fcol == -1)
            return true;
        
        auto npoints = projs[i][fcol]->getNPoint();
        for (int j = fcol + 1; j < nlatent; j++)
        {   if( projs[i][j] != nullptr)
            {
                if ( projs[i][j]->getNPoint() != npoints)
                {
                    messerr("Inconsistency between the IProj Point Numbers.");
                    messerr("Element [%d,%d] should have Point Number = %d  instead of %d.",
                            i,j,npoints, projs[i][j]->getNPoint());
                    return true;
                }
            }
        }
    }

    for (int j = 0; j < nlatent; j++)
    {   
        int frow = findFirstNoNullOnCol(j);
        if (frow == -1)
            return true;

        auto nvertex = projs[frow][j]->getNApex();
        for (int i = frow + 1; i < nvariable; i++)
        {   
            if (projs[i][j] != nullptr)
            {
                if (projs[i][j]->getNApex() != nvertex)
                {
                    messerr("Inconsistency between the IProj Apex Numbers.");
                    messerr("Element [%d,%d] should have Apex Number = %d  instead of %d.",
                            i,j,nvertex,projs[i][j]->getNApex());
                    return true;
                }
            }
        }
    }
    return false;
}

void ProjMulti::_init()
{
    _nvariable = (int)_projs.size();
    _nlatent = (int)_projs[0].size();

    for (int i = 0; i < _nvariable; i++)
    {
        int fcol = findFirstNoNullOnRow(i);
        auto npoints = _projs[i][fcol]->getNPoint();
        _pointNumbers.push_back(npoints);
        _pointNumber += npoints;
    }

    for (int j = 0; j < _nlatent; j++)
    {   
        int frow = findFirstNoNullOnCol(j);
        auto nvertex = _projs[frow][j]->getNApex();
        _apexNumbers.push_back(nvertex);
        _apexNumber += nvertex;
    }


}

ProjMulti::~ProjMulti()
{
     _clear();
}

ProjMulti::ProjMulti(const std::vector<std::vector<const IProj*>> &projs, bool silent)
: _projs(projs)
, _pointNumber(0)
, _apexNumber(0)
, _nlatent(0)
, _nvariable(0)
, _silent(silent)
{   

    if (_checkArg(_projs))
    {   if (_projs.size()!=0)
        {
            messerr("Problem in initialization of ProjMulti.");
        }
        _projs.resize(0);
        return;
    }

    _init();
}

int ProjMulti::_addPoint2mesh(const constvect inv, vect outv) const
{
  vect wms;
  int iadvar = 0;
  for (int i = 0; i < _nlatent; i++)
  {
    int iad     = 0;
    int nvertex = _apexNumbers[i];
    _workmesh.resize(nvertex);
    std::fill(_workmesh.begin(), _workmesh.end(), 0.);
    for (int j = 0; j < _nvariable; j++)
    {
      if (_projs[j][i] != nullptr)
      {
        constvect view(inv.data() + iad, _pointNumbers[j]);
        wms = vect(_workmesh);
        _projs[j][i]->addPoint2mesh(view, wms);
      }
      iad += _pointNumbers[j];
    }

    vect outs(outv.data() + iadvar, _workmesh.size());
    VectorHelper::addInPlace(wms, outs);
    iadvar += _apexNumbers[i];
  }
  return 0;
}

int ProjMulti::_addMesh2point(const constvect inv, vect outv) const
{
  vect ws;
  int iadvar = 0;
  for (int i = 0; i < _nvariable; i++)
  {
    int iad    = 0;
    int npoint = _pointNumbers[i];
    _work.resize(npoint);
    std::fill(_work.begin(), _work.end(), 0.);
    for (int j = 0; j < _nlatent; j++)
    {
      if (_projs[i][j] != nullptr)
      {
        constvect view(inv.data() + iad, _apexNumbers[j]);
        ws = vect(_work);
        _projs[i][j]->addMesh2point(view, ws);
      }
      iad += _apexNumbers[j];
    }
    vect outs(outv.data() + iadvar, _work.size());
    VectorHelper::addInPlace(ws, outs);
    iadvar += _pointNumbers[i];
  }
  return 0;
}

int ProjMulti::getNApex() const
{
  return _apexNumber;
}

int ProjMulti::getNPoint() const
{
  return _pointNumber;
}

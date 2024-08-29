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

#include "LinearOp/ProjMultiMatrix.hpp"
#include "Basic/AStringable.hpp"
#include "LinearOp/IProjMatrix.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "Matrix/MatrixSparse.hpp"

static std::vector<std::vector<const IProjMatrix*>> castToBase(std::vector<std::vector<const ProjMatrix*>> vect)
{
    std::vector<std::vector<const IProjMatrix*>> casted;
    for (auto &e : vect)
    {
        std::vector<const IProjMatrix*> temp;
        for (auto &f: e)
        {
            temp.push_back(static_cast<const IProjMatrix*>(f));
        }
        casted.push_back(temp);
    } 
    return casted;
}

std::vector<std::vector<const ProjMatrix*>> ProjMultiMatrix::create(std::vector<const ProjMatrix*> &vectproj, int nvariable)
{
    int nlatent = vectproj.size();
    std::vector<std::vector<const ProjMatrix*>> result;

    
    for (int i = 0; i < nlatent; i++)
    {
        if (vectproj[i] == nullptr)
        {
            messerr("Projmatrix shouldn't be nullptr.");
            return result;
        }
    }

    int npoint = vectproj[0]->getPointNumber();

    for (int i = 1; i < nlatent; i++)
    {
        if (vectproj[i]->getPointNumber() != npoint)
        {
            messerr("All the ProjMatrix should have the same number of Point.");
            messerr("Element %d has %d Point instead of %d.",i,vectproj[i]->getPointNumber(),npoint);
            return result;
        }
    }
    result.resize(nvariable);
    for (int i = 0; i < nvariable; i++)
    {
        std::vector<const ProjMatrix*> e(nlatent * nvariable,nullptr);
        for (int j = 0; j < nlatent; j++)
        {   
            e[j * nvariable + i] = vectproj[j];
        }
        result[i] = e;
    }
    return result;
}

ProjMultiMatrix::ProjMultiMatrix(const std::vector<std::vector<const ProjMatrix*>> &proj)
: ProjMulti(castToBase(proj))
{
    const VectorInt& pointNumbers = getPointNumbers();
    const VectorInt& apexNumbers  = getApexNumbers();

    MatrixSparse currentrow;
    for (int i = 0; i < getNVariable(); i++)
    {   
        currentrow = MatrixSparse(0,0);
        for (int j = 0; j < getNLatent(); j++)
        {
            if (_projs[i][j] != nullptr)
            {
                MatrixSparse::glueInPlace(&currentrow,(MatrixSparse*)proj[i][j],0,1);
            }
            else 
            {
                auto tempMat = MatrixSparse(pointNumbers[i],apexNumbers[j]);
                MatrixSparse::glueInPlace(&currentrow,&tempMat,0,1);
            }
         
        }
        
        MatrixSparse::glueInPlace((MatrixSparse*)this,&currentrow,1,0);
    }   
}

int  ProjMultiMatrix::_addPoint2mesh(const Eigen::VectorXd& inv,
                                        Eigen::VectorXd& outv) const
{
    addProdMatVecInPlaceToDest(inv, outv,true);
    return 0;
}
int  ProjMultiMatrix::_addMesh2point(const Eigen::VectorXd& inv,
                                        Eigen::VectorXd& outv) const
{
    addProdMatVecInPlaceToDest(inv, outv,false);
    return 0;
}

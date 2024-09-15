#include "Covariances/NoStatFunctional.hpp"
#include "Basic/VectorNumT.hpp"
#include "geoslib_define.h"

NoStatFunctional::NoStatFunctional(const AFunctional* func)
:_func(func)
{

}

NoStatFunctional::~NoStatFunctional()
{

}

void NoStatFunctional::_informField(const VectorVectorDouble& coords,
                                    VectorDouble& tab,
                                    bool verbose) 
{
    DECLARE_UNUSED(verbose)
    int size = (int)coords.size();
    int ndim = (int)coords[0].size();
    VectorDouble vec(ndim);
    for (int icoords = 0; icoords < size; icoords++)
    {
        for (int idim = 0; idim < ndim; idim ++)
        {
            vec[idim] = coords[idim][icoords];
        }
        tab[icoords] =  _func->getFunctionValue(vec);
    }
}

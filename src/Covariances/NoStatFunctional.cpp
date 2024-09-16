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
    int size = (int)coords[0].size();
    int ndim = (int)coords.size();
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

String NoStatFunctional::toString(const AStringFormat* strfmt) const
{
    DECLARE_UNUSED(strfmt)
    std::stringstream sstr;
    if (_func == nullptr) return sstr.str();
    sstr << ANoStat::toString(strfmt);
    AStringFormat sf;
    if (strfmt != nullptr) sf = *strfmt;
    if (sf.getLevel() > 0)
        sstr << "Functional" << std::endl;
    return sstr.str();
}

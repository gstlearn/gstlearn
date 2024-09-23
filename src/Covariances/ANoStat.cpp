#include "Covariances/ANoStat.hpp"
#include "Basic/VectorNumT.hpp"
#include "Mesh/AMesh.hpp"
#include "Db/Db.hpp"
#include "geoslib_define.h"

ANoStat::ANoStat()
{

}

ANoStat::~ANoStat()
{

}

double ANoStat::getValueOnDb(int iech,int icas) const
{
    if (icas == 1)
      return getValueOnDbIn(iech);
    if (icas == 2)
      return getValueOnDbOut(iech);
    return TEST;
}

double ANoStat::getValueOnMesh(int iapex,bool center) const
{
    if (center)
        return getValueOnMeshByMesh(iapex);
    return getValueOnMeshByApex(iapex);
}


bool ANoStat::getValuesOnDb(int icas1, int iech1, double* val1,
                            int icas2, int iech2, double* val2) const
{
    *val1 = getValueOnDb(iech1,icas1);
    *val2 = getValueOnDb(iech2,icas2);    
    if (FFFF(*val1) && FFFF(*val2)) return false;
    
    if (FFFF(*val1)) *val1 = *val2;
    if (FFFF(*val2)) *val2 = *val1;

    return *val1 != *val2;                    
}

double ANoStat::getValueOnDbOut(int iech) const
{
    return _tabdbout[iech];
}

double  ANoStat::getValueOnDbIn(int iech) const
{
    return _tabdbin[iech];
}

double ANoStat::getValueOnMeshByMesh(int imesh) const
{
    return _tabmesh[imesh];

}

String ANoStat::toString(const AStringFormat* strfmt) const
{
    DECLARE_UNUSED(strfmt)
    std::stringstream sstr;
    return sstr.str();
}

double ANoStat::getValueOnMeshByApex(int iapex) const
{
    return  _tabvertices[iapex];
}

void ANoStat::informField(const VectorVectorDouble & coords, VectorDouble& tab, bool verbose)
{
    tab.resize(coords[0].size());
    _informField(coords, tab, verbose);
}


void ANoStat::informMeshByMesh(const AMesh* amesh,bool verbose) 
{
    VectorVectorDouble coords = amesh->getAllCenterCoordinates();
    informField(coords,_tabmesh, verbose);
}

void ANoStat::informMeshByApex(const AMesh* amesh, bool verbose)
{
    VectorVectorDouble coords = amesh->getAllCoordinates();
    informField(coords,_tabvertices, verbose);
}
  
void ANoStat::informDbIn(const Db* dbin, bool verbose)
{
    _informDb(dbin,_tabdbin, verbose);
}

void ANoStat::informDbOout(const Db* dbout, bool verbose)
{
    _informDb(dbout,_tabdbout, verbose);
}

void ANoStat::_informDb(const Db* db, VectorDouble &res, bool verbose)
{
    VectorVectorDouble coords = db->getAllCoordinates(false);
    informField(coords, res, verbose);
}


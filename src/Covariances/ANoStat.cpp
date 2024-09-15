#include "Covariances/ANoStat.hpp"
#include "Basic/VectorNumT.hpp"
#include "Mesh/AMesh.hpp"
#include "Db/Db.hpp"

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
    return getValueOnDbOut(iech);
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

double ANoStat::getValueOnMeshByApex(int iapex) const
{
    return  _tabvertices[iapex];
}

void ANoStat::informMeshByMesh(const AMesh* amesh)
{
    VectorVectorDouble coords = amesh->getAllCenterCoordinates();
    _informField(coords,_tabmesh);
}

void ANoStat::informMeshByApex(const AMesh* amesh)
{
    VectorVectorDouble coords = amesh->getAllCoordinates();
    _informField(coords,_tabvertices);
}
  
void ANoStat::informDbIn(const Db* dbin)
{
    _informDb(dbin,_tabdbin);
}

void ANoStat::informDbOout(const Db* dbout)
{
    _informDb(dbout,_tabdbout);
}

void ANoStat::_informDb(const Db* db, VectorDouble &res)
{
    VectorVectorDouble coords = db->getAllCoordinates(false);
    _informField(coords, res);
}


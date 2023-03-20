/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Basic/AException.hpp"
#include "Basic/String.hpp"
#include "Basic/Utilities.hpp"
#include "Model/ANoStat.hpp"
#include "Model/NoStatFunctional.hpp"
#include "Db/Db.hpp"
#include "Mesh/AMesh.hpp"

#include <math.h>

NoStatFunctional::NoStatFunctional()
    : ANoStat(),
      _func(nullptr)
{
  VectorString code = {"A"};
  (void) addNoStatElems(code);
}

NoStatFunctional::NoStatFunctional(const AFunctional* func)
    : ANoStat(),
      _func(nullptr)
{
  VectorString code = {"A"};
  (void) addNoStatElems(code);
  _func = func;
}

NoStatFunctional::NoStatFunctional(const NoStatFunctional &m)
    : ANoStat(m),
      _func(m._func)
{
}

NoStatFunctional& NoStatFunctional::operator= (const NoStatFunctional &m)
{
  if (this != &m)
  {
    ANoStat::operator=(m);
    _func = m._func;
  }
  return *this;
}

NoStatFunctional::~NoStatFunctional()
{
}


int NoStatFunctional::attachToMesh(const AMesh* mesh, bool verbose) const
{
  if (mesh->getNDim() != 2)
  {
    messerr("This function is only defined in 2-D space");
    return 1;
  }
  return ANoStat::attachToMesh(mesh,verbose);
}

int NoStatFunctional::attachToDb(Db* db, int icas, bool verbose) const
{
  if (db == nullptr) return 0;
  if (db->getNDim() != 2)
  {
    messerr("This function is only defined in 2-D space");
    return 1;
  }
  return ANoStat::attachToDb(db,icas,verbose);
}


/**
 * Returns the value of a non-stationary parameter at a target sample
 * @param igrf  Rank of the GRF
 * @param icov  Rank of the Covariance
 * @param type  Type of non-stationary element (EConsElem)
 * @param iv1   Rank of the first variable (optional)
 * @param iv2   Rank of the second variable (optional)
 * @param icas  Additional identifier (0 for Meshing; 1 for Dbin; 2 for Dbout)
 * @param rank  Rank of the target (in Meshing (0); in Dbin (1) or in Dbout (2)
 * @return
 */
double NoStatFunctional::getValue(int igrf,
                                  int icov,
                                  const EConsElem& type,
                                  int iv1,
                                  int iv2,
                                  int icas,
                                  int rank) const
{
  int ipar = getRank(igrf, icov, type, iv1, iv2);
  return getValueByParam(ipar, icas, rank);
}

/**
 * Return the value of the non-stationary parameter (ipar) at target (rank)
 * @param ipar  Rank of the non-stationary parameter
 * @param icas  Additional identifier
 * @param rank  Rank of the target
 * @return
 */
double NoStatFunctional::getValueByParam(int ipar, int icas, int rank) const
{
  if (ipar != 0)
    my_throw("Invalid rank when searching for Non-stationary parameter");

  // Dispatch

  VectorDouble vec;
  if (icas == 0)
  {

    // From Meshing

    if (_amesh == nullptr) return TEST;
    if (rank < 0 || rank > _amesh->getNApices()) return TEST;
    for (int idim = 0; idim < _amesh->getNDim(); idim++)
      vec.push_back(_amesh->getApexCoor(rank, idim));
  }
  else if (icas == 1)
  {

    // From Dbin

    if (_dbin == nullptr) return TEST;
    if (rank < 0 || rank > _dbin->getSampleNumber()) return TEST;
    for (int idim = 0; idim < _dbin->getNDim(); idim++)
      vec.push_back(_dbin->getCoordinate(rank, idim));
  }
  else if (icas == 2)
  {

    // From Dbout

    if (_dbout == nullptr) return TEST;
    if (rank < 0 || rank > _dbout->getSampleNumber()) return TEST;
    for (int idim = 0; idim < _dbin->getNDim(); idim++)
       vec.push_back(_dbout->getCoordinate(rank, idim));
  }
  else
  {
    my_throw("Invalid argument 'icas'");
  }
  return _func->getFunctionValue(vec);
}

String NoStatFunctional::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  if (_func == nullptr) return sstr.str();

  sstr << ANoStat::toString(strfmt);

  AStringFormat sf;
  if (strfmt != nullptr) sf = *strfmt;
  if (sf.getLevel() > 0)
    sstr << "Functional" << std::endl;
  return sstr.str();
}


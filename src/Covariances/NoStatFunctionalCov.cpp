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
#include "Basic/AException.hpp"
#include "Covariances/ANoStatCov.hpp"
#include "Covariances/NoStatFunctionalCov.hpp"
#include "Db/Db.hpp"
#include "Mesh/AMesh.hpp"

#include <math.h>

NoStatFunctionalCov::NoStatFunctionalCov(const VectorString& code)
    : ANoStatCov(),
      _func(nullptr)
{
  (void) addNoStatElems(code);
}

NoStatFunctionalCov::NoStatFunctionalCov(const AFunctional* func, const VectorString& code)
    : ANoStatCov(),
      _func(nullptr)
{
  (void) addNoStatElems(code);
  _func = func;
}

NoStatFunctionalCov::NoStatFunctionalCov(const NoStatFunctionalCov &m)
    : ANoStatCov(m),
      _func(m._func)
{
}

NoStatFunctionalCov& NoStatFunctionalCov::operator= (const NoStatFunctionalCov &m)
{
  if (this != &m)
  {
    ANoStatCov::operator=(m);
    _func = m._func;
  }
  return *this;
}

NoStatFunctionalCov::~NoStatFunctionalCov()
{
}


int NoStatFunctionalCov::attachToMesh(const AMesh* mesh, bool center,bool verbose) const
{
  if (mesh->getNDim() != 2)
  {
    messerr("This function is only defined in 2-D space");
    return 1;
  }
  return ANoStatCov::attachToMesh(mesh,center,verbose);
}

int NoStatFunctionalCov::attachToDb(Db* db, int icas, bool verbose) const
{
  if (db == nullptr) return 0;
  if (db->getNDim() != 2)
  {
    messerr("This function is only defined in 2-D space");
    return 1;
  }
  return ANoStatCov::attachToDb(db,icas,verbose);
}

/**
 * Returns the value of a non-stationary parameter at a target sample
 * @param type  Type of non-stationary element (EConsElem)
 * @param icas  Additional identifier (0 for Meshing; 1 for Dbin; 2 for Dbout)
 * @param rank  Rank of the target (in Meshing (0); in Dbin (1) or in Dbout (2)
 * @param icov  Rank of the Covariance
 * @param iv1   Rank of the first variable (optional)
 * @param iv2   Rank of the second variable (optional)
 * @param igrf  Rank of the GRF
 * @return
 */
double NoStatFunctionalCov::getValue(const EConsElem &type,
                                  int icas,
                                  int rank,
                                  int iv1,
                                  int iv2) const
{
  if (! _isValid(icas, rank)) return TEST;
  int ipar = getRank(type, iv1, iv2);
  return getValueByParam(ipar, icas, rank);
}

/**
 * Return the value of the non-stationary parameter (ipar) at target (rank)
 * @param ipar  Rank of the non-stationary parameter
 * @param icas  Additional identifier
 * @param rank  Rank of the target
 * @return
 */
double NoStatFunctionalCov::getValueByParam(int ipar, int icas, int rank) const
{
  DECLARE_UNUSED(ipar);
  if (! _isValid(icas, rank)) return TEST;

  // Dispatch

  VectorDouble vec;
  if (icas == 0)
  {
    const AMesh* amesh = _getAMesh();

    // From Meshing

    for (int idim = 0, ndim = amesh->getNDim(); idim < ndim; idim++)
      vec.push_back(amesh->getCenterCoordinate(rank, idim));
  }
  else if (icas == 1)
  {

    // From Dbin

    const Db* dbin = _getDbin();
    int ndim = dbin->getNDim();
    vec.resize(ndim);
    dbin->getCoordinatesPerSampleInPlace(rank, vec);
  }
  else if (icas == 2)
  {

    // From Dbout

    const Db* dbout = _getDbout();
    int ndim = dbout->getNDim();
    vec.resize(ndim);
    dbout->getCoordinatesPerSampleInPlace(rank, vec);
  }
  else
  {
    my_throw("Invalid argument 'icas'");
  }
  return _func->getFunctionValue(vec);
}

String NoStatFunctionalCov::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  if (_func == nullptr) return sstr.str();

  sstr << ANoStatCov::toString(strfmt);

  AStringFormat sf;
  if (strfmt != nullptr) sf = *strfmt;
  if (sf.getLevel() > 0)
    sstr << "Functional" << std::endl;
  return sstr.str();
}


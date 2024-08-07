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
#include "Model/ANoStat.hpp"
#include "Model/NoStatFunctional.hpp"
#include "Db/Db.hpp"
#include "Mesh/AMesh.hpp"

#include <math.h>

NoStatFunctional::NoStatFunctional(const VectorString& code)
    : ANoStat(),
      _func(nullptr)
{
  (void) addNoStatElems(code);
}

NoStatFunctional::NoStatFunctional(const AFunctional* func, const VectorString& code)
    : ANoStat(),
      _func(nullptr)
{
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
 * @param type  Type of non-stationary element (EConsElem)
 * @param icas  Additional identifier (0 for Meshing; 1 for Dbin; 2 for Dbout)
 * @param rank  Rank of the target (in Meshing (0); in Dbin (1) or in Dbout (2)
 * @param icov  Rank of the Covariance
 * @param iv1   Rank of the first variable (optional)
 * @param iv2   Rank of the second variable (optional)
 * @param igrf  Rank of the GRF
 * @return
 */
double NoStatFunctional::getValue(const EConsElem &type,
                                  int icas,
                                  int rank,
                                  int icov,
                                  int iv1,
                                  int iv2,
                                  int igrf) const
{
  if (! _isValid(icas, rank)) return TEST;
  int ipar = getRank(type, icov, iv1, iv2, igrf);
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
  DECLARE_UNUSED(ipar);
  if (! _isValid(icas, rank)) return TEST;

  // Dispatch

  VectorDouble vec;
  if (icas == 0)
  {

    // From Meshing

    for (int idim = 0; idim < _amesh->getNDim(); idim++)
      vec.push_back(_amesh->getCenterCoordinate(rank, idim));
  }
  else if (icas == 1)
  {

    // From Dbin

    int ndim = _dbin->getNDim();
    vec.resize(ndim);
    _dbin->getCoordinatesPerSampleInPlace(rank, vec);
  }
  else if (icas == 2)
  {

    // From Dbout

    int ndim = _dbout->getNDim();
    vec.resize(ndim);
    _dbout->getCoordinatesPerSampleInPlace(rank, vec);
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


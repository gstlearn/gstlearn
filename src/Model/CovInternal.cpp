/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Model/CovInternal.hpp"
#include "Db/Db.hpp"

CovInternal::CovInternal()
    : _icas1(-1),
      _iech1(-1),
      _icas2(-1),
      _iech2(-1),
      _ndim(0),
      _db1(nullptr),
      _db2(nullptr),
      _x1(),
      _x2()
{
}

CovInternal::CovInternal(int icas1,
                         int iech1,
                         int icas2,
                         int iech2,
                         int ndim,
                         const Db* db1,
                         const Db* db2)
    : _icas1(-1),
      _iech1(-1),
      _icas2(-1),
      _iech2(-1),
      _ndim(0),
      _db1(nullptr),
      _db2(nullptr),
      _x1(),
      _x2()
{
  init(icas1, iech1, icas2, iech2, ndim, db1, db2);
}

CovInternal::CovInternal(const CovInternal &m)
    : _icas1(m._icas1),
      _iech1(m._iech1),
      _icas2(m._icas2),
      _iech2(m._iech2),
      _ndim (m._ndim),
      _db1  (m._db1),
      _db2  (m._db2),
      _x1   (m._x1),
      _x2   (m._x2)
{

}

CovInternal& CovInternal::operator=(const CovInternal &m)
{
  if (this != &m)
  {
    _icas1 = m._icas1;
    _iech1 = m._iech1;
    _icas2 = m._icas2;
    _iech2 = m._iech2;
    _ndim  = m._ndim;
    _db1   = m._db1;
    _db2   = m._db2;
    _x1    = m._x1;
    _x2    = m._x2;
  }
  return *this;
}

CovInternal::~CovInternal()
{

}

void CovInternal::init(int icas1,
                       int iech1,
                       int icas2,
                       int iech2,
                       int ndim,
                       const Db* db1,
                       const Db* db2)
{
  _icas1 = icas1;
  _iech1 = iech1;
  _icas2 = icas2;
  _iech2 = iech2;
  _ndim  = ndim;
  _db1   = db1;
  _db2   = db2;
  _calculateCoordinates();
}

void CovInternal::setIech1(int iech1)
{
  _iech1 = iech1;
  _calculateCoordinates();
}

void CovInternal::setIech2(int iech2)
{
  _iech2 = iech2;
  _calculateCoordinates();
}

void CovInternal::_calculateCoordinates()
{
  if (_ndim <= 0) return;

  // First point
  _x1.resize(_ndim);
  if (_icas1 == 1)
  {
    if (_db1 != nullptr && _iech1 >= 0)
    {
      for (int idim = 0; idim<_ndim; idim++)
        _x1[idim] = _db1->getCoordinate(_iech1, idim);
    }
  }
  else
  {
    if (_db2 != nullptr && _iech2 >= 0)
    {
      for (int idim = 0; idim< _ndim; idim++)
        _x1[idim] = _db2->getCoordinate(_iech2, idim);
    }
  }

  // Second point
  _x2.resize(_ndim);
  if (_icas2 == 1)
  {
    if (_db1 != nullptr && _iech1 >= 0)
    {
      for (int idim = 0; idim<_ndim; idim++)
        _x2[idim] = _db1->getCoordinate(_iech1, idim);
    }
  }
  else
  {
    if (_db2 != nullptr && _iech2 >= 0)
    {
      for (int idim = 0; idim< _ndim; idim++)
        _x2[idim] = _db2->getCoordinate(_iech2, idim);
    }
  }
}

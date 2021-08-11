/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
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
      _db2(nullptr)
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
      _db2(nullptr)
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
      _db2  (m._db2)
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
}

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
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"

namespace gstlrn
{ 
class Db;

class GSTLEARN_EXPORT CovInternal
{
public:
  CovInternal();
  CovInternal(Id icas1,
              Id iech1,
              Id icas2,
              Id iech2,
              Id ndim,
              const Db* db1,
              const Db* db2);
  CovInternal(const CovInternal &m);
  CovInternal& operator= (const CovInternal &m);
  virtual ~CovInternal();

  void init(Id icas1,
            Id iech1,
            Id icas2,
            Id iech2,
            Id ndim,
            const Db* db1,
            const Db* db2);

  const Db* getDb1() const { return _db1;   }
  const Db* getDb2() const { return _db2;   }
  Id getIcas1()     const { return _icas1; }
  Id getIcas2()     const { return _icas2; }
  Id getIech1()     const { return _iech1; }
  Id getIech2()     const { return _iech2; }
  Id getNdim()      const { return _ndim;  }
  const VectorDouble& getX1() const { return _x1; }
  const VectorDouble& getX2() const { return _x2; }

  void setDb1(const Db* db1) { _db1   = db1;   }
  void setDb2(const Db* db2) { _db2   = db2;   }
  void setIcas1(Id icas1)   { _icas1 = icas1; }
  void setIcas2(Id icas2)   { _icas2 = icas2; }
  void setIech1(Id iech1);
  void setIech2(Id iech2);
  void setNdim(Id ndim)     { _ndim  = ndim;  }
  void setX1(const VectorDouble& x1) { _x1 = x1; }
  void setX2(const VectorDouble& x2) { _x2 = x2; }

private:
  void _calculateCoordinates();

private:
  Id _icas1;     // Type of Db for first point: 1 for Dbin; 2 for Dbout
  Id _iech1;     // Rank of the first sample within Db1
  Id _icas2;     // Type of Db for second point: 1 for Dbin; 2 for Dbout
  Id _iech2;     // Rank of the second sample within Db2
  Id _ndim;      // Space dimension
  const Db* _db1; // Pointer to the first Db
  const Db* _db2; // Pointer to the second Db
  mutable VectorDouble _x1;
  mutable VectorDouble _x2;
};
}
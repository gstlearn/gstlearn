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
#pragma once

class Db;

class CovInternal
{
public:
  CovInternal();
  CovInternal(int icas1,
              int iech1,
              int icas2,
              int iech2,
              int ndim,
              const Db* db1,
              const Db* db2);
  CovInternal(const CovInternal &m);
  CovInternal& operator= (const CovInternal &m);
  virtual ~CovInternal();

  void init(int icas1,
            int iech1,
            int icas2,
            int iech2,
            int ndim,
            const Db* db1,
            const Db* db2);

  const Db* getDb1() const { return _db1;   }
  const Db* getDb2() const { return _db2;   }
  int getIcas1()     const { return _icas1; }
  int getIcas2()     const { return _icas2; }
  int getIech1()     const { return _iech1; }
  int getIech2()     const { return _iech2; }
  int getNdim()      const { return _ndim;  }

  void setDb1(const Db* db1) { _db1   = db1;   }
  void setDb2(const Db* db2) { _db2   = db2;   }
  void setIcas1(int icas1)   { _icas1 = icas1; }
  void setIcas2(int icas2)   { _icas2 = icas2; }
  void setIech1(int iech1)   { _iech1 = iech1; }
  void setIech2(int iech2)   { _iech2 = iech2; }
  void setNdim(int ndim)     { _ndim  = ndim;  }

private:
  int _icas1;     // Type of Db for first point: 1 for Dbin; 2 for Dbout
  int _iech1;     // Rank of the first sample within Db1
  int _icas2;     // Type of Db for second point: 1 for Dbin; 2 for Dbout
  int _iech2;     // Rank of the second sample within Db2
  int _ndim;      // Space dimension
  const Db* _db1; // Pointer to the first Db
  const Db* _db2; // Pointer to the second Db
};

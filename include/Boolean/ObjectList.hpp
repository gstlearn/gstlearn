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

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Basic/AStringable.hpp"

class AToken;
class Tokens;
class Object;
class DbGrid;
class Db;

class GSTLEARN_EXPORT ObjectList: public AStringable
{
public:
  ObjectList();
  ObjectList(const ObjectList &r);
  ObjectList& operator=(const ObjectList &r);
  virtual ~ObjectList();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int getNObjects(int mode = 0) const;
  void countConditioning(const Db* db,
                          int *nbgrain_arg,
                          int *nbpore_arg,
                          bool verbose);
  int generatePrimary(Db* dbin,
                      DbGrid* dbout,
                      const Tokens* tokens,
                      bool flagStat,
                      double thetaCst,
                      const VectorDouble& dilate = VectorDouble(),
                      int maxiter = 100000);
  int generateSecondary(Db* dbin,
                        DbGrid* dbout,
                        const Tokens* tokens,
                        bool flagStat,
                        double thetaCst,
                        double tmax,
                        const VectorDouble& dilate = VectorDouble(),
                        int maxiter = 100000);
  void projectToGrid(DbGrid* dbout,
                     int iptr_simu,
                     int iptr_rank,
                     int facies);

private:
  int _getRankUncovered(const Db* db, int rank);
  int _getObjectRank(int mode, int rank);
  int _deleteObject(int mode, Db* dbin);
  int _getAverageCount(const DbGrid* dbout,
                       bool flagStat,
                       double thetaCst,
                       const VectorDouble& dilate);

private:
  std::vector<Object*> _objlist;
};

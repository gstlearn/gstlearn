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

#include "Covariances/ANoStat.hpp"
#include "gstlearn_export.hpp"

class Db;

class GSTLEARN_EXPORT NoStatArray : public ANoStat
{
public:
  NoStatArray(){};
  NoStatArray(const Db *dbref,const String& colname);
  NoStatArray(const NoStatArray &m) = delete;
  NoStatArray& operator=(const NoStatArray &m) = delete;
  virtual ~NoStatArray(){};
  String toString(const AStringFormat* strfmt = nullptr) const override;

private :
  void _informField(const VectorVectorDouble& coords,
                    VectorDouble& tab,
                    bool verbose = false) override;

private:
  const Db* _dbNoStat;
  const String _colName;

};

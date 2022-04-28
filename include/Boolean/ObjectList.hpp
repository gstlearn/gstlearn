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
class Object;
class DbGrid;

class GSTLEARN_EXPORT ObjectList: public AStringable
{
public:
  ObjectList(const AToken* atoken);
  ObjectList(const ObjectList &r);
  ObjectList& operator=(const ObjectList &r);
  virtual ~ObjectList();

  /// Interface to AStringable
  String toString(const AStringFormat* strfmt) const;

  void projectToGrid(DbGrid* dbout,
                     int iptr_simu,
                     int iptr_rank,
                     int background,
                     int facies);

private:
  std::vector<Object> _objlist;
};

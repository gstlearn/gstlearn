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
#include "Db/DbStringFormat.hpp"

DbStringFormat::DbStringFormat(int level)
    : AStringFormat(level),
      _params(0),
      _cols(),
      _names(),
      _flagSel(true),
      _mode(1)
{
  if (level == 0)
    _params = FLAG_RESUME;
  else
    _params = FLAG_RESUME | FLAG_VARS;
}

DbStringFormat::DbStringFormat(const DbStringFormat& r)
    : AStringFormat(r),
      _params(r._params),
      _cols(r._cols),
      _names(r._names),
      _flagSel(r._flagSel),
      _mode(r._mode)
{
}

DbStringFormat& DbStringFormat::operator=(const DbStringFormat& r)
{
  if (this != &r)
  {
    AStringFormat::operator=(r);
    _params = r._params;
    _cols = r._cols;
    _names = r._names;
    _flagSel = r._flagSel;
    _mode = r._mode;
  }
  return *this;
}

DbStringFormat::~DbStringFormat()
{
}

bool DbStringFormat::_matchFlag(int flag) const
{
  int reste = _params & flag;
  if (reste > 0)
    return true;
  else
    return false;
}

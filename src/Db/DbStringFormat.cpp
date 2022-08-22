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

DbStringFormat::DbStringFormat(unsigned char params,
                               const VectorString& names,
                               const VectorInt& cols,
                               bool useSel)
    : AStringFormat(1),
      _params(params),
      _cols(cols),
      _names(names),
      _useSel(useSel),
      _mode(1)
{
}

DbStringFormat::DbStringFormat(const DbStringFormat& r)
    : AStringFormat(r),
      _params(r._params),
      _cols(r._cols),
      _names(r._names),
      _useSel(r._useSel),
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
    _useSel = r._useSel;
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

DbStringFormat* DbStringFormat::create(unsigned char params,
                                       const VectorString &names,
                                       const VectorInt &cols,
                                       bool useSel)
{
  return new DbStringFormat(params, names, cols, useSel);
}

DbStringFormat* DbStringFormat::createFromKeys(unsigned char params,
                                               const VectorString &names,
                                               const VectorInt &cols,
                                               bool useSel)
{
  return new DbStringFormat(params, names, cols, useSel);
}

void DbStringFormat::setFlags(bool flag_resume,
                              bool flag_vars,
                              bool flag_extend,
                              bool flag_stats,
                              bool flag_array,
                              bool flag_locator,
                              const VectorString& names,
                              const VectorInt& cols,
                              bool useSel)
{
  _cols = cols;
  _names = names;
  _useSel = useSel;

  _params = 0;
  if (flag_resume)  _params = _params | FLAG_RESUME;
  if (flag_vars)    _params = _params | FLAG_VARS;
  if (flag_extend)  _params = _params | FLAG_EXTEND;
  if (flag_stats)   _params = _params | FLAG_STATS;
  if (flag_array)   _params = _params | FLAG_ARRAY;
  if (flag_locator) _params = _params | FLAG_LOCATOR;
}


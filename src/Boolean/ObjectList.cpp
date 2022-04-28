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
#include "Boolean/ObjectList.hpp"
#include "Boolean/Object.hpp"
#include "Boolean/AToken.hpp"
#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"

#include <math.h>

ObjectList::ObjectList(const AToken* atoken)
    : AStringable(),
      _objlist()
{
}

ObjectList::ObjectList(const ObjectList &r)
    : AStringable(r),
      _objlist(r._objlist)
{
}

ObjectList& ObjectList::operator=(const ObjectList &r)
{
  if (this != &r)
  {
    AStringable::operator =(r);
    _objlist = r._objlist;
  }
  return *this;
}

ObjectList::~ObjectList()
{
}

String ObjectList::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  for (int iobj = 0; iobj < (int) _objlist.size(); iobj++)
  {
    sstr << "Characteristics of the Object: " << iobj+1 << std::endl;
    sstr << _objlist[iobj].toString(strfmt);
  }
  return sstr.str();
}

/*****************************************************************************/
/*!
 **  Project the objects on the output grid
 **
 ** \param[in]  background    Value for the background
 ** \param[in]  facies        Value of the facies assigned
 ** \param[in]  iptr_simu     Pointer for storing the Boolean simulation (or <0)
 ** \param[in]  iptr_rank     Pointer for storing the object ranks (or <0)
 **
 *****************************************************************************/
void ObjectList::projectToGrid(DbGrid* dbout,
                               int iptr_simu,
                               int iptr_rank,
                               int background,
                               int facies)
{

  /* Initialize the output result (not accounting for mask) */

  for (int iech = 0; iech < dbout->getSampleNumber(); iech++)
  {
    if (iptr_simu >= 0) dbout->setArray(iech, iptr_simu, background);
    if (iptr_rank >= 0) dbout->setArray(iech, iptr_rank, TEST);
  }

  /* Loop on the objects */

  for (int iobj = 0; iobj < (int) _objlist.size(); iobj++)
  {
    _objlist[iobj].projectToGrid(dbout, iptr_simu, iptr_rank, facies, iobj+1);
  }
  return;
}


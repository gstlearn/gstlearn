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
#include "Drifts/DriftF.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"
#include "Drifts/EDrift.hpp"

DriftF::DriftF(const CovContext& ctxt)
    : ADriftElem(EDrift::F, ctxt)
{
}

DriftF::DriftF(const DriftF &r)
    : ADriftElem(r)
{
}

DriftF& DriftF::operator=(const DriftF &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftF::~DriftF()
{
}

double DriftF::eval(const Db* db, int iech) const
{
  return db->getExternalDrift(iech,getRankFex());
}

IClonable* DriftF::clone() const
{
  return new DriftF(*this);
}

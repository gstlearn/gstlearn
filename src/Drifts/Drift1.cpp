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
#include "Drifts/Drift1.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"
#include "Drifts/EDrift.hpp"

Drift1::Drift1(const CovContext& ctxt)
    : ADriftElem(EDrift::UC, ctxt)
{
}

Drift1::Drift1(const Drift1 &r)
    : ADriftElem(r)
{
}

Drift1& Drift1::operator=(const Drift1 &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

Drift1::~Drift1()
{
}

double Drift1::eval(const Db* /*db*/, int /*iech*/) const
{
  return 1;
}

ICloneable* Drift1::clone() const
{
  return new Drift1(*this);
}

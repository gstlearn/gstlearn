/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Enum/EDrift.hpp"

#include "Drifts/Drift1.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"

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


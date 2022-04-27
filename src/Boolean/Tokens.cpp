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
#include "../../include/Boolean/Tokens.hpp"

Tokens::Tokens()
    : _tokens()
{
}

Tokens::Tokens(const Tokens &r)
    : _tokens(r._tokens)
{

}

Tokens& Tokens::operator=(const Tokens &r)
{
  if (this != &r)
  {
    _tokens = r._tokens;
  }
  return *this;
}

Tokens::~Tokens()
{
}


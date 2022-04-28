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
#include "Boolean/Object.hpp"
#include "Basic/Vector.hpp"

class AToken;

class GSTLEARN_EXPORT Tokens
{
public:
  Tokens();
  Tokens(const Tokens &r);
  Tokens& operator=(const Tokens &r);
  virtual ~Tokens();

  void addToken(const AToken& token);
  void normalizeProportions();
  Object* generateObject(int ndim) const;

private:
  std::vector<AToken*> _tokens; // List of the Token
};

#pragma once

#include "Basic/AStringable.hpp"

#include "LithoRule/Rule.hpp"


class PGS : public AStringable
{

public:
  PGS();
  PGS(const PGS& pgs);
  PGS& operator=(const PGS &pgs);
  virtual ~PGS();

private:
  Rule _rule;
};



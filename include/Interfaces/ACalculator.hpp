#pragma once

#include "gstlearn_export.hpp"

class GSTLEARN_EXPORT ACalculator
{
  public:
    ACalculator(){};
    virtual ~ACalculator(){};
    virtual void run() = 0;
    virtual bool check() const = 0;
};


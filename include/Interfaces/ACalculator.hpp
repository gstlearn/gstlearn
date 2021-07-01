#ifndef ACALCULATOR_HPP
#define ACALCULATOR_HPP

class ACalculator
{
  public:
    ACalculator(){};
    virtual ~ACalculator(){};
    virtual void run() = 0;
    virtual bool check() const = 0;
};

#endif

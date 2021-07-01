#ifndef CALC_VARIO_EXP_HPP
#define CALC_VARIO_EXP_HPP

#include "Interfaces/ACalculator.hpp"  // ISA
#include "Space/ASpaceObject.hpp"      // ISA
#include "Interfaces/Param.hpp"        // HASA ParamVario
#include "Interfaces/VarioExp.hpp"     // HASA

class Database;
class ASpace;

class CalcVarioExp : public ASpaceObject, public ACalculator
{
  public:
    CalcVarioExp(const ASpace* space = nullptr);
    /// TODO : Not clonable !
    virtual ~CalcVarioExp();
    
    void setParamVario(const ParamVario &p);
    void setInputData(const Database &db);
    Vario* getVario() const; 
    VarioExp getVarioExp() const;

    virtual bool isConsistent(const ASpace* space) const override;
    virtual void run() override;
    bool  check() const override;
    
  private:
    ParamVario  _param;
    Database*   _data;
    VarioExp    _res;
};

#endif

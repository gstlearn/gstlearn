#ifndef VARIO_EXP
#define VARIO_EXP

#include "Space/ASpaceObject.hpp"  // ISA
#include "Interfaces/VarioDir.hpp" // HASA
#include "Interfaces/Param.hpp"    // HASA ParamVario

class Vario;

class VarioExp: public ASpaceObject
{
  public:
    VarioExp(const ASpace* space = nullptr);
    VarioExp(const VarioExp& ref);
    VarioExp& operator=(const VarioExp& ref);
    virtual ~VarioExp();

    virtual bool isConsistent(const ASpace* space) const override;
    
    void fromGeoslib(Vario* vario,const ParamVario& pvario);
    virtual void display_old() const;
    
    int                 getNDir()     const { return _ndir; }
    int                 getNVar()     const { return _nvar; }
    const VectorDouble& getVariance() const { return _variance; }

    VarioDir getVarioDir(int i) const;
    Vario* toGeoslib() const;

  private:
    int                   _ndir;
    int                   _nvar;
    VectorDouble          _variance;
    std::vector<VarioDir> _varioDirs;
    ParamVario            _paramVario;
};

#endif

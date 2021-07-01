#ifndef VARIO_DIR_HPP
#define VARIO_DIR_HPP

#include <iostream>
#include <iomanip>

#include "Interfaces/interface_d.hpp"

#include "Interfaces/VarioValue.hpp"
#include "Interfaces/Param.hpp"
#include "Variogram/Dir.hpp"

#include "Basic/AStringable.hpp"

class VarioDir : public AStringable
{
  public:
    VarioDir(Dir* dir, ParamVarioDir pdir, int nvar);
    VarioDir(const VarioDir& ref);
    VarioDir& operator=(const VarioDir& ref);
    ~VarioDir();

    virtual void display_old() const;
    
    Dir* getDir(int nvar,bool flag_asym) const;
    SpacePoint getNormDir() const;
    int  getNLag() const;
    int  getLag() const;
    VectorDouble getGgs(int i);
    VectorDouble getDists();

  private:
    std::vector<VarioValue> _lagValues;     // Vario values for each lags
    ParamVarioDir           _paramVarioDir; // Definition of the calculation Parameter
};

#endif

#ifndef VARIO_VALUE_HPP
#define VARIO_VALUE_HPP

#include <iostream>        // USE cout
#include <iomanip>         // USE setw

#include "Interfaces/interface_d.hpp" // USE define VectorDouble

#include "Space/SpacePoint.hpp"  // HASA

#include "Basic/AStringable.hpp"

//:WARNING:  dist à la place de _incr
class VarioValue : public AStringable
{
public:
  VarioValue();
  VarioValue(double sw, double dist, VectorDouble gg);
  VarioValue(const VarioValue& ref);
  VarioValue& operator=(const VarioValue& ref);
  virtual ~VarioValue();

  double getSw() const;
  double getDist() const;
  VectorDouble getGg() const;
  virtual void display_old() const;

private:
  double _sw;   //Pondération
  // SpacePoint   _incr;
  double _dist;
  VectorDouble _gg;  //Values of variogrammes
};

#endif

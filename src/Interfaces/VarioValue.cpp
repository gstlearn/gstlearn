#include "Interfaces/VarioValue.hpp"

/*********************************************************************
** Default Constructor
*********************************************************************/
VarioValue::VarioValue(): _sw(0), _dist(0), _gg()
{
}

/*********************************************************************
**  Constructor via value of Vario*
*********************************************************************/
VarioValue::VarioValue(double sw, double dist, VectorDouble gg): _sw(sw), _dist(dist), _gg(gg)
{
}

/*********************************************************************
** Copy Constructor
*********************************************************************/
VarioValue::VarioValue(const VarioValue& ref): _sw(ref._sw), _dist(ref._dist), _gg(ref._gg)
{
}

/*********************************************************************
** Assignment operator
*********************************************************************/
VarioValue& VarioValue::operator=(const VarioValue& ref)
{
  if (this != &ref)
  {
    _sw = ref._sw;
    _dist = ref._dist;
    _gg = ref._gg;
  }
  return(*this);
}

/*********************************************************************
** Destructor
*********************************************************************/
VarioValue::~VarioValue()
{
}

/*********************************************************************
** Getter _sw
*********************************************************************/
double VarioValue::getSw() const
{
  return(_sw);
}

/*********************************************************************
** Getter _dist
*********************************************************************/
double VarioValue::getDist() const
{
  return(_dist);
}

/*********************************************************************
** Getter _gg
*********************************************************************/
VectorDouble VarioValue::getGg() const
{
  return(_gg);
}

void VarioValue::display_old() const
{
  std::cout << std::setw(20) << _sw;
  std::cout << std::setw(20) << _dist;
  for (auto g: _gg)
    std::cout << std::setw(20) << g;
}

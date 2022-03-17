#include "Interfaces/VariableCategorical.hpp"

/**********************************************************************
**Copy Constructor
**********************************************************************/
VariableCategorical::VariableCategorical(const VariableCategorical& ref):AVariableTemplate(ref)
{
  _dico= ref._dico;
}


VariableCategorical::VariableCategorical(const String& name,
                                         const Dictionary& dico)
    : AVariableTemplate(name),
      _dico(dico)
{
}

VariableCategorical::~VariableCategorical()
{
}

/**
 * Clone
 */
VariableCategorical* VariableCategorical::clone() const
{
  return(new VariableCategorical(*this));
}

void VariableCategorical::resize(int n, const Category& val)
{
  if (_dico.hasCategory(val))
  {
    AVariableTemplate::resize(n, val);
  }
  else
  {
    AVariableTemplate::resize(n, _dico.getUndefCategory());
  }
}

VectorDouble VariableCategorical::getValues() const
{
  VectorDouble vals;
  for (const auto& val: _values)
  {
    vals.push_back(val.getValue());
  }
  return(vals);
}

void VariableCategorical::setValue(int i, const Category& cat)
{
  if (_dico.hasCategory(cat))
  {
    AVariableTemplate::setValue(i,cat);
  }
  else
  {
    AVariableTemplate::setValue(i, _dico.getUndefCategory());
  }
} 

void VariableCategorical::setValue(int i, int value)
{
  if (_dico.hasCategory(value))
  {
    AVariableTemplate::setValue(i, _dico.getCategory(value));
  }
  else
  {
    AVariableTemplate::setValue(i, _dico.getUndefCategory());
  }
}

void VariableCategorical::setValue(int i, const String& label)
{
  if (_dico.hasCategory(label))
  {
    AVariableTemplate::setValue(i, _dico.getCategory(label));
  }
  else
  {
    AVariableTemplate::setValue(i, _dico.getUndefCategory());
  }
}

#ifdef _USE_NETCDF
/// TODO: need to Serialize Dictionnary and Categories
netCDF::NcVar VariableCategorical::serialize(netCDF::NcFile& file, std::vector<netCDF::NcDim>& dims) const
{
  netCDF::NcVar res;
  return (res);
}

/// TODO: need to Deserialize Dictionnary and Categories
void VariableCategorical::deserialize(const netCDF::NcFile& file,const netCDF::NcVar& var)
{
}
#endif


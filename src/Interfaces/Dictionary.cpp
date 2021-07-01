#include "Interfaces/Dictionary.hpp"

Dictionary::Dictionary()
{
}

Dictionary::Dictionary(const String& name): _name(name)
{
}

Dictionary::~Dictionary()
{
}

const String& Dictionary::getName() const
{
  return (_name);
}

const Category& Dictionary::getUndefCategory() const
{
  return (_undefCategory);
}

const Category& Dictionary::getCategory(int value) const
{
  for (unsigned int i = 0; i < _categories.size(); i++)
  {
    if (_categories[i].getValue()== value)
    {
      return (_categories[i]);
    }
  }
  return (_undefCategory);
}

const Category& Dictionary::getCategory(const String& label) const
{
  for (unsigned int i = 0; i < _categories.size(); i++)
  {
    if (_categories[i].getLabel() == label)
    {
      return (_categories[i]);
    }
  }
  return (_undefCategory);
}

void Dictionary::addCategory(int value, const String& label)
{
  if ( !hasCategory(value) && !hasCategory(label))
  {
    _categories.push_back(Category(value,label));
  }
}

bool Dictionary::hasCategory(int value) const
{
  for (const auto& cat: _categories)
  {
    if (cat.getValue() == value)
    {
      return (true);
    }
  }
  return (false);
}

bool Dictionary::hasCategory(const String& label) const
{
  for (const auto& cat: _categories)
  {
    if (cat.getLabel() == label)
    {
      return (true);
    }
  }
  return (false);
}

bool Dictionary::hasCategory(const Category& cat) const
{
  if (hasCategory(cat.getValue()))
  {
    return getCategory(cat.getValue()).getLabel() == cat.getLabel();
  }
  return (false);
}

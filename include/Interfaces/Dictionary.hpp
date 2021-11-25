#pragma once

#include "gstlearn_export.hpp"
#include "Interfaces/Category.hpp"
#include "Interfaces/interface_d.hpp"

#include <vector>
#include <string>

class GSTLEARN_EXPORT Dictionary
{
  public:
    Dictionary();
    Dictionary(const std::string&name);
    virtual ~Dictionary();
    
    const String& getName() const;
    const Category& getUndefCategory() const;
    
    const Category& getCategory(int value) const;
    const Category& getCategory(const String& label) const;
    void  addCategory(int value, const String& label);
    bool  hasCategory(int value) const;
    bool  hasCategory(const String& label) const;
    bool  hasCategory(const Category& cat) const;

  private:
    String _name;
    Category _undefCategory;
    std::vector<Category> _categories;
};


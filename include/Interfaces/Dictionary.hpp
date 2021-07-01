#ifndef DICTIONARY_HPP
#define DICTIONARY_HPP

#include <vector>           //USE
#include <string>           //USE

#include "Interfaces/Category.hpp"     //HASA

#include "Interfaces/interface_d.hpp"  //USE

class Dictionary
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

#endif

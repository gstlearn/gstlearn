#ifndef MY_INI_H
#define MY_INI_H
    
#include "geoslib_e.h"

class MyIni {

public:
  MyIni() :
   _filename(),
   _content() {
  }

  ~MyIni() {
    // Write the file
    if (!_filename.empty()) {
      std::fstream file;
      file.open(_filename.c_str(), std::fstream::out);
      std::map<std::string, std::map< std::string, std::string> >::const_iterator itsect = _content.begin();
      while (itsect != _content.end()) {
        file << "[" << itsect->first << "]" << std::endl;
         std::map<std::string, std::string>::const_iterator itkey = itsect->second.begin();
         while (itkey != itsect->second.end()) {
           file << itkey->first << "=" << itkey->second << std::endl;
           itkey++;
         }
         itsect++;
      }
      file.close();
    }
  }

  bool Open(const std::string& filename) {
    // Open the file (create it if not exist)
    std::fstream file;
    _filename = filename;
    file.open(_filename.c_str(), std::fstream::in);
    if (!file.is_open()) {
      file.open(_filename.c_str(), std::fstream::out);
      file.close();
      file.open(_filename.c_str(), std::fstream::in);
    }
    if (!file.is_open()) {
      message("Cannot open %s\n",_filename.c_str());
      return false;
    }

    // Load the file
    std::string line;
    std::string section = "General";
    while (ReadLine(file, line)) {
      if (line.find_first_of("[") == 0 &&
          line.find_first_of("]") == line.size()-1 &&
          line.find_first_of("=") == std::string::npos) {
        section = line.substr(1, line.size()-2);
      }
      else {
        size_t pos = line.find_first_of("=");
        if (pos != std::string::npos) {
          std::string key = Trim(line.substr(0, pos));
          std::string value = Trim(line.substr(pos+1, line.size()-pos-1));
          _content[section][key] = value;
        }
      }
    }

    // Close the file
    file.close();
    return true;
  }

  bool SetValue(const std::string& section,
                const std::string& key,
                const std::string& value) {
    if (_filename.empty())
      return false;
    _content[section][key] = value;
    return true;
  }

  bool GetValue(const std::string& section,
                const std::string& key,
                std::string& value) {
    value.erase();
    if (_content.find(section) != _content.end()) {
      if (_content[section].find(key) != _content[section].end()) {
        value = _content[section][key];
        return true;
      }
    }
    return false;
  }

private:
  bool ReadLine(std::fstream& file, std::string& line) {
    char l[1000];
    if (file.is_open() && !file.eof()) {
      file.getline(l, 1000);
      line = l;
      line = Trim(line);
      return true;
    }
    return false;
  }

  std::string Trim(const std::string& str) {
    std::string s = str;
    // Trim trailing spaces
    if (s.size() > 0) {
      size_t endpos = s.find_last_not_of(" \t");
      s = s.substr(0, endpos+1);
    }

    // Trim leading spaces
    if (s.size() > 0) {
      size_t startpos = s.find_first_not_of(" \t");
      s = s.substr(startpos);
    }
    return s;
  }

private:
  std::string                                                _filename;
  std::map<std::string, std::map<std::string, std::string> > _content;
} ;

#endif // MY_INI_H

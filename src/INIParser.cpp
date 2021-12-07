#include "INIParser.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"

INIParser::INIParser(bool autoSave) :
  ini(),
  FileName(""),
  AutoSave(autoSave)
{
}

INIParser::~INIParser()
{
  if(AutoSave)
    save();
  ini.clear();
}

bool INIParser::load(const std::string &filename)
{
  std::string section, line, key, valeur;
  size_t pos;
  
  FileName = filename;
  ini.clear();
  
  std::ifstream file(filename.c_str());
  if (!file.is_open())
  {
    messerr("Cannot open file %s\n", filename.c_str());
    return false;
  }
  
  while (!file.eof())
  {
    std::getline (file, line);
    
    // Suppress comments
    pos = line.find_first_of(';');
    if(pos != std::string::npos)
      line.erase (line.begin() + pos, line.end());
    else {
      pos = line.find_first_of('#');
      if(pos != std::string::npos)
        line.erase (line.begin() + pos, line.end());
    }
    
    // Continue when the line is not empty
    if(!line.empty())
    {
      // Test if the line corresponds to a section
      pos = line.find_first_of('[');
      if(pos != std::string::npos)
      {
        line.erase(line.begin(), line.begin() + pos+1);
        line.erase(line.begin() + line.find_first_of (']'), line.end());
        section = line;
      }
      else  // Otherwise it corresponds to a keyword
      {
        pos = line.find_first_of('=');

        // If '=' has been found
        if(pos != std::string::npos)
        {
          key = line.substr(0, pos);
          valeur = line.substr(pos+1);
          
          // Suppresses spaces in the keyword
          while(std::string::npos != (pos = key.find_first_of(' ')))
            key.erase(pos);

          // Remove leading and trailing spaces from value
          valeur = trim(valeur);
          // Remove leading and trailing quotes from value
          valeur = trim(valeur, "\"'");
          
          ini[section][key] = valeur;
        }
      }
    }
  }
  
  file.close();
  return (!isEmpty());
}
    
bool INIParser::save(const std::string& filename)
{
  std::map<std::string, std::map<std::string, std::string> >::iterator _it;
  std::map<std::string, std::string>::iterator it;
  
  std::ofstream file;
  
  if(filename == "")
    file.open(FileName.c_str());
  else
    file.open(filename.c_str());
  
  if(!file.is_open())
    return false;
  
  // Write the map in the file
  for(_it=ini.begin(); _it!=ini.end(); ++_it)
  {
    file << "[" << _it->first << "]" << std::endl;
    
    for(it=_it->second.begin(); it!=_it->second.end(); ++it)
      file << it->first << "=" << it->second << std::endl;
  }
  
  file.close();
  return true;
}

template <class T> T INIParser::GetValue(const std::string &section,
                                         const std::string &key,
                                         const T &defaultValue)
{
  std::map<std::string, std::map<std::string, std::string> >::iterator _it=ini.find(section);
  if (_it != ini.end())
  {
    // If the value is found
    std::map<std::string, std::string>::iterator it=_it->second.find(key);
    if (it != _it->second.end())
    {
      // If the keyword is found,
      // convert the value from std::string to the requested type
      T val;
      std::istringstream iss(it->second);
      iss >> val;
      return val;
    }
    else
      return defaultValue;
  }
  else
    return defaultValue;
}

template <class T> void INIParser::SetValue(const std::string &section,
                                            const std::string &key,
                                            const T &Value)
{
  std::ostringstream oss;
  // Value conversion to std::string
  oss << Value;
  // Store the new value in the map
  ini[section][key] = oss.str();
}

// Explicit instanciation for common basic types handled by istream
template bool           INIParser::GetValue<bool>          (const std::string &, const std::string &, const bool &);
template short          INIParser::GetValue<short>         (const std::string &, const std::string &, const short &);
template unsigned short INIParser::GetValue<unsigned short>(const std::string &, const std::string &, const unsigned short &);
template int            INIParser::GetValue<int>           (const std::string &, const std::string &, const int &);
template unsigned int   INIParser::GetValue<unsigned int>  (const std::string &, const std::string &, const unsigned int &);
template long           INIParser::GetValue<long>          (const std::string &, const std::string &, const long &);
template unsigned long  INIParser::GetValue<unsigned long> (const std::string &, const std::string &, const unsigned long &);
template float          INIParser::GetValue<float>         (const std::string &, const std::string &, const float &);
template double         INIParser::GetValue<double>        (const std::string &, const std::string &, const double &);
template long double    INIParser::GetValue<long double>   (const std::string &, const std::string &, const long double &);

template <>
std::string INIParser::GetValue<std::string>(const std::string &section,
                                             const std::string &key,
                                             const std::string &defaultValue)
{
  std::map<std::string, std::map<std::string, std::string> >::iterator _it=ini.find(section);
  if (_it != ini.end())
  {
    // If the value is found
    std::map<std::string, std::string>::iterator it=_it->second.find(key);
    if (it != _it->second.end())
    {
      // If the keyword is found
      // directly returns the value
      return it->second;
    }
    else
      return defaultValue;
  }
  else
    return defaultValue;
}

template void INIParser::SetValue<bool>          (const std::string &, const std::string &, const bool &);
template void INIParser::SetValue<short>         (const std::string &, const std::string &, const short &);
template void INIParser::SetValue<unsigned short>(const std::string &, const std::string &, const unsigned short &);
template void INIParser::SetValue<int>           (const std::string &, const std::string &, const int &);
template void INIParser::SetValue<unsigned int>  (const std::string &, const std::string &, const unsigned int &);
template void INIParser::SetValue<long>          (const std::string &, const std::string &, const long &);
template void INIParser::SetValue<unsigned long> (const std::string &, const std::string &, const unsigned long &);
template void INIParser::SetValue<float>         (const std::string &, const std::string &, const float &);
template void INIParser::SetValue<double>        (const std::string &, const std::string &, const double &);
template void INIParser::SetValue<long double>   (const std::string &, const std::string &, const long double &);
template void INIParser::SetValue<std::string>   (const std::string &, const std::string &, const std::string &);



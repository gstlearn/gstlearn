#pragma once

#include "gstlearn_export.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <map>

/**
 * Class for Managing a Resource File
 *
 * Tis Class is a simple C library offering ini file parsing services.
 * The library is pretty small (less than 1500 lines of C) and robust,
 * and does not depend on any other external library to compile.
 * It is written in ANSI C and should compile on most platforms without difficulty.
 *
 * This Resource File is an ASCII file describing simple parameters
 * (character strings, integers, floating-point values or booleans)
 * in an explicit format, easy to use and modify for users.
 *
 * The resource file is segmented into Sections, declared by the following syntax:
 *    [Section Name]
 *
 * i.e. the section name enclosed in square brackets, alone on a line.
 * Sections names are allowed to contain any character but square brackets or linefeeds.
 *
 * In any section are zero or more variables, declared with the following syntax:
 *     Key = value ; comment
 *
 * The key is any string (possibly containing blanks). The value is any character on the
 * right side of the equal sign. Values can be given enclosed with quotes.
 * If no quotes are present, the value is understood as containing all characters
 * between the first and the last non-blank characters before the comment.
 *
 * The semicolon (;) or Pound sign (#) correspond to the optional delimiter before comments.
 * If there is a comment, it starts from the first character after the delimiter up
 * to the end of the line.
 *
 */
class GSTLEARN_EXPORT INIParser
{
  public:
    // Constructor
    INIParser(bool autoSave=false);
    // Destructor
    ~INIParser();
    // Load the INIparser from a file
    bool load(const std::string &filename);
    // Test if the INI file has data
    bool isEmpty() const  { return ini.empty(); }
    // Get the value and return the default value if the keyword is not found
    template <class T> T GetValue(const std::string &, const std::string &, const T &);
    // Save the value or modify an existing one
    template <class T> void SetValue(const std::string &, const std::string &, const T &);
    // Save the map in the file
    bool save(const std::string &filename="");

  private:
    std::map<std::string, std::map<std::string, std::string> > ini;
    std::string FileName;
    bool AutoSave;
};



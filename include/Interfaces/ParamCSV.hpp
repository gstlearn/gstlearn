#ifndef PARAM_CSV_HPP
#define PARAM_CSV_HPP

#include "Interfaces/interface_d.hpp"

class ParamCSV
{

public:
  ParamCSV(std::string filePath,
           char sep,
           char decimalChar,
           bool useHeader,
           int skipNLines);
  ~ParamCSV();

  std::string getFilePath() const;
  char getSeparatorChar() const;
  char getDecimalChar() const;
  bool getUseHeader() const;
  int getSkipNLines() const;

private:
  std::string _filePath;
  char _separatorChar;
  char _decimalChar;
  bool _useHeader;
  int _skipNLines;
};
#endif

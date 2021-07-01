#ifndef PARAM_CSV_HPP
#define PARAM_CSV_HPP

#include "Interfaces/interface_d.hpp"

class ParamCSV
{

public:
  ParamCSV(std::string filePath,
           std::string sep,
           std::string decimalChar,
           bool useHeader,
           int skipNLines);
  ~ParamCSV();

  std::string getFilePath() const;
  std::string getSeparatorChar() const;
  std::string getDecimalChar() const;
  bool getUseHeader() const;
  int getSkipNLines() const;

private:
  std::string _filePath;
  std::string _separatorChar;
  std::string _decimalChar;
  bool _useHeader;
  int _skipNLines;
};
#endif

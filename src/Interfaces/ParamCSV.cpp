#include "Interfaces/ParamCSV.hpp"

ParamCSV::ParamCSV(std::string filePath,
                   std::string sep,
                   std::string decimalChar,
                   bool useHeader,
                   int skipNLines)
    : _filePath(filePath),
      _separatorChar(sep),
      _decimalChar(decimalChar),
      _useHeader(useHeader),
      _skipNLines(skipNLines)
{
}

ParamCSV::~ParamCSV()
{
}

std::string ParamCSV::getFilePath() const
{
  return (_filePath);
}

std::string ParamCSV::getSeparatorChar() const
{
  return (_separatorChar);
}

std::string ParamCSV::getDecimalChar() const
{
  return (_decimalChar);
}

bool ParamCSV::getUseHeader() const
{
  return (_useHeader);
}

int ParamCSV::getSkipNLines() const
{
  return (_skipNLines);
}

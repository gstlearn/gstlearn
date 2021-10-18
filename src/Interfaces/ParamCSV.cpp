#include "Interfaces/ParamCSV.hpp"

ParamCSV::ParamCSV(std::string filePath,
                   char sep,
                   char decimalChar,
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

char ParamCSV::getSeparatorChar() const
{
  return (_separatorChar);
}

char ParamCSV::getDecimalChar() const
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

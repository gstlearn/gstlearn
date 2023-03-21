/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Basic/CSVformat.hpp"

CSVformat::CSVformat(bool flagHeader,
                     int nSkip,
                     char charSep,
                     char charDec,
                     const String& naString)
  : AStringable(),
    _flagHeader(flagHeader)
  , _nSkip(nSkip)
  , _charSep(charSep)
  , _charDec(charDec)
  , _naString(naString)
{
}

CSVformat::CSVformat(const CSVformat &r)
  : AStringable(r),
    _flagHeader(r._flagHeader)
  , _nSkip(r._nSkip)
  , _charSep(r._charSep)
  , _charDec(r._charDec)
  , _naString(r._naString)
{
}

CSVformat& CSVformat::operator= (const CSVformat &r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _flagHeader = r._flagHeader;
    _nSkip = r._nSkip;
    _charSep = r._charSep;
    _charDec = r._charDec;
    _naString = r._naString;
  }
  return *this;
}

CSVformat::~CSVformat()
{
}

CSVformat* CSVformat::create(bool flagHeader,
                             int nSkip,
                             char charSep,
                             char charDec,
                             const String& naString)
{
  return new CSVformat(flagHeader, nSkip, charSep, charDec, naString);
}

String CSVformat::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << toTitle(1, "CSV Format");

  if (_flagHeader)
    sstr << "The first line contains a Header" << std::endl;
  if (_nSkip > 0)
    sstr << "The first" << _nSkip << "lines should be skipped" << std::endl;

  sstr << "Separator character: '" << _charSep << "'" << std::endl;
  sstr << "Decimal symbol: '" << _charDec << "'" << std::endl;
  sstr << "Missing information string: '" << _naString << "'" << std::endl;
  return sstr.str();
}


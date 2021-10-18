/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "Basic/CSVformat.hpp"

CSVformat::CSVformat(int flagHeader,
                     int nSkip,
                     char charSep,
                     char charDec,
                     const String& naString)
  : _flagHeader(flagHeader)
  , _nSkip(nSkip)
  , _charSep(charSep)
  , _charDec(charDec)
  , _naString(naString)
{
}

CSVformat::CSVformat(const CSVformat &r)
  : _flagHeader(r._flagHeader)
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

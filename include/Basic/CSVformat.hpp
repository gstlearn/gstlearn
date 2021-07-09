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
#pragma once

#include "Basic/String.hpp"

class CSVformat
{
public:
  CSVformat(int flagHeader = true,
            int nSkip = 0,
            const String& charSep = ",",
            const String& charDec = ".",
            const String& naString = "NA");
  CSVformat(const CSVformat &r);
  CSVformat& operator=(const CSVformat &r);
  virtual ~CSVformat();

  const String getCharDec()  const { return _charDec; }
  const String getCharSep()  const { return _charSep; }
  int   getFlagHeader()      const { return _flagHeader; }
  const String getNaString() const { return _naString; }
  int   getNSkip()           const { return _nSkip; }

  void  setFlagHeader(int flagHeader)       { _flagHeader = flagHeader; }
  void  setCharDec(const String& charDec)   { _charDec    = charDec;    }
  void  setCharSep(const String& charSep)   { _charSep    = charSep;    }
  void  setNaString(const String& naString) { _naString   = naString;   }
  void  setNSkip(int nskip)                 { _nSkip      = nskip;       }

private:
  int _flagHeader;
  int _nSkip;
  String _charSep;
  String _charDec;
  String _naString;
};

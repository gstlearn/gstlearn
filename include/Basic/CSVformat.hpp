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

#include "gstlearn_export.hpp"
#include "Basic/String.hpp"

class GSTLEARN_EXPORT CSVformat
{
public:
  CSVformat(int flagHeader = true,
            int nSkip = 0,
            char charSep = ',',
            char charDec = '.',
            const String& naString = "NA");
  CSVformat(const CSVformat &r);
  CSVformat& operator=(const CSVformat &r);
  virtual ~CSVformat();

  char  getCharDec()         const { return _charDec; }
  char  getCharSep()  const { return _charSep; }
  int   getFlagHeader()      const { return _flagHeader; }
  const String getNaString() const { return _naString; }
  int   getNSkip()           const { return _nSkip; }

  void  setFlagHeader(int flagHeader)       { _flagHeader = flagHeader; }
  void  setCharDec(char charDec)            { _charDec    = charDec;    }
  void  setCharSep(char charSep)            { _charSep    = charSep;    }
  void  setNaString(const String& naString) { _naString   = naString;   }
  void  setNSkip(int nskip)                 { _nSkip      = nskip;       }

private:
  int _flagHeader;
  int _nSkip;
  char _charSep;
  char _charDec;
  String _naString;
};

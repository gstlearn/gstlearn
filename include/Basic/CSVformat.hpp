/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AStringable.hpp"

// TODO : Inherits from AParam which inherits from ASerializable, AStringable, IClonable
class GSTLEARN_EXPORT CSVformat: public AStringable
{
public:
  CSVformat(bool flagHeader = true,
            int nSkip = 0,
            char charSep = ',',
            char charDec = '.',
            const String& naString = STRING_NA);
  CSVformat(const CSVformat &r);
  CSVformat& operator=(const CSVformat &r);
  virtual ~CSVformat();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  char  getCharDec()         const { return _charDec; }
  char  getCharSep()         const { return _charSep; }
  bool  getFlagHeader()      const { return _flagHeader; }
  const String getNaString() const { return _naString; }
  int   getNSkip()           const { return _nSkip; }

  void  setFlagHeader(bool flagHeader)      { _flagHeader = flagHeader; }
  void  setCharDec(char charDec)            { _charDec    = charDec;    }
  void  setCharSep(char charSep)            { _charSep    = charSep;    }
  void  setNaString(const String& naString) { _naString   = naString;   }
  void  setNSkip(int nskip)                 { _nSkip      = nskip;      }

  static CSVformat *create(bool flagHeader = true,
                          int nSkip = 0,
                          char charSep = ',',
                          char charDec = '.',
                          const String& naString = STRING_NA);

private:
  bool _flagHeader;
  int _nSkip;
  char _charSep;
  char _charDec;
  String _naString;
};

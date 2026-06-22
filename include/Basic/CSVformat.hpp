/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"
#include "Enum/ECSV.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
  // TODO : Inherits from AParam which inherits from ASerializable, AStringable, IClonable
  class GSTLEARN_EXPORT CSVformat: public AStringable
  {
  public:
    CSVformat(
      bool flagHeader = true,
      Id nSkip = 0,
      char charSep = ',',
      char charDec = '.',
      const String& naString = STRING_NA);
    CSVformat(const CSVformat& r);
    CSVformat& operator=(const CSVformat& r);
    virtual ~CSVformat();

    /// Interface to AStringable
    String toString(const AStringFormat* strfmt = nullptr) const override;

    char getCharDec() const { return _charDec; }

    char getCharSep() const { return _charSep; }

    bool getFlagHeader() const { return _flagHeader; }

    String getNaString() const { return _naString; }

    Id getNSkip() const { return _nSkip; }

    void setFlagHeader(bool flagHeader) { _flagHeader = flagHeader; }

    void setCharDec(char charDec) { _charDec = charDec; }

    void setCharSep(char charSep) { _charSep = charSep; }

    void setNaString(const String& naString) { _naString = naString; }

    void setNSkip(Id nskip) { _nSkip = nskip; }

    static CSVformat* create(
      bool flagHeader = true,
      Id nSkip = 0,
      char charSep = ',',
      char charDec = '.',
      const String& naString = STRING_NA);

    static CSVformat* createStandard(
      const ECSV& csvStyle = ECSV::fromKey("ENGLISH"),
      bool flagHeader = true,
      Id nSkip = 0,
      const String& naString = STRING_NA);

  private:
    bool _flagHeader;
    Id _nSkip;
    char _charSep;
    char _charDec;
    String _naString;
  };
} // namespace gstlrn

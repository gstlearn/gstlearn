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

#include "Basic/AStringable.hpp"
#include "Basic/Vector.hpp"

class GSTLEARN_EXPORT SimuSubstitutionParam: public AStringable
{
public:
  SimuSubstitutionParam();
  SimuSubstitutionParam(const SimuSubstitutionParam &r);
  SimuSubstitutionParam& operator=(const SimuSubstitutionParam &r);
  virtual ~SimuSubstitutionParam();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  const VectorDouble& getColang() const { return _colang; }
  void setColang(const VectorDouble& colang) { _colang = colang; }
  int getColfac() const { return _colfac; }
  void setColfac(int colfac) { _colfac = colfac; }
  double getFactor() const { return _factor; }
  void setFactor(double factor) { _factor = factor; }
  bool isFlagAuto() const { return _flagAuto; }
  void setFlagAuto(bool flagAuto) { _flagAuto = flagAuto; }
  bool isFlagCoding() const { return _flagCoding; }
  void setFlagCoding(bool flagCoding) { _flagCoding = flagCoding; }
  bool isFlagDirect() const { return _flagDirect; }
  void setFlagDirect(bool flagDirect) { _flagDirect = flagDirect; }
  bool isFlagOrient() const { return _flagOrient; }
  void setFlagOrient(bool flagOrient) { _flagOrient = flagOrient; }
  double getIntensity() const { return _intensity; }
  void setIntensity(double intensity) { _intensity = intensity; }
  int getNfacies() const { return _nfacies; }
  void setNfacies(int nfacies) { _nfacies = nfacies; }
  int getNstates() const { return _nstates; }
  void setNstates(int nstates) { _nstates = nstates; }
  const VectorDouble& getTrans() const { return _trans; }
  void setTrans(const VectorDouble& trans) { _trans = trans; }
  const VectorDouble& getVector() const { return _vector; }
  void setVector(const VectorDouble& vector) { _vector = vector; }
  double getColang(int idim) const { return _colang[idim]; }

  bool SimuSubstitutionParam::isValid(bool verbose) const;

private:
  bool _isIrreductibility(bool verbose);

private:
  int _nfacies;
  int _nstates;
  int _colfac;
  bool _flagDirect;
  bool _flagCoding;
  bool _flagOrient;
  bool _flagAuto;
  double _intensity;
  double _factor;
  VectorDouble _vector;
  VectorDouble _colang;
  VectorDouble _trans;
};

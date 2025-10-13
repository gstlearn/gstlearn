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

#include "gstlearn_export.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"

namespace gstlrn
{
class GSTLEARN_EXPORT SimuSubstitutionParam: public AStringable
{
public:
  SimuSubstitutionParam(Id nfacies       = 2,
                        double intensity = 0.1,
                        bool flag_direct = true,
                        bool flag_coding = true,
                        bool flag_orient = false);
  SimuSubstitutionParam(const SimuSubstitutionParam& r);
  SimuSubstitutionParam& operator=(const SimuSubstitutionParam& r);
  virtual ~SimuSubstitutionParam();

  /// Interface to AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  const VectorInt& getColang() const { return _colang; }
  void setColang(const VectorInt& colang) { _colang = colang; }
  Id getColfac() const { return _colfac; }
  void setColfac(Id colfac) { _colfac = colfac; }
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
  Id getNfacies() const { return _nfacies; }
  void setNfacies(Id nfacies) { _nfacies = nfacies; }
  Id getNstates() const { return _nstates; }
  void setNstates(Id nstates) { _nstates = nstates; }
  VectorDouble getTrans() const { return _trans; }
  void setTrans(const VectorDouble& trans) { _trans = trans; }
  const VectorDouble& getVector() const { return _vector; }
  void setVector(const VectorDouble& vector) { _vector = vector; }
  Id getColang(Id idim) const;
  double getVector(Id idim) const { return _vector[idim]; }

  bool isValid(bool verbose = false);
  static void isValidOrientation(VectorDouble& vector, bool verbose = false);
  static void isValidFactor(double* factor, bool verbose = false);

  bool isAngleLocal() const;
  bool isLocal() const;

private:
  bool _isIrreductibility(bool verbose = false);
  bool _isValidTransition(bool verbose = false, double eps = EPSILON3);

private:
  Id _nfacies;
  Id _nstates;
  Id _colfac;
  bool _flagDirect;
  bool _flagCoding;
  bool _flagOrient;
  bool _flagAuto;
  double _intensity;
  double _factor;
  VectorInt _colang;
  VectorDouble _vector;
  VectorDouble _trans;
};
} // namespace gstlrn
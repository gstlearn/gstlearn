/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Enum/ELaw.hpp"
#include "Basic/AStringable.hpp"

// TODO Will be replaced by future class"Law" or "Distribution" which does not
// actually exist
class GSTLEARN_EXPORT ShapeParameter: public AStringable
{
public:
  ShapeParameter(ELaw law = ELaw::fromKey("CONSTANT"), double value = 0.);
  ShapeParameter(const ShapeParameter &r);
  ShapeParameter& operator=(const ShapeParameter &r);
  virtual ~ShapeParameter();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  ELaw getLaw() const { return _law; }
  const VectorDouble& getValarg() const { return _valarg; }
  double getValarg(int iarg) const;
  int getNbValarg() const { return (int) _valarg.size(); }

  void setLaw(ELaw law) { _law = law; }
  void setValarg(int iarg, double value);

  double generateValue() const;

private:
  bool _isValidArgIndex(int iarg) const;

private:
  ELaw _law; /* Type of law */
  VectorDouble _valarg; /* Randomization arguments */
};

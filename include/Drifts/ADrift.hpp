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
#include "Basic/ICloneable.hpp"
#include "Basic/VectorNumT.hpp"

/* Elementary Drift function
 * */

class Db;

class GSTLEARN_EXPORT ADrift : public AStringable, public ICloneable
{
public:
  ADrift();
  ADrift(const ADrift &r);
  ADrift& operator= (const ADrift &r);
  virtual ~ADrift();

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface for daughter classes
  virtual VectorInt getPowers() const { return VectorInt(); }
  virtual int       getRankFex() const { return 0; }

  // ADriftelem Interface
  virtual String getDriftName() const = 0;
  virtual int    getOrderIRF() const = 0;
  virtual int    getOrderIRFIdim(int idim) const = 0;
  virtual double eval(const Db* db,int iech) const = 0;
  virtual int    getDriftNDimMax() const { return 0; }
  virtual bool   isDriftExternal() const { return false; }

};
